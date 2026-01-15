import MDAnalysis as mda
import os
import numpy as np
from scipy.spatial.transform import Rotation
from typing import Dict, Tuple, List
from src.utils import apbs_tool as apbs_tool
from src.utils import tools as tools
from pathlib import Path
from src.classes.handler_functions import HandlerFunction

class TestHandler(HandlerFunction):

    def handle(
            self,
            core_pdb: str,
            target_dir: str,
            handler_options: Dict,
            coarse_grain: bool=False
        ) -> Tuple[Path, Dict[str,List[str]]]:

        rot_angle = handler_options['rotation_angle']
        distance = handler_options['distance']
        pdb_name = os.path.split(core_pdb)[-1]
        pdb_name,_ = os.path.splitext(pdb_name)
        out_pdb_filename = f'{pdb_name}.pdb'
        left_out_pdb_filename = f'left_{pdb_name}.pdb'
        right_out_pdb_filename = f'right_{pdb_name}.pdb'

        #Crea PDB nuevo con proteintas separadas.
        pdb = mda.Universe(core_pdb)
        left = pdb.select_atoms('chainID D or chainID F')
        right = pdb.select_atoms('chainID C or chainID E')
        left_center = left.center_of_mass()
        right_center =  right.center_of_mass()
        center_vector = right_center - left_center
        u_center_vector = center_vector/np.linalg.norm(center_vector)
        foo = np.ones(right.atoms.positions.shape)

        #Solo se separará la proteina right. Left queda fija. Se separa en el eje que conecta CC.MM.
        if distance != 0:
            #Dist entre centros de masa.
            right.atoms.positions = right.atoms.positions - center_vector + distance*u_center_vector*foo

        #Rota el ARN segun el angulo del caso.
        theta = np.deg2rad(rot_angle)
        center = right.atoms.center_of_geometry()
        eje_z = right.atoms.principal_axes(wrap=False)[2]
        coords = right.atoms.positions -  center
        rotated = Rotation.from_rotvec(theta*eje_z).apply(coords)
        right.atoms.positions = rotated + center

        u = mda.Merge(left.atoms,right.atoms)
        u.atoms.write(os.path.join(target_dir, out_pdb_filename)) # type: ignore
        left.atoms.write(os.path.join(target_dir, left_out_pdb_filename))
        right.atoms.write(os.path.join(target_dir, right_out_pdb_filename))
        apbs_tool.generate_pqr_files(target_dir,'/home/chris/Software/apbs-pdb2pqr/pdb2pqr/pdb2pqr.py')
        chains = {'left':['D','F'],'right':['C','E']}

        if coarse_grain:
            charge = handler_options['cg_charges'][0]
            pqr_com = os.path.splitext(out_pdb_filename)[0] + '.pqr'
            left_cg_pos,left_cg_chg = self.coarse_grain_rna({'chain':left,'id':['D','F']},charge)
            right_cg_pos,right_cg_chg = self.coarse_grain_rna({'chain':right,'id':['E','C']},charge)
            cg_pos = {'D':left_cg_pos,'C':right_cg_pos}
            cg_chg = {'D':left_cg_chg,'C':right_cg_chg}

            tools.create_pqr_cg(
                os.path.join(target_dir,pqr_com),
                cg_pos,
                cg_chg
                )

        return Path(out_pdb_filename),chains

    def coarse_grain_rna(self,rna_chain: Dict, charge: float) -> tuple:
        rna_chain_universe = rna_chain['chain']
        rna_chain_id = rna_chain['id']

        #Calcula centros para CG
        base_mean = 5
        total_residues = len(rna_chain_universe.residues)//2
        rna_cg_centers = np.zeros((int(total_residues/base_mean),3))
        rna_cg_charge = np.ones((int(total_residues/base_mean),1))*charge

        for i in range(1,total_residues,base_mean):
            i_n = i + base_mean - 1
            j_0 = 40 - i_n + 1
            j_n = 40 - i + 1
            select_str = f'(chainID {rna_chain_id[0]} and resid {i}:{i_n}) or (chainID {rna_chain_id[1]} and resid {j_0}:{j_n})'
            selection = rna_chain_universe.select_atoms(select_str)
            cg_center = selection.positions.mean(axis=0)
            rna_cg_centers[int((i-1)/base_mean),:] = cg_center

        return rna_cg_centers, rna_cg_charge