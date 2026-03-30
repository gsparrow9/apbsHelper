import os
import numpy as np
import MDAnalysis as mda
from pathlib import Path
from scipy.spatial.transform import Rotation
from src.classes.handler_functions import HandlerFunction
from src.utils import apbs_tool as apbs_tool
from src.utils import tools as tools
from typing import Dict, Tuple, List

def coarse_grain_rna(rna_chain: mda.Universe,rna_chain_id: List[str], charge: float) -> tuple:

        # Calculate CG Centers
        base_mean = 5
        total_residues = len(rna_chain.residues)//2 # type: ignore
        rna_cg_centers = np.zeros((int(total_residues/base_mean),3))
        rna_cg_charge = np.ones((int(total_residues/base_mean),1))*charge

        for i in range(1,total_residues,base_mean):
            i_n = i + base_mean - 1
            j_0 = 40 - i_n + 1
            j_n = 40 - i + 1
            select_str = f'(chainID {rna_chain_id[0]} and resid {i}:{i_n}) or (chainID {rna_chain_id[1]} and resid {j_0}:{j_n})'
            selection = rna_chain.select_atoms(select_str)
            cg_center = selection.positions.mean(axis=0)
            rna_cg_centers[int((i-1)/base_mean),:] = cg_center

            left_extremum = i == 1
            right_extremum = i == (total_residues - (base_mean-1))
            
            if left_extremum and right_extremum:
                rna_cg_charge[int((i-1)/base_mean)] = charge + 2
            elif left_extremum or right_extremum:
                rna_cg_charge[int((i-1)/base_mean)] = charge + 1
            else:
                rna_cg_charge[int((i-1)/base_mean)] = charge
            rna_cg_centers[int((i-1)/base_mean),:] = cg_center
        return rna_cg_centers, rna_cg_charge

class many_to_many(HandlerFunction):
    def handle(
            self,
            core_pdb: Path | str,
            target_dir: Path | str,
            handler_options: Dict,
            coarse_grain: bool = False
    ) -> Tuple[Path, Dict[str, List[str]]]:
        
        # Read handler options. See RNA Handler Notes for angles.
        chains = {'left':['D','F'],'right':['C','E']}
        distance = handler_options['distance']
        omega_angle = handler_options['rotation_omega']
        theta_angle = handler_options['rotation_theta']
        core_pdb = Path(core_pdb)
        target_dir = Path(target_dir)
        complex_pdb = target_dir / core_pdb.name
        left_pdb = complex_pdb.with_stem(f'left_{core_pdb.stem}')
        right_pdb = complex_pdb.with_stem(f'right_{core_pdb.stem}')

        # Write new Left/Right PDB
        complex_rna = mda.Universe(core_pdb)
        left = complex_rna.select_atoms('chainID D or chainID F')
        right = complex_rna.select_atoms('chainID C or chainID E')

        left_com = left.center_of_mass()
        right_com = right.center_of_mass()
        center_vector = right_com - left_com
        u_center_vector = center_vector/np.linalg.norm(center_vector)
        ones_matrix = np.ones(right.atoms.positions.shape)

        # Only move Right chain
        if distance != 0:
             #CoM Distance
             displacement = distance*u_center_vector*ones_matrix
             right.atoms.positions = right.atoms.positions - center_vector + displacement
        
        # RNA Rotations
        omega = np.deg2rad(omega_angle)
        theta = np.deg2rad(theta_angle)
        center = right.atoms.center_of_mass() - center_vector/2
        principal_axes = right.atoms.principal_axes(wrap=False)
        theta_axis = principal_axes[1]
        omega_axis = principal_axes[2]
        coords = right.atoms.positions - center

        # Apply omega rotation
        rotated = Rotation.from_rotvec(omega*omega_axis).apply(coords)

        # Apply theta rotation
        rotated = Rotation.from_rotvec(theta*theta_axis).apply(rotated)
        right.atoms.positions = rotated + center

        # Merge and write pdbs
        u = mda.Merge(left.atoms,right.atoms)
        u.atoms.write(target_dir / core_pdb.name) # type: ignore
        left.atoms.write(target_dir / left_pdb.name)
        right.atoms.write(target_dir / right_pdb.name)

        # Generate pqr files
        apbs_tool.generate_pqr_files(
             target_dir,
             '/home/chris/Software/apbs-pdb2pqr/pdb2pqr/pdb2pqr.py'
        )
        
        if coarse_grain:
            charge = handler_options['cg_charges'][0]
            cg_geometry = handler_options['cg_geometry']
            complex_pqr = complex_pdb.with_suffix('.pqr')

            # Calculate CG centers
            left_cg_pos,left_cg_chg = coarse_grain_rna(
                 rna_chain=left,
                 rna_chain_id=['D','F'],
                 charge=charge
            )
            right_cg_pos,right_cg_chg = coarse_grain_rna(
                 rna_chain=right,
                 rna_chain_id=['E','C'],
                 charge=charge
            )

            cg_pos = {'D':left_cg_pos,'C':right_cg_pos}
            cg_chg = {'D':left_cg_chg,'C':right_cg_chg}

            if cg_geometry == 'FullGeometry': # All-Atom geometry w/ cg without radius
                tools.create_pqr_cg(
                     str(complex_pqr),
                     cg_pos,
                     cg_chg,
                     cg_radius=0.0
                )
            elif cg_geometry == 'SimpleGeometry':
                tools.create_pqr_cg(
                     str(complex_pqr),
                     cg_pos,
                     cg_chg,
                     supress_radius=True,
                     cg_radius=11.185
                )
            else:
                 raise RuntimeError("cg_geometry must be 'FullGeometry' or 'SimpleGeometry'")
        return complex_pdb,chains


class one_to_many(HandlerFunction):
    def handle(
            self,
            core_pdb: Path | str,
            target_dir: Path | str,
            handler_options: Dict,
            coarse_grain: bool = False
    ) -> Tuple[Path, Dict[str, List[str]]]:
        
        # Read handler options. See RNA Handler Notes for angles.
        chains = {'left':['D','F'],'right':['C','E']}
        distance = handler_options['distance']
        omega_angle = handler_options['rotation_omega']
        theta_angle = handler_options['rotation_theta']
        core_pdb = Path(core_pdb)
        target_dir = Path(target_dir)
        complex_pdb = target_dir / core_pdb.name
        left_pdb = complex_pdb.with_stem(f'left_{core_pdb.stem}')
        right_pdb = complex_pdb.with_stem(f'right_{core_pdb.stem}')

        # Read new Left/Right PDB
        complex_rna = mda.Universe(core_pdb)
        left = complex_rna.select_atoms('chainID D or chainID F')
        right = complex_rna.select_atoms('chainID C or chainID E')

        left_com = left.center_of_mass()
        right_com = right.center_of_mass()
        center_vector = right_com - left_com
        u_center_vector = center_vector/np.linalg.norm(center_vector)
        ones_matrix = np.ones(right.atoms.positions.shape)

        # Cut Left Chain to 1-cg
        select_str = f'(chainID D and resid 1:5) or (chainID F and resid 35:40)'
        left = left.select_atoms(select_str)

        # Only move Right chain
        if distance != 0:
             #CoM Distance
             displacement = distance*u_center_vector*ones_matrix
             right.atoms.positions = right.atoms.positions - center_vector + displacement

        # RNA Rotations
        omega = np.deg2rad(omega_angle)
        theta = np.deg2rad(theta_angle)
        center = right.atoms.center_of_mass() - center_vector/2
        principal_axes = right.atoms.principal_axes(wrap=False)
        theta_axis = principal_axes[1]
        omega_axis = principal_axes[2]
        coords = right.atoms.positions - center
        
        # Apply omega rotation
        rotated = Rotation.from_rotvec(omega*omega_axis).apply(coords)

        # Apply theta rotation
        rotated = Rotation.from_rotvec(theta*theta_axis).apply(rotated)
        right.atoms.positions = rotated + center

        # Merge and write pdbs
        u = mda.Merge(left.atoms,right.atoms)
        u.atoms.write(target_dir / core_pdb.name) # type: ignore
        left.atoms.write(target_dir / left_pdb.name)
        right.atoms.write(target_dir / right_pdb.name)

        # Generate pqr files
        apbs_tool.generate_pqr_files(
             target_dir,
             '/home/chris/Software/apbs-pdb2pqr/pdb2pqr/pdb2pqr.py'
        )
        
        if coarse_grain:
            charge = handler_options['cg_charges'][0]
            cg_geometry = handler_options['cg_geometry']
            complex_pqr = complex_pdb.with_suffix('.pqr')

            # Calculate CG centers
            left_cg_pos,left_cg_chg = coarse_grain_rna(
                 rna_chain=left,
                 rna_chain_id=['D','F'],
                 charge=charge
            )
            right_cg_pos,right_cg_chg = coarse_grain_rna(
                 rna_chain=right,
                 rna_chain_id=['E','C'],
                 charge=charge
            )

            cg_pos = {'D':left_cg_pos,'C':right_cg_pos}
            cg_chg = {'D':left_cg_chg,'C':right_cg_chg}

            if cg_geometry == 'FullGeometry': # All-Atom geometry w/ cg without radius
                tools.create_pqr_cg(
                     str(complex_pqr),
                     cg_pos,
                     cg_chg,
                     cg_radius=0.0
                )
            elif cg_geometry == 'SimpleGeometry':
                tools.create_pqr_cg(
                     str(complex_pqr),
                     cg_pos,
                     cg_chg,
                     supress_radius=True,
                     cg_radius=11.185
                )
            else:
                 raise RuntimeError("cg_geometry must be 'FullGeometry' or 'SimpleGeometry'")
        return complex_pdb,chains