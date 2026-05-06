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

          # Define new Left/Right PDB
          complex_rna = mda.Universe(core_pdb)
          left = complex_rna.select_atoms('chainID D or chainID F')
          right = complex_rna.select_atoms('chainID C or chainID E')

          right_centroid = right.atoms.positions.mean(axis=0)
          left_centroid = left.atoms.positions.mean(axis=0)

          # Get phosphate backbone axes and define local coordinate system
          def get_right_select_str(chain_length):
               second_resid = 40 - (chain_length-1)
               start_str = f'(chainID E and resid {chain_length-1}:{chain_length})'
               start_str += f' or (chainID C and resid {second_resid}:{second_resid+1}) and name P'
               return start_str

          def get_left_select_str(chain_length):
               second_resid = 40 - (chain_length-1)
               start_str = f'(chainID D and resid {chain_length-1}:{chain_length})'
               start_str += f' or (chainID F and resid {second_resid}:{second_resid+1}) and name P'
               return start_str

          right_chain_length = len(right.residues)//2
          left_chain_length = len(left.residues)//2
          right_start_str = get_right_select_str(right_chain_length)
          right_end_str = '(chainID E and resid 1:2) or (chainID C and resid 39:40) and name P'
          left_start_str = get_left_select_str(left_chain_length)
          left_end_str = '(chainID D and resid 1:2) or (chainID F and resid 39:40) and name P'

          right_backbone_axis = (right.select_atoms(right_end_str).
                                 positions.mean(axis=0) - 
                                 right.select_atoms(right_start_str).
                                 positions.mean(axis=0)
                                 )
          left_backbone_axis = (left.select_atoms(left_end_str).
                                 positions.mean(axis=0) - 
                                 left.select_atoms(left_start_str).
                                 positions.mean(axis=0)
                                 )
          right_backbone_vector = right_backbone_axis/np.linalg.norm(right_backbone_axis)
          right_x_axis = np.linalg.cross([0,1,0],right_backbone_vector)
          right_x_axis /= np.linalg.norm(right_x_axis)
          right_y_axis = np.linalg.cross(right_backbone_vector,right_x_axis)
          right_y_axis /= np.linalg.norm(right_y_axis)
          right_rotation = np.array([right_x_axis,right_y_axis,right_backbone_vector])

          left_backbone_vector = left_backbone_axis/np.linalg.norm(left_backbone_axis)
          left_x_axis = np.linalg.cross([0,1,0],left_backbone_vector)
          left_x_axis /= np.linalg.norm(left_x_axis)
          left_y_axis = np.linalg.cross(left_backbone_vector,left_x_axis)
          left_y_axis /= np.linalg.norm(left_y_axis)
          left_rotation = np.array([left_x_axis,left_y_axis,left_backbone_vector])

          # RNA orient local CS to global CS
          right_coords = right.atoms.positions - right_centroid
          right_coords = right_coords @ right_rotation.T + right_centroid

          left_coords = left.atoms.positions - left_centroid
          left_coords = left_coords @ left_rotation.T + left_centroid
          right.atoms.positions = right_coords
          left.atoms.positions = left_coords

          # Align centroids on the ZY plane
          right_centroid = right.atoms.positions.mean(axis=0)
          left_centroid = left.atoms.positions.mean(axis=0)

          right.atoms.positions = (right.atoms.positions - 
                                   np.array([0,right_centroid[1],right_centroid[-1]]))
          left.atoms.positions = (left.atoms.positions - 
                                   np.array([0,left_centroid[1],left_centroid[-1]]))

          # RNA Rotations
          x_axis = np.array([1,0,0])
          omega_axis = np.array([0,0,1])
          theta_axis = x_axis
          right_centroid = right.positions.mean(axis=0)
          right_centered_pos = right.atoms.positions - right_centroid
          omega = np.deg2rad(omega_angle)
          theta = np.deg2rad(theta_angle)

          # Apply omega rotation
          rotated = Rotation.from_rotvec(omega*omega_axis).apply(right_centered_pos)

          # Apply theta rotation
          rotated = Rotation.from_rotvec(theta*theta_axis).apply(rotated)
          right.atoms.positions = rotated + right_centroid

          # Only move Right chain
          ones_matrix = np.ones(right.atoms.positions.shape)
          left_centroid = left.positions.mean(axis=0)
          right_centroid = right.positions.mean(axis=0)
          centroid_vector = right_centroid - left_centroid
          if distance != 0:
               displacement = distance*x_axis*ones_matrix
               right.atoms.positions = right.atoms.positions + displacement - centroid_vector

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

          # Define new Left/Right PDB
          complex_rna = mda.Universe(core_pdb)
          left = complex_rna.select_atoms('chainID D or chainID F')
          right = complex_rna.select_atoms('chainID C or chainID E')

          # Cut Left Chain to 1-cg bulk sphere
          select_str = f'(chainID D and resid 1:5) or (chainID F and resid 36:40)'
          left = left.select_atoms(select_str)

          right_centroid = right.atoms.positions.mean(axis=0)
          left_centroid = left.atoms.positions.mean(axis=0)

          # Get phosphate backbone axes and define local coordinate system
          def get_right_select_str(chain_length):
               second_resid = 40 - (chain_length-1)
               start_str = f'(chainID E and resid {chain_length-1}:{chain_length})'
               start_str += f' or (chainID C and resid {second_resid}:{second_resid+1}) and name P'
               return start_str

          right_chain_length = len(right.residues)//2
          right_start_str = get_right_select_str(right_chain_length)
          right_end_str = '(chainID E and resid 1:2) or (chainID C and resid 39:40) and name P'
          left_start_str = '(chainID F and resid 36:37) or (chainID D and resid 4:5) and name P'
          left_end_str = '(chainID F and resid 39:40) or (chainID D and resid 1:2) and name P'

          right_backbone_axis = (right.select_atoms(right_end_str).
                                 positions.mean(axis=0) - 
                                 right.select_atoms(right_start_str).
                                 positions.mean(axis=0)
                                 )
          left_backbone_axis = (left.select_atoms(left_end_str).
                                 positions.mean(axis=0) - 
                                 left.select_atoms(left_start_str).
                                 positions.mean(axis=0)
                                 )
          right_backbone_vector = right_backbone_axis/np.linalg.norm(right_backbone_axis)
          right_x_axis = np.linalg.cross([0,1,0],right_backbone_vector)
          right_x_axis /= np.linalg.norm(right_x_axis)
          right_y_axis = np.linalg.cross(right_backbone_vector,right_x_axis)
          right_y_axis /= np.linalg.norm(right_y_axis)
          right_rotation = np.array([right_x_axis,right_y_axis,right_backbone_vector])

          left_backbone_vector = left_backbone_axis/np.linalg.norm(left_backbone_axis)
          left_x_axis = np.linalg.cross([0,1,0],left_backbone_vector)
          left_x_axis /= np.linalg.norm(left_x_axis)
          left_y_axis = np.linalg.cross(left_backbone_vector,left_x_axis)
          left_y_axis /= np.linalg.norm(left_y_axis)
          left_rotation = np.array([left_x_axis,left_y_axis,left_backbone_vector])

          # RNA orient local CS to global CS
          right_coords = right.atoms.positions - right_centroid
          right_coords = right_coords @ right_rotation.T + right_centroid

          left_coords = left.atoms.positions - left_centroid
          left_coords = left_coords @ left_rotation.T + left_centroid
          right.atoms.positions = right_coords
          left.atoms.positions = left_coords

          # Align centroids on the ZY plane
          right_centroid = right.atoms.positions.mean(axis=0)
          left_centroid = left.atoms.positions.mean(axis=0)

          right.atoms.positions = (right.atoms.positions - 
                                   np.array([0,right_centroid[1],right_centroid[-1]]))
          left.atoms.positions = (left.atoms.positions - 
                                   np.array([0,left_centroid[1],left_centroid[-1]]))

          # RNA Rotations
          x_axis = np.array([1,0,0])
          omega_axis = np.array([0,0,1])
          theta_axis = x_axis
          right_centroid = right.positions.mean(axis=0)
          right_centered_pos = right.atoms.positions - right_centroid
          omega = np.deg2rad(omega_angle)
          theta = np.deg2rad(theta_angle)

          # Apply omega rotation
          rotated = Rotation.from_rotvec(omega*omega_axis).apply(right_centered_pos)

          # Apply theta rotation
          rotated = Rotation.from_rotvec(theta*theta_axis).apply(rotated)
          right.atoms.positions = rotated + right_centroid

          # Only move Right chain
          ones_matrix = np.ones(right.atoms.positions.shape)
          left_centroid = left.positions.mean(axis=0)
          right_centroid = right.positions.mean(axis=0)
          centroid_vector = right_centroid - left_centroid
          if distance != 0:
               displacement = distance*x_axis*ones_matrix
               right.atoms.positions = right.atoms.positions + displacement - centroid_vector

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