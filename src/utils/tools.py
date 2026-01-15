import numpy as np
import MDAnalysis as mda
from scipy.spatial import KDTree
import pandas as pd
from typing import Dict, List, Tuple
from pathlib import Path
import subprocess,os,glob,time
from src.utils import apbs_tool as apbs_tool

def def_mesh(pqr_file: str, spacing: float) -> List[List]: #APBS Tool

    #Define mesh size
    length_x = 0.
    length_y = 0.
    length_z = 0.

    pqr_data = mda.Universe(pqr_file)
    pos = pqr_data.coord.positions
    xmin,xmax = np.min(pos[:,0]), np.max(pos[:,0])
    ymin,ymax = np.min(pos[:,1]), np.max(pos[:,1])
    zmin,zmax = np.min(pos[:,2]), np.max(pos[:,2])

    length_x = max(length_x, xmax-xmin)
    length_y = max(length_y, ymax-ymin)
    length_z = max(length_z, zmax-zmin)

    tol = 15 #Ideal clearance 10 - 20 A.
    length_x += 2*tol
    length_y += 2*tol
    length_z += 2*tol
    
    possible_dime = np.arange(1,36,1)*2**(4+1)+1
    dx = length_x/(possible_dime-1)
    dy = length_y/(possible_dime-1)
    dz = length_z/(possible_dime-1)

    index_dime = np.argwhere(dx<spacing)[0][0]
    dime_x = possible_dime[index_dime]
    index_dime = np.argwhere(dy<spacing)[0][0]
    dime_y = possible_dime[index_dime]
    index_dime = np.argwhere(dz<spacing)[0][0]
    dime_z = possible_dime[index_dime]

    clearance_x = (dime_x*spacing - (length_x - 2*tol))/2
    clearance_y = (dime_y*spacing - (length_y - 2*tol))/2
    clearance_z = (dime_z*spacing - (length_z - 2*tol))/2

    dime = [dime_x,dime_y,dime_z]
    bbox = [xmax-xmin,ymax-ymin,zmax-zmin]
    clearance = [clearance_x,clearance_y,clearance_z]

    return [dime,bbox,clearance] # type: ignore

def gen_apbs(pqr_file: str, spacing: float, keep_dx: bool=False, linear: bool=False) -> Tuple[List,List]:
    in_file =""" 
    read
        mol pqr {{PQR_FILE}}
    end
    elec name solv 
        mg-manual
        dime {{DIME_X}} {{DIME_Y}} {{DIME_Z}}
        grid {{SPACING}} {{SPACING}} {{SPACING}}
        gcent mol 1 
        mol 1 
        {{PB_MODE}}
        bcfl mdh 
        ion charge 1 conc 0.150 radius 2.0
        ion charge -1 conc 0.150 radius 2.0 
        pdie 2.0000
        sdie 78.5400
        srfm smol 
        chgm spl2 
        sdens 10.00 
        srad 1.4
        swin 0.30
        temp 298.15
        calcenergy total
        calcforce no
        {{WRITE_SOLV_DX}}
    end

    elec name vacc 
        mg-manual
        dime {{DIME_X}} {{DIME_Y}} {{DIME_Z}}
        grid {{SPACING}} {{SPACING}} {{SPACING}}
        gcent mol 1 
        mol 1
        {{PB_MODE}}
        bcfl mdh
        pdie 2.0000
        sdie 2.0000
        srfm smol
        chgm spl2
        sdens 10.00
        srad 1.4
        swin 0.30
        temp 298.15
        calcenergy total
        calcforce no
        {{WRITE_VAC_DX}}
    end

    print elecEnergy solv - vacc end
    quit"""

    working_dir, pqr_filename = os.path.split(pqr_file)
    pqr_name,_ = os.path.splitext(pqr_filename)
    dim,bbox,clearance = def_mesh(pqr_file, spacing)

    in_file = in_file.replace("{{DIME_X}}",str(dim[0]))
    in_file = in_file.replace("{{DIME_Y}}",str(dim[1])) # type: ignore
    in_file = in_file.replace("{{DIME_Z}}",str(dim[2])) # type: ignore
    in_file = in_file.replace("{{SPACING}}",str(spacing))
    in_file = in_file.replace("{{PQR_FILE}}",pqr_filename)
    
    if linear:
        in_file = in_file.replace("{{PB_MODE}}",'lpbe')
    else:
        in_file = in_file.replace("{{PB_MODE}}",'npbe')

    if keep_dx:
        write_dx = "write pot dx {{DX_LOC}}"
        dx_solv_path = os.path.join(pqr_name+'_solv')
        dx_vacc_path = os.path.join(pqr_name+'_vacc')
        solv_dx = write_dx.replace("{{DX_LOC}}",dx_solv_path)
        vacc_dx = write_dx.replace("{{DX_LOC}}",dx_vacc_path)
        in_file = in_file.replace("{{WRITE_VAC_DX}}",vacc_dx)
        in_file = in_file.replace("{{WRITE_SOLV_DX}}",solv_dx)
    else:
        write_dx = "#write pot dx {{DX_LOC}}"
        dx_solv_path = os.path.join(pqr_name+'_solv')
        dx_vacc_path = os.path.join(pqr_name+'_vacc')
        solv_dx = write_dx.replace("{{DX_LOC}}",dx_solv_path)
        vacc_dx = write_dx.replace("{{DX_LOC}}",dx_vacc_path)
        in_file = in_file.replace("{{WRITE_VAC_DX}}",vacc_dx)
        in_file = in_file.replace("{{WRITE_SOLV_DX}}",solv_dx)

    with open(os.path.join(working_dir,pqr_name + '.in'),'w') as f:
        f.write(in_file)

    bbox[0] = float(f'{bbox[0]:.2f}')
    bbox[1] = float(f'{bbox[1]:.2f}')
    bbox[2] = float(f'{bbox[2]:.2f}')

    clearance[0] = float(f'{clearance[0]:.2f}')
    clearance[1] = float(f'{clearance[1]:.2f}')
    clearance[2] = float(f'{clearance[2]:.2f}')
    
    return bbox,clearance # type: ignore

def create_xyzr(pqr_file: str) -> str:
    working_dir, pqr_filename = os.path.split(pqr_file)
    basename,_ = os.path.splitext(pqr_filename)
    xyzr_file = open(os.path.join(working_dir, basename + '.xyzr'),'w')
    #subprocess.run(f"awk '/^ATOM/ {{print $6,$7,$8,$10}}' {pqr_file}",shell=True,stdout=xyzr_file)
    subprocess.run(
        f"""awk '/^ATOM/ {{
            x = substr($0,33,8);
            y = substr($0,42,8);
            z = substr($0,51,8);
            radius = substr($0,68,6);
            printf "%s %s %s %s\\n", x, y, z, radius
        }}' {pqr_filename}""",
        shell=True,
        stdout=xyzr_file,
        cwd=working_dir
    ) 
    
    xyzr_file.close()
    return os.path.join(working_dir,pqr_filename)+'.xyzr'

def run_nanoshaper(xyzr_file: str):
    working_dir, xyzr_filename = os.path.split(xyzr_file)
    xyzr_name,_ = os.path.splitext(xyzr_filename)
    config = """Grid_scale = 2.0
    Grid_perfil = 90.0
    XYZR_FileName = {{XYZR_FILE}}
    Build_epsilon_maps  = false
    Build_status_map = false
    Surface = skin
    Smooth_Mesh = true
    Number_thread = 32 #default value is 32
    Skin_Surface_Parameter = 0.45
    Blobbyness = -2.5 
    Surface_File_Name = {{OFF_FILE}}
    Cavity_Detection_Filling = false
    Conditional_Volume_Filling_Value = 11.4
    Probe_Radius = 1.4
    Accurate_Triangulation = true
    Triangulation = true
    Check_duplicated_vertices = true
    Save_Status_map = false
    Save_PovRay = false
    Max_mesh_auxiliary_grid_size = 100
    Max_mesh_patches_per_auxiliary_grid_cell = 250;
    Max_mesh_auxiliary_grid_2d_size = 100
    Max_mesh_patches_per_auxiliary_grid_2d_cell = 250 
    Max_ses_patches_auxiliary_grid_size = 100
    Max_ses_patches_per_auxiliary_grid_cell = 400
    Max_ses_patches_auxiliary_grid_2d_size = 50
    Max_ses_patches_per_auxiliary_grid_2d_cell = 400
    Max_skin_patches_auxiliary_grid_size = 100
    Max_skin_patches_per_auxiliary_grid_cell = 400
    Max_skin_patches_auxiliary_grid_2d_size = 150
    Max_skin_patches_per_auxiliary_grid_2d_cell = 200"""

    config = config.replace('{{XYZR_FILE}}',xyzr_filename)
    config = config.replace('{{OFF_FILE}}',xyzr_name+'.off')
    config_file = os.path.join(working_dir,'conf.prm')
    with open(config_file,'w') as conf:
        conf.write(config)

    subprocess.run(['/home/chris/Software/nanoshaper/NanoShaper','conf.prm'],
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL,
                   cwd=working_dir)

    with open(os.path.join(working_dir,xyzr_name) + '.off','w') as out:
        subprocess.run(['awk',r'NR > 4 && NF == 3 {print $1,$2,$3}','triangulatedSurf.off'],
                       stdout=out,
                       cwd=working_dir)

    polygon_data = np.loadtxt(os.path.join(working_dir, xyzr_name + '.off'), delimiter=' ')

    #Cleanup
    if os.path.exists(os.path.join(working_dir,'triangleAreas.txt')):
        os.remove(os.path.join(working_dir,'triangleAreas.txt'))
        os.remove(os.path.join(working_dir,'triangulatedSurf.off'))
        os.remove(os.path.join(working_dir,'stderror.txt'))
        os.remove(os.path.join(working_dir,'conf.prm'))
        os.remove(os.path.join(working_dir, xyzr_name + '.off'))
    return polygon_data

def get_wall_distance(xyzr_file_left, xyzr_file_right) -> float:
    left = run_nanoshaper(xyzr_file_left)
    right = run_nanoshaper(xyzr_file_right)
    right = KDTree(right)
    wall_dist = right.query(left,k=1,workers=5)[0].min()
    return wall_dist

def calc_coulomb(pqr_file: str|Path, epsilon_r:float=2) -> float:
    """
	Calcula la energía de Coulomb para un archivo PQR.
    """
    e = 1.60217663E-19 #C
    epsilon_0 = 8.8541878188E-12 #F/m 
    NA = 6.02214076E23
    pqr = mda.Universe(pqr_file)
    charges = pqr.atoms.charges  # type: ignore 
    coords = pqr.atoms.positions # type: ignore
    r_ij = np.linalg.norm(coords[:,None,:]-coords[None,:,:],axis=-1)*1E-10 #m
    r_ij = np.where(r_ij==0,np.nan,r_ij)
    #np.fill_diagonal(r_ij,np.nan)
    q_ij = np.outer(charges,charges)*e*e
    a_ij = q_ij/(4*np.pi*epsilon_0*epsilon_r*r_ij)
    a_ij = np.nan_to_num(a_ij,nan=0)
    col_energy = np.triu(a_ij,k=1).sum() #J
    col_energy = col_energy*NA/1000 #kJ/mol
    return col_energy


def create_pqr_cg(
        pqr_file: str,
        cg_pos: Dict,
        cg_charge: Dict,
        cg_radius: float=11.185
    ) -> int:
  """
  Funcion que espera un pqr de ARN, un arreglo de posiciones de cg_pos y una carga.
  Inserta el CG en el PQR y genera el CG PQR con la carga de parametrizacion.
  """

  working_dir, pqr_filename = os.path.split(pqr_file)
  pqr_name,_ = os.path.splitext(pqr_filename)

  base_line = 'ATOM   {atom_id}  H41  RC  {chain_id}  {resd_id}     {pos_x}  {pos_y}  {pos_z}  {charge} {radius}'

  new_lines = []
  n = 0
  for chain in cg_pos.keys():
    insert_line = []
    pos = cg_pos[chain]
    chg = cg_charge[chain]
    for i in range(pos.shape[0]):
        line = (
            f"ATOM   "              # 1-7
            f"{90000+n:>5}"         # 8-12
            f" "                    # 13
            f"{'H41':>4}"           # 14-17
            f" "                    # 18
            f"{'RC':>3}"            # 19-21
            f"  "                   # 22-23
            f"{chain:<1}"           # 24
            f" "                    # 25
            f"{i:>3}"               # 26-28
            f"    "                 # 29-32
            f"{pos[i,0]:>8.3f}"     # 33-40
            f" "                    # 41
            f"{pos[i,1]:>8.3f}"     # 42-49
            f" "                    # 50
            f"{pos[i,2]:>8.3f}"     # 51-58
            f" "                    # 59
            f"{chg[i,0]:>7.4f}"     # 60-66
            f" "                    # 67
            f"{cg_radius:>7.4f}"    # 68-73
        ).ljust(73)[:73]
        insert_line.append(line)
        n += 1



        # insert_line.append(base_line.format(
        #     resd_id=i,
        #     atom_id=i+10000,
        #     chain_id=chain,
        #     pos_x=f'{pos[i,0]:.3f}',
        #     pos_y=f'{pos[i,1]:.3f}',
        #     pos_z=f'{pos[i,2]:.3f}',
        #     charge=f'{chg[i,0]:.3f}',
        #     radius=f'{cg_radius:.3f}'
        # ))
    insert_line = '\n'.join(insert_line)
    new_lines.append(insert_line)     

  #Eliminar todas las cargas existentes
  temp_file = open(os.path.join(working_dir,'temp.pqr'),'w')
  awk_command = '''
    $1 == "ATOM" {
    $0 = substr($0, 1, 59) " 0.0000" substr($0, 67)
    }
    {print}
    '''
  subprocess.run(['awk',awk_command, pqr_filename], stdout=temp_file, cwd=working_dir)
  temp_file.close()
  subprocess.run(['mv','temp.pqr', pqr_filename],stdout=subprocess.DEVNULL, cwd=working_dir)

  #Crea PQR del caso CG
  with open(pqr_file,'r') as f:
      lines = f.readlines()
  with open(pqr_file,'a') as f:
      for cg_set in new_lines:
        f.writelines(cg_set + '\n') # type: ignore
  return 0


def get_energies(pqr_file: str, left_iso_pqr: str, right_iso_pqr) -> Dict[str,float]:

    pqr_file_path,pqr_filename = os.path.split(pqr_file)
    prot_basename,_ = os.path.splitext(pqr_filename)
    apbs_log = os.path.join(pqr_file_path, prot_basename + '.log')

    left_iso_pqr_path, left_iso_pqr_filename = os.path.split(left_iso_pqr)
    left_iso_basename,_ = os.path.splitext(left_iso_pqr_filename)
    left_apbs_log = os.path.join(left_iso_pqr_path,left_iso_basename + '.log')

    right_iso_pqr_path, right_iso_pqr_filename = os.path.split(right_iso_pqr)
    rigth_iso_basename,_ = os.path.splitext(right_iso_pqr_filename)
    right_apbs_log = os.path.join(right_iso_pqr_path,rigth_iso_basename + '.log')

    g_col_iso_left = calc_coulomb(left_iso_pqr) #kJ/mol
    g_col_iso_right = calc_coulomb(right_iso_pqr) #kJ/mol
    g_col_com = calc_coulomb(pqr_file) #kJ/mol
    g_solv_iso_left = apbs_tool.extract_energy_from_log(left_apbs_log) #kJ/mol
    g_solv_iso_right = apbs_tool.extract_energy_from_log(right_apbs_log) #kJ/mol
    g_solv_com = apbs_tool.extract_energy_from_log(apbs_log) #kJ/mol
    
    dd_g_col = g_col_com - g_col_iso_left - g_col_iso_right
    dd_g_solv = g_solv_com - g_solv_iso_left - g_solv_iso_right
    g_binding = dd_g_col + dd_g_solv

    energies = {'G Col Com':g_col_com,
        'G Col Iso L':g_col_iso_left,
        'G Col Iso R':g_col_iso_right,
        'G Solv Com':g_solv_com,
        'G Solv Iso L':g_solv_iso_left,
        'G Solv Iso R':g_solv_iso_right,
        'DD Col G':dd_g_col,
        'DD Solv G':dd_g_solv,
        'G Binding':g_binding
    }
                            
    return energies

def gen_isolated_pqr(pqr_file: str, chains: List[str]) -> int:
    """
	Function that generates the isolated pqr file by masking the charges and radius of chain ids not in chains.

    Args:
	    pqr_file: Target pqr file. Needs to be a path. 
        chains: List of chain id's to keep.
        sim_type: Identifies if pqr is 'left' or 'right'.
	"""
 
    working_dir, pqr_filename = os.path.split(pqr_file)
    pqr_basename,_ = os.path.splitext(pqr_filename)

    matching_chains = '$5=='
    for i, chain in enumerate(chains):
        if i == (len(chains)-1):
            matching_chains += f'"{chain}"'
        else:
            matching_chains += f'"{chain}"||$5=='
    awk_command = '''
    $1 == "ATOM" && !({chains_pattern}) {
    $0 = substr($0, 1, 59) " 0.0000 0.0000"
    }
    {print}
    '''
    awk_command = awk_command.replace('{chains_pattern}',matching_chains)
    awk_args = [
        'awk',
        awk_command,
        pqr_filename]

    isolated_pqr = open(os.path.join(working_dir,'temp.pqr'),'w')
    subprocess.run(
        awk_args,
        cwd=working_dir,
        stdout=isolated_pqr)
    subprocess.run(
        ['mv','temp.pqr',pqr_filename],
        cwd=working_dir,
        stdout=subprocess.DEVNULL
    )
    return 0
