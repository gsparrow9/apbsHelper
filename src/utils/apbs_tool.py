import os
import subprocess
import MDAnalysis as mda
from gridData import Grid
import numpy as np
import matplotlib.pyplot as plt

class DXFile:
    def __init__(self, dx_filename):
        self.dx_filename = dx_filename
        self.data = self.load_dx_file(dx_filename)
    
    def load_dx_file(self, dx_filename):
        """
        Load the data from a .dx file.
        """
        data = Grid(dx_filename)
        
        return data
    
    def visualize(self, cut_position='middle', axis='x', vmin=None, vmax=None):
        """
        Visualize a 2D slice of the electrostatic potential.
        
        Parameters:
            cut_position (str or float): Position along the chosen axis to slice. 
                                         Options: 'middle', 'min', 'max', or a specific float value.
            axis (str): Axis along which to slice ('x', 'y', or 'z').
        """
        if axis == 'x':
            if cut_position == 'middle':
                slice_index = np.shape(self.data.grid)[0] // 2
            elif cut_position == 'min':
                slice_index = np.unravel_index(np.argmin(self.data.grid, axis=None), np.shape(self.data.grid))[0]
            elif cut_position == 'max':
                slice_index = np.unravel_index(np.argmax(self.data.grid, axis=None), np.shape(self.data.grid))[0]
            else:
                slice_index = int((cut_position - self.data.origin[0]) / self.data.delta[0])
            slice_data = self.data.grid[slice_index, :, :]
            x_data = self.data.midpoints[1]
            y_data = self.data.midpoints[2]
            xlabel, ylabel = 'y', 'z'
        
        elif axis == 'y':
            if cut_position == 'middle':
                slice_index = np.shape(self.data.grid)[1] // 2
            elif cut_position == 'min':
                slice_index = np.unravel_index(np.argmin(self.data.grid, axis=None), np.shape(self.data.grid))[1]
            elif cut_position == 'max':
                slice_index = np.unravel_index(np.argmax(self.data.grid, axis=None), np.shape(self.data.grid))[1]
            else:
                slice_index = int((cut_position - self.data.origin[1]) / self.data.delta[1])
            slice_data = self.data.grid[:, slice_index, :]
            x_data = self.data.midpoints[0]
            y_data = self.data.midpoints[2]
            xlabel, ylabel = 'x', 'z'
        
        elif axis == 'z':
            if cut_position == 'middle':
                slice_index = np.shape(self.data.grid)[2] // 2
            elif cut_position == 'min':
                slice_index = np.unravel_index(np.argmin(self.data.grid, axis=None), np.shape(self.data.grid))[2]
            elif cut_position == 'max':
                slice_index = np.unravel_index(np.argmax(self.data.grid, axis=None), np.shape(self.data.grid))[2]
            else:
                slice_index = int((cut_position - self.data.origin[2]) / self.data.delta[2])           
            slice_data = self.data.grid[:, :, slice_index]
            x_data = self.data.midpoints[0]
            y_data = self.data.midpoints[1]
            xlabel, ylabel = 'x', 'y'
        
        else:
            raise ValueError("Invalid axis; choose 'x', 'y', or 'z'")
        
        X,Y = np.meshgrid(x_data, y_data, indexing="ij")
        plt.scatter(X, Y, c=slice_data, cmap="coolwarm", vmin=vmin, vmax=vmax)
        plt.colorbar(label="Potential (kT/e)")
        plt.xlabel(f'{xlabel} axis')
        plt.ylabel(f'{ylabel} axis')
        plt.title(f'Slice at {axis} = {cut_position}')
        plt.show()

    def __sub__(self, other):
        """
        Overload the subtraction operator to calculate the difference between two DX files.
        """
        if not isinstance(other, DXFile):
            raise TypeError("Subtraction is only supported between DXFile instances.")
        
        if np.shape(self.data.grid) != np.shape(other.data.grid) or not np.allclose(self.data.delta, other.data.delta):
            raise ValueError("DX files must have the same grid shape and spacing to be subtracted.")
        
        # Perform element-wise subtraction of the data grids
        diff_data = self.data.grid - other.data.grid
        
        # Return a new DXFile instance with the difference
        diff_dx = DXFile.__new__(DXFile)  # Create an uninitialized instance of DXFile

        diff_dx.data = Grid(diff_data, edges=self.data.edges)
        
        return diff_dx

    def __repr__(self):
        return f"<DXFile: {self.dx_filename} with grid shape {self.data.shape}>" # type: ignore

def generate_pqr_files(pdb_folder, pdb2pqr_path):
    """
    Generate PQR files from PDB files using pdb2pqr.

    Parameters:
        pdb_folder (str): Path to the folder containing PDB files.
        pdb2pqr_path (str): Path to the pdb2pqr executable.
    """
    # Define the replacement mappings for nucleotide names
    nucleotide_mapping = {
        "GUA": "RG",
        "CYT": "RC",
        "ADE": "RA",
        "URA": "RU"
    }
    
    for pdb_file in os.listdir(pdb_folder):
        if pdb_file.endswith(".pdb"):
            pdb_filepath = os.path.join(pdb_folder, pdb_file)
            pqr_filepath = pdb_filepath.replace(".pdb", ".pqr")
            
            # Check if PQR file already exists
            if not os.path.exists(pqr_filepath):
                # Replace nucleotide names in the PDB file
                with open(pdb_filepath, 'r') as file:
                    pdb_data = file.read()
                for old, new in nucleotide_mapping.items():
                    pdb_data = pdb_data.replace(old, new)
                with open(pdb_filepath, 'w') as file:
                    file.write(pdb_data)
                
                # Run pdb2pqr to generate PQR file
                subprocess.run(
                    [pdb2pqr_path, pdb_filepath, pqr_filepath, "--ff=amber","--whitespace", "--chain"],
                    stderr=subprocess.DEVNULL,
                    stdout=subprocess.DEVNULL
                )

def modify_and_run_apbs(pdb_folder, apbs_path, ionic_concentration=0.150, perm_solvent=78.54, perm_solute=2.0, spacing=0.3):
    """
    Modify the base input file for APBS and run APBS for each PQR file.

    Parameters:
        pdb_folder (str): Path to the folder containing PDB and PQR files.
        apbs_path (str): Path to the APBS executable.
        ionic_concentration (float): Ionic concentration for the APBS calculation.
        perm_solvent (float): Permittivity of the solvent.
        perm_solute (float): Permittivity of the solute.
        spacing (float): Spacing between grid points in Angstrom.
    """

    input_file_text = """
    read
    mol pqr {{PQR_FILE}}
    end
    elec name pb_run
        mg-auto
        dime {{DIME_X}} {{DIME_Y}} {{DIME_Z}} 
        cglen {{CGLEN_X}} {{CGLEN_Y}} {{CGLEN_Z}}
        fglen {{FGLEN_X}} {{FGLEN_Y}} {{FGLEN_Z}}
        cgcent mol 1
        fgcent mol 1
        mol 1
        {{PB_TYPE}}
        bcfl mdh
        pdie {{PERM_SOLUTE}}
        sdie {{PERM_SOLVENT}}
        ion charge 1 conc {{IONIC_CONCENTRATION}} radius 0
        ion charge -1 conc {{IONIC_CONCENTRATION}} radius 0
        srfm smol
        chgm spl2
        sdens 10.00
        srad 1.40
        swin 0.30
        temp 298.15
        calcenergy total
        calcforce no
        write pot dx {{OUTPUT_FILENAME_DX}}
    end
    print 
        elecEnergy pb_run
    end
    """

    # define mesh size
    length_x = 0.
    length_y = 0.
    length_z = 0.
    for pqr_file in os.listdir(pdb_folder):
        if pqr_file.endswith(".pqr"):
            pqr_filepath = os.path.join(pdb_folder, pqr_file)
            pqr_data = mda.Universe(pqr_filepath)
            pos = pqr_data.coord.positions
            xmin,xmax = np.min(pos[:,0]), np.max(pos[:,0])
            ymin,ymax = np.min(pos[:,1]), np.max(pos[:,1])
            zmin,zmax = np.min(pos[:,2]), np.max(pos[:,2])

            length_x = max(length_x, xmax-xmin)
            length_y = max(length_y, ymax-ymin)
            length_z = max(length_z, zmax-zmin)

    length_x *= 1.9 # length of mesh is 1.9x larger than molecule
    length_y *= 1.9 # length of mesh is 1.9x larger than molecule
    length_z *= 1.9 # length of mesh is 1.9x larger than molecule
    
    possible_dime = np.array([65, 97, 129, 161, 193, 257, 513])
    
    dx = length_x/(possible_dime-1)
    dy = length_y/(possible_dime-1)
    dz = length_z/(possible_dime-1)

    index_dime = np.argwhere(dx<spacing)[0][0]
    dime_x = possible_dime[index_dime]
    index_dime = np.argwhere(dy<spacing)[0][0]
    dime_y = possible_dime[index_dime]
    index_dime = np.argwhere(dz<spacing)[0][0]
    dime_z = possible_dime[index_dime]

    energy_file = open(pdb_folder+"/energy.csv","w")
    energy_file.write("pdb_file,pqr_file,energy_linear_kJ,energy_nonlinear_kJ,energy_vacuum_kJ,solv_energy_kJ,dx_file_linear_kT,dx_file_nonlinear_kT,dx_file_vacuum_kT\n")
    for pqr_file in os.listdir(pdb_folder):
        if pqr_file.endswith(".pqr"):
            pqr_filepath = os.path.join(pdb_folder, pqr_file)
            pdb_basename = pqr_file.replace(".pqr", "")
            
            # Modify the base input file for each case
            energy_cache = []
            for case, pb_type in zip(['linear_solvated', 'nonlinear_solvated', 'linear_vacuum'], ['lpbe', 'npbe', 'lpbe']):
                case_input_file = os.path.join(pdb_folder, f"{pdb_basename}_{case}.in")
               
                input_data = input_file_text
                
                input_data = input_data.replace("{{PQR_FILE}}", pqr_filepath)
                input_data = input_data.replace("{{DIME_X}}", str(dime_x))
                input_data = input_data.replace("{{DIME_Y}}", str(dime_y))
                input_data = input_data.replace("{{DIME_Z}}", str(dime_z))
                input_data = input_data.replace("{{CGLEN_X}}", str(length_x))
                input_data = input_data.replace("{{CGLEN_Y}}", str(length_y))
                input_data = input_data.replace("{{CGLEN_Z}}", str(length_z))
                input_data = input_data.replace("{{FGLEN_X}}", str(length_x/1.12)) # just 12% smaller
                input_data = input_data.replace("{{FGLEN_Y}}", str(length_y/1.12))
                input_data = input_data.replace("{{FGLEN_Z}}", str(length_z/1.12))
                input_data = input_data.replace("{{PB_TYPE}}", pb_type)
                input_data = input_data.replace("{{IONIC_CONCENTRATION}}", str(ionic_concentration))
                input_data = input_data.replace("{{PERM_SOLVENT}}", str(perm_solvent))
                input_data = input_data.replace("{{PERM_SOLUTE}}", str(perm_solute))
                input_data = input_data.replace("{{OUTPUT_FILENAME_DX}}", str(case_input_file.replace(".in","")))
                
                with open(case_input_file, 'w') as file:
                    file.write(input_data)
                
                # Run APBS
                log_filename = pdb_basename + ".log"
                with open(log_filename, "w") as log_file:
                    subprocess.run([apbs_path, case_input_file], stdout=log_file)
                log_file.close()
                print(f"Ran APBS for case: {case} - {case_input_file}")

                energy = extract_energy_from_log(log_filename)
                energy_cache.append(energy)
                os.remove(log_filename)
                os.remove(case_input_file)

            energy_file.write(f"{pdb_basename}.pdb,{pqr_file},{energy_cache[0]},{energy_cache[1]},{energy_cache[2]},{energy_cache[1]-energy_cache[2]},{pdb_folder}{pdb_basename}_linear_solvated.dx,{pdb_folder}{pdb_basename}_nonlinear_solvated.dx,{pdb_folder}{pdb_basename}_linear_vacuum.dx\n")
    energy_file.close()
            

def extract_energy_from_log(log_filename):
    """
    Extract the electrostatic energy from an APBS log file.

    Parameters:
        log_filename (str): Path to the log file.

    Returns:
        float: The total electrostatic energy extracted from the log file.
    """
    with open(log_filename, "r") as f:
        for line in f:
            if "Global net ELEC" in line:
                # The energy is typically the last number on the line
                energy_value = float(line.split()[-2])
                return energy_value

    raise ValueError(f"Energy not found in {log_filename}")



if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Electrostatic potential calculation using APBS for RNA molecules.")
    parser.add_argument("--pdb_folder", type=str, default=".", help="Folder containing the PDB files (default is current directory).")
    parser.add_argument("--pdb2pqr_path", type=str, default="/home/chris/Software/apbs-pdb2pqr/pdb2pqr/pdb2pqr.py" ,help="Path to the pdb2pqr executable.")
    parser.add_argument("--apbs_path", type=str, default="apbs", help="Path to the APBS executable.")
    parser.add_argument("--ionic_concentration", type=float, default=0.15, help="Ionic concentration for the APBS calculation.")
    parser.add_argument("--eps_solvent", type=float, default=78.54, help="Permittivity of the solvent.")
    parser.add_argument("--eps_solute", type=float, default=2.0, help="Permittivity of the solute.")
    parser.add_argument("--spacing", type=float, default=0.2, help="Mesh spacing for APBS.")
    
    args = parser.parse_args()
    
    generate_pqr_files(args.pdb_folder, args.pdb2pqr_path)
    modify_and_run_apbs(args.pdb_folder, args.apbs_path, args.ionic_concentration, args.eps_solvent, args.eps_solute, args.spacing)
