# Class Definitions

## ConfigHandler

Functions:
	- Ingest yaml config file.
	- Verify yaml integrity for all needed fields and files. Raise error otherwise.
	- Import pdb_handler function.
	- Create CaseHandler objects.
	- Create top level folder structure.
	- Configure all cases.
	- Run all cases.
	- Dry run option. (Config only)
	- Export all results in csv.

Variables:
	- yaml config file.
	- Dictionary of cases. {case_name:CaseHandler}
	- dry_run flag.
	- results dictionary.
	- configured/ready_to_run flag
	- finished_run flag
	- top level path

Input Variables:
	- yaml file path.

## CaseHandler

Functions:
	- Create all SimulationRunner objects for CG and AllAtom (full) simulations.
	- Run simulations in parallel.
	- Configure all simulations.
	- Create wall distance csv in case root path.

Variables:
	- finished_run flag.
	- configured/ready_to_run flag.
	- Dictionary of all atom case simulations. {sim_name:SimulationRunner}
	- Dictionary of CG case simulations. {sim_name:SimulationRunner}
	- Dictionary of all atom isolated simulatinos. {sim_name:SimulationRunner}
	- Dictionary of CG isolated simulations. {sim_name:SimulationRunner}
	- case_results dictionary. {sim_name:SimulationRunner.results[Dict]}

Input Variables:
	- Case root path.
	- Case simulation options dictionary.
	- Case PDB Handler Options dictionary (optional).
	- dry_run: bool (optional). 

## SimulationHandler

Functions:
	- Configure simulation.
	- Run simulation.
	- Recover simulation.
	- Store results in json in sim_path.
	- Creates SimulationOutput objects. (Dataclass with simulation results)

Variables:
	- Simulation files dictionary.
	- Simulation root path.
	- Simulation options dictionary:
		Name: str
		keep_dx: bool
		pb_type: bool
		core_pdb_path: str
		type: str (allAtom or coarseGrain)
		isolated: bool
		isolated_type: str (left or right) (only if isolated True)
		Parent Case Name (optional)
	- PDB Handler options dictionary:
		Handler Parameters Dictionary
	- Results
	- Configured flag
	- Finished flag
		
Input Variables:
	- Simulation path
	- Simulation Options
	- Handler Options
	- Isolated Left SimulationResult object (Only if Complex Simulation)
	- Isolated Right SimulationResult object (Only if Complex Simulation)

## Handler Functions

Input Variables:
	- Core PDB
	- Target Path
	- Handler Options

Functions:
	- Read Core PDB and modify  according to handler options.
	- Create main PDB and PQR files.
	- Create LEFT and RIGHT PDB and PQR files.

Return:
	- Filename of created main PDB
	- Left/Right Chains dictionary of PDB: {'left':['A','B'],'right':['D','E']
