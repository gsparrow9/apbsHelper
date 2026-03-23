from src.utils.setup_logs import setup_log
from src.utils import tools as tools                
from src.utils import apbs_tool as apbs_tool         
from src.classes.handler_functions import HandlerFunction
from typing import Optional, Dict, Tuple, Literal
from dataclasses import dataclass
from pathlib import Path
import subprocess
import logging
import json
import time

@dataclass(frozen=True)
class SimulationResult:
    pqr_file: Path
    sim_path: Path
    finished_ok: bool                            #Sim status. Ok or Error
    yukawa: bool
    apbs_log: Optional[Path] = None
    misc: Optional[Dict] = None
    is_isolated: Optional[bool] = False          #Isolated or Complex
    left_iso_energies: Optional[Dict] = None     #Needed if is_isolated
    right_iso_energies: Optional[Dict] = None    #Needed if is_isolated
    wall_distance: Optional[float] = None        #Needed if complex

    def __post_init__(self):
        if not self.yukawa:
            if self.apbs_log is None:
                raise ValueError('Need apbs log when yukawa=False')
        if self.is_isolated:
            if self.left_iso_energies is not None or self.right_iso_energies is not None:
                raise ValueError('left_iso_sim/right_iso_sim must be None when is_isolated=True')
        else:
            if self.left_iso_energies is None or self.right_iso_energies is None:
                raise ValueError('Need left_iso_sim and right_iso_sim when is_isolated=False')
            if self.wall_distance is None:
                raise ValueError('Need wall_distance when is_isolated=False')                        

    @property
    def energies(self) -> Dict[str, float]:
        if not self.finished_ok:
            return {}

        if self.is_isolated:
            if self.yukawa:
                _energies = {
                    'solv':tools.calc_yukawa(self.pqr_file),
                    'col':tools.calc_coulomb(self.pqr_file)
                }
            else:
                _energies = {
                    'solv':apbs_tool.extract_energy_from_log(self.apbs_log),
                    'col':tools.calc_coulomb(self.pqr_file)
                }
            return _energies
        else:
            if self.yukawa:
                solvation = tools.calc_yukawa(self.pqr_file)
            else:
                solvation = apbs_tool.extract_energy_from_log(self.apbs_log)
            coulomb = tools.calc_coulomb(self.pqr_file)
            dd_g_solv = solvation - self.left_iso_energies['solv'] - self.right_iso_energies['solv'] # type: ignore
            dd_g_col = coulomb - self.left_iso_energies['col'] - self.right_iso_energies['col'] # type: ignore

            if self.yukawa:
                _energies = {
                    'solv':solvation,
                    'col':coulomb,
                    'solv_left':self.left_iso_energies['solv'], # type: ignore
                    'col_left':self.left_iso_energies['col'], # type: ignore
                    'solv_right':self.right_iso_energies['solv'], # type: ignore
                    'col_right':self.right_iso_energies['col'], # type: ignore
                    'dd_g_solv':dd_g_solv,
                    'dd_g_col':dd_g_col,
                    'binding':dd_g_solv
                }
            else:
                _energies = {
                    'solv':solvation,
                    'col':coulomb,
                    'solv_left':self.left_iso_energies['solv'], # type: ignore
                    'col_left':self.left_iso_energies['col'], # type: ignore
                    'solv_right':self.right_iso_energies['solv'], # type: ignore
                    'col_right':self.right_iso_energies['col'], # type: ignore
                    'dd_g_solv':dd_g_solv,
                    'dd_g_col':dd_g_col,
                    'binding':dd_g_solv + dd_g_col
                }
            return _energies

    def to_json(self,filename: str|Path) -> None:
        out_dict = {}
        out_dict['sim_path'] = str(self.sim_path)
        out_dict['pqr_file'] = str(self.pqr_file)
        out_dict['apbs_log'] = str(self.apbs_log) if self.apbs_log is not None else None
        out_dict['finished_ok'] = self.finished_ok
        out_dict['is_isolated'] = self.is_isolated
        out_dict['is_yukawa'] = self.yukawa

        if self.misc is not None: 
            out_dict['misc'] = self.misc
        else:
            out_dict['misc'] = None
        
        out_dict['energies'] = self.energies
        
        if not self.is_isolated:
            out_dict['wall_distance'] = self.wall_distance

        with open(filename,'w') as f:
            json.dump(out_dict,f)
        return None


class SimulationHandler:

    def __init__(
            self,
            sim_path: str|Path,
            mesh_size: float,
            sim_name: str,
            core_pdb_path: str|Path,
            sim_type: Literal['allAtom','coarseGrained','yukawaCoarseGrained','yukawaAllAtom'],
            isolated: bool,
            handler_func: HandlerFunction,
            handler_options: Dict,
            keep_dx: bool=False,
            linear_pb: bool=False,
            log_level: int=logging.INFO,
            iso_left_sim: Optional["SimulationHandler"]=None,
            iso_right_sim: Optional["SimulationHandler"]=None,
            isolated_type: Literal['left','right']|None=None
    ) -> None:
        setup_log(log_level)
        self.logger = logging.getLogger(self.__class__.__name__)
        self.core_pdb = Path(core_pdb_path)
        self.sim_path = Path(sim_path)

        self.sim_files = {}
        self.sim_files['pdb'] = self.sim_path / Path(self.core_pdb.name)
        self.sim_files['pqr'] = self.sim_files['pdb'].with_suffix('.pqr')
        self.sim_files['json'] = self.sim_files['pdb'].with_suffix('.json')
        self.sim_files['apbs_log'] = self.sim_files['pdb'].with_suffix('.log')
        self.sim_files['apbs_in'] = self.sim_files['pdb'].with_suffix('.in')

        self.mesh_size = mesh_size
        self.handler_func = handler_func       
        self.handler_options = handler_options
        self.keep_dx = keep_dx
        self.linear_pb = linear_pb
        self.sim_name = sim_name
        self.sim_type = sim_type
        self.isolated = isolated
        self.sim_configured = False
        self._finished_ok = False
        self.sim_exists = self.sim_files['json'].exists()

        if self.sim_exists:
            self.result = self._recover_result(self.sim_files['json'])

        if isolated:
            if isolated_type is None:
                raise ValueError('isolated_type is needed when isolated=True')
            elif isolated_type != 'left' and isolated_type != 'right':
                raise ValueError("isolated_type must be 'left' or 'right'")
            self.isolated_type = isolated_type.lower()
            self.iso_right_sim = None
            self.iso_left_sim = None
        else:
            if iso_left_sim is None:
                raise ValueError('iso_left_result is needed when isolated=False')
            if iso_right_sim is None:
                raise ValueError('iso_right_result is needed when isolated=False')
            self.isolated_type = None
            self.iso_right_sim = iso_right_sim
            self.iso_left_sim = iso_left_sim
            self.sim_files['left_pqr'] = self.sim_path / Path('left_' + self.sim_files['pqr'].stem + '.pqr')
            self.sim_files['right_pqr'] = self.sim_path / Path('right_' + self.sim_files['pqr'].stem + '.pqr')
            self.sim_files['left_xyzr'] = self.sim_path / Path('left_' + self.sim_files['pqr'].stem + '.xyzr')
            self.sim_files['right_xyzr'] = self.sim_path / Path('right_' + self.sim_files['pqr'].stem + '.xyzr')

        if (sim_type != 'allAtom' and sim_type != 'coarseGrained'
            and sim_type != 'yukawaCoarseGrained' and sim_type != 'yukawaAllAtom'):
            raise ValueError("sim_type must be 'allAtom', 'coarseGrained', 'yukawaAllAtom' or 'yukawaCoarseGrained'")

        self.logger.debug(f'Created {self}')
        pass

    def __repr__(self) -> str:
        if self.isolated:
            return f'<SimulationHandler(name={self.sim_name},type=isolated_{self.isolated_type})>'
        else:
            return f'<SimulationHandler(name={self.sim_name}),type=complex)>'

    def configure_simulation(self) -> None:
        
        if self.sim_configured or self.sim_exists:
            return None

        if self.isolated:
            self.logger.info(f'Configuring simulation {self.isolated_type} {self.sim_name} - {self.sim_path}')
        else:
            self.logger.info(f'Configuring simulation {self.sim_name} - {self.sim_path}')
        self.sim_path.mkdir(exist_ok=True, parents=True)

        #Generating PDB file with handler_func
        try:
            if self.sim_type == 'allAtom' or self.sim_type == 'yukawaAllAtom':
                out_pdb,chains = self.handler_func.handle(
                    self.core_pdb,
                    self.sim_path,
                    handler_options=self.handler_options,
                    coarse_grain=False
                ) 
            elif self.sim_type == 'coarseGrained' or self.sim_type == 'yukawaCoarseGrained':
                out_pdb,chains = self.handler_func.handle(
                    self.core_pdb,
                    self.sim_path,
                    handler_options=self.handler_options,
                    coarse_grain=True
                )        
        except Exception as e:
            self.logger.error(f'Encountered error in pdb_handler. {e}')
            raise
        
        #Generating isolated pqr file
        try:
            if self.isolated:
                self.logger.info(f'Generating isolated pqr for {self.sim_name} - {self.sim_path}')
                if self.isolated_type == 'left':
                    tools.gen_isolated_pqr(self.sim_files['pqr'],chains['left']) # type: ignore
                elif self.isolated_type == 'right':
                    tools.gen_isolated_pqr(self.sim_files['pqr'],chains['right']) # type: ignore
                else:
                    raise RuntimeError
        except RuntimeError:
            self.logger.error('Wrong sim_type')
        except Exception as e:
            self.logger.error(f'Encountered an error in generating isolated pqr. {e}')
            raise

        #Getting wall distance with Nanoshaper
        try:
            if not self.isolated:
                self.logger.info(f'Getting wall distance for {self.sim_name} with Nanoshaper')
                tools.create_xyzr(self.sim_files['left_pqr'])
                tools.create_xyzr(self.sim_files['right_pqr'])
                self.wall_distance = tools.get_wall_distance(
                    self.sim_files['left_xyzr'],
                    self.sim_files['right_xyzr']
                )
        except Exception as e:
            self.logger.error(f'Encountered an error in getting wall distance for {self.sim_name}. {e}')
            raise

        if self.sim_type != 'yukawaCoarseGrained' or self.sim_type != 'yukawaAllAtom':
            #Generating apbs config file
            self.logger.info(f'Generating APBS config file for {self.sim_name}')
            if self.isolated:
                self.bbox,self.clearance = tools.gen_apbs(
                    self.sim_files['pqr'],
                    self.mesh_size,
                    keep_dx=self.keep_dx,
                    linear=self.linear_pb
                )
            else:
                self.bbox,self.clearance = tools.gen_apbs(
                    self.sim_files['pqr'],
                    self.mesh_size,
                    keep_dx=False,
                    linear=self.linear_pb
                )
            self.logger.info(f'{self.sim_name} BBOX: {self.bbox} - Clearance: {self.clearance}')
        self.logger.info(f'Finished configuring {self.sim_name} - {self.sim_path}')
        self.sim_configured = True
        return None

    def run_simulation(self) -> None:
        if self.isolated:
            self.logger.info(f'Running simulation {self.isolated_type} {self.sim_name} - {self.sim_path}')
        else:
            self.logger.info(f'Running simulation {self.sim_name} - {self.sim_path}')
        if ((not self.isolated)
            and (not self.iso_right_sim._finished_ok or not self.iso_left_sim._finished_ok) # type: ignore
        ):
            self.logger.error(f'Isolated simulations for {self.sim_name} have errors')
            raise RuntimeError('Isolated simulations for complex simulation have not run')

        if not self.sim_configured:
            self.configure_simulation()

        try:

            if 'yukawa' in self.sim_type.lower():        
                misc = {'sim_name':self.sim_name,
                        'wall_clock':None,
                        'status':'ok',
                        'sim_type':self.sim_type,
                        'linear_pb':self.linear_pb}
            else:
                _,wall_clock = self._execute_apbs()
                misc = {'sim_name':self.sim_name,
                        'wall_clock':wall_clock,
                        'status':'ok',
                        'sim_type':self.sim_type,
                        'linear_pb':self.linear_pb}
            
            if self.handler_options:
                    misc.update(self.handler_options)

            if self.isolated:
                if 'yukawa' in self.sim_type.lower():
                    self.result = SimulationResult(
                        pqr_file=self.sim_files['pqr'],
                        sim_path=self.sim_path,
                        finished_ok=True,
                        yukawa=True,
                        misc=misc,
                        is_isolated=self.isolated
                    )
                    self._finished_ok = True
                else:
                    self.result = SimulationResult(
                        pqr_file=self.sim_files['pqr'],
                        sim_path=self.sim_path,
                        finished_ok=True,
                        yukawa=False,
                        apbs_log=self.sim_files['apbs_log'],
                        misc=misc,
                        is_isolated=self.isolated
                    )
            else:
                if 'yukawa' in self.sim_type.lower():
                    self.result = SimulationResult(
                        pqr_file=self.sim_files['pqr'],
                        sim_path=self.sim_path,
                        finished_ok=True,
                        yukawa=True,
                        misc=misc,
                        is_isolated=self.isolated,
                        wall_distance=self.wall_distance,
                        left_iso_energies=self.iso_left_sim.result.energies, # type: ignore
                        right_iso_energies=self.iso_right_sim.result.energies # type: ignore
                    )
                    self._finished_ok = True
                else:
                    self.result = SimulationResult(
                        pqr_file=self.sim_files['pqr'],
                        sim_path=self.sim_path,
                        finished_ok=True,
                        yukawa=False,
                        apbs_log=self.sim_files['apbs_log'],
                        misc=misc,
                        is_isolated=self.isolated,
                        wall_distance=self.wall_distance,
                        left_iso_energies=self.iso_left_sim.result.energies, # type: ignore
                        right_iso_energies=self.iso_right_sim.result.energies # type: ignore
                    )

        except Exception as e:
            self.logger.error(f'Error in executing apbs for {self.sim_name}. {e}')
            misc = {'sim_name':self.sim_name,'status':'failed'}

            if self.handler_options:
                misc.update(self.handler_options)

            self.result = SimulationResult(
                pqr_file=self.sim_files['pqr'],
                sim_path=self.sim_path,
                finished_ok=False,
                yukawa=True,
                misc=misc,
                is_isolated=True,
            )
            self.logger.info(f'Simulation {self.sim_name} - ERROR')
            raise
        self.logger.info(f'Simulation {self.sim_name} - OK')
        self.result.to_json(self.sim_files['json'])
        return None

    def _recover_result(self,json_file: str|Path) -> SimulationResult:
        self.logger.info(f'Attempting to recover {self.sim_name} - {self.sim_path}')
        with open(json_file,'r') as f:
            dump = json.load(f)
        pqr_file = dump['pqr_file']
        apbs_log = dump['apbs_log']
        sim_path = dump['sim_path']
        finished_ok = dump['finished_ok']
        is_isolated = dump['is_isolated']
        misc = dump['misc']
        energies = dump['energies']
        self._finished_ok = finished_ok

        if is_isolated:
            if self.sim_type == 'yukawaCoarseGrained' or self.sim_type == 'yukawaAllAtom':
                result = SimulationResult(
                    pqr_file=pqr_file,
                    sim_path=sim_path,
                    finished_ok=True,
                    yukawa=True,
                    misc=misc,
                    is_isolated=self.isolated
                )
            else:
                result = SimulationResult(
                    pqr_file=pqr_file,
                    sim_path=sim_path,
                    finished_ok=True,
                    yukawa=False,
                    apbs_log=apbs_log,
                    misc=misc,
                    is_isolated=self.isolated
                )
        else:
            left_energies = {'solv':energies['solv_left'],'col':energies['col_left']}
            right_energies = {'solv':energies['solv_right'],'col':energies['col_right']}
            if self.sim_type == 'yukawaCoarseGrained' or self.sim_type == 'yukawaAllAtom':
                result = SimulationResult(
                    pqr_file=pqr_file,
                    sim_path=sim_path,
                    finished_ok=True,
                    yukawa=True,
                    misc=misc,
                    is_isolated=self.isolated,
                    left_iso_energies=left_energies, # type: ignore
                    right_iso_energies=right_energies # type: ignore
                )
            else:
                result = SimulationResult(
                    pqr_file=pqr_file,
                    sim_path=sim_path,
                    finished_ok=True,
                    yukawa=False,
                    apbs_log=apbs_log,
                    misc=misc,
                    is_isolated=self.isolated,
                    left_iso_energies=left_energies, # type: ignore
                    right_iso_energies=right_energies # type: ignore
                )
            self.wall_distance = dump['wall_distance']
        
        return result

    def _execute_apbs(self) -> Tuple:
        if self.isolated:
            self.logger.info(f'Executing APBS for simulation {self.isolated_type} {self.sim_name}')
        else:
            self.logger.info(f'Executing APBS for simulation {self.sim_name}')
        run_apbs = ['apbs',self.sim_files['apbs_in'].name]
        
        start_time = time.time()
        apbs_log = open(self.sim_files['apbs_log'],'w')
        self.logger.debug(f'Executing APBS for {self.sim_name}')
        result = subprocess.run(
            run_apbs,
            stdout=apbs_log,
            stderr=subprocess.DEVNULL,
            cwd=self.sim_path
        )
        apbs_log.close()
        wall_clock = time.time() - start_time
        if result.returncode != 0:
            raise RuntimeError(f'Error in running APBS. Return code: {result.returncode}')
        self._finished_ok = True
        return result.returncode, wall_clock # type: ignore