from src.classes.simulation_handler import SimulationHandler
from src.classes.handler_functions import HandlerFunction
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, Tuple, Optional, List
from src.utils.setup_logs import setup_log
from string import Template
from pathlib import Path
import pandas as pd
import logging

class CaseHandler:
    def __init__(
            self,
            case_name: str,
            case_path: str|Path,
            core_pdb_path: str|Path,
            handler_func: HandlerFunction,
            mesh_size: float,
            coarse_grain: bool=False,
            keep_complex_dx: bool=False,
            max_threads: int=1,
            linear_pb: bool=False,
            log_level:int=logging.INFO,
            handler_options: Optional[Dict]=None,
            sim_name_format: Optional[str]=None
    ) -> None:
        setup_log(log_level)
        self.logger = logging.getLogger(self.__class__.__name__)
        self.case_name = case_name
        self.case_path = Path(case_path)
        self.core_pdb_path = Path(core_pdb_path)
        self.mesh_size = mesh_size
        self.finished_run = False
        self.ready_to_run = False
        self.left_isolated_sims = {}
        self.right_isolated_sims = {}
        self.complex_sims = {}
        self.coarse_grain = coarse_grain
        self.keep_dx = keep_complex_dx
        self.max_threads = max_threads
        self.linear_pb = linear_pb
        self.handler_func = handler_func

        if handler_options is None:
            self.handler_options = {}
            self.logger.warning('Recieved empty handler options. Are you sure?')
        else:
            self.handler_options = handler_options

        if sim_name_format is None:
            sim_name_format = 'N${N}_${variant}'

        self._create_simulation_objects(self.handler_options,coarse_grained=coarse_grain)
        self.logger.debug(f'Created {self}')
        pass

    def __repr__(self) -> str:
        return f'<CaseHandler(name={self.case_name},path={self.case_path})>'

    def _create_simulation_objects(
            self,
            handler_options: Dict,
            sim_name_format: str='N${N}_${variant}',
            coarse_grained:bool=False
    ) -> Tuple[Dict,Dict,Dict]:

        left_isolated_runners = {}
        right_isolated_runners = {}
        complex_runners = {}

        def _format_variants(x):
            parts = []
            for key,value in x.items():
                value = str(value)
                parts.append(f'{key}_{value}')
            return '_'.join(parts)
        
        values = list(handler_options['variants'].values())
        keys = list(handler_options['variants'].keys())
        for v in handler_options['variants'].values():
            assert isinstance(v,list), f'Expected list, got {type(v)}: {v}'

        for i,variant in enumerate(zip(*values)):
            variant = dict(zip(keys,variant))
            temp_handler_options = {}
            temp_handler_options.update(self.handler_options['base'])
            temp_handler_options.update(variant)
            fmt = Template(sim_name_format)
            values = {
                'N':i,
                'variant':_format_variants(variant),
            }
            sim_name = fmt.safe_substitute(values)

            complex_sim_path = self.case_path / 'complex' / sim_name
            iso_left_sim_path = self.case_path / 'isolatedLeft' / sim_name
            iso_right_sim_path = self.case_path / 'isolatedRight' / sim_name

            if coarse_grained:
                sim_type = 'coarseGrained'
            else:
                sim_type = 'allAtom'

            left_isolated_sim = SimulationHandler(
                iso_left_sim_path,
                self.mesh_size,
                sim_name,
                self.core_pdb_path,
                sim_type=sim_type,
                isolated=True,
                handler_func=self.handler_func,
                handler_options=temp_handler_options,
                keep_dx=self.keep_dx,
                linear_pb=self.linear_pb,
                isolated_type='left'
            )

            right_isolated_sim = SimulationHandler(
                iso_right_sim_path,
                self.mesh_size,
                sim_name,
                self.core_pdb_path,
                sim_type=sim_type,
                isolated=True,
                handler_func=self.handler_func,
                handler_options=temp_handler_options,
                keep_dx=self.keep_dx,
                linear_pb=self.linear_pb,
                isolated_type='right'
            )

            complex_sim = SimulationHandler(
                complex_sim_path,
                self.mesh_size,
                sim_name,
                self.core_pdb_path,
                sim_type=sim_type,
                isolated=False,
                handler_func=self.handler_func,
                handler_options=temp_handler_options,
                keep_dx=self.keep_dx,
                linear_pb=self.linear_pb,
                iso_left_sim=left_isolated_sim,
                iso_right_sim=right_isolated_sim
            )

            left_isolated_runners[sim_name] = left_isolated_sim
            right_isolated_runners[sim_name] = right_isolated_sim
            complex_runners[sim_name] = complex_sim

        self.left_isolated_sims = left_isolated_runners
        self.right_isolated_sims = right_isolated_runners
        self.complex_sims = complex_runners
        return left_isolated_runners,right_isolated_runners,complex_runners

    def configure_case_sims(self,save_wall_distances: bool=True) -> None:
        self.logger.info(f'Configuring simulations for case: {self.case_name}')
        self.case_path.mkdir(exist_ok=True, parents=True)
        isolated_runners = []
        isolated_runners.extend(self.left_isolated_sims.values())
        isolated_runners.extend(self.right_isolated_sims.values())

        if self.keep_dx:
            self.logger.info(f'Keeping dx files for complex simulations')

        self.logger.info(f'Configuring isolated simulations for case: {self.case_name}')
        isolated_to_run = [runner for runner in isolated_runners if not runner.sim_exists]
        self._run_parallel_simulations(isolated_to_run,config_only=True)

        complex_runners = []
        complex_runners.extend(self.complex_sims.values())
        self.logger.info(f'Configuring complex simulations for case: {self.case_name}')
        complex_to_run = [runner for runner in complex_runners if not runner.sim_exists]
        self._run_parallel_simulations(complex_to_run,config_only=True)

        if save_wall_distances:
            self.wall_distances.to_csv(self.case_path / f'{self.case_name}_wall_dst.csv',index=False)
        self.ready_to_run = True

        return None

    def _run_parallel_simulations(
            self,
            runners: List[SimulationHandler],
            config_only: bool=False
        ) -> None:

        if not runners:
            return None

        with ThreadPoolExecutor(max_workers=self.max_threads) as executor:

            #Creates dictionary of runners to execute.
            if config_only:
                future_to_runner = {executor.submit(runner.configure_simulation):
                                    runner for runner in runners}
                self.logger.debug(f'Configuring {len(future_to_runner)} simulations')
            else:
                future_to_runner = {executor.submit(runner.run_simulation):
                                    runner for runner in runners}
                self.logger.debug(f'Running {len(future_to_runner)} simulations')
            
            for future in as_completed(future_to_runner):
                runner = future_to_runner[future]
                try:
                    future.result()
                except Exception as e:
                    self.logger.error(f'Encountered error while running {runner.sim_name} - {runner.sim_path}. {e}')

            return None

    def run_case_sims(self) -> None:
        if not self.ready_to_run:
            self.configure_case_sims()
        self.logger.info(f'Running simulations for case: {self.case_name}')
        isolated_runners = []
        isolated_runners.extend(self.left_isolated_sims.values())
        isolated_runners.extend(self.right_isolated_sims.values())

        if self.keep_dx:
            self.logger.info(f'Keeping dx files for complex simulations')

        self.logger.info(f'Running isolated simulations for case: {self.case_name}')
        isolated_to_run = [runner for runner in isolated_runners if not runner.sim_exists]
        self._run_parallel_simulations(isolated_to_run,config_only=False)

        complex_runners = []
        complex_runners.extend(self.complex_sims.values())
        self.logger.info(f'Running complex simulations for case: {self.case_name}')
        complex_to_run = [runner for runner in complex_runners if not runner.sim_exists]
        self._run_parallel_simulations(complex_to_run,config_only=False)
        return None
    
    @property
    def case_results(self) -> pd.DataFrame:

        if not self.complex_sims:
            return pd.DataFrame()
        
        results = {
            'core_pdb':[],
            'case_name':[],
            'sim_name':[],
            'wall_dst':[],
            'status':[]
        }

        energy_keys = [
                'solv',
                'col',
                'solv_left',
                'col_left',
                'solv_right',
                'col_right',
                'dd_g_solv',
                'dd_g_col',
                'binding'
            ]

        #Add handler base options as columns in DataFrame
        if self.handler_options['base']:
            base_columns = []
            for column in self.handler_options['base'].keys():
                base_columns.append(column)
                results[column] = []

        #Add handler variants as columns in DataFrame
        if self.handler_options['variants']:
            variant_columns = []
            for column in self.handler_options['variants'].keys():
                variant_columns.append(column)
                results[column] = []
                
        for key in energy_keys:
            results[key] = []

        results['case_path'] = [] 
        results['sim_path'] = []
        results['iso_left_path'] = []
        results['iso_right_path'] = []

        for name, runner in self.complex_sims.items():
            result = runner.result
            results['core_pdb'].append(self.core_pdb_path.name)
            results['case_name'].append(self.case_name)
            results['sim_name'].append(name)
            results['wall_dst'].append(runner.wall_distance)
            results['case_path'].append(self.case_path)
            results['sim_path'].append(result.sim_path)
            results['iso_left_path'].append(runner.iso_left_sim.sim_path)
            results['iso_right_path'].append(runner.iso_right_sim.sim_path)

            if result.finished_ok:
                results['status'].append('OK')
            else:
                results['status'].append('ERROR')

            if self.handler_options['base']:
                for column in base_columns: # type: ignore
                    results[column].append(result.misc[column])
            if self.handler_options['variants']: 
                for column in variant_columns: # type: ignore
                    results[column].append(result.misc[column])
            
            for key in energy_keys:
                if result.finished_ok:
                    results[key].append(result.energies[key])
                else:
                    results[key].append(0)

        results_df = pd.DataFrame(results)
        return results_df

    @property
    def wall_distances(self) -> pd.DataFrame:

        if not self.complex_sims:
            return pd.DataFrame()

        wall_distances = {
            'sim_name':[],
            'wall_dst':[]
        }

        #Add handler base options as columns in DataFrame
        if self.handler_options['base']:
            base_columns = []
            for column in self.handler_options['base'].keys():
                base_columns.append(column)
                wall_distances[column] = []

        #Add handler variants as columns in DataFrame
        if self.handler_options['variants']:
            variant_columns = []
            for column in self.handler_options['variants'].keys():
                variant_columns.append(column)
                wall_distances[column] = []

        wall_distances['sim_path'] = []

        for name, runner in self.complex_sims.items():
            wall_distances['sim_name'].append(name)
            wall_distances['wall_dst'].append(runner.wall_distance)
            wall_distances['sim_path'].append(runner.sim_path)

            if self.handler_options:
                for column in base_columns: # type: ignore
                    wall_distances[column].append(runner.handler_options[column])
                for column in variant_columns: # type: ignore
                    wall_distances[column].append(runner.handler_options[column])
        
        wall_dst_df = pd.DataFrame(wall_distances)
        return wall_dst_df 