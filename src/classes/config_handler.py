from src.classes.handler_functions import HandlerFunction
from src.classes.case_handler import CaseHandler
from src.utils.setup_logs import setup_log
from typing import Dict, List
from pathlib import Path
import pandas as pd
import importlib
import logging
import yaml

class ConfigHandler:

    def __init__(
            self,
            yaml_file: Path|str,
            log_level=logging.INFO
    ) -> None:
        yaml_file = Path(yaml_file)
        self.config = self._ingest_config(yaml_file)
        if 'log_level' in self.config.keys():
            if self.config['log_level'] == 'info':
                log_level = logging.INFO
            elif self.config['log_level'] == 'debug':
                log_level = logging.DEBUG
            else:
                raise ValueError("Unknown log level in config. Must be 'info' or 'debug'")
        setup_log(log_level)
        self.logger = logging.getLogger(self.__class__.__name__)
        self.config_name = self.config['name']
        self.pdb_handler_func = self._import_pdb_handler()
        self.ran_ok = False
        self.configured_ok = False
        self.root_path = Path('data',self.config_name)
        self.core_pdb_path = Path(self.config['pdb_directory'],
                                  self.config['core_pdb_filename'])
        self._create_top_level()
        self.cases = self._create_cases()
        if not self.core_pdb_path.exists():
            raise ValueError('Missing core pdb file. Check path')
        pass

    def __repr__(self) -> str:
        return f'<ConfigHandler(name={self.config_name}, cases={len(self.cases)})'

    def _ingest_config(self, yaml_file: Path) -> Dict:
        with open(yaml_file,'r') as f:
            config = yaml.safe_load(f)
        self.config_ok = self._verify_config(config)
        return config

    
    def _import_pdb_handler(self) -> HandlerFunction:
        module_name, class_name = self.config['pdb_handler'].rsplit('.',1)
        module = importlib.import_module(module_name)
        user_class = getattr(module,class_name)

        if not issubclass(user_class,HandlerFunction):
            raise TypeError(f'{class_name} is not a HandlerFunction subclass')
        return user_class() # type: ignore

    def _create_cases(self) -> Dict[str,CaseHandler]:
        self.logger.info(f'Creating cases for config {self.config_name}')
        get_options = lambda key,opt: opt[key] if key in opt.keys() else {}
        global_handler_opts = get_options('global_handler_options',self.config)
        global_sim_opts = get_options('global_sim_options',self.config)

        cases = {}
        for case_config in self.config['cases']:
            case_handler_opts = get_options('handler_options',case_config) 
            case_sim_opts = get_options('sim_options',case_config)

            if case_config['name'] not in self.config['cases_to_run']:
                continue
            handler_options = self._merge_config(global_handler_opts,case_handler_opts)
            sim_options = self._merge_config(global_sim_opts,case_sim_opts)

            case = CaseHandler(
                case_name=case_config['name'],
                case_path=Path(self.root_path / case_config['name']),
                core_pdb_path=self.core_pdb_path,
                handler_func=self.pdb_handler_func,
                mesh_size=sim_options['mesh_size'],
                coarse_grain=sim_options['coarse_grain'],
                keep_complex_dx=sim_options['keep_complex_dx'],
                max_threads=self.config['max_threads'],
                linear_pb=sim_options['linear_pb'],
                handler_options=handler_options.copy()
            )
            self.logger.debug(f'Created case {case.case_name} for config {self.config_name}')
            cases[case.case_name] = case
        return cases

    def _create_top_level(self) -> None:
        mesh_path = self.root_path
        self.logger.info(f'Creating top level directory for config {self.config_name}')
        mesh_path.mkdir(exist_ok=True,parents=True)
        return None

    @staticmethod 
    def _verify_config(config_yaml: Dict) -> bool:
        is_ok = True

        needed_keys = [
            'name',
            'config_only',
            'max_threads',
            'cases_to_run',
            'pdb_handler',
            'pdb_directory',
            'core_pdb_filename',
            'cases'
        ]
        missing_keys = []
        for key in needed_keys:
            if key not in config_yaml.keys():
                missing_keys.append(key)
        
        if  missing_keys:
            is_ok = False
            raise RuntimeError(f'Missing config options: {missing_keys}')

        return is_ok
    
    @staticmethod
    def _merge_config(global_config: Dict, case_config: Dict) -> Dict:
        result = global_config.copy()
        for key, value in case_config.items():
            if (key in result
                and isinstance(result[key],dict)
                and isinstance(value,dict)
            ):
                result[key] = ConfigHandler._merge_config(result[key],value)
            else:
                result[key] = value
        return result

    def configure_all_cases(self) -> None:

        for name,case in self.cases.items():
            self.logger.info(f'Configuring case {name}')
            case.configure_case_sims(save_wall_distances=True)
        self.configured_ok = True
        return None

    def run_all_cases(self) -> None:

        for name, case in self.cases.items():
            self.logger.info(f'Running case {name}')
            case.run_case_sims()
        self.configured_ok = True
        self.ran_ok = True
        return None

    @property
    def results(self) -> pd.DataFrame:

        if not self.ran_ok:
            return pd.DataFrame()

        results_df = []
        for case in self.cases.values():
            results_df.append(case.case_results)

        results_df = pd.concat(results_df,ignore_index=True,sort=False)
        return results_df