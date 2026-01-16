"""
apbsHelper.py

Python tool to aid in the running of multiple APBS simulations.
See README for configuration details.

Usage:
    apbsHelper YAML_FILE [YAML_FILE ...]

Arguments:
    YAML_FILE     Paths to config yaml files
"""
from src.classes.config_handler import ConfigHandler
from src.utils.setup_logs import setup_log
from pathlib import Path
from datetime import date
import pandas as pd
import random
import sys

RESULTS_PATH = Path('results')
RESULTS_PATH.mkdir(exist_ok=True,parents=True)

def generate_uid(length: int = 7) -> str:
    return ''.join(random.choices('0123456789abcdef', k=length))

def main() -> None:
    if len(sys.argv) < 2:
        print(__doc__.strip()) # type: ignore
        sys.exit(1)

    meshes_to_run = []
    for yaml in sys.argv[1:]:
        mesh = ConfigHandler(yaml)
        meshes_to_run.append(mesh)

    if len(meshes_to_run) == 1:
        for mesh in meshes_to_run:
            if mesh.config['config_only']:
                mesh.configure_all_cases()
            else:
                uid = generate_uid(8)
                out_csv = RESULTS_PATH / f'{date.today().strftime('%Y-%m-%d')}_{uid}_{mesh.config_name}.csv'
                mesh.run_all_cases()
                mesh.results.to_csv(out_csv,index=False)
    else:
        all_results = [] 
        uid = generate_uid(8)
        for mesh in meshes_to_run:
            if mesh.config['config_only']:
                mesh.configure_all_cases()
            else:
                out_csv = RESULTS_PATH / f'{date.today().strftime('%Y-%m-%d')}_{uid}_{mesh.config_name}.csv'
                mesh.run_all_cases()
                result = mesh.results
                all_results.append(result)
                result.to_csv(out_csv,index=False)
        if all_results:
            all_results_df = pd.concat(all_results)
            out_csv = RESULTS_PATH / f'{date.today().strftime('%Y-%m-%d')}_{uid}_allResults.csv'
            all_results_df.to_csv(out_csv,index=False)

if __name__ == "__main__":
    main()