from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Tuple

class HandlerFunction(ABC):

    @abstractmethod
    def handle(
        self,
        core_pdb: Path,
        target_dir: Path,
        handler_options: Dict,
        coarse_grain: bool=False
    ) -> Tuple[Path, Dict[str,List[str]]]:
        
        pass