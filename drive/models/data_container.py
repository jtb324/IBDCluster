from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Protocol

from .networks import Network_Interface


class Data_Interface(Protocol):
    """Protocol defining what attributes the DataHolder Interface needs to have"""

    networks: List[Network_Interface]
    output_path: Path
    carriers: Dict[str, Dict[str, List[str]]]
    phenotype_descriptions: Dict[str, Dict[str, str]]


@dataclass
class Data:
    """main class to hold the data from the network analysis and the different pvalues"""

    networks: List[Network_Interface]
    output_path: Path
    carriers: Dict[str, Dict[str, List[str]]]
    phenotype_descriptions: Dict[str, Dict[str, str]]
