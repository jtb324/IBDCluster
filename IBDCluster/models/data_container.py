from dataclasses import dataclass, field
from typing import Dict, Tuple, List
from models import Network
import pandas as pd


@dataclass
class DataHolder:
    networks_dict: Dict[Tuple[str, int], List[Network]]
    affected_inds: Dict[float, List[str]]
    phenotype_table: pd.DataFrame
    phenotype_cols: List[str]
    ibd_program: str
    phenotype_percentages: Dict[str, float] = field(default_factory=dict)
