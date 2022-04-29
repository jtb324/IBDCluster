from dataclasses import dataclass, field
from typing import Dict, List
from models import Network
import pandas as pd


@dataclass
class DataHolder:
    gene_name: str
    chromosome: int
    networks_list: List[Network]
    affected_inds: Dict[float, List[str]]
    phenotype_table: pd.DataFrame
    phenotype_cols: List[str]
    ibd_program: str
    phenotype_description: Dict[str, str] = None
    phenotype_percentages: Dict[str, float] = field(default_factory=dict)
    network_pvalues: Dict[int, Dict[str, float]] = field(default_factory=dict)
