from dataclasses import dataclass, field
from models import Network
import pandas as pd


@dataclass
class DataHolder:
    gene_name: str
    chromosome: int
    networks_list: list[Network]
    affected_inds: dict[float, list[str]]
    phenotype_table: pd.DataFrame
    phenotype_cols: list[str]
    ibd_program: str
    phenotype_description: dict[str, str] = None
    phenotype_percentages: dict[str, float] = field(default_factory=dict)
    network_pvalues: dict[int, dict[str, float]] = field(default_factory=dict)
