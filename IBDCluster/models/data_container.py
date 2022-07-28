from dataclasses import dataclass, field
from models import Network
import pandas as pd


@dataclass
class DataHolder:
    gene_name: str
    chromosome: int
    networks: Network
    affected_inds: dict[float, list[str]]
    phenotype_prevalance: dict[str, float]
    phenotype_cols: list[str]
    ibd_program: str
    phenotype_description: None | dict[str, str] = None
    network_pvalues: dict[int, dict[str, float]] = field(default_factory=dict)
