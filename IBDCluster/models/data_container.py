from dataclasses import dataclass


@dataclass
class DataHolder:
    gene_name: str
    chromosome: int
    affected_inds: dict[str, list[str]]
    phenotype_prevalence: dict[str, float]
    phenotype_cols: list[str]
    ibd_program: str
    phenotype_description: None | dict[str, str] = None
