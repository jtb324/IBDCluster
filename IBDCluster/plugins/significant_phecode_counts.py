# pluging that will be used to form an output that has the
# phecode, the # of networks that pass bonferroni, the number
# of networks that pass 10^-5, and the phecode description
from dataclasses import dataclass, field
from typing import Protocol, Any
import pandas as pd
from models import Network
from tqdm import tqdm
import os


@dataclass
class DataHolder(Protocol):
    gene_name: str
    chromosome: int
    networks_list: list[Network]
    affected_inds: dict[float, list[str]]
    phenotype_table: pd.DataFrame
    phenotype_cols: list[str]
    phenotype_description: dict[str, str] = None
    network_pvalues: dict[int, dict[str, float]] = field(default_factory=dict)


@dataclass
class PhecodeCounter:
    """Class object that will handle the counting of each phecode and writting the counts out for all networks"""

    name: str = "Significant Phecode Counter"

    @staticmethod
    def _get_min_phecode(pvalue_dict: dict[str, float]) -> dict[str, list | float]:
        """Method that will return the the networks that have the min

        Parameters
        ----------
        pvalue_dict : dict[str, float]
            dictionary that has all the phecodes for each
            network. The keys are the phecodes and the values
            are the pvalues

        Returns
        -------
        dict[str, list | float]
            returns a dictionary where the key is 'phecodes'
            or 'min_pvalue' and the valus are a list of
            phecodes or the pvalue float
        """

        min_pvalue: float = min(pvalue_dict)

        # if the min_pvalue is 1 then the network had no
        # affected individuals and the dictionary returns an
        # empty list
        if min_pvalue == 1:
            return {"phecodes": [], "min_pvalue": 1}
        # filtering down the dictionary to all the phecodes
        # that could have the same
        min_pvalue_dict = {
            phecode: value
            for phecode, value in pvalue_dict.items()
            if value == min_pvalue
        }

        return {"phecodes": list(min_pvalue_dict.keys()), "min_pvalue": min_pvalue}

    @staticmethod
    def _generate_counts_structure(
        phecode_desc: dict[str, str]
    ) -> dict[str, dict[str, str | int]]:
        """Method that will create the data structure that will have hte phecode counts

        Parameters
        ----------
        phecode_desc : dict[str, str]
            dictionary where the keys are the phecodes and the values are the phecode descriptions

        Returns
        -------
        dict[str, dict[str | int]]
            returns a dictioary where the key is the phecode
            with an inner dictionary that has the values
            "count_bonferroni", "count_under", "desc"
        """
        counts_structure = {}

        for phecode, desc in phecode_desc.items():
            counts_structure[phecode] = {
                "count_above_bonf": 0,
                "count_above_minus5": 0,
                "desc": desc,
            }

        return counts_structure

    def analyze(self, **kwargs) -> dict[str, Any]:
        """Function that will determine the phecode counts for all the networks"""

        data: DataHolder = kwargs["data"]
        output_path: str = kwargs["output"]

        # iterating over each network to get the different pvalues
        for network_dict in tqdm(
            data.network_pvalues,
            desc="Networks that phecode counts have been determined",
        ):
            # getting a dictionary with the most significant
            # phecodes and the most significant pvalue
            min_phecodes = self._get_min_phecode(network_dict)

            # create a function that takes the phecode
            # descriptions and then creates a dictionary for
            # the counts
            self._generate_counts_structure(data.phenotype_description)

        # returning an object with the list of strings, the
        # path, and the gene name
        return {
            "output": ...,
            "path": os.path.join(output_path, data.gene_name),
            "gene": data.gene_name,
        }
