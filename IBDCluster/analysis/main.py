import os
import pandas as pd
from typing import Protocol, Dict, Tuple, List, Set
from collections import namedtuple

# named tuple that will contain info about the carriers in a specific network as well the percentage, and the IIDs
Carriers_Comp = namedtuple("Carrier_Comp", ["ind_in_network", "percentage", "IIDs"])


class Write_Object(Protocol):
    """Interface indicating that the writer object needs to have the methods write to file"""

    def _form_phenotype_header(self) -> str:
        ...

    def write(self, **kwargs) -> None:
        ...


class Writer(Protocol):
    def set_writer(self, writer: Write_Object) -> None:
        ...

    def write_to_file(self, output: str, networks_dict: Dict) -> None:
        ...


def _generate_carrier_list(carriers_matrix: pd.DataFrame) -> Dict[str, List[str]]:
    """Function that will take the carriers_pheno_matrix and generate a dictionary that has the list of carriers for each phenotype"""
    return_dict = {}

    # iterating over each phenotype which starts with the
    # second column
    for column in carriers_matrix[carriers_matrix.columns[1:]]:

        filtered_matrix: pd.DataFrame = carriers_matrix[carriers_matrix[column] == 1][
            "grids", column
        ]

        return_dict[column] = filtered_matrix.grids.values.tolist()

    return return_dict


def determine_in_networks(
    carriers_list: Dict[str, List[str]], iids: Set[str]
) -> Dict[str, Carriers_Comp]:
    """Method that will determine how many carriers for a specific phenotype are in a network and then create a Carriers_Comp tuple with the info"""
    ...


def analyze(
    gene_info: Tuple[str, str],
    network_info: Dict[int, Dict],
    carriers_pheno_matrix: pd.DataFrame,
    writer: Writer,
):
    """Main function from the analyze module that will determine the pvalues and the number of haplotypes, and individuals and will write this all to file.

    Parameters

    gene_info : Tuple[str, str]
        Tuple that contains the gene name and the chromosome number

    network_info : Dict[int, Dict]
        Dictionary that has the network id as the outer key,
        and the inner dictionary has information about the
        individuals in the network, the haplotypes, and the
        pairs

    carriers_pheno_matrix : pd.DataFrame
        dataframe matrix that has the carrier status for each
        individual for each phenotype

    writer : Writer
        object that will be used to write output to both the
        networks.txt file as well as the allpairs.txt file.
        This object follows the Writer protocol above
    """

    # generate list of carriers for each phenotype
    carriers_lists: Dict[int, List[str]] = _generate_carrier_list(carriers_pheno_matrix)

    # we will iterate over each network and basically
    # determine the number of carriers for each grid as well
    # as the percent
    for network_id, info in network_info.items():

        ...
