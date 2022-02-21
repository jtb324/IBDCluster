from typing import Protocol, Dict, Tuple, List, Set
from collections import namedtuple
import pandas as pd
from scipy.stats import binom
from models import PairWriter, NetworkWriter
import analysis
from tqdm import tqdm
import log

logger = log.get_logger(__name__)
# named tuple that will contain info about the carriers in a specific network
# as well the percentage, and the IIDs
class CarriersInfo(
    namedtuple("Carrier_Comp", ["ind_in_network", "percentage", "IIDs", "pvalue"])
):
    """This is an extension of the class tuple so that I can overwrite the __str__ method"""
    def __str__(self):
        return f"Carrier Info Object - Individuals in Network: {self.ind_in_network}, Percentages: {self.percentage}, IID List: {self.IIDs}, pvalue: {self.pvalue}"


class WriteObject(Protocol):
    """Interface indicating that the writer object needs to have the methods write to file"""

    def _form_phenotype_header(self) -> str:
        ...

    def write(self, **kwargs) -> None:
        ...


class Writer(Protocol):
    """Interface that defines that these objects need a write_to_file method
    and a set_writer method"""

    def set_writer(self, writer: WriteObject) -> None:
        ...

    def write_to_file(self, output: str, networks_dict: Dict) -> None:
        ...


def _generate_carrier_list(carriers_matrix: pd.DataFrame) -> Dict[str, List[str]]:
    """Function that will take the carriers_pheno_matrix and generate a dictionary that has the list of carriers for each phenotype"""
    return_dict = {}

    # iterating over each phenotype which starts with the
    # second column

    for column in carriers_matrix.columns[1:]:

        filtered_matrix: pd.DataFrame = carriers_matrix[carriers_matrix[column] == 1][
            ["grids", column]
        ]

        return_dict[column] = filtered_matrix.grids.values.tolist()

    logger.debug(f"identified carriers for {len(return_dict.keys())} phenotypes")

    return return_dict


def determine_pvalue(
    phenotype: str,
    percentage_pop_phenotypes: Dict[str, int],
    carriers_count: int,
    network_size: int,
) -> float:
    """Function that will determine the pvalue for each network

    Returns

    float
        Returns the calculated pvalue
    """
    # the probability is 1 if the carrier count is zero because it is chances of finding
    # 0 or higher which is everyone
    if carriers_count == 0:
        return 1
    else:
        prob: float = 1 - binom.cdf(
            carriers_count - 1, network_size, percentage_pop_phenotypes[phenotype]
        )

        return prob


# FIXME: This needs to be refactored. The function is doing many different things
def determine_in_networks(
    carriers_list: Dict[str, List[str]],
    info_dict: Dict[int, Dict],
    percentage_pop_phenotypes: Dict[str, float],
) -> None:
    """Function that will determine information about how many carriers are in each
    network, the percentage, the IIDs of the carriers in the network, and the pvalue for
    the network. These values are stored in the namedtuple CarriersInfo, and then
    finally stored in the phenotype key of the inner info_dict dictionary

    Parameters

    carriers_list : dict[str, List[str]]
        Dictionaary that has all the carriers in list for each phenotype of interest

    info_dict : Dict[int, Dict]
        dictionary where the keys are network ids and the inner dictionary has information about the network ind different keys
    """

    for _, info in tqdm(info_dict.items(), desc="Analyzing networks: "):

        iids_in_network: Set[str] = info["IIDs"]

        phenotype_info: Dict[str, CarriersInfo] = {}

        for phenotype, carriers in carriers_list.items():

            carriers_in_network: List[str] = [
                iid for iid in iids_in_network if iid in carriers
            ]

            num_carriers_in_network: int = len(carriers_in_network)

            carrier_info: CarriersInfo = CarriersInfo(
                num_carriers_in_network,
                num_carriers_in_network / len(iids_in_network),
                carriers_in_network,
                determine_pvalue(
                    phenotype,
                    percentage_pop_phenotypes,
                    num_carriers_in_network,
                    len(iids_in_network),
                ),
            )

            phenotype_info[phenotype] = carrier_info

            logger.debug(carrier_info)

        info["phenotype"] = phenotype_info


def _determine_pair_carrier_status(
    phenotype_list: List[str], carrier_list: Dict[str, List[str]], pair_id: str
) -> str:
    """Function that will iterate over the phenotype list and get a list of carriers for each one and check if the pair id is in each list

    Parameters

    phenotype_list : List[str]
        List of phenotypes. This gives us a consistent order for the columns in
        the file

    carrier_list : Dict[str, List[str]]
        Dictionary that has the phenotypes as keys and list of IIDs that carry that phenotype as values

    pair_id : str
        id of the pair of interest

    Returns

    str
        returns a string with either 0's or 1's for whether the individual is a carrier or not. These values will be tab separated
    """
    # iterrating over each phenotype to get the carriers
    carrier_str: str = ""

    for phenotype in phenotype_list:

        carriers: List[str] = carrier_list[phenotype]

        if pair_id in carriers:
            logger.debug(f"id {pair_id} is in the list {', '.join(carriers)}")
            carrier_str += "1\t"

        else:
            logger.debug(f"id {pair_id} is not in the list {', '.join(carriers)}")
            carrier_str += "0\t"

    carrier_str = carrier_str.strip("\t")

    return carrier_str


def analyze_pair_carrier_status(
    phenotype_list: List[str],
    carrier_list: Dict[str, List[str]],
    info_dict: Dict[int, Dict],
) -> None:
    """Function that will iterate through each pair in a network and add the carrier status for each id for each phenotype

    Parameters

    phenotype_list : List[str]
        List of phenotypes. This gives us a consistent order for the columns in
        the file

    carriers_list : dict[str, List[str]]
        Dictionaary that has all the carriers in list for each phenotype of interest

    info_dict : Dict[int, Dict]
        dictionary where the keys are network ids and the inner dictionary has information about the network ind different keys
    """
    for _, info in info_dict.items():

        # getting all of the pair objects and iterate over them
        pair_list = info["pairs"]

        for pair in pair_list:

            carrier_str_1: str = _determine_pair_carrier_status(
                phenotype_list, carrier_list, pair.pair_1
            )

            carrier_str_2: str = _determine_pair_carrier_status(
                phenotype_list, carrier_list, pair.pair_2
            )

            pair.carrier_str_1 = carrier_str_1

            pair.carrier_str_2 = carrier_str_2


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

    # getting the percentage of carriers for each phenotype per the input population
    percent_carriers_in_pop: Dict[str, float] = analysis.get_percentages(
        carriers_pheno_matrix
    )
    # recording this percentage to keep track of it. CONSIDER JUST
    # ADDING THIS AS ANOTHER CLASS WITHIN THE WRITER
    analysis.write_to_file(percent_carriers_in_pop, writer.output)

    phenotype_list: List[str] = list(carriers_lists.keys())
    # we will iterate over each network and basically
    # determine the number of carriers for each grid as well
    # as the percent
    determine_in_networks(carriers_lists, network_info, percent_carriers_in_pop)

    # creating the _networks.txt file using the writer object
    writer.set_writer(NetworkWriter(gene_info[0], gene_info[1], phenotype_list))

    writer.write_to_file(network_info)

    # adding information to the pair about the carrier string
    analyze_pair_carrier_status(phenotype_list, carriers_lists, network_info)

    # creating the allpairs.txt file
    writer.set_writer(PairWriter(gene_info[0], gene_info[1], phenotype_list))

    writer.write_to_file(network_info)
