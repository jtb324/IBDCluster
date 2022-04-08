from typing import Protocol, Dict, Tuple, List
from collections import namedtuple
import pandas as pd
from scipy.stats import binom
from tqdm import tqdm
from models import PairWriter, NetworkWriter, Network
import analysis
import log


logger = log.get_logger(__name__)

# named tuple that will contain info about the carriers in a specific network
# as well the percentage, and the IIDs
class CarriersInfo(
    namedtuple(
        "Carrier_Comp",
        ["ind_in_network", "percentage", "IIDs", "pvalue", "network_len"],
    )
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
        """Method that determine the writer for the appropriate output"""
        ...

    def write_to_file(
        self,
        networks_list: List[Network] = None,
    ) -> None:
        """Method that will call the write method of the writer object"""
        ...


def determine_pvalue(
    phenotype: str,
    percentage_pop_phenotypes: Dict[str, int],
    carriers_count: int,
    network_size: int,
) -> str:
    """Function that will determine the pvalue for each network

    Returns

    float
        Returns the calculated pvalue
    """
    # the probability is 1 if the carrier count is zero because it is chances of finding
    # 0 or higher which is everyone
    if carriers_count == 0:
        logger.debug(f"carrier count = 0 therefore pvalue for {phenotype} = 1")
        return "1"

    prob = str(
        1
        - binom.cdf(
            carriers_count - 1, network_size, percentage_pop_phenotypes[phenotype]
        )
    )

    logger.debug(f"pvalue for {phenotype} = {prob}")
    return prob


def determine_in_networks(
    carriers_list: Dict[str, List[str]],
    network_list: List[Network],
    percentage_pop_phenotypes: Dict[str, float],
) -> None:
    """Function that will determine information about how many carriers are in each
    network, the percentage, the IIDs of the carriers in the network, and the pvalue for
    the network. These values are stored in the namedtuple CarriersInfo, and then
    finally stored in the phenotype attribute of the network object

    Parameters

    carriers_list : dict[str, List[str]]
        Dictionaary that has all the carriers in list for each phenotype of interest

    network_list : List[Network]
        List of Network objects which have attributes for iids, pairs, and haplotypes
    """

    for network in tqdm(network_list, desc="Analyzing networks: "):

        phenotype_info: Dict[str, CarriersInfo] = {}

        for phenotype, carriers in carriers_list.items():

            carriers_in_network: List[str] = [
                iid for iid in network.iids if iid in carriers
            ]

            num_carriers_in_network: int = len(carriers_in_network)

            carrier_info: CarriersInfo = CarriersInfo(
                num_carriers_in_network,
                num_carriers_in_network / len(network.iids),
                carriers_in_network,
                determine_pvalue(
                    phenotype,
                    percentage_pop_phenotypes,
                    num_carriers_in_network,
                    len(network.iids),
                ),
                len(network.iids),
            )

            phenotype_info[phenotype] = carrier_info

            logger.debug(carrier_info)

        network.phenotype = phenotype_info


def analyze(
    gene_info: Tuple[str, str],
    network_list: List[Network],
    carriers_pheno_matrix: pd.DataFrame,
    writer: Writer,
    carriers_lists: Dict[float, List[str]],
) -> None:
    """Main function from the analyze module that will determine the
    pvalues and the number of haplotypes, and individuals and will
    write this all to file.

    Parameters

    gene_info : Tuple[str, str]
        Tuple that contains the gene name and the chromosome number

    network_list : List[Network]
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

    carriers_list : Dict[float, List[str]]
        Dictionary where the keys are the phecodes and the values
        are list of IIDs affected with that phecode
    """

    phenotype_list: List[str] = carriers_pheno_matrix.columns[1:]

    # getting the percentage of carriers for each phenotype per the input population
    percent_carriers_in_pop: Dict[str, float] = analysis.get_percentages(
        carriers_pheno_matrix
    )

    # recording this percentage to keep track of it. CONSIDER JUST
    # ADDING THIS AS ANOTHER CLASS WITHIN THE WRITER
    analysis.write_to_file(percent_carriers_in_pop, writer.output)

    # we will iterate over each network and basically
    # determine the number of carriers for each grid as well
    # as the percent
    determine_in_networks(carriers_lists, network_list, percent_carriers_in_pop)

    # creating the _networks.txt file using the writer object
    writer.set_writer(NetworkWriter(gene_info[0], gene_info[1], phenotype_list))

    writer.write_to_file(networks_list=network_list)

    # creating the allpairs.txt file
    writer.set_writer(PairWriter(gene_info[0], gene_info[1], phenotype_list))

    writer.write_to_file(networks_list=network_list)
