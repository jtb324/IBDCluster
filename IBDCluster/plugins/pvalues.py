from typing import Protocol
from models import Network
from dataclasses import dataclass, field
from numpy import float64
from scipy.stats import binomtest
from plugins import factory_register
import log

logger = log.get_logger(__name__)


@dataclass
class Network(Protocol):
    """
    General Interface to define what the Network object should
    look like
    """

    gene_name: str
    gene_chr: str
    network_id: int
    pairs: list = field(default_factory=list)
    iids: set[str] = field(default_factory=set)
    haplotypes: set[str] = field(default_factory=set)
    pvalues: None | str = None


@dataclass
class DataHolder(Protocol):
    affected_inds: dict[float, list[str]]
    phenotype_prevalence: dict[str, float]
    phenotype_description: dict[str, str] = None


@dataclass
class Pvalues:
    """Class that is responsible for determining the pvalues for each network"""

    name: str = "Pvalue plugin"

    @staticmethod
    def _determine_pvalue(
        phenotype: str,
        phenotype_percent: int,
        carriers_count: int,
        network_size: int,
    ) -> float64:
        """Function that will determine the pvalue for each network

        Returns

        float
            Returns the calculated pvalue
        """
        # the probability is 1 if the carrier count is zero because it is chances of finding
        # 0 or higher which is everyone
        if carriers_count == 0:
            logger.debug(f"carrier count = 0 therefore pvalue for {phenotype} = 1")
            return 1

        result = binomtest(carriers_count - 1, network_size, phenotype_percent)

        pvalue = result.pvalue

        logger.debug(f"pvalue for {phenotype} = {pvalue}")

        return pvalue

    def _determine_pvalues(
        self,
        carriers_list: dict[str, list[str]],
        network: Network,
        phenotype_percentages: dict[str, float],
    ) -> tuple[str, str]:
        """Function that will determine information about how many carriers are in each
        network, the percentage, the IIDs of the carriers in the network, and use this to calculate the pvalue for the network. The function keeps track of the smallest non-zero pvalue and returns it or NA

        Parameters
        ----------
        carriers_list : dict[str, list[str]]
            Dictionary that has all the carriers in list for each phenotype of interest

        network : Network
            Network objectattributes for iids, pairs, and haplotypes

        phenotype_percentages : dict[str, float]
            dictionary where the keys are phecode strings and the values are the phecode
            frequencies in the population

        Returns
        -------
        tuple(str, str)
            returns a tuple where the first element is the string of all pvalues for all
            phecodes and the second values it a string of the minimum phecode and pvalue
        """
        # dictionary that will contain the phecodes as keys
        # and the pvalues as floats
        pvalue_dictionary: dict[str, float] = {}

        output_str = ""

        cur_min_pvalue = 1

        cur_min_phecode = ""

        # iterating over each phenotype
        for phenotype, phenotype_freq in phenotype_percentages.items():
            # getting the list of iids in our population that carry the phenotype
            # POTENTIAL BUG: POTENTIALLY the program could fail here. It may be good ot add some error catching
            carriers = carriers_list[phenotype]
            # getting a list of iids in the network that are a carrier
            carriers_in_network: list[str] = [
                iid for iid in network.iids if iid in carriers
            ]

            num_carriers_in_network: int = len(carriers_in_network)
            # "ind_in_network", "percentage", "IIDs", "pvalue", "network_len"
            # we want to keep this value incase it could be added back to the program
            _percentage_in_network: float = num_carriers_in_network / len(network.iids)

            _network_size = len(network.iids)

            # calling the sub function that determines the pvalue
            pvalue: float = self._determine_pvalue(
                phenotype,
                phenotype_freq,
                num_carriers_in_network,
                len(network.iids),
            )

            # Next two lines create the string and then concats it to the output_str
            phenotype_str = f"{num_carriers_in_network}\t{pvalue}\t"

            output_str += phenotype_str

            # logging the string in debug mode. This logs the individual phenotype string not the total output for size
            logger.debug(
                f"network_id {network.network_id}: phenotype_str - {output_str}"
            )

            # Now we will see if the phecode is lower then the cur_min_pvalue. If it is then
            # we will change the cur_min_pvalue and we will update the cur_min_phecode
            if pvalue < cur_min_pvalue and pvalue != 0:
                cur_min_pvalue = pvalue

                cur_min_phecode = phenotype
        # if a minimum phecode is identified then we need to create a string, otherwise we
        # use N/A's
        if cur_min_phecode:
            min_phecode_str = f"{cur_min_pvalue}\t{cur_min_phecode}"
        else:
            min_phecode_str = "N/A\tN/A"

        # remove the trailing tab space
        output_str = output_str.rstrip("\t")
        # return the pvalue_output string first and either a tuple of N/As or the min pvalue/min_phecode
        return output_str + "\n", min_phecode_str

    @staticmethod
    def _get_descriptions(
        phecode_description: dict[str, dict[str, str]], min_phecode: str
    ) -> str:
        """Method to get the description for the minimum phecode

        Parameters
        ----------
        phecode_descriptions : dict[str, dict[str, str]]
            dictionary with descriptions of each phecode

        min_phecode : str
            minimum phecode string

        Returns
        -------
        str
            returns a string that has the phecode description
        """

        logger.debug(f"min phecode: {min_phecode}")
        # getting the inner dictionary if the key exists, otherwise getting
        # an empty dictionary
        desc_dict = phecode_description.get(min_phecode, {})
        # getting the phenotype string if key exists,
        # otherwise returns an empty string
        return desc_dict.get("phenotype", "N/A")

    def analyze(self, **kwargs) -> None:

        # this is the DataHolder model. We will use the networks, the affected_inds, and the phenotype_prevalances attribute
        data: DataHolder = kwargs["data"]
        network: Network = kwargs["network"]

        # Determining the pvalua and the tuple
        pvalue_str, min_pvalue_str = self._determine_pvalues(
            data.affected_inds,
            network,
            data.phenotype_prevalence,
        )

        min_phecode_description = self._get_descriptions(
            data.phenotype_description, min_pvalue_str.split("\t")[0]
        )

        network.pvalues = "\t".join(
            [min_pvalue_str, min_phecode_description, pvalue_str]
        )


def initialize() -> None:
    factory_register("pvalues", Pvalues)
