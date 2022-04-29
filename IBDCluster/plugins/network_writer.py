from typing import List, Tuple, Dict, Protocol, Any, Set
from dataclasses import dataclass, field
import log
import pandas as pd
from scipy.stats import binom
import os
from tqdm import tqdm
from plugins import factory_register
import pathlib

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
    pairs: List = field(default_factory=list)
    iids: Set[str] = field(default_factory=set)
    haplotypes: Set[str] = field(default_factory=set)


@dataclass
class DataHolder(Protocol):
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


@dataclass
class NetworkWriter:
    """Class that is responsible for creating the _networks.txt file from the information provided"""

    name: str = "NetworkWriter plugin"

    @staticmethod
    def _form_header(phenotype_columns) -> str:
        """Method that will form the phenotype section of the header string. Need to
        append the words ind_in_network and pvalue to the phenotype name."""
        # pulling out all of the phenotype names from the carriers matrix

        column_list: List[str] = []

        for column in phenotype_columns:

            column_list.extend(
                [
                    column + ending
                    for ending in [
                        "_ind_in_network",
                        "_pvalue",
                    ]
                ]
            )

        return "\t".join(column_list)

    @staticmethod
    def _determine_pvalue(
        phenotype: str,
        phenotype_percent: int,
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
            logger.debug(f"carrier count = 0 therefore pvalue for {phenotype} = 1")
            return 1

        prob = 1 - binom.cdf(carriers_count - 1, network_size, phenotype_percent)

        logger.debug(f"pvalue for {phenotype} = {prob}")

        return prob

    @staticmethod
    def _check_min_pvalue(phenotype_pvalues: Dict[str, float]) -> Tuple[str, str]:
        """
        Function that will determine the smallest pvalue and return the corresponding phecode

        phenotype_pvalues : Dict[str, float]

            dictionary where the phecode strings are keys and the values are the
            pvalues as floats

        Returns

        Tuple[str, str]

            returns a tuple where the first value is the pvalue and the
            second is the phecode. If the minimum phecode is 1 (meaning that none of the phecodes had carriers) then the program returns N/A for both spots.
        """

        min_pvalue, min_phecode = min(
            zip(phenotype_pvalues.values(), phenotype_pvalues.keys())
        )

        if min_pvalue == 1:
            return "N/A", " N/A"
        else:
            return str(min_pvalue), min_phecode

    def _determine_pvalues(
        self,
        carriers_list: Dict[str, List[str]],
        network: Network,
        phenotype_list: List[str],
        percentage_pop_phenotypes: Dict[str, float],
    ) -> str:
        """Function that will determine information about how many carriers are in each
        network, the percentage, the IIDs of the carriers in the network, and use this to calculate the pvalue for the network. The function keeps track of the smallest non-zero pvalue and returns it or NA

        Parameters

        carriers_list : dict[str, List[str]]
            Dictionaary that has all the carriers in list for each phenotype of interest

        network : Network
            Network objectattributes for iids, pairs, and haplotypes

        Returns

        str
            returns a string that has the number of carriers and the pvalue for each phenotype
        """
        # dictionary that will contain the phecodes as keys
        # and the pvalues as floats
        pvalue_dictionary: Dict[str, float] = {}

        output_str = ""

        # iterating over each phenotype
        for phenotype in phenotype_list:
            # getting the list of iids in our population that carry the phenotype
            carriers = carriers_list[phenotype]
            # getting a list of iids in the network that are a carrier
            carriers_in_network: List[str] = [
                iid for iid in network.iids if iid in carriers
            ]

            num_carriers_in_network: int = len(carriers_in_network)
            # "ind_in_network", "percentage", "IIDs", "pvalue", "network_len"
            # we want to keep this value incase it could be added to the program
            _percentage_in_network: float = num_carriers_in_network / len(network.iids)

            _network_size = len(network.iids)

            # calling the sub function that determines the pvalue
            pvalue: float = self._determine_pvalue(
                phenotype,
                percentage_pop_phenotypes[phenotype],
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

            # updating the pvalue_dictionary for this phenotype
            pvalue_dictionary[phenotype] = pvalue

        # remove the trailing tab space
        output_str = output_str.rstrip("\t")
        # return the pvalue_output string first and either a tuple of N/As or the min pvalue/min_phecode
        return output_str + "\n", pvalue_dictionary

    @staticmethod
    def get_descriptions(
        phecode_desc: Dict[str, Dict[str, str]], min_phecode: str
    ) -> str:
        """Method to get the description for the minimum phecode

        Parameters
        ----------
        phecode_desc : Dict[str, Dict[str, str]]
            dictionary with descriptions of each phecode

        min_phecode : str
            minimum phecode string

        Returns
        -------
        str
            returns a string that has the phecode description
        """
        if phecode_desc:
            # getting the inner dictionary if the key exists, otherwise getting
            # an empty dictionary
            desc_dict = phecode_desc.get(min_phecode, {})
            # getting the phenotype string if key exists,
            # otherwise returns an empty string
            return desc_dict.get("phenotype", "N/A")

        return "N/A"

    def analyze(self, **kwargs) -> Dict[str, Any]:
        """main function of the plugin. It needs to determine the pvalue"""

        data: DataHolder = kwargs["data"]
        output_path = kwargs["output"]

        # This will be a list of strings that has the output for each network
        networks_analysis_list: List[str] = []

        for network in tqdm(
            data.networks_list, desc="networks with calculated pvalues: "
        ):

            # adding a key for gene id and the network id to the data.network_pvalues
            data.network_pvalues[data.gene_name] = {}
            # string that has the network information such as the
            # network_id, ibd_program, the gene it is for and the
            # chromosome number
            networks: str = f"{data.ibd_program}\t{data.gene_name}\t{network.network_id}\t{data.chromosome}"
            # string that has the number of individuals in the
            # network as well as the the number of haplotypes
            counts: str = f"{len(network.iids)}\t{len(network.haplotypes)}"
            # string that has the list of GRID IIDs and the haplotype phases
            iids: str = f"{', '.join(network.iids)}\t{', '.join(network.haplotypes)}"
            # Determining the pvalua and the tuple
            pvalue_str, phenotype_pvalue_dict = self._determine_pvalues(
                data.affected_inds,
                network,
                data.phenotype_cols,
                data.phenotype_percentages,
            )

            # getting a string that has the phecode and the minimum pvalue
            # for the network
            min_pvalue, min_phecode = self._check_min_pvalue(phenotype_pvalue_dict)

            phecode_desc = self.get_descriptions(
                data.phenotype_description, min_phecode
            )

            networks_analysis_list.append(
                f"{networks}\t{counts}\t{iids}\t{min_pvalue}\t{min_phecode}\t{phecode_desc}\t{pvalue_str}"
            )

            # adding the pvalue_dictionary to the network_pvalues attribute of the dataHolder
            data.network_pvalues[data.gene_name][
                network.network_id
            ] = phenotype_pvalue_dict

        # returning an object with the list of strings, the
        # path, and the gene name
        return {
            "output": networks_analysis_list,
            "path": os.path.join(output_path, data.gene_name),
            "gene": data.gene_name,
        }

    def write(self, **kwargs) -> None:
        """Method to write the output to a networks.txt file"""
        data = kwargs["input_data"]
        ibd_program = kwargs["ibd_program"]
        phenotype_list: List[str] = kwargs["phenotype_list"]

        gene_name = data["gene"]
        # forming the correct output path based on the gene name
        gene_output = data["path"]

        pathlib.Path(gene_output).mkdir(parents=True, exist_ok=True)

        output_file_name = os.path.join(
            gene_output,
            "".join([ibd_program, "_", gene_name, "_networks.txt"]),
        )

        logger.debug(f"Information written to a networks.txt at: {output_file_name}")
        # Opening the file and writing the head to it and then each network
        with open(
            output_file_name,
            "w+",
            encoding="utf-8",
        ) as output_file:
            output_file.write(
                f"program\tgene\tnetwork_id\tchromosome\tIIDs_count\thaplotypes_count\tIIDs\thaplotypes\tmin_pvalue\tmin_pvalue_phecode\tmin_pvalue_desc\t{self._form_header(phenotype_list)}\n"
            )
            # iterating over each network and writing the values to file
            for network in tqdm(data["output"], desc="Networks written to file: "):

                # if debug mode is choosen then it will write the output string to a file/console
                logger.debug(network)

                output_file.write(network)


def initialize() -> None:
    factory_register("network_writer", NetworkWriter)
