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
    networks_dict: Dict[Tuple[str, int], List[Network]]
    affected_inds: Dict[float, List[str]]
    phenotype_table: pd.DataFrame
    phenotype_cols: List[str]
    ibd_program: str
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
    def _check_min_pvalue(
        min_pvalue: float, cal_pvalue: float, phenotype: str
    ) -> Tuple:
        """
        Function that will compare the min_pvalue to the calculate pvalue to see which is smaller
        """
        if cal_pvalue < min_pvalue and cal_pvalue != 0:
            min_pvalue = cal_pvalue
            return min_pvalue, phenotype
        else:
            return min_pvalue, "N/A"

    def _determine_pvalues(
        self,
        carriers_list: Dict[str, List[str]],
        network: Network,
        phenotype_list: List[str],
        percentage_pop_phenotypes: Dict[str, float],
    ) -> Dict[str, Any]:
        """Function that will determine information about how many carriers are in each
        network, the percentage, the IIDs of the carriers in the network, and use this to calculate the pvalue for the network. The function keeps track of the smallest non-zero pvalue and returns it or NA

        Parameters

        carriers_list : dict[str, List[str]]
            Dictionaary that has all the carriers in list for each phenotype of interest

        network : Network
            Network objectattributes for iids, pairs, and haplotypes

        Returns

        Dict[str, Any]
            returns a Dictionary that has the pvalue string,
            the min pvalue tuple, and the pvalue_dictionary as
            values
        """
        # dictionary that will contain the phecodes as keys
        # and the pvalues as floats
        pvalue_dictionary: Dict[str, float] = {}

        output_str = ""
        # create place holders for the min_phecode and the minimum pvalue
        min_pvalue: int = 1
        min_phecode = "N/A"

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

            min_pvalue, min_phecode = self._check_min_pvalue(
                min_pvalue, pvalue, phenotype
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
        return {
            "pvalue_str": output_str + "\n",
            "min_pvalue_str": ("N/A", "N/A")
            if min_pvalue == 1
            else (
                str(min_pvalue),
                min_phecode,
            ),
            "pvalue_dict": pvalue_dictionary,
        }

    def analyze(self, **kwargs) -> Tuple[int, Any]:
        """main function of the plugin. It needs to determine the pvalue"""

        data: DataHolder = kwargs["data"]
        output_path = kwargs["output"]

        # creating a dictionary that will have the pvalues for each network that way we can add them to teh data container
        # This iwll be a list of strings that has the output for each network
        output_dict: Dict[str, Dict[str, Any]] = {}
        # iterating over each gene
        for gene_info, network_list in data.networks_dict.items():

            networks_analysis_list: List[str] = []

            for network in tqdm(
                network_list, desc="networks with calculated pvalues: "
            ):
                # adding a key for gene id and the network id to the data.network_pvalues
                data.network_pvalues[gene_info[0]] = {}
                # string that has the network information such as the
                # network_id, ibd_program, the gene it is for and the
                # chromosome number
                networks: str = f"{data.ibd_program}\t{network.gene_name}\t{network.network_id}\t{network.gene_chr}"

                # string that has the number of individuals in the
                # network as well as the the number of haplotypes
                counts: str = f"{len(network.iids)}\t{len(network.haplotypes)}"
                # string that has the list of GRID IIDs and the haplotype phases
                iids: str = (
                    f"{', '.join(network.iids)}\t{', '.join(network.haplotypes)}"
                )
                # Determining the pvalua and the tuple
                pvalue_output: Dict[str, Any] = self._determine_pvalues(
                    data.affected_inds,
                    network,
                    data.phenotype_cols,
                    data.phenotype_percentages,
                )
                # pulling out the min_pvalue_tuple and the pvalue_str from the pvalue_output dictionary
                min_pvalue_tuple = pvalue_output["min_pvalue_str"]

                pvalue_str = pvalue_output["pvalue_str"]

                min_pvalue_str = f"{min_pvalue_tuple[0]}\t{min_pvalue_tuple[1]}"

                networks_analysis_list.append(
                    f"{networks}\t{counts}\t{iids}\t{min_pvalue_str}\t{pvalue_str}\n"
                )

                # adding the pvalue_dictionary to the network_pvalues attribute of the dataHolder
                data.network_pvalues[gene_info[0]][network.network_id] = pvalue_output[
                    "pvalue_dict"
                ]

            output_dict[gene_info[0]] = {
                "output": networks_analysis_list,
                "path": os.path.join(output_path, gene_info[0]),
            }

        return output_dict

    def write(self, **kwargs) -> None:
        """Method to write the output to a networks.txt file"""
        data = kwargs["input_data"]
        ibd_program = kwargs["ibd_program"]
        phenotype_list: List[str] = kwargs["phenotype_list"]

        # Iterate over every gene so that we can write to different files.
        for gene_name in data.keys():
            # forming the correct output path based on the gene name
            gene_output = data[gene_name]["path"]

            pathlib.Path(gene_output).mkdir(parents=True, exist_ok=True)

            output_file_name = os.path.join(
                gene_output,
                "".join([ibd_program, "_", gene_name, "_networks.txt"]),
            )

            logger.debug(
                f"Information written to a networks.txt at: {output_file_name}"
            )
            # Opening the file and writing the head to it and then each network
            with open(
                output_file_name,
                "w+",
                encoding="utf-8",
            ) as output_file:
                output_file.write(
                    f"program\tgene\tnetwork_id\tchromosome\tIIDs_count\thaplotypes_count\tIIDs\thaplotypes\tmin_pvalue\tmin_pvalue_phecode\t{self._form_header(phenotype_list)}\n"
                )
                # iterating over each network and writing the values to file
                for network in tqdm(
                    data[gene_name]["output"], desc="Networks written to file: "
                ):

                    # if debug mode is choosen then it will write the output string to a file/console
                    logger.debug(network)

                    output_file.write(network)


def initialize() -> None:
    factory_register("network_writer", NetworkWriter)
