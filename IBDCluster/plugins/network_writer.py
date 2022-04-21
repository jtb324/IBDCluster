from typing import List, Tuple, Dict, Protocol, Any, Set
from dataclasses import dataclass, field
import log
import pandas as pd
from scipy.stats import binom
from collections import namedtuple
import os
from tqdm import tqdm
from plugins import factory_register

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

        prob = str(1 - binom.cdf(carriers_count - 1, network_size, phenotype_percent))

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
    ) -> Tuple[str, Any]:
        """Function that will determine information about how many carriers are in each
        network, the percentage, the IIDs of the carriers in the network, and use this to calculate the pvalue for the network. The function keeps track of the smallest non-zero pvalue and returns it or NA

        Parameters

        carriers_list : dict[str, List[str]]
            Dictionaary that has all the carriers in list for each phenotype of interest

        network : Network
            Network objectattributes for iids, pairs, and haplotypes

        Returns

        str
            returns a string that has the number of carriers first and then the pvalue. These values are tab separated and end with a newline
        """

        output_str = ""
        # crea
        min_pvalue = 1
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
            phenotype_str = f"{num_carriers_in_network}\t{pvalue}"

            output_str += phenotype_str
            # logging the string in debug mode. This logs the individual phenotype string not the total output for size
            logger.debug(
                f"network_id {network.network_id}: phenotype_str - {phenotype_str}"
            )
        # return the pvalue_output string first and either a tuple of N/As or the min pvalue/min_phecode
        return output_str + "\n", ("N/A", "N/A") if min_pvalue == 1 else (
            min_pvalue,
            min_phecode,
        )

    def analyze(self, **kwargs) -> Tuple[int, Any]:
        """main function of the plugin. It needs to determine the pvalue"""

        data: DataHolder = kwargs["data_container"]
        output_path = kwargs["output"]

        # This iwll be a list of strings that has the output for each network
        output_dict: Dict[str, Dict[str, Any]] = {}
        # iterating over each gene
        for gene_info, network_list in data.networks_dict.items():

            networks_analysis_list: List[str] = []

            for network in tqdm(
                network_list, desc="networks with calculated pvalues: "
            ):
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
                pvalue_str, min_pvalue_tuple = self._determine_pvalues(
                    data.affected_inds,
                    network,
                    data.phenotype_cols,
                    data.phenotype_percentages,
                )

                min_pvalue_str = f"{min_pvalue_tuple[0]}\t{min_pvalue_tuple[1]}"

                networks_analysis_list.append(
                    f"{networks}\t{counts}\t{iids}\t{min_pvalue_str}\t{pvalue_str}\n"
                )

            output_dict[gene_info[0]] = {
                "output": networks_analysis_list,
                "path": os.path.join(output_path, gene_info[0]),
            }
            # TODO: need to return this values as well as the output path.

    def write(self, **kwargs) -> None:
        """Method to write the output to a networks.txt file"""
        data = kwargs["input_data"]
        ibd_program = kwargs["ibd_program"]
        phenotype_list: List[str] = kwargs["phenotype_list"]

        # Iterate over every gene so that we can write to different files.
        for gene_name in data.keys():
            # forming the correct output path based on the gene name
            output_file_name = os.path.join(
                data[gene_name]["path"],
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
                for network in tqdm(data["output"], desc="Networks written to file: "):

                    # if debug mode is choosen then it will write the output string to a file/console
                    logger.debug(network)

                    output_file.write(network)


def initialize() -> None:
    factory_register("network_writer", NetworkWriter)
