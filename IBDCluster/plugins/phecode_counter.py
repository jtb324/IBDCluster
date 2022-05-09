# pluging that will be used to form an output that has the
# phecode, the # of networks that pass bonferroni, the number
# of networks that pass 10^-5, and the phecode description
from dataclasses import dataclass, field
from typing import Protocol, Any
import pandas as pd
from models import Network
from tqdm import tqdm
from plugins import factory_register
import os
import log
import pathlib

logger = log.get_logger(__name__)


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

    name: str = "Significant Phecode Counter plugin"

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

        min_pvalue: float = min(pvalue_dict.values())

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
            "count_above_bonf", "count_under", "desc"
        """
        logger.debug("Create a data structure for the phecode counts data structure")

        counts_structure = {}

        for phecode, desc in phecode_desc.items():
            counts_structure[phecode] = {
                "count_above_bonf": 0,
                "count_above_minus5": 0,
                "desc": desc["phenotype"],
            }

        return counts_structure

    @staticmethod
    def _count_significant_networks(
        min_phecodes: dict[str, list | float],
        phecode_count_struct: dict[str, dict[str, str | int]],
        BONF_CORR: float,
        ARBITARY_THRES: float,
    ) -> None:
        """Method that will update the phecode_count_struct for the min phecode

        Parameters
        ----------
        min_phecodes : dict[str, list | float]
            dictionary where the key is 'phecodes'
            or 'min_pvalue' and the values are a list of
            phecodes or the pvalue float

        phecode_count_struct : dict[str, dict[str, str | int]]
            dictioary where the key is the phecode
            with an inner dictionary that has the values
            "count_above_bonf", "count_above_minus5", "desc"
        """

        pvalue = min_phecodes["min_pvalue"]

        # if the pvalue is significant then you iterate over
        # the phecodes and adjust the count in the
        # phecode_count_struct
        if pvalue < BONF_CORR:

            for phecode in min_phecodes["phecodes"]:

                phecode_count_struct[phecode]["count_bonferroni"] += 1

        if pvalue < ARBITARY_THRES:

            for phecode in min_phecodes["phecodes"]:

                phecode_count_struct[phecode]["count_above_minus5"] += 1

    @classmethod
    def _convert_counts_to_str(
        cls, phecode_counts: dict[str, dict[str, str | int]]
    ) -> list[str]:
        """Method that will combine all the values of the phecode_dict into a string

        Parameters
        ----------
        phecode_counts: dict[str, dict[str, str | int]]
            dictioary where the key is the phecode
            with an inner dictionary that has the values
            "count_above_bonf", "count_above_minus5", "desc"

        Returns
        -------
        list[str]
            returns a list where each value is the network
            string
        """
        logger.debug("converted the phecode counts to a string")
        output_strings = []
        # iterating over each phecode to get the counts and form then into a stirng
        for phecode, counts_dict in phecode_counts.items():
            output_strings.append(
                f"{phecode}\t{counts_dict['count_above_bonf']}\t{counts_dict['count_above_minus5']}\t{counts_dict['desc']}\n"
            )

        return output_strings

    def analyze(self, **kwargs) -> dict[str, Any]:
        """Function that will determine the phecode counts for all the networks"""

        data: DataHolder = kwargs["data"]
        output_path: str = kwargs["output"]

        BONFERRONI_THRESHOLD: float = 1 / len(data.networks_list)

        ARBITARY_THRESHOLD: float = 1e-5
        # create a function that takes the phecode
        # descriptions and then creates a dictionary for
        # the counts
        all_networks_phecode_counts = self._generate_counts_structure(
            data.phenotype_description
        )

        # getting a phecode counts dictionary for all the
        # networks with more than 3 individuals
        networks_3inds_phecode_counts = self._generate_counts_structure(
            data.phenotype_description
        )
        # iterating over each network to get the different pvalues
        for network in tqdm(
            data.networks_list,
            desc="Networks that phecode counts have been determined",
        ):
            # getting a dictionary with the most significant
            # phecodes and the most significant pvalue
            min_phecodes_dict = self._get_min_phecode(network.pvalues)

            # match the length of the network iids. If it is
            # two then you count the phecodes in both the
            # all_networks_phecode_counts and the
            # networks_3indx_phecode_counts dict
            match len(network.iids):
                case 2:

                    logger.debug(
                        f"Identified that the network {network.network_id} was of size 2"
                    )

                    self._count_significant_networks(
                        min_phecodes_dict,
                        networks_3inds_phecode_counts,
                        BONFERRONI_THRESHOLD,
                        ARBITARY_THRESHOLD,
                    )
                    self._count_significant_networks(
                        min_phecodes_dict,
                        all_networks_phecode_counts,
                        BONFERRONI_THRESHOLD,
                        ARBITARY_THRESHOLD,
                    )
                case _:
                    logger.debug(
                        f"Identified that the network {network.network_id} was not of size 2"
                    )

                    self._count_significant_networks(
                        min_phecodes_dict,
                        networks_3inds_phecode_counts,
                        BONFERRONI_THRESHOLD,
                        ARBITARY_THRESHOLD,
                    )

        # turning the phecode counts into a list of strings
        all_networks_counts = self._convert_counts_to_str(all_networks_phecode_counts)

        greater_2inds_counts = self._convert_counts_to_str(
            networks_3inds_phecode_counts
        )

        # returning an object with the list of strings, the
        # path, and the gene name
        return {
            "output": {
                "all_networks": all_networks_counts,
                "greater_2inds_networks": greater_2inds_counts,
            },
            "path": os.path.join(output_path, data.gene_name),
            "gene": data.gene_name,
        }

    def write(self, **kwargs) -> None:
        """Method to write the output to a networks.txt file"""

        data = kwargs["input_data"]

        # pulling out the list of strings for all the phecode
        # counts
        all_networks = data["output"]["all_networks"]
        filtered_networks = data["output"]["greater_2inds_networks"]

        gene_name = data["gene"]
        # forming the correct output path based on the gene name
        file_output = data["path"]

        pathlib.Path(file_output).mkdir(parents=True, exist_ok=True)

        all_networks_output = os.path.join(
            file_output,
            "".join([gene_name, "_all_networks_phecode_counts.txt"]),
        )

        filtered_networks_output = os.path.join(
            file_output,
            "".join([gene_name, "_greater_than_2inds_phecode_counts.txt"]),
        )

        logger.info(
            f"Information about the phecode counts written to a file at: {', '.join([all_networks_output, filtered_networks_output])}"
        )
        # Opening the file and writing the phecode counts for
        # all the networks to the file
        with open(
            all_networks_output,
            "w+",
            encoding="utf-8",
        ) as output_file:
            # writing a header string
            output_file.write(
                f"phecode\tnetworks_passed_bonf\tnetworks_passed_10e-5\tphecode_description\n"
            )
            # iterating over each network and writing the
            # values to file
            for count_str in tqdm(all_networks, desc="Networks written to file: "):

                # if debug mode is choosen then it will write
                # the output string to a file/console
                logger.debug(count_str)

                output_file.write(count_str)
        # Opening the file and writing the phecode counts for
        # all the networks to the file
        with open(
            filtered_networks_output,
            "w+",
            encoding="utf-8",
        ) as output_file:
            # writing a header string
            output_file.write(
                f"phecode\tnetworks_passed_bonf\tnetworks_passed_10e-5\tphecode_description\n"
            )
            # iterating over each network and writing the
            # values to file
            for count_str in tqdm(filtered_networks, desc="Networks written to file: "):

                # if debug mode is choosen then it will write
                # the output string to a file/console
                logger.debug(count_str)

                output_file.write(count_str)


def initialize() -> None:
    factory_register("phecode_counter", PhecodeCounter)
