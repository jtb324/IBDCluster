from dataclasses import dataclass
from typing import List

from drive.factory import factory_register
from drive.log import CustomLogger
from drive.models import Data_Interface, Network_Interface

logger = CustomLogger.get_logger(__name__)


@dataclass
class NetworkWriter:
    """Class that is responsible for creating the *_networks.
    txt file from the information provided"""

    name: str = "NetworkWriter plugin"

    @staticmethod
    def _form_header(phenotypes: List[str]) -> str:
        """Method that will form the header line for the output file

        Parameters
        ----------
        phenotypes : List[str]
            list of all the phenotypes from the Data.carriers attribute

        Returns
        -------
        str
            returns the header string for the file
        """

        # making a string for the initial first few columns
        header_str = "clstID\tn.total\tn.haplotype\ttrue.positive.n\ttrue.positive\tfalst.postive\tIDs\tID.haplotype"  # noqa: E501

        if not phenotypes:
            return header_str + "\n"
        else:
            # We need to add columns for the min pvalue descriptions
            header_str += "\tmin_pvalue\tmin_phenotype\tmin_phenotype_description"
            # for each phenotype we are going to create 4 columns for the number
            # of cases in the network, The case ids in the network the number of
            # excluded individuals in the network, and the pvalue for the phenotype
            for column in phenotypes:
                header_str += f"\t{column + '_case_count_in_network'}\t{column + 'cases_in_network'}\t{column + '_excluded_in_network'}\t{column + '_pvalue'}"  # noqa: E501

            return header_str + "\n"

    @staticmethod
    def _create_network_info_str(
        network: Network_Interface, phenotypes: List[str]
    ) -> str:
        """create a string that has all the information from the network
        such as cluster id, member count, pvalues, etc...

        Parameters
        ----------
        network : Network_Interface
            network object that hsas all the per network
            information from the clusters such as cluster
            ids, networks, pvalues, etc...

        phenotypes : List[str]
            list of the phenotypes provided the program.

        Returns
        -------
        str
            returns a string formatted for the output file
        """
        # fill in the initial few columns of the output string
        output_str = f"{network.clst_id}\t{len(network.members)}\t{len(network.haplotypes)}\t{network.true_positive_count}\t{network.true_positive_percent:.4f}\t{network.false_negative_count}\t{','.join(network.members)}\t{','.join(network.haplotypes)}"  # noqa: E501

        if not phenotypes:
            return output_str + "\n"
        else:
            output_str += f"\t{network.min_pvalue_str}"

            for phenotype in phenotypes:
                output_str += f"\t{network.pvalues[phenotype]}"

            return output_str + "\n"

    def analyze(self, **kwargs) -> None:
        """main function of the plugin that will create the
        output path and then use helper functions to write
        information to a file"""

        data: Data_Interface = kwargs["data"]

        # creating the full output path for the output file
        network_file_output = data.output_path.parent / (
            data.output_path.name + ".drive_networks.txt"
        )  # noqa: E501

        logger.debug(
            f"The output in the network_writer plugin is being written to: {network_file_output}"  # noqa: E501
        )

        # we are going to pull out the phenotypes into a list so that we
        # are guarenteed to maintain order as we are creating the rows
        phenotypes = list(data.carriers.keys())

        with open(network_file_output, "w", encoding="utf-8") as networks_output:
            header_str = NetworkWriter._form_header(phenotypes)
            # iterate over each network and pull out the appropriate
            # information into strings
            _ = networks_output.write(header_str)

            for network in data.networks:
                network_info_str = NetworkWriter._create_network_info_str(
                    network, phenotypes
                )

                networks_output.write(network_info_str)


def initialize() -> None:
    factory_register("network_writer", NetworkWriter)
