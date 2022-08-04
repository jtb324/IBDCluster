from typing import Protocol
from dataclasses import dataclass, field
import log
import pandas as pd
import os
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
    pairs: list = field(default_factory=list)
    iids: set[str] = field(default_factory=set)
    haplotypes: set[str] = field(default_factory=set)


@dataclass
class DataHolder(Protocol):
    gene_name: str
    chromosome: int
    networks: Network
    affected_inds: dict[float, list[str]]
    phenotype_table: pd.DataFrame
    phenotype_cols: list[str]
    ibd_program: str
    phenotype_description: dict[str, str] = None
    phenotype_percentages: dict[str, float] = field(default_factory=dict)
    network_pvalues: dict[int, dict[str, float]] = field(default_factory=dict)


@dataclass
class NetworkWriter:
    """Class that is responsible for creating the _networks.txt file from the information provided"""

    name: str = "NetworkWriter plugin"

    @staticmethod
    def _form_header(phenotype_columns) -> str:
        """Method that will form the phenotype section of the header string. Need to
        append the words ind_in_network and pvalue to the phenotype name."""
        # pulling out all of the phenotype names from the carriers matrix

        column_list: list[str] = []

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
    def _check_min_pvalue(phenotype_pvalues: dict[str, float]) -> tuple[str, str]:
        """
        Function that will determine the smallest pvalue and return the corresponding phecode

        phenotype_pvalues : Dict[str, float]

            dictionary where the phecode strings are keys and the values are the
            pvalues as floats

        Returns

        tuple[str, str]

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

    def analyze(self, **kwargs) -> None:
        """main function of the plugin. It needs to determine the pvalue"""

        data: DataHolder = kwargs["data"]
        network: Network = kwargs["network"]
        output_path = kwargs["output"]

        # string that has the network information such as the
        # network_id, ibd_program, the gene it is for and the
        # chromosome number
        networks: str = f"{data.ibd_program}\t{data.gene_name}\t{network.network_id}\t{data.chromosome}"

        # string that has the number of individuals in the
        # network as well as the the number of haplotypes
        counts: str = f"{len(network.iids)}\t{len(network.haplotypes)}"
        # string that has the list of GRID IIDs and the haplotype phases
        iids: str = f"{', '.join(network.iids)}\t{', '.join(network.haplotypes)}"

        output_str = "\t".join([networks, counts, iids, network.pvalues])

        self._write(
            output_str,
            output_path,
            data.ibd_program,
            data.gene_name,
            data.phenotype_cols,
        )
        # # returning an object with the list of strings, the
        # # path, and the gene name
        # return {
        #     "output": networks_analysis_list,
        #     "path": os.path.join(output_path, data.gene_name),
        #     "gene": data.gene_name,
        # }

    def _write(
        self,
        output_str: str,
        output_path: str,
        ibd_program: str,
        gene_name: str,
        phenotype_list: list[str],
    ) -> None:
        """Method to write the output to a networks.txt file
        Parameters
        ----------
        output_str : str
            This is the str created that has all the information
            for each row of the networks.txt file

        output_path : str
            path to write the output to. This is different then
            the output path that the user provides because the
            gene name has been appended to the end of it

        ibd_program : str
            IBD program used to detect segments

        gene_name : str
            name of the gene that is being used as a locus

        phenotype_list : list[str]
            list of phecodes to form each column
        """

        # making sure that the output path is creating
        pathlib.Path(output_path).mkdir(parents=True, exist_ok=True)

        # appending the output file name to the output_path
        output_file_name = os.path.join(
            output_path,
            "".join([ibd_program, "_", gene_name, "_networks.txt"]),
        )

        logger.debug(f"Information written to a networks.txt at: {output_file_name}")
        # Opening the file and writing the head to it and then each network. We are
        # going to open in append mode since we are doing this network by network
        with open(
            output_file_name,
            "a+",
            encoding="utf-8",
        ) as output_file:

            if os.path.getsize(output_file_name) == 0:
                output_file.write(
                    f"program\tgene\tnetwork_id\tchromosome\tIIDs_count\thaplotypes_count\tIIDs\thaplotypes\tmin_pvalue\tmin_pvalue_phecode\tmin_phecode_desc\t{self._form_header(phenotype_list)}\n"
                )

            # if debug mode is choosen then it will write the output string to a file/console
            logger.debug(output_str)

            output_file.write(output_str)


def initialize() -> None:
    factory_register("network_writer", NetworkWriter)
