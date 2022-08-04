"""
Plugin that iterates over the pairs for each network and will write that information to an allpairs.txt file 
"""
import os
import pathlib
from dataclasses import dataclass, field
from typing import Protocol

import log
from models import DataHolder
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
    pairs: list = field(default_factory=list)
    iids: set[str] = field(default_factory=set)
    haplotypes: set[str] = field(default_factory=set)


@dataclass
class AllpairWriter:
    """Class object that will handle the string formatting and output for the allpairs.txt file"""

    name: str = "PairWriter plugin"

    def analyze(self, **kwargs) -> None:
        """
        Main function that will create a dictionary of strings that will be used when writting the file
        """
        data: DataHolder = kwargs["data"]
        network: Network = kwargs["network"]
        output_path: str = kwargs["output"]

        # upacking the pairs for each network
        pairs = network.pairs
        logger.debug(f"Network id: {network.network_id}, number of pairs: {len(pairs)}")
        for pair in pairs:

            # creating a string that has all the information
            output_str = f"{data.ibd_program}\t{network.network_id}\t{pair.form_id_str()}\t{pair.chromosome}\t{data.gene_name}\t{pair.form_segment_info_str()}\n"

            self._write(
                output_str,
                output_path,
                data.gene_name,
            )

    def _write(
        self,
        output_str: str,
        output_path: str,
        gene_name: str,
    ) -> None:
        """Method to write the output to an allpairs.txt file
        Parameters
        ----------
        output_str : str
            This is the str created that has all the information
            for each row of the allpairs.txt file

        output_path : str
            path to write the output to. This is different then
            the output path that the user provides because the
            gene name has been appended to the end of it

        ibd_program : str
            IBD program used to detect segments

        gene_name : str
            name of the gene that is being used as a locus
        """
        # making sure that the output path exists
        pathlib.Path(output_path).mkdir(parents=True, exist_ok=True)

        # full filepath to write the output to
        output_file_name = os.path.join(
            output_path, "".join(["IBD_", gene_name, "_allpairs.txt"])
        )

        # opening the file and then writting the information from each
        # pair to that file. A file will be created for each gene
        with open(
            output_file_name,
            "a+",
            encoding="utf-8",
        ) as output_file:

            if os.path.getsize(output_file_name) == 0:
                output_file.write(
                    "program\tnetwork_id\tpair_1\tpair_2\tphase_1\tphase_2\tchromosome\tgene_name\tstart\tend\tlength\n"
                )

            logger.debug(f"allpair_output_str: {output_str}")
            output_file.write(output_str)


def initialize() -> None:
    factory_register("allpair_writer", AllpairWriter)
