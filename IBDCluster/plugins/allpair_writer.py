from dataclasses import dataclass, field
from typing import List, Set, Protocol, Dict, Any, Tuple
import pandas as pd
import log
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
class AllpairWriter:
    """Class object that will handle the string formatting and output for the allpairs.txt file"""

    name: str = "PairWriter plugin"

    @staticmethod
    def _form_header(phenotype_list: List[str]) -> str:
        """Method that will form the phenotype section of the header string"""

        # Appending the word Pair_1/Pair_2 to the column label and then joining them
        # into a string
        logger.debug(
            "Creating a string for all the phecode statuses to the allpairs.txt file"
        )
        header_list = []

        for phecode in phenotype_list:

            header_list.extend(
                [
                    "".join([str(phecode), "_Pair_1_status"]),
                    "".join([str(phecode), "_Pair_2_status"]),
                ]
            )

        return "\t".join(header_list)

    def analyze(self, **kwargs) -> Dict[str, Any]:
        """
        Main function that will create a dictionary of strings that will be used when writting the file
        """
        data: DataHolder = kwargs["data_container"]
        output_path: str = kwargs["output"]

        # This iwll be a list of strings that has the output for each network
        output_dict: Dict[str, Dict[str, Any]] = {}

        for gene_info, network_list in data.networks_dict.items():

            pair_analysis_list: List[str] = []

            for network in tqdm(
                network_list, desc="Networks with pairs written to file: "
            ):
                # upacking the pairs for each network
                pairs = network.pairs

                for pair in pairs:

                    output_str = f"{data.ibd_program}\t{network.network_id}\t{pair.form_id_str()}\t{pair.chromosome}\t{gene_info[0]}\t{pair.form_affected_string()}\t{pair.form_segment_info_str()}\n"

                    pair_analysis_list.append(output_str)

            output_dict[gene_info[0]] = {
                "output": pair_analysis_list,
                "path": os.path.join(output_path, gene_info[0]),
            }

        return output_dict

    def write(self, **kwargs) -> None:
        """Method to write the output to an allpairs.txt file"""
        # get the necessary information from the kwargs
        data = kwargs["input_data"]

        for gene_name in data.keys():
            # full filepath to write the output to
            output_file_name = os.path.join(
                data[gene_name]["path"], "".join(["IBD_", gene_name, "_allpairs.txt"])
            )

            logger.info(f"Writing the allpairs.txt file to: {output_file_name}")
            # opening the file and then writting the information from each
            # pair to that file. A file will be created for each gene
            with open(
                output_file_name,
                "w",
                encoding="utf-8",
            ) as output_file:

                output_file.write(
                    f"program\tnetwork_id\tpair_1\tpair_2\tphase_1\tphase_2\tchromosome\tgene_name\t{self._form_header()}\tstart\tend\tlength\n"
                )
                for pair in tqdm(
                    data[gene_name]["output"], desc="Networks written to file: "
                ):
                    output_file.write(pair)


def initialize() -> None:
    factory_register("allpair_writer", AllpairWriter)
