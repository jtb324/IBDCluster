from dataclasses import dataclass
from typing import Protocol, Dict, List
import os
from collections import namedtuple
from tqdm import tqdm

CarriersInfo = namedtuple(
    "Carrier_Comp", ["ind_in_network", "percentage", "IIDs", "pvalue"]
)


class WriteObject(Protocol):
    """Interface indicating that the writer object needs to have the methods write to file"""

    def _form_phenotype_header(self) -> str:
        """method used to form a dynamic column header for different phenotypes. That way the columns can scale with different number of phenotypes that are added"""
        ...

    def write(self, **kwargs) -> None:
        """Method used to write out to a file"""
        ...


@dataclass
class Writer:
    """This will be the overall writer class that is responsible for write information to a file"""

    output: str
    ibd_program: str

    def set_writer(self, writer: WriteObject) -> None:
        """Method that will set whether we are writing for networks or pairs"""
        self.writer: WriteObject = writer

    def write_to_file(self, networks_dict: Dict) -> None:
        """Method that will call the write method of the self.writer class"""
        self.writer.write(
            output=self.output, network_info=networks_dict, program=self.ibd_program
        )


@dataclass
class PairWriter:
    """Class object that will handle the string formatting and output for the allpairs.txt file"""

    gene_name: str
    chromosome: str
    phenotype_list: List[str]

    def _form_phenotype_header(self) -> str:
        """Method that will form the phenotype section of the header string"""

        # Appending the word Pair_1/Pair_2 to the column label and then joining them
        # into a string
        pair_1_section: str = "\t".join(
            ["Pair_1_" + label for label in self.phenotype_list]
        )
        pair_2_section: str = "\t".join(
            ["Pair_2_" + label for label in self.phenotype_list]
        )

        return "\t".join([pair_1_section, pair_2_section])

    def write(self, **kwargs) -> None:
        """Method to write the output to an allpairs.txt file"""
        # get the necessary information from the kwargs
        output_dir: str = kwargs["output"]
        networks_info: Dict = kwargs["network_info"]
        ibd_program: str = kwargs["program"]

        # opening the file and then writting the information from each
        # pair to that file. A file will be created for each gene
        with open(
            os.path.join(
                output_dir, "".join(["IBD_", self.gene_name, "_allpairs.txt"])
            ),
            "w",
            encoding="utf-8",
        ) as output_file:

            output_file.write(
                f"program\tnetwork_id\tpair_1\tpair_2\tphase_1\tphase_2\tchromosome\tgene_name\t{self._form_phenotype_header()}\tstart\tend\tlength\n"
            )
            for network_id, info in tqdm(networks_info.items(), desc="Networks written to file: "):

                pairs = info["pairs"]

                for pair in pairs:

                    output_file.write(
                        f"{ibd_program}\t{network_id}\t{pair.form_id_str()}\t{pair.chromosome}\t{self.gene_name}\t{pair.carrier_str_1}\t{pair.carrier_str_2}\t{pair.form_segment_info_str()}"
                    )


@dataclass
class NetworkWriter:
    """Class that is responsible for creating the _networks.txt file from the information provided"""

    gene_name: str
    chromosome: str
    carriers_columns: List[str]

    def _form_phenotype_header(self) -> str:
        """Method that will form the phenotype section of the header string. Need to
        append the words ind_in_network, percentage_in_net, carrier_IIDs, and pvalue to
        the phenotype name."""

        column_list: List[str] = []
        # pulling out all of the phenotype names from the carriers matrix
        for column in self.carriers_columns:

            column_list.extend(
                [
                    column + ending
                    for ending in [
                        "_ind_in_network",
                        "_percentage_in_net",
                        "_carrier_IIDs",
                        "_pvalue",
                    ]
                ]
            )

        return "\t".join(column_list)

    def _form_analysis_string(self, analysis_dict: Dict[str, CarriersInfo]) -> str:
        """Method that will form a string for each phenotype. We need to use the self.carriers_columns so that the values line up with the columns

        Parameters

        analysis_dict : Dict[str, CarriersInfo]
            Dictionary that has the phenotypes as keys and a namedtuple, Carrier_Comp as the key with info about the carrier count/percentage in network, the IIDs in network and the pvalue

        Returns

        str
            returns a tab separated string that has all of the values
        """

        output_str: str = ""

        # iteratign over every phenotype in the carrier_columns list.
        # Adds all of the CarrierInfo attributes to these string for each phenotype.
        # These are all tab separated
        for phenotype in self.carriers_columns:

            analysis_obj: CarriersInfo = analysis_dict[phenotype]

            if analysis_obj.ind_in_network != 0:
                output_str += f"{analysis_obj.ind_in_network}\t{analysis_obj.percentage}\t{', '.join(analysis_obj.IIDs)}\t{analysis_obj.pvalue}\t"
            else:
                output_str += f"{analysis_obj.ind_in_network}\t{analysis_obj.percentage}\t{'N/A'}\t{analysis_obj.pvalue}\t"

        # strips of the final tab space and then replaces it with a newline
        output_str.rstrip("\t")

        output_str += "\n"

        return output_str

    def write(self, **kwargs) -> None:
        """Method to write the output to a networks.txt file"""
        output_dir: str = kwargs["output"]
        networks_info: Dict = kwargs["network_info"]
        ibd_program: str = kwargs["program"]

        with open(
            os.path.join(
                output_dir, "".join([ibd_program, "_", self.gene_name, "_networks.txt"])
            ),
            "w+",
            encoding="utf-8",
        ) as output_file:
            output_file.write(
                f"network_id\tprogram\tgene\tchromosome\tIIDs_count\thaplotypes_count\tIIDs\thaplotypes\t{self._form_phenotype_header()}\n"
            )

            for network_id, info in tqdm(networks_info.items(), desc="Networks written to file: "):

                # string that has the network information such as the network_id, ibd_program, the gene it is for and the chromosome number
                networks: str = (
                    f"{network_id}\t{ibd_program}\t{self.gene_name}\t{self.chromosome}"
                )

                # string that has the number of individuals in the network as well as the the number of haplotypes
                counts: str = f"{len(info['IIDs'])}\t{len(info['haplotypes'])}"
                # string that has the list of GRID IIDs and the haplotype phases
                iids: str = (
                    f"{', '.join(info['IIDs'])}\t{', '.join(info['haplotypes'])}"
                )

                # Getting the phenotype analysis infofrom the object
                analysis_obj: Dict[str, CarriersInfo] = info["phenotype"]

                analysis_str: str = f"{self._form_analysis_string(analysis_obj)}"

                output_file.write(f"{networks}\t{counts}\t{iids}\t{analysis_str}")
