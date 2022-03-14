from dataclasses import dataclass
from typing import Protocol, Dict, List, Tuple, Optional
import os
from collections import namedtuple
from tqdm import tqdm
import log
from models import Network

logger = log.get_logger(__name__)


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
    writer: Optional[WriteObject] = None

    def set_writer(self, writer: WriteObject) -> None:
        """Method that will set whether we are writing for networks or pairs"""
        self.writer: WriteObject = writer

    def write_to_file(self, networks_dict: Dict) -> None:
        """Method that will call the write method of the self.writer class"""
        logger.info(f"writing the output to the directory: {self.output}")
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
        logger.debug(
            "Creating a string for all the phecode statuses to the allpairs.txt file"
        )
        header_list = []

        for phecode in self.phenotype_list:

            header_list.extend(
                [
                    "".join([str(phecode), "_Pair_1_status"]),
                    "".join([str(phecode), "Pair_2_status"]),
                ]
            )

        return "\t".join(header_list)

    def write(self, **kwargs) -> None:
        """Method to write the output to an allpairs.txt file"""
        # get the necessary information from the kwargs
        output_dir: str = kwargs["output"]
        network_list: List[Network] = kwargs["network_info"]
        ibd_program: str = kwargs["program"]

        # full filepath to write the output to
        output_file_name = os.path.join(
            output_dir, "".join(["IBD_", self.gene_name, "_allpairs.txt"])
        )

        logger.debug(f"Writing the allpairs.txt file to: {output_file_name}")
        # opening the file and then writting the information from each
        # pair to that file. A file will be created for each gene
        with open(
            output_file_name,
            "w",
            encoding="utf-8",
        ) as output_file:

            output_file.write(
                f"program\tnetwork_id\tpair_1\tpair_2\tphase_1\tphase_2\tchromosome\tgene_name\t{self._form_phenotype_header()}\tstart\tend\tlength\n"
            )
            for network in tqdm(network_list, desc="Networks written to file: "):

                pairs = network.pairs

                for pair in pairs:
                    output_str = f"{ibd_program}\t{network.network_id}\t{pair.form_id_str()}\t{pair.chromosome}\t{self.gene_name}\t{pair.form_affected_string()}\t{pair.form_segment_info_str()}"

                    logger.debug(output_str)

                    output_file.write(output_str)


@dataclass
class NetworkWriter:
    """Class that is responsible for creating the _networks.txt file from the information provided"""

    gene_name: str
    chromosome: str
    carriers_columns: List[str]

    def _form_phenotype_header(self) -> str:
        """Method that will form the phenotype section of the header string. Need to
        append the words ind_in_network and pvalue to the phenotype name."""
        # pulling out all of the phenotype names from the carriers matrix

        column_list: List[str] = []

        for column in self.carriers_columns:

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

    def _find_min_phecode(self, analysis_dict) -> Tuple[float, str]:
        """Function that will identify the lowest phecode and will write that and the phecode name to a column

        Parameters

        analysis_dict : Dict[str, CarriersInfo]
            Dictionary that has the phenotypes as keys and a namedtuple, Carrier_Comp as the key with info about the carrier count/percentage in network, the IIDs in network and the pvalue

        Returns

        Tuple[float, str]
            returns a tuple where the first element is the pvalue and
            the second element is the phecode
        """
        pvalue_dict: Dict[str, float] = {}

        for phenotype in self.carriers_columns:

            carrier_obj: CarriersInfo = analysis_dict[phenotype]

            pvalue_dict[phenotype] = carrier_obj.pvalue
        # retu
        pvalue_tuple = min(zip(pvalue_dict.values(), pvalue_dict.keys()))

        return str(pvalue_tuple[0]), pvalue_tuple[1]

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

            carrier_obj: CarriersInfo = analysis_dict[phenotype]

            output_str += f"{carrier_obj.ind_in_network}\t{carrier_obj.pvalue}\t"

        # strips of the final tab space and then replaces it with a newline
        output_str.rstrip("\t")

        output_str += "\n"

        return output_str

    def write(self, **kwargs) -> None:
        """Method to write the output to a networks.txt file"""
        network_list: List[Network] = kwargs["network_info"]
        ibd_program: str = kwargs["program"]

        output_file_name = os.path.join(
            kwargs["output"],
            "".join([ibd_program, "_", self.gene_name, "_networks.txt"]),
        )

        logger.debug(f"Information written to a networks.txt at: {output_file_name}")

        with open(
            output_file_name,
            "w+",
            encoding="utf-8",
        ) as output_file:
            output_file.write(
                f"program\tgene\tnetwork_id\tchromosome\tIIDs_count\thaplotypes_count\tIIDs\thaplotypes\tmin_pvalue\tmin_pvalue_phecode\t{self._form_phenotype_header()}\n"
            )

            for network in tqdm(network_list, desc="Networks written to file: "):
                # string that has the network information such as the
                # network_id, ibd_program, the gene it is for and the
                # chromosome number
                networks: str = f"{ibd_program}\t{self.gene_name}\t{network.network_id}\t{self.chromosome}"

                # string that has the number of individuals in the
                # network as well as the the number of haplotypes
                counts: str = f"{len(network.iids)}\t{len(network.haplotypes)}"
                # string that has the list of GRID IIDs and the haplotype phases
                iids: str = (
                    f"{', '.join(network.iids)}\t{', '.join(network.haplotypes)}"
                )

                # Getting the phenotype analysis infofrom the object
                analysis_obj: Dict[str, CarriersInfo] = network.phenotype

                # Creating a string that has the information about the
                # minimum pvalue and corresponding phecode for the network
                min_pvalue_str: str = "\t".join(self._find_min_phecode(analysis_obj))

                # creating a string that has the carrier count and the pvalues for each
                # phecode
                analysis_str: str = f"{self._form_analysis_string(analysis_obj)}"

                # if debug mode is choosen then it will write the output string to a file/console
                logger.debug(
                    f"{networks}\t{counts}\t{iids}\t{min_pvalue_str}\t{analysis_str}"
                )

                output_file.write(
                    f"{networks}\t{counts}\t{iids}\t{min_pvalue_str}\t{analysis_str}"
                )
