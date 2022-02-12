from dataclasses import dataclass
from typing import Protocol, Dict, List
import os
import pandas as pd


class Write_Object(Protocol):
    """Interface indicating that the writer object needs to have the methods write to file"""

    def _form_phenotype_header(self) -> str:
        ...

    def write(self, **kwargs) -> None:
        ...


@dataclass
class Writer:
    """This will be the overall writer class that is responsible for write information to a file"""

    ibd_program: str

    def set_writer(self, writer: Write_Object) -> None:
        self.writer: Write_Object = writer

    def write_to_file(self, output: str, networks_dict: Dict) -> None:
        self.writer.write(
            output=output, network_info=networks_dict, program=self.ibd_program
        )


@dataclass
class Pair_Writer:
    gene_name: str
    chromosome: str
    carriers_matrix: pd.DataFrame

    def check_carriers(self, pair_id: str) -> int:
        """Function to check if pair 1 or pair 2 is a carrier. If not then it returns a 0, if it is then it returns a 1

        pair_id : str
            id of the pair of interest"""

        return int(pair_id in self.carriers_matrix)

    def _form_phenotype_header(self) -> str:
        """Method that will form the phenotype section of the header string"""
        # pulling out all of the phenotype names from the carriers matrix
        columns: List[str] = list(self.carriers_matrix.columns)
        # Appending the word Pair_1/Pair_2 to the column label and then joining them
        # into a string
        pair_1_section: str = "\t".join(["Pair_1_" + label for label in columns])
        pair_2_section: str = "\t".join(["Pair_2_" + label for label in columns])

        return "\t".join([pair_1_section, pair_2_section])

    def _determine_carrier_status(self, pair_id: str) -> str:
        """Method that will determine if the pair id is a carrier for any of the phenotypes listed

        Parameters

        pair_id : str
            IID for each individual

        Returns

        str
            returns a string that has tab separated 0's or 1's whether the individual is
            a carrier or not
        """
        pair_df: pd.DataFrame = self.carriers_matrix[
            self.carriers_matrix["grids"] == pair_id
        ]

        if not pair_id.empty:
            carrier_status: List[int] = pair_df[pair_df.columns[1:]].values.tolist()

            return "\t".join([str(status) for status in carrier_status])

        else:
            return "\t".join(["N/A"] * len(self.carriers_matrix.columns[1:]))

    def write(self, **kwargs) -> None:
        """Method to write the output to an allpairs.txt file"""
        # get the necessary information from the kwargs
        output_dir: str = kwargs["output"]
        networks_info: Dict = kwargs["network_info"]
        ibd_program: str = kwargs["program"]

        # opening the file and then writting the information from each
        # pair to that file. A file will be created for each gene
        with open(
            os.path.join(output_dir, "IBD_", self.gene_name, "_allpairs.txt"), "w"
        ) as output_file:

            output_file.write(
                f"program\tnetwork_id\tpair_1\tpair_2\tchromosome\tgene_name\tphase_1\tphase_2\t{self._form_phenotype_header()}\tstart\tend\tlength\n"
            )
            for network_id, info in networks_info.items():

                pairs = info["pairs"]

                for pair in pairs:

                    output_file.write(
                        pair.format_info(
                            self.gene_name,
                            self.check_carriers(pair.pair_1),
                            self.check_carriers(pair.pair_2),
                            network_id,
                            ibd_program,
                        )
                    )


@dataclass
class Network_Writer:
    gene_name: str
    chromosome: str
    carriers_df: str

    def _form_phenotype_header(self) -> str:
        """Method that will form the phenotype section of the header string"""
        # pulling out all of the phenotype names from the carriers matrix
        columns: List[str] = list(self.carriers_matrix.columns)

        return "\t".join(columns)

    def write(self, **kwargs) -> None:
        """Method to write the output to a networks.txt file"""
        output_dir: str = kwargs["output"]
        networks_info: Dict = kwargs["network_info"]
        ibd_program: str = kwargs["program"]
        pvalue_dict: Dict[int, float] = kwargs["pvalues"]

        with open(
            os.path.join(output_dir, self.gene_name + "_networks.txt"), "w+"
        ) as output_file:
            output_file.write(
                "network_id\tprogram\tgene\tchromosome\tIIDs in Network\tHaplotypes in Network\tIIDs\thaplotypes\n"
            )

            for network_id, info in networks_info.items():

                # string that has the network information such as the network_id, ibd_program, the gene it is for and the chromosome number
                networks: str = (
                    f"{network_id}\t{ibd_program}\t{self.gene_name}\t{self.chromosome}"
                )

                # string that has the number of individuals in the network as well as the the number of haplotypes
                counts: str = f"{len(info['IIDs'])}\t{len(info['haplotypes'])}"
                # string that has the list of GRID IIDs and the haplotype phases
                iids: str = (
                    f"{', '.join(info['IIDs'])}\t{', '.join(info['haplotypes'])}\n"
                )

                output_file.write(f"{networks}\t{counts}\t{iids}")
