from dataclasses import dataclass
from typing import Protocol, Dict, List
import os


class Write_Object(Protocol):
    """Interface indicating that the writer object needs to have the methods write to file"""

    def write(self, **kwargs) -> None:
        pass


# class protocol that makes sure that the indices object has these methods
class File_Info(Protocol):
    def gather_files(self) -> None:
        ...

    def find_files(self, chr_num: str) -> str:
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
    carriers: List[str]

    def check_carriers(self, pair_id: str) -> int:
        """Function to check if pair 1 or pair 2 is a carrier. If not then it returns a 0, if it is then it returns a 1
        
        pair_id : str
            id of the pair of interest"""

        return int(pair_id in self.carriers)

    def write(self, **kwargs) -> None:
        """Method to write the output to an allpairs.txt file"""
        # get the necessary information from the kwargs
        output_dir: str = kwargs["output"]
        networks_info: Dict = kwargs["network_info"]
        ibd_program: str = kwargs["program"]
        
        # opening the file and then writting the information from each 
        # pair to that file. A file will be created for each gene
        with open(os.path.join(output_dir, "IBD_", self.gene_name, "_allpairs.txt"), "w") as output_file:

            output_file.write(
                "program\tnetwork_id\tpair_1\tpair_2\tchromosome\tgene_name\tphase_1\tphase_2\tstart\tend\tlength\n"
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

    def write(self, **kwargs) -> None:
        """Method to write the output to a networks.txt file"""
        output_dir: str = kwargs["output"]
        networks_info: Dict = kwargs["network_info"]
        ibd_program: str = kwargs["program"]

        with open(
            os.path.join(output_dir, self.gene_name + "_networks.txt"), "w+"
        ) as output_file:
            output_file.write("network_id\tprogram\tgene\tchromosome\tIIDs\thaplotypes\n")
            for network_id, info in networks_info.items():
                output_file.write(f"{network_id}\t{ibd_program}\t{self.gene_name}\t{self.chromosome}\t{', '.join(info['IIDs'])}\t{', '.join(info['haplotypes'])}\n")
                


