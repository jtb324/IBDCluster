from dataclasses import dataclass
from typing import Protocol, Dict

class Write_Object(Protocol):
    """Interface indicating that the writer object needs to have the methods write to file"""
    def write(self, **kwargs) -> None:
        pass

@dataclass
class Writer:
    """This will be the overall writer class that is responsible for write information to a file"""

    def set_writer(self, writer: Write_Object) -> None:
        self.writer: Write_Object = writer
    
    def write_to_file(self, output: str, networks_dict: Dict) -> None:
        self.writer.write(
            output=output,
            network_info=networks_dict
            )

class Pair_Writer:
    def write(self, **kwargs) -> None:
        """Method to write the output to an allpairs.txt file"""
        output_dir: str = kwargs["output"]
        networks_info: Dict = kwargs["network_info"]

class Network_Writer:
    def write(self, **kwargs) -> None:
        """Method to write the output to a networks.txt file"""
        output_dir: str = kwargs["output"]
        networks_info: Dict = kwargs["network_info"]

        for network_id, info in networks_info.items():
            pass

