from dataclasses import dataclass
from typing import List

@dataclass
class Pairs:
    pair_1: str
    phase_1: str
    pair_2: str 
    phase_2: str
    chromosome: int
    segment_start: int
    segment_end: int
    length: float  

    def format_info(self) -> str:
        """Function that will formation the information so that it can easily be written to a file"""
        return f"{self.pair_1}\t{self.pair_2}\t{self.phase_1}\t{self.phase_2}\t{self.chromosome}\t{self.segment_start}\t{self.segment_end}\t{self.length}"
    
    def check_confirmed_status(self, carrier_list: List[str]) -> int:
        """Function that will return a 1 if the pair_2 is a carrier and pair_1 is a carrier. Otherwise it returns a 0
        
        Parameters
        
        carrier_list : List[str]
            list that has the iid of all individuals who are identified as carriers
        """
        pass
