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

    def format_info(
        self,
        gene: str,
        network_id: int,
        ibd_program: str,
    ) -> str:
        """Function that will formation the information so that it can easily be written to the allpairs.txt file"""
        return f"{ibd_program}\t{network_id}{self.pair_1}\t{self.pair_2}\t{self.chromosome}\t{gene}\t{self.phase_1}\t{self.phase_2}\t{self.segment_start}\t{self.segment_end}\t{self.length}\n"

