from dataclasses import dataclass
import pandas as pd


@dataclass
class Pairs:
    """Object that has the data for each pair and contains a method
    to format the output string for the allpairs.txt file using this data"""

    pair_1: str
    phase_1: str
    pair_2: str
    phase_2: str
    chromosome: int
    segment_start: int
    segment_end: int
    length: float
    affected_statuses: pd.Series

    def form_affected_string(self) -> str:
        """Method that will buld a string of affected values from
        the affected_status series attribute"""
        return "\t".join(self.affected_statuses.astype(str))

    def form_id_str(self) -> str:
        """Method that will return the pair ids and phases in a formated string"""
        return f"{self.pair_1}\t{self.pair_2}\t{self.phase_1}\t{self.phase_2}"

    def form_segment_info_str(self) -> str:
        """Method that will return a string with all of the segment info"""
        return f"{self.segment_start}\t{self.segment_end}\t{self.length}\n"
