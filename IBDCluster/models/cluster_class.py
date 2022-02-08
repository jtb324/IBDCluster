from dataclasses import dataclass
import pandas as pd
from typing import Dict, List
from models import Pairs
import os

@dataclass
class Cluster:
    ibd_file: str
    network_id: int = 1

    def load_file(self, start: int, end: int, start_indx: int, end_indx:int) -> None:
        """Method filters the ibd file based on location and loads this into memory as a dataframe 
        attribute of the class
        
        Parameters
        
        start : int
            integer that describes the start position of the gene of interest
            
        end : int
            integer that desribes the end position of the gene of interest
        """
        self.ibd_df: pd.DataFrame = pd.DataFrame()

        for chunk in pd.read_csv(self.ibd_file, sep="\t", header=None, chunksize=1000000):
            
            filtered_chunk: pd.DataFrame = chunk[((chunk[start_indx] <= int(start)) & (chunk[end_indx] >= int(start))) | ((chunk[start_indx] >= int(start)) & (chunk[end_indx] <= int(end))) | ((chunk[start_indx] <= int(end)) & (chunk[end_indx] >= int(end)))]

            if not filtered_chunk.empty:
                self.ibd_df = pd.concat(self.ibd_df, filtered_chunk)

        if os.environ.get("verbose", "False") == "True":
            print(f"Found {self.ibd_df.shape[0]} individuals that overlap the given gene region")
            
    def filter_for_haplotype(self, phase_1_indx: int, phase_2_indx: int) -> None:
        """Method that will filter the dataframe for only those haplotypes phases for pair 1 and pair 2 that match
        
        Parameters
        
        phase_1_indx : int
            integer that tells which column in the file has the phasing information for the first individual
            
        phase_2_indx : int
            integer that tells whichh columb in the file has the phasing information for the second individual
        """

        self.ibd_df = self.ibd_df[self.ibd_df[phase_1_indx] == self.ibd_df[phase_2_indx]]

        if os.environ.get("verbose", "False") == "True":
            print(f"Found {self.ibd_df.shape[0]} pairs that have matching haplotypes")

    def find_networks(self) -> Dict[int: Dict]:
        """Method that will go through the dataframe and will identify networks.
        
        Parameters
        
        
        Returns
        Dict[int: Dict]
            returns a dictionary with the following structure:
            {network_id: {pairs: List[Pair_objects], in_network: List[str], haplotypes_list: List[str]}}. Other data information will be added to this file
        """
