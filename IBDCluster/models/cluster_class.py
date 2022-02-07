from dataclasses import dataclass
from email import header
import pandas as pd

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