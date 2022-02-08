from dataclasses import dataclass, field
import pandas as pd
from typing import Dict, List, Set
from models import Pairs, File_Info
import os

@dataclass
class Cluster:
    gene_name: str
    ibd_file: str
    individuals_in_network: Set[str] = field(default_factory=set)
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

    @staticmethod
    def _find_all_grids(dataframe: pd.DataFrame, ind_1_indx: int, ind_2_indx: int) -> List[str]:
        """Method that will take the dataframe that is filtered for the location and the haplotype and return 
        a list of all unique individuals in that dataframe
        
        Parameters
        
        dataframe : pd.DataFrame
            dataframe that has the pairs of individuals who share an ibd segment. This dataframe is filtered for a 
            specific loci and matching haplotype phase.

        ind_1_indx : int
            integer for the column that has the individual 1 IIDs

        ind_2_indx : int
            integer for the column that has the individual 2 IIDs
            
        Returns
        
        List[str]
            returns a list of unique individuals in the dataframe
        """
        return list(set(dataframe[ind_1_indx].values.tolist()+dataframe[ind_2_indx].values.tolist()))

    @staticmethod
    def _determine_pairs(ibd_row: pd.Series, indices: File_Info) -> Pairs:
        """Method that will take each row of the dataframe and convert it into a pair object"""

        return Pairs(
            ibd_row[indices.id1_indx],
            ibd_row[indices.id1_phase_indx],
            ibd_row[indices.id2_indx],
            ibd_row[indices.id2_phase_indx],
            ibd_row[indices.chr_indx],
            ibd_row[indices.str_indx],
            ibd_row[indices.end_indx],
            ibd_row[indices.cM_indx]
            )
    
    @staticmethod
    def _determine_iids_in_network(dataframe: pd.DataFrame, indices: File_Info) -> List[str]:
        """Staticmethod that will find all the unique individuals in the dataframe and return that list
        
        Parameters
        
        dataframe : pd.DataFrame
            dataframe that has all the ibd sharing for pairs of individuals
            
        indices : Union[Hapibd_Info | Ilash_Info]
            class object that has the indices for either hapibd files or ilash files
        
        Returns

        List[str]
            returns a list of unique iids that are in the dataframe
        """
        return list(set(dataframe[indices.id1_indx].values.tolist()+dataframe[indices.id2_indx].values.tolist()))

    @staticmethod
    def _gather_haplotypes(dataframe: pd.DataFrame, indices: File_Info) -> List[str]:
        """Static Method that will determine the unique haplotypes within the pairs
        
        Parameters
        
        dataframe : pd.DataFrame
            dataframe that has all the ibd sharing for pairs of individuals
            
        indices : Union[Hapibd_Info | Ilash_Info]
            class object that has the indices for either hapibd files or ilash files
        
        Returns

        List[str]
            returns a list of unique haplotypes that are in the dataframe
        """
        pair_1_haplotypes: List[str] = list(set(list(dataframe[indices.id1_indx] + "." + dataframe[indices.id1_phase_indx])))

        pair_2_haplotypes: List[str] = list(set(list(dataframe[indices.id2_indx] + "." + dataframe[indices.id2_phase_indx])))

        return list(set(pair_1_haplotypes+pair_2_haplotypes))

    def find_networks(self, indices: File_Info) -> Dict[int: Dict]:
        """Method that will go through the dataframe and will identify networks.
        
        Parameters
        
        ind_1_indx : int
            integer for the column that has the individual 1 IIDs

        ind_2_indx : int
            integer for the column that has the individual 2 IIDs

        Returns

        Dict[int: Dict]
            returns a dictionary with the following structure:
            {network_id: {pairs: List[Pair_objects], in_network: List[str], haplotypes_list: List[str]}}. Other data information will be added to this file
        """
        # First need to get a list of all the unique individuals in the network
        iid_list: List[str] = self._find_all_grids(self.ibd_df, indices.id1_indx, indices.id2_indx)

        # creating a dictionary that will return information of interest
        network_info: Dict[int: Dict] = {}

        # iterate over each iid in the original dataframe
        for ind in iid_list:

            # creating a key in the dictionary for the network id
            network_info[self.network_id] = {}

            # if this iid has already been associated with a network then we need to skip it. If not then we can get the network connected to it
            if ind not in self.individuals_in_network:
                
                filtered_df: pd.DataFrame = self.ibd_df[(self.ibd_df[indices.id1_indx] == ind) | (self.ibd_df[indices.id2_indx] == ind)]

                # forming a list of Pair objects that will be added to the network_info_dict in a pairs key
                network_info[self.network_id]["pairs"] = list(filtered_df.apply(lambda row: self._determine_pairs(row, indices), axis=1))

                # need to create a list of people in the network
                network_info[self.network_id]["in_network"] = self._determine_iids_in_network(filtered_df, indices)

                # now need to get a list of haplotypes 




