from dataclasses import dataclass, field
from hashlib import new
import pandas as pd
from typing import Dict, List, Set, Protocol, Tuple
from .pairs import Pairs
import os

# class protocol that makes sure that the indices object has these methods
class File_Info(Protocol):
    def gather_files(self) -> None:
        ...
    def find_files(self, chr_num: str) -> str:
        ...

@dataclass
class Cluster:
    gene_name: str
    ibd_file: str
    pairs_found: List[Tuple] = field(default_factory=list) # These tuples will have be formatted (pair1, pair2) and (pair2, pair1) so that we can match based on any order
    individuals_in_network: Set[str] = field(default_factory=set) # This attribute will keep track of all individuals that are in any network
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
                self.ibd_df = pd.concat([self.ibd_df, filtered_chunk])

        if os.environ.get("verbose", "False") == "True":
            individ: List[str] = list(set(self.ibd_df[0].values.tolist() + self.ibd_df[2].values.tolist()))

            print(f"Found {len(individ)} individuals that overlap the given gene region")
            
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
    def _determine_iids_in_network(dataframe: pd.DataFrame, indices: File_Info) -> Set[str]:
        """Staticmethod that will find all the unique individuals in the dataframe and return that list
        
        Parameters
        
        dataframe : pd.DataFrame
            dataframe that has all the ibd sharing for pairs of individuals
            
        indices : Union[Hapibd_Info | Ilash_Info]
            class object that has the indices for either hapibd files or ilash files
        
        Returns

        Set[str]
            returns a list of unique iids that are in the dataframe
        """
        return set(dataframe[indices.id1_indx].values.tolist()+dataframe[indices.id2_indx].values.tolist())

    @staticmethod
    def _gather_haplotypes(dataframe: pd.DataFrame, indices: File_Info) -> Set[str]:
        """Static Method that will determine the unique haplotypes within the pairs
        
        Parameters
        
        dataframe : pd.DataFrame
            dataframe that has all the ibd sharing for pairs of individuals
            
        indices : Union[Hapibd_Info | Ilash_Info]
            class object that has the indices for either hapibd files or ilash files
        
        Returns

        Set[str]
            returns a set of unique haplotypes that are in the dataframe
        """
        pair_1_haplotypes: List[str] = list(set(list(dataframe[indices.id1_indx] + "." + dataframe[indices.id1_phase_indx].astype('str'))))

        pair_2_haplotypes: List[str] = list(set(list(dataframe[indices.id2_indx] + "." + dataframe[indices.id2_phase_indx].astype('str'))))

        return set(pair_1_haplotypes+pair_2_haplotypes)
    
    def _add_pair_tuples(self, id1: List[str], id2: List[str]) -> None:
        """Method that will add each of the pairs found to the class attribute 'pairs_found'. This atrribute provides a list of tuples that I can later check against to make sure that we are not repeating pair objects
        
        id1 : List[str]
            List of all the pair1 iids
            
        id2 : List[str]
            List of all the pair2 iids
        """
        self.pairs_found.extend(list(zip(id1, id2))+list(zip(id2, id1)))

    def _find_secondary_connections(self, new_individuals: List[str], exclusion: Set[str], indices: File_Info, network_dict: Dict) -> None:
        """Function that will find the secondary connections within the graph
        
        Parameters
        
        new_individuals : List[str]
            List of iids that are not the exclusion iid. The list will be all 
            individuals who share with the iid
        
        exclusion : Set[str]
            This is a set of iids that we want to keep out of the secondary cluster because we 
            seeded from this iid so we have these connections
        
        indices : File_Info
            object that has all the indices values for the correct column in the hapibd file
        
        network_dict : Dict[int: Dict]
            returns a dictionary with the following structure:
            {network_id: {pairs: List[Pair_objects], in_network: List[str], haplotypes_list: List[str]}}. Other data information will be added to this file
        """
        
        # getting all the connections based on the new_individuals grids
        secondary_connections: pd.DataFrame = self.ibd_df[(self.ibd_df[indices.id1_indx].isin(new_individuals)) | (self.ibd_df[indices.id2_indx].isin(new_individuals))]
        # print(f"secondary carriers: {secondary_connections.shape[0]}")
        # excluding the exclusion individual from the secondary connections because we 
        # already have these connections
        excluded_df: pd.DataFrame = secondary_connections[~(secondary_connections[indices.id1_indx].isin(exclusion)) & ~(secondary_connections[indices.id2_indx].isin(exclusion))]
        # print(f"excluded carriers: {excluded_df.shape[0]}")
        if not excluded_df.empty:
            # getting a list of all new individuals 
            new_connections: List[str] = list(set(excluded_df[indices.id1_indx].values.tolist() + excluded_df[indices.id2_indx].values.tolist())) 
            # print([iid for iid in new_connections if iid in exclusion])
            # This will extend all the list of pairs for the network to include the new 
            # pairs 
            network_dict[self.network_id]["pairs"].extend(list(excluded_df.apply(lambda row: self._determine_pairs(row, indices), axis=1)))
            
            # need to create a list of people in the network
            network_dict[self.network_id]["in_network"].update(self._determine_iids_in_network(excluded_df, indices))

            # now need to get a list of haplotypes 
            network_dict[self.network_id]["haplotypes"].update(self._gather_haplotypes(excluded_df, indices))
            # This adds all the grids that are in the network to the class attribute that stores all the individuals that are found in any network
            self.individuals_in_network.update(network_dict[self.network_id]["in_network"])
            # updating the exclusion set so that the individuals that we seeded of of this iteration are added to it
            exclusion.update(new_individuals)
            # now need to recursively find more connections
            self._find_secondary_connections(
                [iid for iid in new_connections if iid not in new_individuals],
                exclusion,
                indices,
                network_dict
                )

    def find_networks(self, indices: File_Info) -> Dict[int, Dict]:
        """Method that will go through the dataframe and will identify networks.
        
        Parameters
        
        indices : File_Info
            object that has all the indices values for the correct column in the hapibd file

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

            # if this iid has already been associated with a network then we need to skip it. If not then we can get the network connected to it
            if ind not in self.individuals_in_network:
                
                # creating a key in the dictionary for the network id
                network_info[self.network_id] = {}

                filtered_df: pd.DataFrame = self.ibd_df[(self.ibd_df[indices.id1_indx] == ind) | (self.ibd_df[indices.id2_indx] == ind)]
                # print(f"initial cluster: {filtered_df.shape[0]}")
                # forming a list of Pair objects that will be added to the network_info_dict in a pairs key
                network_info[self.network_id]["pairs"] = list(filtered_df.apply(lambda row: self._determine_pairs(row, indices), axis=1))

                # need to create a list of people in the network
                network_info[self.network_id]["in_network"] = self._determine_iids_in_network(filtered_df, indices)

                # now need to get a list of haplotypes 
                network_info[self.network_id]["haplotypes"] = self._gather_haplotypes(filtered_df, indices)
                
                # This adds all the grids that are in the network to the class attribute that stores all the individuals that are found in any network
                self.individuals_in_network.update(network_info[self.network_id]["in_network"])
                
                # Updating the pairs_found attribute for the pairs that are connected
                # self._add_pair_tuples(filtered_df[indices.id1_indx].values.tolist() , filtered_df[indices.id2_indx].values.tolist())
                # finding the secondary connections. These are connections to all 
                # individuals in the self.individuals_in_network set that are not the 
                # seeding individual
                self._find_secondary_connections(
                    [iid for iid in network_info[self.network_id]["in_network"] if iid != ind], 
                    set([ind]), 
                    indices,
                    network_info
                    )
                
                print(f"number of iids in network {self.network_id}: {len(network_info[self.network_id]['in_network'])}")
                print(len(self.individuals_in_network))
                self.network_id += 1


        # we can reset the network_id for the next cluster
        self.network_id = 1
        print(network_info)
        return network_info
        



