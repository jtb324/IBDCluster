from dataclasses import dataclass, field
import pandas as pd
from tqdm import tqdm
from typing import Dict, List, Set, Protocol, Optional
import log
import os
from numpy import where

from .pairs import Pairs


logger = log.get_logger(__name__)

# class protocol that makes sure that the indices object has these methods
class FileInfo(Protocol):
    """Protocol that enforces these two methods for the
    FileInfo object"""

    id1_indx: int
    ind1_with_phase: int
    id2_indx: int
    ind2_with_phase: int
    chr_indx: int
    str_indx: int
    end_indx: int

    def set_program_indices(self, program_name: str) -> None:
        """Method that will set the proper ibd_program indices"""
        ...

    @staticmethod
    def find_file(chr_num: str, file_list: List[str]) -> str:
        """Method that will find the proper ibd files"""
        ...


@dataclass
class Network:
    """This class is going to be responsible for the clustering of each network"""

    gene_name: str
    gene_chr: str
    network_id: int
    pairs: List[Pairs] = field(default_factory=list)
    iids: Set[str] = field(default_factory=set)
    haplotypes: Set[str] = field(default_factory=set)

    def filter_for_seed(
        self,
        ibd_df: pd.DataFrame,
        ind_seed: List[str],
        indices: FileInfo,
        exclusion: Set[str] = None,
    ) -> pd.DataFrame:
        """Method to filter the ibd_df for the first individual. This gets the first level new_connections"""
        filtered_df: pd.DataFrame = ibd_df[
            (ibd_df[indices.ind1_with_phase].isin(ind_seed))
            | (ibd_df[indices.ind2_with_phase].isin(ind_seed))
        ]

        # in the secondary connections
        if exclusion:
            logger.debug(f"excluding {', '.join(exclusion)}")

            filtered_df: pd.DataFrame = filtered_df[
                ~(filtered_df[indices.ind1_with_phase].isin(exclusion))
                & ~(filtered_df[indices.ind2_with_phase].isin(exclusion))
            ]
        logger.debug(
            f"found {filtered_df.shape[0]} pairs that contain the individuals: {', '.join(ind_seed)}"
        )

        return filtered_df

    @staticmethod
    def _determine_pairs(ibd_row: pd.Series, indices: FileInfo) -> Pairs:
        """Method that will take each row of the dataframe and convert it into a pair object"""

        carrier_cols = [col for col in ibd_row.keys() if "status" in str(col)]

        affected_values: pd.Series = ibd_row[carrier_cols]

        return Pairs(
            ibd_row[indices.id1_indx],
            ibd_row[indices.ind1_with_phase],
            ibd_row[indices.id2_indx],
            ibd_row[indices.ind2_with_phase],
            ibd_row[indices.chr_indx],
            ibd_row[indices.str_indx],
            ibd_row[indices.end_indx],
            ibd_row[indices.program_indices.cM_indx],
            affected_values,
        )

    @staticmethod
    def gather_grids(
        dataframe: pd.DataFrame, pair_1_indx: int, pair_2_indx: int
    ) -> Set[str]:
        """Staticmethod that will find all the unique values in two columns that the user passed in

        Parameters

        dataframe : pd.DataFrame
            dataframe that has all the ibd sharing for pairs of individuals

        indices : Union[Hapibd_Info | Ilash_Info]
            class object that has the indices for either hapibd files or ilash files

        Returns

        Set[str]
            returns a list of unique iids that are in the dataframe
        """
        return set(
            dataframe[pair_1_indx].values.tolist()
            + dataframe[pair_2_indx].values.tolist()
        )

    def update(self, ibd_df: pd.DataFrame, indices: FileInfo) -> None:
        """Method to update the pair attribute, the iids, and the haplotypes"""

        logger.debug(
            f"updating the pairs, iids, and haplotype attributes for the network {self.network_id}"
        )

        self.pairs.extend(
            list(ibd_df.apply(lambda row: self._determine_pairs(row, indices), axis=1))
        )

        # updating the iids attribute with what is in the ibd_df
        self.iids.update(self.gather_grids(ibd_df, indices.id1_indx, indices.id2_indx))

        # now need to get a list of haplotypes
        self.haplotypes.update(
            ibd_df[indices.ind1_with_phase].values.tolist()
            + ibd_df[indices.ind2_with_phase].values.tolist()
        )


@dataclass
class Cluster:
    """Class object that will handle preparing the data to be clustered"""

    ibd_file: str
    ibd_program: str
    count: int = 0  # this is a counter that is used in testing to speed the process
    ibd_df: Optional[pd.DataFrame] = field(default_factory=pd.DataFrame)

    def load_file(self, start: int, end: int, start_indx: int, end_indx: int) -> None:
        """Method filters the ibd file based on location and loads this into memory as a dataframe
        attribute of the class

        Parameters

        start : int
            integer that describes the start position of the gene of interest

        end : int
            integer that desribes the end position of the gene of interest
        """
        logger.debug(
            f"Gathering shared ibd segments that overlap the gene region from {start} to {end} using the file {self.ibd_file}"
        )

        if self.ibd_program == "hapibd":
            name_cols = [0, 1, 2, 3, 4, 5, 6, 7]
        else:
            name_cols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        for chunk in pd.read_csv(
            self.ibd_file, sep="\t", header=None, chunksize=1000000, names=name_cols
        ):

            filtered_chunk: pd.DataFrame = chunk[
                ((chunk[start_indx] <= int(start)) & (chunk[end_indx] >= int(start)))
                | ((chunk[start_indx] >= int(start)) & (chunk[end_indx] <= int(end)))
                | ((chunk[start_indx] <= int(end)) & (chunk[end_indx] >= int(end)))
            ]

            if not filtered_chunk.empty:
                self.ibd_df = pd.concat([self.ibd_df, filtered_chunk])

        logger.info(f"identified {self.ibd_df.shape[0]} pairs within the gene region")

    def add_carrier_status(
        self,
        carriers: Dict[float, List[str]],
        pair_1_indx: int,
        pair_2_indx: int,
        phecode_list: List[float],
    ) -> None:
        """Method that will determine the carrier status for each pair

        Parameters

        carrier_list : Dict[float, List[str]]
            dictionary where the keys are the phecodes and the values are list of individuals who are carriers for a specific phecode
        """
        logger.debug(
            f"Adding a carrier status of either 1 or 0 for all {len(carriers)} phecodes"
        )
        status_dict = {}

        for phecode in phecode_list:
            carrier_list = carriers[phecode]

            status_dict["_".join([str(phecode), "pair_1_status"])] = where(
                self.ibd_df[pair_1_indx].isin(carrier_list), 1, 0
            )

            status_dict["_".join([str(phecode), "pair_2_status"])] = where(
                self.ibd_df[pair_2_indx].isin(carrier_list), 1, 0
            )

        self.ibd_df = pd.concat(
            [
                self.ibd_df.reset_index(drop=True),
                pd.DataFrame.from_dict(status_dict).reset_index(drop=True),
            ],
            axis=1,
        )

        # release the memory from the status_dict
        del status_dict

    def filter_cm_threshold(self, cM_threshold: int, len_index: int) -> None:
        """Method that will filter the self.ibd_df to only individuals larger than the specified threshold. This should be run after the load_file method.

        Parameters

        cM_threshold : int
            threshold to filter the file lengths on

        len_index : index
            column index for the hapibd or ilash file to find the
            lengths of each segment for
        """
        logger.info(
            f"Filtering the dataframe of shared segments to greater than or equal to {cM_threshold}cM"
        )
        logger.info(self.ibd_df)
        self.ibd_df = self.ibd_df[self.ibd_df[len_index] >= cM_threshold]

    def _find_all_grids(self, indices) -> List[str]:
        """Method that will take the dataframe that is filtered for the location and the haplotype and return
        a list of all unique individuals in that dataframe

        Parameters

        dataframe : pd.DataFrame
            dataframe that has the pairs of individuals who share an ibd segment. This dataframe is filtered for a
            specific loci and matching haplotype phase.

        indices : FileInfo
            object that has all the indices for the file of interest

        Returns

        List[str]
            returns a list of unique individuals in the dataframe
        """

        self.ibd_df["ind_1"] = (
            self.ibd_df[indices.id1_indx]
            + "."
            + self.ibd_df[indices.id1_phase_indx].astype(str)
        )

        self.ibd_df["ind_2"] = (
            self.ibd_df[indices.id2_indx]
            + "."
            + self.ibd_df[indices.id2_phase_indx].astype(str)
        )

        # adding these new categories as indices to the indices object
        indices.ind1_with_phase = "ind_1"
        indices.ind2_with_phase = "ind_2"

        grids: List[str] = list(
            set(
                self.ibd_df["ind_1"].values.tolist()
                + self.ibd_df["ind_2"].values.tolist()
            )
        )
        logger.debug(f"Found {len(grids)} unique individual that will be clustered")

        return grids

    def _find_secondary_connections(
        self,
        new_individuals: List[str],
        exclusion: Set[str],
        indices: FileInfo,
        network: Network,
        inds_in_network: Set[str],
    ) -> None:
        """Function that will find the secondary connections within the graph

        Parameters

        new_individuals : List[str]
            List of iids that are not the exclusion iid. The list will be
            all individuals who share with the iid

        exclusion : Set[str]
            This is a set of iids that we want to keep out of the secondary
            cluster because we seeded from this iid so we have these
            connections

        indices : FileInfo
            object that has all the indices values for the correct column
            in the hapibd file

        network : Network
            class that has attributes for the pairs, iids, and haplotypes


        """
        # filtering the dataframe with the new individuals
        second_filter: pd.DataFrame = network.filter_for_seed(
            self.ibd_df, new_individuals, indices, exclusion
        )

        if not second_filter.empty:
            # getting a list of all new individuals
            new_connections: List[str] = list(
                network.gather_grids(
                    second_filter, indices.ind1_with_phase, indices.ind2_with_phase
                )
            )

            network.update(second_filter, indices)

            # need to create a list of people in the network
            inds_in_network.update(new_connections)

            # updating the exclusion set so that the individuals that we seeded of of this iteration are added to it
            exclusion.update(new_individuals)

            # now need to recursively find more connections
            self._find_secondary_connections(
                [iid for iid in new_connections if iid not in new_individuals],
                exclusion,
                indices,
                network,
                inds_in_network,
            )

    def find_networks(
        self, gene_name: str, chromosome: str, indices: FileInfo
    ) -> List[Network]:
        """Method that will go through the dataframe and will identify networks.

        Parameters

        indices : FileInfo
            object that has all the indices values for the correct column in the hapibd file

        Returns

        Dict[int: Dict]
            returns a dictionary with the following structure:
            {network_id: {pairs: List[Pair_objects], in_network: List[str], haplotypes_list: List[str]}}. Other data information will be added to this file
        """
        # First need to get a list of all the unique individuals in the ibd_file
        iid_list: List[str] = self._find_all_grids(indices)

        # determining values that will be used in the function
        network_id = 1
        inds_in_network = set()

        # making a list that will be returned that has the Network objects as values
        network_list = []

        # count = 1
        # iterate over each iid in the original dataframe
        # creating a progress bar
        for ind in tqdm(iid_list, desc="pairs in clusters"):
            # if this iid has already been associated with a network then we need to skip it. If not then we can get the network connected to it
            if ind not in inds_in_network:

                network_obj = Network(gene_name, chromosome, network_id)

                # filtering the network object for the
                # connections to the first seed
                filtered_df = network_obj.filter_for_seed(self.ibd_df, [ind], indices)

                network_obj.update(filtered_df, indices)

                # need to create a list of people in the network. This will be individuals with the phase
                inds_in_network.update(
                    network_obj.gather_grids(
                        filtered_df, indices.ind1_with_phase, indices.ind2_with_phase
                    )
                )

                # finding the secondary connections. These are connections to all
                # individuals in the self.individuals_in_network set that are not the
                # seeding individual
                logger.debug(
                    f"Finding secondary connections to {ind} in network {network_id}"
                )

                self._find_secondary_connections(
                    [iid for iid in network_obj.haplotypes if iid != ind],
                    set([ind]),
                    indices,
                    network_obj,
                    inds_in_network,
                )

                logger.debug(
                    f"number of iids in network {network_id}: {len(network_obj.iids)}"
                )

                network_id += 1

                # adding the full network to the list
                network_list.append(network_obj)

            # if the program is being run in debug mode then it will only this loop four times. This gives enough information
            # to see if the program is behaving properly
            if int(os.environ.get("program_loglevel")) == 10:
                if self.count == 3:
                    break
                self.count += 1

        return network_list
