"""
Contains the Network and Cluster class that are used during the 
clustering and network creation steps in the cluster/main.py file
"""
import collections
from dataclasses import dataclass, field
from typing import Protocol

import log
import pandas as pd
from .pairs import Pairs

logger = log.get_logger(__name__)

# class protocol that makes sure that the indices object has these methods
class FileInfo(Protocol):
    """Protocol that enforces the set_program_indices method for the
    FileInfo object"""

    id1_indx: int
    id1_phase_indx: int
    id2_indx: int
    id2_phase_indx: int
    chr_indx: int
    str_indx: int
    end_indx: int
    cM_indx: None | int

    def set_program_indices(self, program_name: str) -> None:
        """Method that will set the proper ibd_program indices"""


class Network(Protocol):
    """Protocol that will define what a Network looks like"""

    def filter_for_seed(
        self,
        ibd_df: pd.DataFrame,
        ind_seed: list[str],
        indices: FileInfo,
        exclusion: set[str] | None = None,
    ) -> pd.DataFrame:
        """Method to filter the ibd_df for the first individual. This gets the first level new_connections"""

    def _determine_pairs(ibd_row: pd.Series, indices: FileInfo) -> Pairs:
        """Method that will take each row of the dataframe and convert it into a pair object"""

    @staticmethod
    def gather_grids(
        dataframe: pd.DataFrame, pair_1_indx: int, pair_2_indx: int
    ) -> set[str]:
        """Staticmethod that will find all the unique values in two columns that the user passed in"""

    def update(self, ibd_df: pd.DataFrame, indices: FileInfo) -> None:
        """Method to update the pair attribute, the iids, and the haplotypes"""


@dataclass
class Cluster:
    """Class object that will handle preparing the data to be clustered"""

    ibd_file: str
    ibd_program: str
    indices: FileInfo
    count: int = 0  # this is a counter that is used in testing to speed up the process
    ibd_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    network_id: str = 1  # this is a id that the cluster object will use when it updates each network in the find networks function. This will be increased by 1 for each network
    inds_in_network: set[str] = field(
        default_factory=set
    )  # This attribute will be used to keep track of the individuals that are in any network

    def load_file(self, start: int, end: int, cM_threshold: int) -> None:
        """Method filters the ibd file based on location and loads this into memory as a dataframe
        attribute of the class

        Parameters

        start : int
            integer that describes the start position of the gene of interest

        end : int
            integer that desribes the end position of the gene of interest

        cM_threshold : int
            integer that gives the cM threshold for smallest ibd
            segment to allow
        """
        logger.debug(
            f"Gathering shared ibd segments that overlap the gene region from {start} to {end} using the file {self.ibd_file}"
        )
        cols = [
            self.indices.id1_indx,
            self.indices.id1_phase_indx,
            self.indices.id2_indx,
            self.indices.id2_phase_indx,
            self.indices.chr_indx,
            self.indices.str_indx,
            self.indices.end_indx,
            self.indices.cM_indx,
        ]

        for chunk in pd.read_csv(
            self.ibd_file, sep="\t", header=None, chunksize=1000000, usecols=cols
        ):

            filtered_chunk: pd.DataFrame = chunk[
                (
                    (
                        (chunk[self.indices.str_indx] <= int(start))
                        & (chunk[self.indices.end_indx] >= int(start))
                    )
                    | (
                        (chunk[self.indices.str_indx] >= int(start))
                        & (chunk[self.indices.end_indx] <= int(end))
                    )
                    | (
                        (chunk[self.indices.str_indx] <= int(end))
                        & (chunk[self.indices.end_indx] >= int(end))
                    )
                )
                & (chunk[self.indices.cM_indx] >= cM_threshold)
            ]

            if not filtered_chunk.empty:
                self.ibd_df = pd.concat([self.ibd_df, filtered_chunk])

        logger.info(f"identified {self.ibd_df.shape[0]} pairs within the gene region")

    def filter_connections(self, connection_threshold: int) -> None:
        """Function that will filter out individuals that have fewer
        connections than the threshold.

        Parameters
        ----------
        connection_threshold : int
            number of connections that an individual is required to have to
            be considered for clustering.
        """
        logger.debug(
            f"filtering out individuals who have fewer than {connection_threshold} connections"
        )
        # pull out all of the iids that are in pairs. This list should include
        # duplicates
        pairs = self.ibd_df.ind_1.values.tolist() + self.ibd_df.ind_2.values.tolist()
        # use the collections. Counter to form a dictionary type object where keys
        # are the grids and the values are how many times that grid shows up
        counts = collections.Counter(pairs)
        # iterating over each item to get the iiids that have fewer connections
        # than the threshold

        for key, value in counts.items():
            if value < connection_threshold:
                self.inds_in_network.add(key)

        logger.info(
            f"{len(self.inds_in_network)} individuals that have fewer connections than the threshold of {connection_threshold}"
        )
        logger.debug(
            f"Individuals that failed to meet the connection threshold: {self.inds_in_network}"
        )

        # Filter out these individuals from the original dataframe that way
        # we don't include them in the networks
        self.ibd_df = self.ibd_df[
            ~(self.ibd_df.ind_1.isin(self.inds_in_network))
            & ~(self.ibd_df.ind_2.isin(self.inds_in_network))
        ]

    def get_ids(self) -> list[str]:
        """Method that will return the list of individuals iids with the
        phase information

        Returns
        -------
        list[str]
            returns a list of unique individuals in the dataframe
        """
        grids = list(
            set(
                self.ibd_df["ind_1"].values.tolist()
                + self.ibd_df["ind_2"].values.tolist()
            )
        )

        logger.info(
            f"{len(grids)} unique haplotypes identified that will be clustered into networks"
        )

        return grids

    def create_unique_ids(self, indices: FileInfo) -> None:
        """Method that will take the dataframe that is and create a
        unique id that combines the phase information

        Parameters
        ----------
        indices : FileInfo
            object that has all the indices for the file of interest
        """
        logger.debug(
            "creating new unique ids that incorporate the haplotype phasing information"
        )

        if self.ibd_program.lower() == "hapibd":
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

        # This else statment will be used if we were to try to run ilash or any other ibd program because I haven't implemented the method yet
        else:
            raise NotImplementedError

    def _find_secondary_connections(
        self,
        new_individuals: list[str],
        exclusion: set[str],
        indices: FileInfo,
        network: Network,
        inds_in_network: set[str],
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
            new_connections: list[str] = list(
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

    def construct_network(self, ind: str, network_obj: Network) -> None | str:
        """Method that will go through the dataframe and will identify networks.

        Parameters
        ----------
        indices : FileInfo
            object that has all the indices values for the correct column in the hapibd file

        Returns
        -------
        str
            Returns either 'None' when a network is successfully found or it returns a string message if the individual is already in a network
        """

        # iterate over each iid in the original dataframe
        # creating a progress bar
        # if this iid has already been associated with a network then we need to skip it. If not then we can get the network connected to it

        if ind not in self.inds_in_network:

            # filtering the network object for the
            # connections to the first seed
            filtered_df = network_obj.filter_for_seed(self.ibd_df, [ind], self.indices)

            network_obj.update(filtered_df, self.indices)

            # need to create a list of people in the network. This will be individuals with the phase
            self.inds_in_network.update(
                network_obj.gather_grids(
                    filtered_df,
                    self.indices.ind1_with_phase,  # pylint: disable=no-member
                    self.indices.ind2_with_phase,  # pylint: disable=no-member
                )
            )

            # finding the secondary connections. These are connections to all
            # individuals in the self.individuals_in_network set that are not the
            # seeding individual
            logger.debug(
                "Finding secondary connections to %s in network %d",
                ind,
                self.network_id,
            )

            self._find_secondary_connections(
                [iid for iid in network_obj.haplotypes if iid != ind],
                set([ind]),
                self.indices,
                network_obj,
                self.inds_in_network,
            )

            logger.debug(
                f"number of iids in network {self.network_id}: {len(network_obj.iids)}"
            )

            self.network_id += 1
        else:
            return "Individual already in network"
