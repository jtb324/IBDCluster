from dataclasses import dataclass, field
import pandas as pd
from typing import Protocol
import log
import models
from .pairs import Pairs


logger = log.get_logger(__name__)

# class protocol that makes sure that the indices object has these methods
class FileInfo(Protocol):
    """Protocol that enforces these two methods for the
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
        ...


@dataclass
class Network:
    """This class is going to be responsible for the clustering of each network"""

    gene_name: str
    gene_chr: str
    network_id: int
    pairs: list[Pairs] = field(default_factory=list)
    iids: set[str] = field(default_factory=set)
    haplotypes: set[str] = field(default_factory=set)
    connections: dict[str, int] = field(default_factory=dict)
    pvalues: None | str = None

    def filter_for_seed(
        self,
        ibd_df: pd.DataFrame,
        ind_seed: list[str],
        indices: FileInfo,
        exclusion: set[str] = None,
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
            ibd_row[indices.cM_indx],
        )

    @staticmethod
    def gather_grids(
        dataframe: pd.DataFrame, pair_1_indx: int, pair_2_indx: int
    ) -> set[str]:
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
    indices: models.FileInfo
    count: int = 0  # this is a counter that is used in testing to speed up the process
    ibd_df: None | pd.DataFrame = field(default_factory=pd.DataFrame)
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

    # def filter_cm_threshold(self, cM_threshold: int, len_index: int) -> None:
    #     """Method that will filter the self.ibd_df to only individuals larger than the specified threshold. This should be run after the load_file method.

    #     Parameters

    #     cM_threshold : int
    #         threshold to filter the file lengths on

    #     len_index : index
    #         column index for the hapibd or ilash file to find the
    #         lengths of each segment for
    #     """
    #     logger.info(
    #         f"Filtering the dataframe of shared segments to greater than or equal to {cM_threshold}cM"
    #     )

    #     self.ibd_df = self.ibd_df[self.ibd_df[len_index] >= cM_threshold]

    # def map_indices
    def find_all_grids(self, indices: FileInfo) -> list[str]:
        """Method that will take the dataframe that is filtered for the location and the haplotype and return
        a list of all unique individuals in that dataframe

        Parameters
        ----------
        indices : FileInfo
            object that has all the indices for the file of interest

        Returns
        -------
        list[str]
            returns a list of unique individuals in the dataframe
        """
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

            grids: list[str] = list(
                set(
                    self.ibd_df["ind_1"].values.tolist()
                    + self.ibd_df["ind_2"].values.tolist()
                )
            )

        # This else statment will be used if we were to try to run ilash because the grids in it do not have to be formatted
        else:
            raise NotImplementedError

        logger.info(f"Found {len(grids)} unique individual that will be clustered")

        return grids

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
            # if network_obj.network_id == 107:
            #     logger.warning(f"individual: {ind}")
            # filtering the network object for the
            # connections to the first seed
            filtered_df = network_obj.filter_for_seed(self.ibd_df, [ind], self.indices)

            network_obj.update(filtered_df, self.indices)

            # need to create a list of people in the network. This will be individuals with the phase
            self.inds_in_network.update(
                network_obj.gather_grids(
                    filtered_df,
                    self.indices.ind1_with_phase,
                    self.indices.ind2_with_phase,
                )
            )

            # finding the secondary connections. These are connections to all
            # individuals in the self.individuals_in_network set that are not the
            # seeding individual
            logger.debug(
                f"Finding secondary connections to {ind} in network {self.network_id}"
            )

            self._find_secondary_connections(
                [iid for iid in network_obj.haplotypes if iid != ind],
                set([ind]),
                self.indices,
                network_obj,
                self.inds_in_network,
            )

            logger.info(
                f"number of iids in network {self.network_id}: {len(network_obj.iids)}"
            )

            self.network_id += 1
        else:
            return "duplicate found"
