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
        exclusion: set[str] | None = None,
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
        """
        Method that will take each row of the dataframe and convert it into a pair object

        Parameters
        ----------
        ibd_row : pd.Series
            pandas series that has all the information in each row
            of the ibd file

        indices : FileInfo
            object that has the indices to pull out the correct
            information from the ibd_row for each ibd program

        Returns
        -------
        Pairs
            Returns a pairs object that has all that information as
            attributes
        """

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
