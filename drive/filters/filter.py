import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, TypeVar

from log import CustomLogger
from models import FileIndices, Genes
from pandas import DataFrame, concat, read_csv

logger = CustomLogger.get_logger(__name__)

# we are going to create two exception class for the vertex

T = TypeVar("T", bound="IbdFilter")


@dataclass
class IbdFilter:
    ibd_file: Iterator[DataFrame]
    indices: FileIndices
    target_gene: Genes
    ibd_vs: DataFrame = field(default_factory=DataFrame)
    ibd_pd: DataFrame = field(default_factory=DataFrame)
    hapid_map: Dict[str, int] = field(default_factory=dict)
    all_haplotypes: List[str] = field(default_factory=list)
    haplotype_id: int = 0

    @classmethod
    def load_file(
        cls,
        ibd_file: Path,
        indices: FileIndices,
        target_gene: Genes,
    ) -> T:
        """Factory method that returns the IBDFilter model
        This method makes sure that the ibd file exists

        Parameters
        ----------
        ibd_file : Path
            Path object containing the filepath for the ibd
            file from hapibd, iLASH, etc...

        indices: FileIndices
            Object that has all the indices for the necessary
            columns in the ibd file. This object also has a
            get_haplotype_ids method that will form the
            haplotype ibd.

        target_gene : Genes
            namedtuple that has attributes for the
            chromosome, the gene start position, and the
            gene end position.

        Returns
        -------
        IbdFilter
            returns an initialized IbdFilter object

        Raises
        ------
        FileNotFoundError
            raises an error if the file doesn't exist
        """
        logger.verbose(f"Reading in the ibd input file at {ibd_file}")

        if not ibd_file.is_file():
            raise FileNotFoundError(f"The file, {ibd_file}, was not found")

        input_file_chunks = read_csv(ibd_file, sep="\t", header=None, chunksize=100_100)

        return cls(input_file_chunks, indices, target_gene)

    def _generate_map(self, chunk_data: DataFrame) -> None:
        """Method that will generate the dictionary that maps hapibd to integers

        Parameters
        ----------
        data_chunk : pd.DataFrame
            chunk of the ibdfile. The size of this chunk is determined by the
            chunksize argument to pd.read_csv. This value is currently set to 100,000.
        """
        haplotypes = chunk_data.values.ravel()

        logger.verbose(f"identified {len(haplotypes)} haplotypes.")

        # iterate over each haplotype and add it to the dictionary if the value is not present
        for value in haplotypes:
            key_value = self.hapid_map.setdefault(value, self.haplotype_id)
            if key_value == self.haplotype_id:
                self.haplotype_id += 1

    def _map_grids(self, data_chunk: DataFrame) -> None:
        """Function responsible for creating two new columns that
        have the haplotype id numbers.

        data_chunk : pd.DataFrame
            chunk of the ibdfile. The size of this chunk is
            determined by the chunksize argument to
            pd.read_csv. This value is currently set to 100,000.
        """
        # we are going to map the haplotype id to integers in a new
        # column. this needs to be done for both hapid1 an hapid2
        data_chunk.loc[:, "idnum1"] = data_chunk["hapid1"].map(self.hapid_map)

        data_chunk.loc[:, "idnum2"] = data_chunk["hapid2"].map(self.hapid_map)

    def _filter(self, data_chunk: DataFrame, min_cm: int) -> DataFrame:
        """Method that will filter the ibd file on four conditions: Chromosome number is the same, segment start position is <= target start position, segment end position is >= to the start position, and the size of the segment is >= to the minimum centimorgan threshold.

        Parameters
        ----------
        data_chunk : pd.DataFrame
            chunk of the ibdfile. The size of this chunk is
            determined by the chunksize argument to
            pd.read_csv. This value is currently set to 100,000.

        min_cm : int
            centimorgan threshold

        Returns
        -------
        pd.DataFrame
            returns the filtered dataframe

        Raises
        ------
        ValueError
            raises a ValueError if the target chromosome number is not
            found within the provided IBD file. This situation will
            lead to a error later in the program which is why the
            exception is raised. It is assumed to be due the user
            providing the incorrect file by accident
        """
        # we are going to first make sure that the ibd file is for the
        # right chromosome. If the target_gene chromosome number is not
        # found in the file then a ValueError is raised.
        if self.target_gene.chr not in data_chunk[self.indices.chr_indx].values:
            error_msg = f"Expected the value of the chromosome column in the ibd file to be {self.target_gene.chr}. This value was not found in the column. Please ensure that you selected the proper IBD file for chromosome {self.target_gene.chr} before rerunning DRIVE."

            logger.critical(error_msg)

            raise ValueError(error_msg)
        # We are going to filter the data and then make a copy
        # of it to return so that we don't get the
        # SettingWithCopyWarning
        return data_chunk[
            (data_chunk[self.indices.chr_indx] == self.target_gene.chr)
            & (data_chunk[self.indices.str_indx] <= self.target_gene.start)
            & (data_chunk[self.indices.end_indx] >= self.target_gene.end)
            & (data_chunk[self.indices.cM_indx] >= min_cm)
        ].copy()

    def _remove_dups(self, data: DataFrame) -> DataFrame:
        """Filters out rows where the haplotype ids are the
        same

        Parameters
        ----------
        data_chunk : pd.DataFrame
            chunk of the ibdfile. The size of this chunk is
            determined by the chunksize argument to
            pd.read_csv. This value is currently set to 100,000.

        Returns
        -------
        DataFrame
            returns a dataframe where we make sure that hapibd1 does != hapibd2

        Raises
        ------
        KeyError
            raises a key error if hapid1 or hapid2 are not columns in the dataframe
        """
        try:
            return data[data["hapid1"] != data["hapid2"]]
        except KeyError as e:
            logger.critical(
                f"Expected the keys hapid1 and hapid2 to be in the dataframe. Instead the only keys were: {', '.join(data.columns)}"
            )
            raise e(
                f"Expected the keys hapid1 and hapid2 to be in the dataframe. Instead the only keys were: {', '.join(data.columns)}"
            )

    def _generate_vertices(self, data_chunk: DataFrame) -> None:
        """Method that will generate the vertices dataframe which just has the columns idnum, hapID, and IID

        Parameters
        ----------
        data_chunk : pd.DataFrame
            chunk of the ibdfile. The size of this chunk is
            determined by the chunksize argument to
            pd.read_csv. This value is currently set to 100,000.
        """
        id1_df = data_chunk[["idnum1", "hapid1", self.indices.id1_indx]].rename(
            columns={"idnum1": "idnum", "hapid1": "hapID", 0: "IID"}
        )

        self.ibd_vs = concat([self.ibd_vs, id1_df])

        id2_df = data_chunk[["idnum2", "hapid2", self.indices.id2_indx]].rename(
            columns={"idnum2": "idnum", "hapid2": "hapID", 2: "IID"}
        )

        self.ibd_vs = concat([self.ibd_vs, id2_df])

    def _filter_for_cohort(
        self, chunk: DataFrame, cohort_ids: Optional[List[str]] = None
    ) -> DataFrame:
        """filter cohort chunk to individuals in the cohort
        list

        Parameters
        ----------
        chunk : DataFrame
            chunk of pandas dataframe that has information about the shared
            pairwise IBD segment.

        cohort_ids : List[str]
            Lists of ids that make up the cohort. The ibd_file will be filtered to only this list.

        Returns
        -------
        DataFrame
            returns the filtered pandas dataframe
        """
        # if no cohort ids were provided then we just return the chunk, otherwise we
        # filter the dataframe for where id1 and id2 are in the cohort
        if not cohort_ids:
            return chunk
        else:
            return chunk[
                (chunk[self.indices.id1_indx].isin(cohort_ids))
                & (chunk[self.indices.id2_indx].isin(cohort_ids))
            ]

    def _check_for_no_shared_segments(ibd_pd: DataFrame, ibd_vs: DataFrame) -> None:
        """Check to ensure that there were shared IBD segments
        found based on the input conditions

        Parameters
        ----------
        ibd_pd : DataFrame
            dataframe that has all the pairwise shared
            segments that satisfy the input conditions

        ibd_vs : DataFrame
            dataframe that has the information about each
            haplotype and individual that are in ibd_pd

        Raises
        ------
        ValueError
            raises a ValueError if the dataframes are empty
            because they cannot be empty or other steps of the
            program will fail.
        """
        if ibd_pd.empty:
            logger.critical(
                "There were no shared IBD segments that satisfied the provided input conditions"
            )
            raise ValueError(
                "There were no shared IBD segments that satisfied the provided input conditions. Please check to make sure that there is no typo in the target region or that the correct ibd input file was provided."
            )
        elif ibd_vs.empty:
            logger.critical(
                "There were no haplotype ids identified in the filtering step."
            )
            raise ValueError(
                "There were no haplotype ids identified in the filtering step."
            )
        else:
            logger.debug(
                f"Identified {ibd_pd.shape[0]} shared ibd segments with {ibd_vs.shape[0]} unique haplotypes"
            )

    def _check_empty_dataframes(self) -> None:
        """Check if the provided dataframe is empty. If it is then it
        raises an error and exits the program"""
        if self.ibd_pd.empty:
            logger.info(
                "No individuals from the analysis cohort share an IBD segment across the provided target region"
            )
            sys.exit(0)

    def preprocess(self, min_centimorgan: int, cohort_ids: Optional[List[str]] = None):
        """Method that will filter the ibd file.

        Parameters
        ----------
        min_centimorgan : int
            Minimum segment threshold that is used to filter
            the ibd file. Program only keeps segments that
            are greater than or equal to the threshold.

        cohort_ids : List[str]
            Lists of ids that make up the cohort. The ibd_file will be filtered to only this list.
        """
        for chunk in self.ibd_file:
            cohort_restricted_chunk = self._filter_for_cohort(chunk, cohort_ids)

            if cohort_restricted_chunk.empty:
                continue

            size_filtered_chunk = self._filter(cohort_restricted_chunk, min_centimorgan)

            if not size_filtered_chunk.empty:
                # We have to add two column with the haplotype ids
                self.indices.get_haplotype_id(
                    size_filtered_chunk,
                    self.indices.id1_indx,
                    self.indices.hap1_indx,
                    "hapid1",
                )

                self.indices.get_haplotype_id(
                    size_filtered_chunk,
                    self.indices.id2_indx,
                    self.indices.hap2_indx,
                    "hapid2",
                )
                # We then need to make sure that there are no duplicates in the dataframe
                removed_dups = self._remove_dups(size_filtered_chunk)

                # We need to update the mappings for the grids
                self._generate_map(removed_dups[["hapid1", "hapid2"]])

                self._map_grids(removed_dups)
                # concat the filtered dataframe with the ibd_pd attribute to get all the edges in the graph
                self.ibd_pd = concat([self.ibd_pd, removed_dups])

                self._generate_vertices(removed_dups)

        self._check_empty_dataframes()

        self.ibd_pd.reset_index(drop=True, inplace=True)
        self.ibd_vs = self.ibd_vs.drop_duplicates().sort_values(by="idnum")