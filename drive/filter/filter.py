import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List

import log
from models import FileIndices, Genes
from pandas import DataFrame, concat, read_csv
from rich.progress import Progress


@dataclass
class IbdFilter:
    ibd_file: Iterator[DataFrame]
    indices: FileIndices
    target_gene: Genes
    chunk_count: int
    ibd_vs: DataFrame = field(default_factory=DataFrame)
    ibd_pd: DataFrame = field(default_factory=DataFrame)
    hapid_map: Dict[str, int] = field(default_factory=dict)
    all_haplotypes: List[str] = field(default_factory=list)
    haplotype_id: int = 0
    logger: logging.Logger = log.get_logger(__name__)

    @staticmethod
    def _determine_chunk_count(ibd_file: Path) -> int:
        """Staticmethod that will determine how many chunks the dataframe will be split into

        Parameters
        ----------
        ibd_file : Path
            Path object to the ibd file from hapIBD, iLASH, etc...

        Returns
        -------
        int
            returns the number of chunks that the dataframe will be broken into
        """
        input_file_chunks = read_csv(ibd_file, sep="\t", header=None, chunksize=100_100)

        return sum([1 for _ in input_file_chunks])

    @classmethod
    def load_file(
        cls,
        ibd_file: Path,
        indices: FileIndices,
        target_gene: Genes,
    ):
        """Factory method that can returns the cluster model.
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
        cluster
            returns an initialized cluster object

        Raises
        ------
        FileNotFoundError
            raises an error if the file doesn't exist
        """
        if not ibd_file.is_file():
            raise FileNotFoundError(f"The file, {ibd_file}, was not found")

        chunk_count = IbdFilter._determine_chunk_count(ibd_file)

        input_file_chunks = read_csv(ibd_file, sep="\t", header=None, chunksize=100_100)

        return cls(input_file_chunks, indices, target_gene, chunk_count)

    def _generate_map(self, chunk_data: DataFrame) -> None:
        """Method that will generate the dictionary that maps hapibd to integers

        Parameters
        ----------
        data_chunk : pd.DataFrame
            chunk of the ibdfile. The size of this chunk is determined by the
            chunksize argument to pd.read_csv. This value is currently set to 100,000.
        """
        haplotypes = chunk_data.values.ravel()

        self.logger.debug(f"identified {len(haplotypes)} haplotypes.")
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
        self.logger.debug("Mapping the haplotype ids to integers")
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
        """

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
        """Filters out rows where the haplotype ids are the same

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
            raises a key error if hapid1 and hapid2 are not columns in the dataframe
        """
        try:
            return data[data["hapid1"] != data["hapid2"]]
        except KeyError as e:
            self.logger.critical(
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

    def preprocess(self, min_centimorgan: int):
        """Method that will filter the ibd file.

        Parameters
        ----------
        min_centimorgan : int
            Minimum segment threshold that is used to filter the ibd file. Program only
            keeps segments that are greater than or equal to the threshold
        """
        with Progress(transient=True) as progress_bar:
            for chunk in progress_bar.track(
                self.ibd_file,
                description="Filtering ibd file: ",
                total=self.chunk_count,
            ):
                filtered_chunk = self._filter(chunk, min_centimorgan)

                if not filtered_chunk.empty:
                    # We have to add two column with the haplotype ids
                    self.indices.get_haplotype_id(
                        filtered_chunk,
                        self.indices.id1_indx,
                        self.indices.hap1_indx,
                        "hapid1",
                    )

                    self.indices.get_haplotype_id(
                        filtered_chunk,
                        self.indices.id2_indx,
                        self.indices.hap2_indx,
                        "hapid2",
                    )
                    # We then need to make sure that the
                    removed_dups = self._remove_dups(filtered_chunk)

                    # We need to update the mappings for the grids
                    self._generate_map(removed_dups[["hapid1", "hapid2"]])

                    self._map_grids(removed_dups)
                    # concat the filtered dataframe with the ibd_pd attribute to get all the edges in the graph
                    self.ibd_pd = concat([self.ibd_pd, removed_dups])

                    self._generate_vertices(removed_dups)

            self.ibd_vs = self.ibd_vs.drop_duplicates()
