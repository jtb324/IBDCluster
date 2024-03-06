"""Module to faciliate the user in parsing the phenotype file by incorporating multiple 
ecodings, separators, and by handling multiple errors."""

import gzip
from logging import Logger
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, TypeVar, Union
import pandas as pd
import numpy as np

from drive.log import CustomLogger

logger: Logger = CustomLogger.get_logger(__name__)

# creating a type annotation for the PhenotypeFileParser class
T = TypeVar("T", bound="PhenotypeFileParser")

PhenotypeInfo = TypeVar("PhenotypeInfo", bound=dict)


class PhenotypeFileParser:
    """Parser used to read in the phenotype file. This will allow use to account for
    different delimiters in files as well as catch errors."""

    def __init__(
        self, filepath: Union[Path, str], phenotype_name: Optional[str] = None
    ) -> None:
        """Initialize the PhenotypeFileParser class.

        Parameters
        ----------
        filepath : Path | str
            filepath to the phenotype file that has case control status for individuals

        phenotype_name : str
            Phenotype name that can be used specify a specific column in a
            phenotype matrix if the user only wants ot focus on 1 phenotype.

        Raises
        ------
        FileNotFoundError
        """
        self.specific_phenotype: str = phenotype_name
        # we are going to make sure the filepath variable is a
        # PosixPath
        filepath = Path(filepath)

        # now we are going to try to create an attribute for
        # the input file
        if not filepath.exists():
            raise FileNotFoundError(f"The file {filepath} was not found")
        else:
            self.file: Path = filepath

    def __enter__(self) -> T:
        """Open the input file. Method determines the appropriate file type and open
        the file. Method is called automatically by the context manager.

        Raises
        ------
        pd.errors.ParserError
            Raised if the program encounters any error while trying to read in the
            phenotype matrix using pd.read_csv
        """
        try:
            pheno_df = pd.read_csv(
                self.file,
                sep="\t",
                na_values=["na", "n/a", "-1", "-1.0", " ", "", "NA", "N/A"],
            ).fillna("-1")

        except pd.errors.ParserError as e:
            logger.critical(e)
            logger.critical(
                f"Encountered the following error while trying to read in the phenotype matrix: {self.file}"  # noqa: E501
            )
            raise pd.errors.ParserError(
                f"Encountered the following error while trying to read in the phenotype matrix: {self.file}"  # noqa: E501
            )

        self.phenotype_df = pheno_df

        logger.verbose(
            f"Reading in {self.phenotype_df.shape[1] - 1} phecodes for {self.phenotype_df.shape[0]} individuals"
        )  # noqa: E501

        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Close the resource when it is not longer being used. Used by the context manager."""
        logger.verbose("Finished reading in individuals from the phenotype file")

    def _process_matrix(
        self,
    ) -> Tuple[Dict[str, Dict[str, Set[str]]], List[str]]:
        """Function that will generate a dictionary where the keys are
        phenotypes and the values are lists of the cases/exclusions/controls


        Returns
        -------
        Tuple[Dict[str, Dict[str, Set[str]]], List[str]]
            returns a tuple with three elements. The first element is a
            dictionary where the keys are phenotypes. Values are
            dictionaries where the keys are 'cases' or 'controls' or
            'excluded' and values are list of ids. The second element is
            a dictionary that maps the index of the phenotype in the
            header line to the phenotype name. The third element is the
            separator string
        """

        # we need to create the dictionary that will have the case/control/exclusion
        # ids for each phecode
        phenotype_dict = {phecode: None for phecode in self.phenotype_df.columns[1:]}

        # we ultimately are going to just want to filter this grids series
        grids = self.phenotype_df.iloc[:, 0]

        # We will pull out values for each phenotype to determine the cases/controls/
        # exclusions
        for phecode_name, phenotyping_status in self.phenotype_df.iloc[:, 1:].items():
            cases = set(
                grids[phenotyping_status[phenotyping_status == 1].index].unique()
            )

            controls = set(
                grids[phenotyping_status[phenotyping_status == 0].index].unique()
            )

            exclusions = set(
                grids[phenotyping_status[phenotyping_status == -1].index].unique()
            )

            phenotype_dict[phecode_name] = {
                "cases": cases,
                "controls": controls,
                "excluded": exclusions,
            }

        return phenotype_dict, grids.values.tolist()

    def parse_cases_and_controls(
        self,
    ) -> Tuple[Dict[str, Dict[str, Set[str]]], Dict[int, str]]:
        # ) -> PhenotypeInfo:
        """Generate a list for cases, controls, and excluded individuals.

        Returns
        -------
        Tuple[Dict[str, Dict[str, Set[str]]], List[str]]
            returns a tuple where the first element is a dictionary where
            the keys are the phenotypes and the values are dictionary of
            the case/controls/excluded individuals lists. The second
            element is a list of all grids from the file to be used as a
            cohort
        """

        phenotyping_dictionary, cohort_ids = self._process_matrix()

        return phenotyping_dictionary, cohort_ids
