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

    # @staticmethod
    # def _check_separator(line: str) -> str:
    #     """Determine what the separator for the file should be.

    #     Parameters
    #     ----------
    #     line : str
    #         this is the first or second line of the phenotype file
    #         depending on whether or not the file has a header.

    #     Returns
    #     -------
    #     str
    #         returns the separator as a string. At the moment, comma
    #         separated, tab separated, and pipe separated are supported.

    #     Raises
    #     ------
    #     ValueError
    #         raises a value error if the method doesn't identify a
    #         supported separator.
    #     """
    #     if len(line.split(",")) > 1:
    #         return ","
    #     elif len(line.split("\t")) > 1:
    #         return "\t"
    #     elif len(line.split("|")) > 1:
    #         return "|"
    #     else:
    #         raise ValueError(
    #             "The was no appropriate separator found for the file. Currently DRIVE supports: ',' or '\t' or '|'"  # noqa: E501
    #         )

    def _determine_status(
        self,
        line: list[str],
        phenotype_dict: Dict[str, Dict[str, Set[str]]],
        phenotype_indx: Dict[int, str],
    ) -> None:
        """Add the individual to the appropriate case/ control/exclusion list.

        Parameters
        ----------
        line: List[str]
            list of individuals status for each phenotype in the file

        line_indx : int
            index of the line in the phenotype file that is currently being processed.
            This index is being used to map the individual id to a unique number and
            therefore the index will be added to the numpy arrays.

        phenotype_dict : Dict[str, Dict[str, Set[str]]]

            returns a tuple with three elements. The first element is a
            dictionary where the keys are phenotypes. Values are
            dictionaries where the keys are 'cases' or 'controls' or
            'excluded' and values are list of ids. The second element is
            a dictionary that maps the index of the phenotype in the
            header line to the phenotype name. The third element is the
            separator string

        phenotype_indx : Dict[int, str]
            Dictionary that maps the index that the phecode is in the file to
            the PheCode label
        """
        # pull out the correct grid id
        grid_id = line[0]

        # we need to keep track of the total list of grids used in the
        # analysis
        self.individuals.append(grid_id)

        # go through each value in the file
        for indx, value in enumerate(line[1:]):
            # If the user only wants to specify a specific column then there is a
            # possible of the indx number not being in the phenotype indx. We can
            # use the walrus operator to make sure that we only attempt to
            # determine the case control status for the phenotype if .get function
            # returns a value and not none
            if phenotype_mapping := phenotype_indx.get(indx):
                if value == "1" or value == "1.0":
                    phenotype_dict[phenotype_mapping]["cases"].add(grid_id)
                elif value == "0" or value == "0.0":
                    phenotype_dict[phenotype_mapping]["controls"].add(grid_id)
                # we are going to excluded on values na, n/a, -1, -1.
                # 0, "", " " to try to catch different values
                elif value.lower() in ["na", "n/a", "-1", "-1.0", " ", ""]:
                    phenotype_dict[phenotype_mapping]["excluded"].add(grid_id)
                else:
                    logger.warning(
                        f"The status for individual, {grid_id}, was not recognized. The status found in the file was {value} for phenotype {phenotype_mapping}. This individual will be added to the exclusion list but it is recommended that the user checks to ensure that this is not a typo in the phenotype file."  # noqa: E501
                    )
                    phenotype_dict[phenotype_mapping]["excluded"].add(grid_id)

        logger.debug(f"Phenotype dictionary after adding individuals: {phenotype_dict}")

    def _determine_cases_and_controls(pheno_col: pd.Series, grid: pd.Series) -> None:
        """determine cases/controls/exclusions for each phenotype"""

    def _process_matrix(
        self,
    ) -> Tuple[Dict[str, Dict[str, Set[str]]], List[str]]:
        """Function that will generate a dictionary where the keys are
        phenotypes and the values are lists of the cases/exclusions/controls

        Parameters
        ----------
        header_line : str
            first line from the phenotype file. This line has all the PheCodes
            in the file

        Returns
        -------
        Tuple[Dict[str, Dict[str, Set[str]]], Dict[int, str], str]
            returns a tuple with three elements. The first element is a
            dictionary where the keys are phenotypes. Values are
            dictionaries where the keys are 'cases' or 'controls' or
            'excluded' and values are list of ids. The second element is
            a dictionary that maps the index of the phenotype in the
            header line to the phenotype name. The third element is the
            separator string

        Raises
        ------
        ValueError
            This function will raise a Value error if the value 'grid' or 'grids'
            (uppercase or lowercase) is not found in the header line. The function
            can also raise a value error if the user tries to specify a phenotype
            name that is not in the header line.
        """

        # we need to create the dictionary that has the case and control values
        # for each phecode

        phenotype_dict = {phecode: None for phecode in self.phenotype_df.columns[1:]}

        grids = self.phenotype_df.iloc[:, 0]

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
        # logger.debug(f"Identified the separator, {separator}, in the file: {self.file}")

        # # raise an error if there is no header line, otherwise determine all the
        # # phenotypes
        # if "grid" not in header_line.lower() and "grids" not in header_line.lower():
        #     error_msg = "Expected the first line of the phenotype file to have a header line with a column called grid or grids."  # noqa: E501

        #     logger.critical(error_msg)

        #     raise ValueError(error_msg)
        # else:
        #     split_line_phenotypes = header_line.strip("\n").split(separator)[1:]

        # # creating a dictionary that will map an index position to a
        # # phecode
        # phenotype_indx = {}
        # # creating a dictionary to keep track of who are cases and
        # # controls
        # phenotype_dict = {}

        # # we are going to preallocate numpy arrays so we need to determine the nubmer of
        # # phenotypes in the provided file
        # phenotype_count_in_file = len(split_line_phenotypes)

        # # build a dictionary for each phenotype of the cases and the controls.
        # # If the user has specified a phenotype, we will only add that phenotype
        # # and we will add its corresponding index to the phenotype_indx dictionary
        # # This block will catch a value error if the user provides a phenotype name
        # # that is not in the file

        # logger.debug(
        #     f"{phenotype_count_in_file} phenotypes were provided in the file header. The phenotype dictionary will contain numpy arrays preallocated to this size."
        # )  # noqa: E501

        # if self.specific_phenotype:
        #     try:
        #         indx = split_line_phenotypes.index(self.specific_phenotype)
        #     except ValueError:
        #         logger.critical(
        #             f"The value {self.specific_phenotype} was not found in one of the phenotype column files. Please make sure you spelled the phenotype name the exact same as it is in the phenotype file."  # noqa: E501
        #         )
        #         raise ValueError(
        #             "The value {self.specific_phenotype} was not found in the phenotype file."  # noqa: E501
        #         )
        #     else:
        #         phenotype_indx[indx] = self.specific_phenotype
        #         phenotype_dict[self.specific_phenotype] = {
        #             "cases": set(),
        #             "controls": set(),
        #             "excluded": set(),
        #         }
        # else:
        #     for indx, phenotype in enumerate(split_line_phenotypes):
        #         phenotype_indx[indx] = phenotype
        #         phenotype_dict[phenotype] = {
        #             "cases": set(),
        #             "controls": set(),
        #             "excluded": set(),
        #         }

        # logger.debug(f"Phenotype index dictionary:\n {phenotype_indx}")
        # logger.debug(
        #     f"Phenotype counts dictionary has {len(phenotype_dict.keys())} PheCodes as keys after creation."  # noqa: E501
        # )
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

        # (
        #     phenotype_dict,
        #     phenotype_indx_mappings,
        #     separator,
        # ) = self._create_phenotype_dictionary(self.opened_file.readline())

        # logger.info(
        #     "Reading in cases, controls, and excluded individuals from the provided phenotype file."
        # )  # noqa: E501

        phenotyping_dictionary, cohort_ids = self._process_matrix()
        # for indx, line in enumerate(self.opened_file, start=1):
        #     # we need to first check if there is a header row
        #     logger.debug(f"Reading in phecodes for individual, {indx}, in the file")

        #     split_line = line.strip("\n").split(separator)

        #     # self._record_individual(split_line, indx)

        #     self._determine_status(
        #         split_line, phenotype_dict, phenotype_indx_mappings
        #     )

        return phenotyping_dictionary, cohort_ids
