"""Module to faciliate the user in parsing the phenotype file by incorporating multiple 
ecodings, separators, and by handling multiple errors."""

import gzip
from logging import Logger
from pathlib import Path
from typing import Dict, List, Set, Tuple, TypeVar, Union

from drive.log import CustomLogger

logger: Logger = CustomLogger.get_logger(__name__)

# creating a type annotation for the PhenotypeFileParser class
T = TypeVar("T", bound="PhenotypeFileParser")


class PhenotypeFileParser:
    """Parser used to read in the phenotype file. This will allow use to account for
    different delimiters in files as well as catch errors."""

    def __init__(self, filepath: Union[Path, str]) -> None:
        """Initialize the PhenotypeFileParser class.

        Parameters
        ----------
        filepath : Path | str
            filepath to the phenotype file that has case control status for individuals

        Raises
        ------
        FileNotFoundError
        """
        self.individuals: List[str] = []
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
        OSError
            Raised if the program encounters an OSError while opening the file
        """
        suffix = self.file.suffix

        try:
            if suffix == ".gz":
                file = gzip.open(self.file, "rt")
            else:
                file = open(self.file, "r", encoding="utf-8")
        except OSError as e:
            logger.critical(e)
            logger.critical(
                f"Encountered the following error while trying to open the file: {self.file}"  # noqa: E501
            )
            raise OSError(
                f"Encountered the following error while trying to open the file: {self.file}"  # noqa: E501
            )

        self.opened_file = file

        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Close the resource when it is not longer being used. Used by the context
        manager."""
        self.opened_file.close()

    @staticmethod
    def _check_separator(line: str) -> str:
        """Determine what the separator for the file should be.

        Parameters
        ----------
        line : str
            this is the first or second line of the phenotype file
            depending on whether or not the file has a header.

        Returns
        -------
        str
            returns the separator as a string. At the moment, comma
            separated, tab separated, and pipe separated are supported.

        Raises
        ------
        ValueError
            raises a value error if the method doesn't identify a
            supported separator.
        """
        if len(line.split(",")) > 1:
            return ","
        elif len(line.split("\t")) > 1:
            return "\t"
        elif len(line.split("|")) > 1:
            return "|"
        else:
            raise ValueError(
                "The was no appropriate separator found for the file. Currently DRIVE supports: ',' or '\t' or '|'"  # noqa: E501
            )

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

        phenotype_dict : Dict[str, Dict[str, Set[str]]]

            returns a tuple with three elements. The first element is a
            dictionary where the keys are phenotypes. Values are
            dictionaries where the keys are 'cases' or 'controls' or
            'excluded' and values are list of ids. The second element is
            a dictionary that maps the index of the phenotype in the
            header line to the phenotype name. The third element is the
            separator string
        """
        # pull out the correct grid id
        grid_id = line[0]
        # we need to keep track of the total list of grids used in the
        # analysis
        self.individuals.append(grid_id)

        # go through each value in the file
        for indx, value in enumerate(line[1:]):
            phenotype_mapping = phenotype_indx.get(indx)

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
                phenotype_dict[phenotype_mapping]["excluded"].append(grid_id)

    def _create_phenotype_dictionary(
        self,
        header_line: str,
    ) -> Tuple[Dict[str, Dict[str, Set[str]]], Dict[int, str], str]:
        """Function that will generate a dictionary where the keys are
        phenotypes and the values list of cases/exclusions/controls

        Parameters
        ----------
        header_line : str
            line from the phenotype file

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
        """

        # determining what the appropriate separator should be
        separator = PhenotypeFileParser._check_separator(header_line)

        logger.debug(f"Identified the separator, {separator}, in the file: {self.file}")

        # raise an error if there is no header line, otherwise determine all the
        # phenotypes
        if "grid" not in header_line.lower() and "grids" not in header_line.lower():
            error_msg = "Expected the first line of the phenotype file to have a header line with a column called grid or grids."  # noqa: E501

            logger.critical(error_msg)

            raise ValueError(error_msg)
        else:
            split_line_phenotypes = header_line.strip("\n").split(separator)[1:]

        # creating a dictionary that will map an index position to a
        # phecode
        phenotype_indx = {}
        # creating a dictionary to keep track of who are cases and
        # controls
        phenotype_dict = {}

        # build each dictionary
        for indx, phenotype in enumerate(split_line_phenotypes):
            phenotype_indx[indx] = phenotype
            phenotype_dict[phenotype] = {
                "cases": set(),
                "controls": set(),
                "excluded": set(),
            }

        logger.debug(f"Phenotype index dictionary:\n {phenotype_indx}")
        logger.debug(f"Phenotype counts dictionary: \n {phenotype_dict}")
        return phenotype_dict, phenotype_indx, separator

    def parse_cases_and_controls(
        self,
    ) -> Tuple[Dict[str, Dict[str, Set[str]]], List[str]]:
        """Generate a list for cases, controls, and excluded individuals.

        Returns
        -------
        Tuple[Dict[str, Dict[str, List[str]]], List[str]]
            returns a tuple where the first element is a dictionary where
            the keys are the phenotypes and the values are dictionary of
            the case/controls/excluded individuals lists. The second
            element is a list of all grids from the file to be used as a
            cohort
        """

        (
            phenotype_dict,
            phenotype_indx_mappings,
            separator,
        ) = self._create_phenotype_dictionary(self.opened_file.readline())

        for line in self.opened_file:
            # we need to first check if there is a header row
            split_line = line.strip("\n").split(separator)

            self._determine_status(split_line, phenotype_dict, phenotype_indx_mappings)

        return phenotype_dict, self.individuals
