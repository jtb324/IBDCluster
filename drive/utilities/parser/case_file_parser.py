"""Module to faciliate the user in parsing the phenotype file by incorporating multiple ecodings, separators, and by handling multiple errors."""

import gzip
from pathlib import Path
from typing import Dict, List, Tuple, Union, TypeVar
from logging import Logger
import log

logger: Logger = log.get_logger(__name__)

# creating a type annotation for the PhenotypeFileParser class
T = TypeVar("T", bound="PhenotypeFileParser")


class PhenotypeFileParser:
    """Parser used to read in the phenotype file. This will allow use to account for different delimiters in files as well as catch errors."""

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
        # creating list to keep track of cases and controls
        self.case_list: List[str] = []
        self.control_list: List[str] = []
        self.exclusion_list: List[str] = []
        # now we are going to try to create an attribute for the input file
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"The file {filepath} was not found")
        else:
            self.file: Path = filepath

    def __enter__(self) -> T:
        """Open the input file. Method determines the appropriate file type and open the file. Method is called automatically by the context manager.

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
                f"Encountered the following error while trying to open the file: {self.file}"
            )
            raise OSError(
                f"Encountered the following error while trying to open the file: {self.file}"
            )

        self.opened_file = file

        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Close the resource when it is not longer being used. Used by the context manager."""
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
        if len(line.split(",")) == 2:
            return ","
        elif len(line.split("\t")) == 2:
            return "\t"
        elif len(line.split("|")) == 2:
            return "|"
        else:
            raise ValueError(
                "The was no appropriate separator found for the file. Currently DRIVE supports: ',' or '\t' or '|'"
            )

    def _generate_case_control(self, status_list: List[str]) -> None:
        """Add the individual to the appropriate case/ control/exclusion list.

        Parameters
        ----------
        status_list : List[str]
            list of strings where the first element is the individual id
            and the second element is the status coded as either 0, 1, N/
            A, NA, or -1.
        """
        print(status_list)
        if status_list[1] == "1":
            self.case_list.append(status_list[0])
        elif status_list[1] == "0":
            self.control_list.append(status_list[0])
        elif status_list[1].lower() in ["na", "n/a", "-1"]:
            self.exclusion_list.append(status_list[0])
        else:
            logger.warning(
                f"The status for individual, {status_list[0]}, was not recognized. The status found in the file was {status_list[1]}. This individual will be added to the exclusion list but it is recommended that the user checks to ensure that this is not a typo in the phenotype file."
            )
            self.exclusion_list.append(status_list)

    @staticmethod
    def _create_phenotype_dictionary(
        header_line: str,
    ) -> Tuple[Dict[str, Dict[str, List[str]]], str]:
        """Function that will generate a dictionary where the keys are
        phenotypes and the values list of cases/exclusions/controls

        Parameters
        ----------
        header_line : str
            line from the phenotype file

        Returns
        -------
        Tuple[Dict[str, Dict[str, List[str]]], Dict[int, str] str]
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

        if not "grid" in header_line.lower() or not "grids" in header_line.lower():
            raise ValueError(
                "Expected the first line of the phenotype file to have a header line with a column called grid or grids."
            )
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
            phenotype_dict[phenotype] = {"cases": [], "controls": [], "excluded": []}

        return phenotype_dict, phenotype_indx, separator

    def parse_cases_and_controls(self) -> Tuple[List[str], List[str], List[str]]:
        """Generate a list for cases, controls, and excluded individuals.

        Returns
        -------
        Tuple[List[str], List[str], List[str]]
            returns a tuple where the first element is a list of case
            ids, second element is a list of control ids, and the final
            element is a list of excluded individuals.
        """
        separator = ""

        phenotype_dict, separator = PhenotypeFileParser._create_phenotype_dictionary(
            self.opened_file.readline()
        )

        for line in self.opened_file:
            # we need to first check if there is a header row
            split_line = line.strip("\n").split(separator)

            # Now we can split the file using that separator
            split_line = line.strip("\n").split(separator)

            self._generate_case_control(split_line)

        return self.case_list, self.control_list, self.exclusion_list
