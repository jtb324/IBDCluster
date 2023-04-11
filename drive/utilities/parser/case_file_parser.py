import gzip
from pathlib import Path
from typing import Union, TypeVar
from logging import Logger
import log

logger: Logger = log.get_logger(__name__)

# creating a type annotation for the PhenotypeFileParser class
T = TypeVar("T", bound="PhenotypeFileParser")


class PhenotypeFileParser:
    def __init__(self, filepath: Union[Path, str]) -> None:
        """Parser that will be used to read in the phenotype file.
        This will allow use to account for different delimiters in
        files as well as catch errors

        Parameters
        ----------
        filepath : Path | str
            filepath to the phenotype file that has case control status for individuals

        Raises
        ------
        FileNotFoundError
        """
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"The file {filepath} was not found")
        else:
            self.file: Path = filepath

    def __enter__(self) -> T:
        """__enter__ method that will be automatically called when the
        parser is used with the 'with' context manager. Method will
        determine the appropriate file type and open the file

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
        """method that will be called be the context manager to close the resource when it is not longer being used"""
        self.opened_file.close()
