from dataclasses import dataclass, field
import os
from typing import List, Optional, Protocol
from glob import glob
import log

logger = log.get_logger(__name__)


class ProgramIndices(Protocol):
    """Interface indicating that the class has to have a gather files method"""

    file_dir: Optional[str]
    cM_indx: int

    def gather_files(self) -> List[str]:
        """Method to gather a list of output files for the correct ibd program"""
        ...


@dataclass
class HapibdInfo(ProgramIndices):
    """Class that has all of the indices as well as file paths for the hapibd files."""

    file_dir: Optional[str] = None
    cM_indx: int = 7

    def gather_files(self) -> List[str]:
        """Function that will gather together all the files in the file_dir that have the ext and will store them in a list"""
        # saving the initial directory to a variable
        start_dir = os.getcwd()

        file_list = []
        # changing the directory to where the files are stored
        os.chdir(self.file_dir)

        for file in glob(os.path.join("*.ibd.gz")):

            file: str = os.path.join(self.file_dir, file)

            file_list.append(file)
        # changing back to the initial directory
        os.chdir(start_dir)

        logger.debug(
            f"Found {len(file_list)} with the extension '.ibd.gz' in the directory {self.file_dir}"
        )

        return file_list


@dataclass
class IlashInfo(ProgramIndices):
    """Class that has all of the indices as well as file paths for the ilash files."""

    file_dir: Optional[str] = None
    cM_indx: int = 9

    def gather_files(self) -> None:
        """Function that will gather together all the files in the file_dir that have the ext and will store them in a list"""

        file_list = []

        for file in glob(os.path.join(self.file_dir, "*.match.gz")):

            file: str = os.path.join(self.file_dir, file)

            file_list.append(file)

        logger.debug(
            f"Found {len(file_list)} with the extension '.match.gz' in the directory {self.file_dir}"
        )

        return file_list


@dataclass
class FileInfo:
    """The child classes will have attributes such as:
    'id1_indx': 0,
    'id2_indx': 2,
    'chr_indx': 4,
    'str_indx': 5,
    'end_indx': 6
    """

    id1_indx: int = 0
    id1_phase_indx: int = 1
    id2_indx: int = 2
    id2_phase_indx: int = 3
    chr_indx: int = 4
    str_indx: int = 5
    end_indx: int = 6
    ibd_files: List[Optional[str]] = field(default_factory=list)
    program_indices: Optional[ProgramIndices] = None

    def set_program_indices(self, program_name: str) -> None:
        """Method that will set the program indices as either Hapibd_Info or Ilash_Info

        Parameters

        program_name : str
            This will be either ilash or hapibd

        filepath : str
            string to the directory with ibd files
        """
        if program_name == "hapibd":

            self.program_indices = HapibdInfo(os.getenv("HAPIBD_PATH"))
        else:
            self.program_indices = IlashInfo(os.getenv("ILASH_PATH"))

    @staticmethod
    def find_file(chr_num: str, file_list: List[str]) -> str:
        """Method that will find the appropriate file based on the correct chromosome number

        Parameters

        chr_num : str
            identifier for the chromosome number this wilk take the format 'chrX' or 'chrXX' depending on
            whether the X is < 10 or >= 10. If the user names the files chr09 instead of chr9 then the program will fail with a critical message

        file_list : List[str]
            list of filepaths to the ibd files of interest

        Returns

        str
            returns the file as a string
        """
        try:
            ibd_file: str = [
                file for file in file_list if "".join([chr_num, "."]) in file
            ][0]
        except IndexError as traceback:
            logger.critical(
                f"Attempts to find chromosome, {chr_num} in file list, {', '.join(file_list)}, resulted in the following index out of range error:"
            )
            logger.critical(traceback, exc_info=1)

        logger.debug(
            f"Found the ibd file, {ibd_file}, that matches the chromosome, {chr_num}"
        )

        return ibd_file
