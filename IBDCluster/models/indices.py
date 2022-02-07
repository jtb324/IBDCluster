from abc import ABC, abstractmethod
from dataclasses import dataclass
import os
from typing import List
from glob import glob

@dataclass
class File_Info(ABC):
    """The child classes will have attributes such as:
            'id1_indx': 0,
            'id2_indx': 2,
            'chr_indx': 4,
            'str_indx': 5,
            'end_indx': 6
    """
    @abstractmethod
    def load_indices(self):
        pass
    
    @abstractmethod
    def gather_files(self):
        pass

    @abstractmethod
    def find_file(self):
        """Function that will be used to find a specfic ibd from the correct directory"""
        pass

@dataclass
class Hapibd_Info(File_Info):
    """Class that has all of the indices as well as file paths for the hapibd files."""
    file_dir: str = os.environ.get("hapibd_files", None)
    ext: str = "*.ibd.gz"

    def load_indices(self):
        self.id1_indx: int = 0
        self.id1_phase_indx: int = 1
        self.id2_indx: int = 2
        self.id2_phase_indx: int = 3
        self.chr_indx: int = 4
        self.str_indx: int = 5
        self.end_indx: int = 6
        self.cM_indx: int = 7
    
    def gather_files(self):
        """Function that will gather together all the files in the file_dir that have the ext and will store them in a list"""
        self.ibd_files: List[str] = []
        for file in glob(os.path.join(self.file_dir, self.ext)):

            file: str = os.path.join(self.file_dir, file)

            self.ibd_files.append(file)
    
    def find_file(self, chr_num: str) -> str:
        """Method that will find the appropriate file based on the correct chromosome number
        
        Parameters
        
        chr_num : str
            identifier for the chromosome number this wilk take the format 'chrX' or 'chrXX' depending on 
            whether the X is < 10 or >= 10

        Returns

        str
            returns the file as a string
        """ 
        ibd_file: str = [file for file in self.ibd_files if "".join([chr_num, "."]) in file][0]

        return ibd_file

@dataclass
class Ilash_Info(File_Info):
    """Class that has all of the indices as well as file paths for the ilash files."""
    file_dir: str = os.environ.get("ilash_files", None)
    ext: str = "*.match.gz"

    def load_indices(self):
        self.id1_indx: int = 0
        self.id1_phase_indx: int = 1
        self.id2_indx: int = 2
        self.id2_phase_indx: int = 3
        self.chr_indx: int = 4
        self.str_indx: int = 5
        self.end_indx: int = 6
        self.cM_indx: int = 9

    def gather_files(self):
        """Function that will gather together all the files in the file_dir that have the ext and will store them in a list"""
        self.ibd_files: List[str] = []
        for file in glob(os.path.join(self.file_dir, self.ext)):

            file: str = os.path.join(self.file_dir, file)

            self.ibd_files.append(file)

    def find_file(self, chr_num: str) -> str:
        """Method that will find the appropriate file based on the correct chromosome number
        
        Parameters
        
        chr_num : str
            identifier for the chromosome number this wilk take the format 'chrX' or 'chrXX' depending on 
            whether the X is < 10 or >= 10

        Returns

        str
            returns the file as a string
        """ 
        ibd_file: str = [file for file in self.ibd_files if "".join([chr_num, "."]) in file][0]

        return ibd_file
    
        