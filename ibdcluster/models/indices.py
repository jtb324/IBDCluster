from dataclasses import dataclass, field

# from typing import List
import log

logger = log.get_logger(__name__)


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
    cM_indx: None | int = None
    ibd_files: list[None | str] = field(default_factory=list)
    # program_indices: None | ProgramIndices = None

    def set_cM_indx(self, program_name: str) -> None:
        """Method that will set the index for the cM indx

        Parameters
        ----------
        program_name : str
            This will be either ilash or hapibd

        filepath : str
            string to the directory with ibd files

        Raises
        ------
        NotImplementedError
            raises NotImplementedError if the user provides a ibd program that is not supported
        """
        match program_name:
            case "hapibd":
                self.cM_indx = 7
            case "ilash":
                self.cM_indx = 9
            case _:
                logger.critical(
                    "ibd_program {program_name} not supported. Please use either 'hap-IBD' or 'iLASH' (all lowercase)."
                )
                raise NotImplementedError
