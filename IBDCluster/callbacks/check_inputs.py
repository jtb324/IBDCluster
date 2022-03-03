from pathlib import Path


class IncorrectFileType(Exception):
    """Error that will be raised if the gene info file is not a text file"""

    def __init__(self, suffix: str, message: str) -> None:
        self.suffix: str = suffix
        self.message: str = message
        super().__init__(message)


class UnSupportedIBDProgram(Exception):
    """Error that will be raised if the user doesn't provide a supported ibd program"""

    def __init__(self, program: str, message: str) -> None:
        self.program: str = program
        self.message: str = message
        super().__init__(message)

class UnsupportedLogLevel(Exception):
    """Error that will be raised if the user provides an incorrect loglevel"""
    def __init__(self, loglevel: str) -> None:
        self.message: str = f"The provided loglevel, {loglevel} is not supportted. Supported values are verbose, debug, warning (Case Insensitive)"
        super().__init__(self.message)

def check_ibd_program(program: str) -> str:
    """Callback function that will check if the ibd program provided is hapibd or ilash

    Parameters

    program : str
        string that has the value either 'hapibd' or 'ilash'

    Returns

    str
        returns the program if an error is not raise
    """
    if program.lower() not in ["hapibd", "ilash"]:
        raise UnSupportedIBDProgram(
            program,
            f"The provided ibd program {program} is not supported. Supported values are 'hapibd' and 'ilash'.",
        )

    return program.lower()


def check_gene_file(gene_filepath: str) -> str:
    """Function that will check to make sure the gene info file exist or else it will raise an error.

    Parameters

    gene_filepath : str
        filepath to a text file that has information about the gene/genes of interests

    Returns

    str
        returns the filepath
    """
    file_path = Path(gene_filepath)

    if not file_path.exists():
        raise FileNotFoundError
    if gene_filepath[-4:] != ".txt":
        raise IncorrectFileType(
            gene_filepath[-4:],
            "The filetype provided for the gene info file is incorrect. Please provided a tab delimited text file",
        )

    return gene_filepath

def check_loglevel(loglevel: str) -> str:
    """Function that will check to make sure the loglevel is either info, warning, or debug
    
    Parameters
    
    loglevel : str
        parameter passed by the user to indicate what level of logging they want
    
    Returns
    
    str
        returns the lower case log level"""

    if loglevel.lower() not in ["verbose", "debug", "warning"]:
        raise UnsupportedLogLevel(loglevel)

    return loglevel.lower()
