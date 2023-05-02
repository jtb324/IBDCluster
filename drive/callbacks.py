from pathlib import Path


def check_input_exists(ibd_input_file: Path) -> Path:
    """Callback that will check that the input ibd file exists

    Parameters
    ----------
    ibd_input_file : Path
        Path object to the input ibd file. This file should be gzipped

    Returns
    -------
    Path
        returns the Path object if it exists

    Raises
    ------
    FileNotFoundError
        If the ibd input file does not exist then the program will
        immediately raise a FileNotFoundError
    """
    if ibd_input_file.exists():
        return ibd_input_file
    else:
        raise FileNotFoundError(f"The file, {ibd_input_file}, was not found")
