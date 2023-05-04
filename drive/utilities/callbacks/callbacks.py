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


def check_json_path(json_path: Path) -> Path:
    """Callback function that creates the json path string. If the user provides a value then it uses the user provided value else it creates the path to the default file
    Parameters
    ----------
    json_path : str
        path to the json config file or an empty string
    Returns
    -------
    str
        returns the string to the file
    """

    if json_path:
        return json_path
    else:
        root_dir = Path(__file__).parent.parent.parent.parent

        config_path = root_dir / "config.json"

        if not config_path.exists():
            raise FileNotFoundError(
                f"Expected the user to either pass a configuration file path or for the config.json file to be present in the program root directory at {config_path}."
            )

        return config_path
