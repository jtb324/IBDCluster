from pathlib import Path
import typer
import os
import toml
import re

__version__ = "1.2.1"


class IncorrectFileType(Exception):
    """Error that will be raised if the gene info file is not a text file"""

    def __init__(self, suffix: str, message: str) -> None:
        self.suffix: str = suffix
        self.message: str = message
        super().__init__(message)


class IncorrectGeneFileFormat(Exception):
    """Error that will be raised if there are formatting errors in
    the gene info file."""

    def __init__(self, incorrect_value: str, message: str) -> None:
        self.incorrect_value: str = incorrect_value
        self.message: str = message
        super().__init__(message)


def check_gene_pos_str(gene_pos_str: str) -> str:
    """Callback function that will make sure that the user inputs a gene_pos_str of for X:Start-End. If not then it raises a value error

    Parameters
    ----------
    gene_pos_str : str
        string that has the region of interest defined by chromo:start_position-end_position. An example is 10:1234-1234.

    Returns
    -------
    str
        returns the gene_pos_str

    Raises
    ------
    ValueError
        raises a value error if the user inputs a incorrectly formatted string
    """

    if len(re.split(":|-", gene_pos_str)) != 3:
        raise ValueError(
            f"Expected the gene position string to be formatted like chromosome:start_position-end_position. Instead it was formatted as {gene_pos_str}"
        )
    else:
        return gene_pos_str


def check_json_path(json_path: str) -> str:
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

        program_dir: str = "/".join(os.path.realpath(__file__).split("/")[:-3])

        return "/".join([program_dir, "config.json"])


def check_env_path(env_path: str) -> str:
    """Callback function that creates the .env path string. If the user provides a value then it just returns that otherwise it creates the path to the default .env file.

    Parameters
    ----------
    env_path : str
        path to the .env file or an empty string

    Returns
    -------
    str
        returns the string to the file
    """
    if not env_path:

        filepath: str = "/".join(os.path.realpath(__file__).split("/")[:-3] + [".env"])

        return filepath

    else:
        return env_path


def display_version(value: bool):
    """callback function that displays the version number of the program and then terminates the program"""
    if value:
        # typer.echo(f"{__file__}")

        toml_filepath = "/".join(
            os.path.realpath(__file__).split("/")[:-3] + ["pyproject.toml"]
        )

        version = toml.load(toml_filepath)["tool"]["poetry"]["version"]

        typer.echo(f"IBDCluster - v{version}")
        raise typer.Exit(code=1)
