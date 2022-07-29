from pathlib import Path
import typer
import os

VERSION = "1.0.3"


class IncorrectFileType(Exception):
    """Error that will be raised if the gene info file is not a text file"""

    def __init__(self, suffix: str, message: str) -> None:
        self.suffix: str = suffix
        self.message: str = message
        super().__init__(message)


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
    if value:
        typer.echo(f"IBDCluster - v{VERSION}")
        raise typer.Exit(code=1)
