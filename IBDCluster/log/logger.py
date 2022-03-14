import logging
from typing import Dict
import os

level_dict: Dict[str, int] = {
    "verbose": logging.INFO,
    "debug": logging.DEBUG,
    "warning": logging.WARNING,
}


def get_loglevel(loglevel: str):
    """Function that will return a log level based on the input"""

    return level_dict[loglevel]


def configure(
    logger: logging.Logger,
    output: str,
    filename: str = "IBDCluster.log",
    loglevel: str = "warning",
    to_console: bool = False,
) -> None:
    """Function that will configure the level of logging"""

    filename = os.path.join(output, filename)

    logger.setLevel(level_dict.get(loglevel, logging.WARNING))

    file_formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    stream_formatter = logging.Formatter("%(message)s")
    # program defaults to log to a file called IBDCluster.log in the
    # output directory
    fh = logging.FileHandler(filename, mode="w")
    fh.setFormatter(file_formatter)
    logger.addHandler(fh)

    # If the user selects to also log to console then the program will
    # log information to the stderr
    if to_console:
        sh = logging.StreamHandler()
        sh.setFormatter(stream_formatter)
        logger.addHandler(sh)


def get_logger(module_name: str, main_name: str = "__main__") -> logging.Logger:
    """Function that will be responsible for getting the logger for modules"""
    return logging.getLogger(main_name).getChild(module_name)


def create_logger(
    logger_name: str = "__main__",
) -> logging.Logger:
    """function that will get the correct logger for the program

    Parameters

    loglevel : str
        logging level that the user wants to use. The default level is INFO

    Returns

    logging.Logger
    """

    logger = logging.getLogger(logger_name)

    return logger
