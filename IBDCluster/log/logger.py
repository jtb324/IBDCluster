import logging
import enum
from typing import Dict

level_dict = {
    "info": logging.INFO,
    "debug": logging.INFO
}

def create_logger(loglevel: str) -> logging.Logger:
    """function that will get the correct logger for the program
    
    Parameters
    
    loglevel : str
        logging level that the user wants to use. The default level is INFO
    
    Returns

    logging.Logger
    """
    logger = logging.getLogger("__main__")

    return logger

