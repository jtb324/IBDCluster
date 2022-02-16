import logging
import enum
from typing import Dict

level_dict = {
    "info": logging.INFO,
<<<<<<< HEAD
    "debug": logging.INFO
=======
    "debug": logging.DEBUG
>>>>>>> c26b34fc14f148bb1629760a1ba8eb857a2b94ee
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

