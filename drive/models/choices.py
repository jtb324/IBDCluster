from enum import Enum


class LogLevel(str, Enum):
    """Enum used to define the options for the log level in the cli"""

    INFO = 0
    DEBUG = 1


class FormatTypes(str, Enum):
    HAPIBD = "hapibd"
    ILASH = "ilash"
    GERMLINE = "germline"
    RAPID = "rapid"


class OverlapOptions(str, Enum):
    """Enum defining options for how the use wants to define
    overlap with the gene region of interest. Values can be contains or overlaps"""

    CONTAINS = "contains"
    OVERLAPS = "overlaps"
