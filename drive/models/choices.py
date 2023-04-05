from enum import Enum


class LogLevel(str, Enum):
    """Enum used to define the options for the log level in the cli"""

    WARNING = "warning"
    VERBOSE = "verbose"
    DEBUG = "debug"


class FormatTypes(str, Enum):
    HAPIBD = "hapibd"
    ILASH = "ilash"
    GERMLINE = "germline"
    RAPID = "rapid"
