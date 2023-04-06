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
