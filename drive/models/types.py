# File will have different Interfaces used in the program as well as different class for types. This will simplify the imports and allow me to reduce coupling between modules

from collections import namedtuple
from dataclasses import dataclass, field
from typing import Protocol

from pandas import DataFrame

# namedtuple that will contain information about the gene being run
Genes = namedtuple("Genes", ["chr", "start", "end"])


# interface for the filter object
@dataclass
class Filter(Protocol):
    ibd_vs: DataFrame = field(default_factory=DataFrame)
    ibd_pd: DataFrame = field(default_factory=DataFrame)
