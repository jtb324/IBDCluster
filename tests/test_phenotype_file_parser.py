import pytest
import sys

sys.path.append("./drive")

from drive.utilities.parser import PhenotypeFileParser


@pytest.mark.unit
def test_file_not_found() -> None:
    """Unit test to make sure that the parser s raising a FileNotFoundError if the phenotype file doesn't exist"""
    ...
    with pytest.raises(FileNotFoundError):
        with PhenotypeFileParser("./NotRealFile.txt") as pheno_file:
            ...
