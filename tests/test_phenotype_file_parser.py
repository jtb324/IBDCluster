import pytest
import sys

sys.path.append("./drive")

from drive.utilities.parser import PhenotypeFileParser


@pytest.mark.unit
def test_file_not_found() -> None:
    """Test that the parsers raises a FileNotFoundError if the phenotype file doesn't exist."""
    with pytest.raises(FileNotFoundError):
        with PhenotypeFileParser("./NotRealFile.txt") as pheno_file:
            ...


@pytest.mark.unit
@pytest.mark.parametrize(
    "test_str,expected", [("ID1\t1", "\t"), ("ID1,1", ","), ("ID1|0", "|")]
)
def test_check_separator(test_str: str, expected: str) -> None:
    """Test that the method _check_separator to see if it returns proper separator."""

    assert PhenotypeFileParser._check_separator(test_str) == expected


@pytest.mark.unit
@pytest.mark.parametrize("test_str", [("ID1\n1")])
def test_check_separator_error(test_str: str) -> None:
    """Test that the method _check_separator raises a ValueError if it cannot identify the separator."""
    with pytest.raises(ValueError):
        PhenotypeFileParser._check_separator(test_str)


@pytest.mark.unit
def test_case_control_exclusion_counts() -> None:
    """Check that the PhenotypeFileParser identifies the correct case/control/exlcusion counts."""
    with PhenotypeFileParser("./tests/test_inputs/test_phenotype_file.txt") as parser:
        cases, controls, exclusions = parser.parse_cases_and_controls()

        expected_case_counts = 8
        expected_controls_counts = 84
        expected_exclusion_counts = 8

        error_list = []

        if len(cases) != expected_case_counts:
            error_list.append(
                f"Expected parser to identify {expected_case_counts} cases. Instead {len(cases)} cases were identified."
            )
        if len(controls) != expected_controls_counts:
            error_list.append(
                f"Expected parser to identify {expected_controls_counts} cases. Instead {len(controls)} cases were identified"
            )
        if len(exclusions) != expected_exclusion_counts:
            error_list.append(
                f"Expected parser to identify {expected_exclusion_counts} cases. Instead {len(exclusions)} cases were identified"
            )

        assert not error_list, "errors occured:\n{}".format("\n".join(error_list))
