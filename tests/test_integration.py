from typer.testing import CliRunner
import pytest
import sys
import pandas as pd

sys.path.append("../IBDCluster")
from IBDCluster import app

runner = CliRunner()

# integration test using the raw cftr_and_vw_test_phecode_matrix.txt
# affected statuses
@pytest.mark.integtest
def test_sucessful_run():
    result = runner.invoke(
        app,
        [
            "-o",
            "./",
            "-e",
            "./test.env",
            "-c",
            "./test_data/cftr_and_vw_test_phecode_matrix.txt",
            "-g",
            "./test_data/gene_info.txt",
            "--cM",
            "5",
            "-l",
            "debug",
        ],
    )

    assert result.exit_code == 0


@pytest.mark.integtest
def test_percent_carriers_output():
    """This is going to try to test the percent carriers output file to make sure it is formatted properly and giving the right output"""
    # list to append the error messages into

    errors = []

    percentage_df = pd.read_csv("percent_carriers_in_population.txt", sep="\t")

    # first we are going to make sure that the percentage_df has the
    # right column
    if not all(
        [
            col in ["phenotype", "percentage_in_population"]
            for col in percentage_df.columns
        ]
    ):
        error = f"output file percent_carriers_in_population.txt does not have the proper columns. Expected 'phenotype', 'percentage_in_population'. Instead found {', '.join(percentage_df.columns)}"

        errors.append(error)

    if any(percentage_df.percentage_in_population.values):
        error = f"expected the percentages in the percentage_df determined from cftr_and_vw_test_phecode_matrix.txt to both be 0. Instead these values were {', '.join(percentage_df.percentage_in_population.values)}"

        errors.append(error)

    if not all(
        [
            phenotype in [499.0, 286.11]
            for phenotype in percentage_df.phenotype.values.tolist()
        ]
    ):
        error = f"Expected the two phenotypes in the output from cftr_and_vw_test_phecode_matrix.txt to be 499.0, and 286.11. Instead found {', '.join([str(phenotype) for phenotype in percentage_df.phenotype.values.tolist()])}"

        errors.append(error)

    # assert no error message has been registered, else print messages
    assert not errors, "errors occured:\n{}".format("\n".join(errors))


@pytest.mark.integtest
def test_allpairs():
    pass
