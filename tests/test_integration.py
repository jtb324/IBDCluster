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
def test_sucessful_run_with_desc():
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
            "verbose",
            "-j",
            "./config.json",
            "-d",
            "./phecode_descriptions.txt",
        ],
    )

    assert result.exit_code == 0


# integration test using the raw cftr_and_vw_test_phecode_matrix.txt
# affected statuses
@pytest.mark.integtest
def test_sucessful_run_no_desc():
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
            "verbose",
            "-j",
            "./config.json",
        ],
    )

    assert result.exit_code == 0


@pytest.mark.integtest
def test_percent_carriers_output():
    """This is going to try to test the percent carriers output file to make sure it is formatted properly and giving the right output"""
    # list to append the error messages into

    errors = []

    percentage_df = pd.read_csv("./percent_carriers_in_population.txt", sep="\t")

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
        print(percentage_df.percentage_in_population.values)
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
    """This test is going to check the output of the allpairs.txt file."""
    errors = []

    allpairs_df = pd.read_csv("./TEST_GENE/IBD_TEST_GENE_allpairs.txt", sep="\t")

    # first going to make sure there are the right number of columns
    if len(allpairs_df.columns) != 15:
        error = f"Expected the IBD_TEST_GENE_allpairs.txt file to have 15 columns. Instead there were {len(allpairs_df.columns)} columns"

        errors.append(error)

    # next checking to make sure the proper columns are identified
    col_list = [
        "program",
        "network_id",
        "pair_1",
        "pair_2",
        "phase_1",
        "phase_2",
        "chromosome",
        "gene_name",
        "286.11_Pair_1_status",
        "286.11_Pair_2_status",
        "499_Pair_1_status",
        "499_Pair_2_status",
        "start",
        "end",
        "length",
    ]

    dif_cols = [col for col in allpairs_df.columns if col not in col_list]

    if dif_cols:

        error = f"Expected the IBD_TEST_GENE_allpairs.txt file to have the columns {', '.join(col_list)}. Found different columns {', '.join(dif_cols)}"

        errors.append(error)

    # checking to make sure alll 17 networks were found
    if allpairs_df.network_id.max() != 17:

        error = f"Expected the allpairs.txt file to contain 17 unique network ids. Instead only {allpairs_df.network_id.max()} were found"

        errors.append(error)

    # checking the value of the network ids are the same for this person and that they were not assigned to multiple networks
    network_ids = list(
        set(
            allpairs_df[
                allpairs_df.phase_1 == "patient_id_100.2"
            ].network_id.values.tolist()
        )
    )

    if len(network_ids) != 1:
        error = f"expected the individual patient_id_100.2 to be in 1 network but instead the individual was found in networks {', '.join(network_ids)}"

        errors.append(error)

    if len(network_ids) == 1:

        network_id = network_ids[0]

        filtered_df = allpairs_df[allpairs_df.network_id == network_id]

        if filtered_df.shape[0] != 8:

            error = f"expected network {network_id} to have 8 connections. Instead found {filtered_df.shape[0]} connections"

            errors.append(error)

    # checking to make sure the 286.11_Pair_2_status is all zero like it should be
    if not all(
        [val == 0 for val in allpairs_df["286.11_Pair_2_status"].values.tolist()]
    ):

        error = f"expected all of the values in the 286.11_Pair_2_status column to be zero instead the values were {', '.join(allpairs_df['286.11_Pair_2_status'].value_counts().values())}"

        errors.append(error)

    # assert no error message has been registered, else print messages
    assert not errors, "errors occured:\n{}".format("\n".join(errors))
