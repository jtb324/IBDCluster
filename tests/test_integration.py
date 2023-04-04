from typer.testing import CliRunner
import pytest
import sys

sys.path.append("./drive")

from drive.drive import app

runner = CliRunner()

# @pytest.mark.integtest
# def test_drive_full_run():
#     assert 1==1


@pytest.mark.integtest
def test_drive_full_run():

    result = runner.invoke(
        app,
        [
            "-i",
            "./tests/test_inputs/biovu_longQT_EUR_chr21.ibd.gz",
            "-f",
            "hapIBD",
            "-t",
            "22:35818986-35884508",
            "-o",
            "./tests/test_output/integration_test_results",
            "-m",
            "3",
        ],
    )
    # result = runner.invoke(
    #     app,
    #     [
    #         "--help"
    #     ],
    #     color=True
    # )
    # print(result)
    assert result.exit_code == 0
