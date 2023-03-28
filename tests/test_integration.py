from typer.testing import CliRunner
import pytest
import sys

sys.path.append("../drive")

from DRIVE import app

runner = CliRunner()


@pytest.mark.integtest
def test_drive_full_run():
    result = runner.invoke(
        app,
        [
            "-i",
            "",
            "-f",
            "hapIBD",
            "-t",
            "22:35818986-35884508",
            "-o",
            "./test_output/integration_test_results",
            "-m",
            "3",
        ],
    )
    assert result.exit_code == 0
