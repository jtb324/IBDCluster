from typer.testing import CliRunner
import pytest
import sys

sys.path.append("../IBDCluster")
from IBDCluster import app

runner = CliRunner()


@pytest.mark.integtest
def test_integration():
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
