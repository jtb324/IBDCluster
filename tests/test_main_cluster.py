import sys
from collections import namedtuple
import logging
import pytest

sys.path.append("../IBDCluster")
from cluster.main import load_gene_info  # pylint: disable="wrong-import-position"


class TestLociLoader:
    @pytest.fixture
    def loci_info(self, request):
        """Fixture that will generate the gene info object"""
        return load_gene_info(request.param[0], request.param[1])

    @pytest.mark.parametrize(
        "loci_info", [("./test_data/gene_info.txt", False)], indirect=True
    )
    def test_load_loci_info_no_sliding_window(self, loci_info) -> None:
        """Function to make sure that the return object is of type Genes"""

        gene_generator = loci_info

        for gene_tuple in gene_generator:
            assert all(
                field in gene_tuple._asdict().keys()
                for field in ["name", "chr", "start", "end"]
            ), "expected the gene_tuple to have the fields 'name', 'chr', 'start', 'end'"

    @pytest.mark.parametrize(
        "loci_info", [("./test_data/gene_info.txt", True)], indirect=True
    )
    def test_values_of_sliding_windows_formed(self, loci_info) -> None:
        """Method to check how many steps are made in the sliding window list"""
        gene_generator = loci_info

        first_window = next(gene_generator)

        errors = []
        if first_window.name != "TEST_GENE_2769662-2770662":
            errors.append(
                f"Expected the first window to have the name: TEST_GENE_2769662-2770662. Instead it had the name: {first_window.name}"
            )

        if first_window.chr != 2:
            errors.append(
                f"Expected the first window to have the chromosome number: 2. Instead it had the chromosome number: {first_window.chr}"
            )

        if first_window.start != 2769662:
            errors.append(
                f"Expected the first window to have the start position: 2769662. Instead it had the start position: {first_window.start}"
            )

        if first_window.end != 2770662:
            errors.append(
                f"Expected the first window to have the end position: 2770662. Instead it had the end position: {first_window.end}"
            )

        assert not errors, "errors occured:\n{}".format("\n".join(errors))

    @pytest.mark.parametrize(
        "loci_info", [("./test_data/gene_info.txt", True)], indirect=True
    )
    def test_number_of_sliding_windows_formed(self, loci_info) -> None:
        gene_generator = loci_info

        loci_of_interest_list = list(gene_generator)
        assert (
            len(loci_of_interest_list) == 397
        ), f"Expected 397 sliding window steps to be formed. Instead {len(loci_of_interest_list)} were formed"


# def test_loading_loci_info_sliding_window() -> None:
#     """Unit test that makes sure that the program can properly create the sliding window. It calls the cluster.main.load_gene_info function"""
