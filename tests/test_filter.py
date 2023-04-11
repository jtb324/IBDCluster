from collections import namedtuple
from pathlib import Path
import pandas as pd
import pytest
import sys

sys.path.append("./drive")

from drive.filters import IbdFilter
from drive.models.generate_indices import HapIBD

hapibd = HapIBD()

Genes = namedtuple("Genes", ["chr", "start", "end"])


@pytest.mark.unit
def test_ibd_file_not_found() -> None:
    """Unit test that will make sure a FileNotFound Exception is raised if the ibd file doesn't exist"""
    fake_gene = Genes(12, 1235, 1235)
    with pytest.raises(FileNotFoundError):
        IbdFilter.load_file(Path("FileDoesNotExist.txt"), hapibd, fake_gene)


@pytest.mark.kcne1
def test_filter_results() -> None:
    """Unit test that will make sure that the Filter class produces the proper results"""

    true_ibdpd = pd.read_csv("tests/test_true_results/ibdpd_KNCE1.txt", sep="\t")

    kcne1_gene = Genes(21, 35818986, 35884508)

    filter_obj = IbdFilter.load_file(
        Path("tests/test_inputs/biovu_longQT_EUR_chr21.ibd.gz"), hapibd, kcne1_gene
    )

    kcne1_min_cm_thres = 3

    filter_obj.preprocess(kcne1_min_cm_thres)

    error_list = []

    # Make sure that the code identifies the correct number of IBD segments as before
    if filter_obj.ibd_pd.shape[0] != true_ibdpd.shape[0]:
        error_list.append(
            f"The filtered dataframe from the filter_obj is not the same length as the file test ibdpd_KCNE1.txt. Expected the file to have {true_ibdpd.shape[0]} instead there were {filter_obj.ibd_pd.shape[0]}"
        )

    # Make sure that the program identifies the same number of unique haplotypes as before
    max_id_new = filter_obj.ibd_pd[["idnum1", "idnum2"]].values.max()
    max_id_old = true_ibdpd[["idnum1", "idnum2"]].values.max()

    if max_id_new != max_id_old:
        error_list.append(
            f"The maximum haplotype id assigned in the new version is not the same as before. Expected value to be {max_id_old}, instead it was {max_id_new}"
        )

    assert not error_list, "errors occured:\n{}".format("\n".join(error_list))
