import sys
from collections import namedtuple
import logging

sys.path.append("../IBDCluster")
from cluster.main import load_gene_info


def test_loading_gene_info_return() -> None:
    """Function to make sure that the return object is of type Genes"""
    filepath = "./test_data/HK2_gene_info_1_31_22.txt"

    gene_generator = load_gene_info(filepath)

    gene_tuple = list(gene_generator)[0]

    print(gene_tuple._asdict())
    assert all(
        field in gene_tuple._asdict().keys()
        for field in ["name", "chr", "start", "end"]
    ), "expected the gene_tuple to have the fields 'name', 'chr', 'start', 'end'"
