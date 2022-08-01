from typing import Generator
from collections import namedtuple
import models
import log
import pandas as pd
from tqdm import tqdm
import os

Genes = namedtuple("Genes", ["name", "chr", "start", "end"])
# getting the logger object
logger = log.get_logger(__name__)


def load_gene_info(filepath: str) -> Generator[Genes, None, None]:
    """Function that will load in the information for each gene. This function will return a generator

    filepath : str
        filepath to a file that has the information for the genes of interest

    Returns

    Generator
        returns a generator of namedtuples that has the gene information
    """

    with open(filepath, "r", encoding="utf-8") as gene_input:
        logger.debug(
            f"Loaded in the gene information from the file, {filepath}, into a generator"
        )
        for line in gene_input:
            split_line: list[str] = line.split()

            gene_tuple: Genes = Genes(*split_line)

            yield gene_tuple


def _identify_carriers_indx(
    carriers_series: pd.Series, return_dict: dict[str, list[str]]
) -> None:
    """Function that will insert the key, value pair of the phenotype and the carriers into the dictionary

    Parameters
    ----------
    carriers_series : pd.Series
        row of dataframe that is used in the apply function in
        the generate_carrier_list function

    return_dict : dict[str, list[str]]
        dictionary where the phecodes will be the keys and the values will be a list of indices indicating which grids are carriers
    """
    return_dict[carriers_series.name] = list(
        carriers_series[carriers_series == 1].values.tolist()
    )


def generate_carrier_dict(carriers_matrix: pd.DataFrame) -> dict[str, list[str]]:
    """Function that will take the carriers_pheno_matrix and generate a dictionary that has the list of indices for each carrier

    Parameters
    ----------
    carriers_matrix : pd.DataFrame
        dataframe where the columns are the phecodes and have 0's or 1's for whether or not they have the phecodes

    Returns
    -------
    dict[str, list[str]]
        dictionary where the keys are phecodes and the values are list of integers
    """
    return_dict = {}

    # iterating over each phenotype which starts with the
    # second column
    carriers_matrix.T.apply(lambda x: _identify_carriers_indx(x, return_dict), axis=1)

    return return_dict


def find_clusters(
    ibd_program: str,
    gene: Genes,
    cm_threshold: int,
    ibd_file: str,
) -> Generator[models.Network, None, None]:
    """Main function that will handle the clustering into networks

    Parameters
    ----------
    ibd_program : str
        ibd program used to identify pairwise ibd. This value should either be hap-IBD or iLASH

    gene : Gene
        This is the Gene namedtuple that has things like the name, chromosome, start, and
        end point

    cm_threshold : int
        This is the minimum threshold to filter the ibd segments to

    carriers : dict[float, list[str]]
        dictionary that has the phecodes as keys and list of grid IDs that are affected by the phecode as values

    ibd_file : str
        filepath to the ibd file from hap-IBD or iLASH

    Returns
    -------
    Generator[models.Network]
        returns a generator of the models.Network. A generator is used to avoid massive memory consumption
    """

    # Next two lines create an object with the shared indices for each
    # ibd program. Then it loads the proper unique indices for the correct
    # program
    indices = models.FileInfo()

    # adding the correct index for the ibd program
    logger.debug(f"adding the cM index based on the ibd program: {ibd_program}")

    indices.set_cM_indx(ibd_program)

    logger.debug(f"gene named tuple: {gene}")
    logger.info(f"finding clusters for the gene: {gene.name}")

    cluster_model: models.Cluster = models.Cluster(ibd_file, ibd_program, indices)

    # loading in all the dataframe for the genetic locus
    cluster_model.load_file(gene.start, gene.end, cm_threshold)

    logger.info(f"finished loading in the file for the dataframe")
    # # filtering the dataframe to >= specific centimorgan threshold
    # cluster_model.filter_cm_threshold(cm_threshold, indices.cM_indx)

    all_grids: list[str] = cluster_model.find_all_grids(indices)

    for ind in tqdm(all_grids, desc="pairs in clusters: "):

        network_obj = models.Network(gene.name, gene.chr, cluster_model.network_id)
        # function will return error string if the individual is already in a network
        err = cluster_model.construct_network(ind, network_obj)

        # if the program is being run in debug mode then it
        # will only this loop how ever many times the user
        # wants. This gives enough information to see if
        # the program is behaving properly
        if int(os.environ.get("program_loglevel")) == 10:
            if cluster_model.network_id == int(os.environ.get("debug_iterations")):
                break

        if err:
            continue

        yield network_obj
