"""
Module that is responsible for running the clustering and 
initializing the network and cluster models
"""
import os
import re
from collections import namedtuple
from typing import Generator
from itertools import chain

import log
import models
import pandas as pd
from tqdm import tqdm

Genes = namedtuple("Genes", ["name", "chr", "start", "end"])
# getting the logger object
logger = log.get_logger(__name__)


def load_gene_info(
    chromo_pos_str: str, gene_name: str, sliding_window: bool
) -> list[Genes]:
    """Function that will load in the information for each gene. The function
    can also can create a list of sliding windows for the target region. This
    will be done for 1 MB at a time.

    Parameters
    ----------
    chromo_pos_str : str
        string that has the region of interest in bbase pairs. This string
        will look like 10:1234-1234 where the first number is the chromosome
        number, then the start position, and then the end position of the
        region of interest.

    gene_name : str
        name of the region of interest

    sliding_window : bool
        boolean value indicating if the user wishes to use a sliding window
        approach for clustering. If they do then it will break everything down
        into 1 MB regions

    Returns
    -------
    list[Genes]
        returns a list of Genes namedtuples that have information about the region
        of interest
    """
    gene_regions_list = []
    # we are going to try to split this string if the user
    chr_num, start_pos, end_pos = re.split(":|-", chromo_pos_str)
    if sliding_window:
        logger.info(
            f"creating a sliding window every 1MB from {start_pos} to {end_pos}"
        )
        # We cain first creawte an
        window_iterator = chain(range(start_pos, end_pos, 1_000_000), [end_pos])
        # we are going to pull the first value out of the iterator
        first_pos = next(window_iterator)

        for pos in window_iterator:
            # we need to build a name for each window region
            name = f"{gene_name}_{first_pos}-{pos}"
            region_of_interest = Genes(name, chr_num, start_pos, end_pos)
            logger.debug(
                f"Identified the following region of interest: {region_of_interest}"
            )
            first_pos = pos
            gene_regions_list.append(region_of_interest)
    else:

        region_of_interest = Genes(gene_name, chr_num, start_pos, end_pos)
        logger.debug(
            f"Identified the following region of interest: {region_of_interest}"
        )
        gene_regions_list.append(region_of_interest)

    return gene_regions_list


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

    return_dict[carriers_series.name] = carriers_series[
        carriers_series == 1
    ].index.tolist()


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

    new_indx_df = carriers_matrix.set_index("grids")

    # iterating over each phenotype which starts with the
    # second column
    new_indx_df.apply(lambda x: _identify_carriers_indx(x, return_dict), axis=0)

    return return_dict


def find_clusters(
    ibd_program: str,
    gene: Genes,
    cm_threshold: int,
    ibd_file: str,
    connections_threshold: int,
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

    ibd_file : str
        filepath to the ibd file from hap-IBD or iLASH

    connections_threshold : int
        integer value of the minimum number of connections an individual needs to have to be used
        in the analysis. Individuals who don't meet this threshold will not be considered for the clustering

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

    logger.info(f"finding clusters for the gene: {gene.name}")

    cluster_model: models.Cluster = models.Cluster(ibd_file, ibd_program, indices)

    # loading in all the dataframe for the genetic locus
    cluster_model.load_file(gene.start, gene.end, cm_threshold)

    logger.info("finished loading in the file for the dataframe")
    # # filtering the dataframe to >= specific centimorgan threshold

    cluster_model.create_unique_ids(indices)

    # filtering for individuals who have fewer connections than we want
    cluster_model.filter_connections(connections_threshold)

    all_grids = cluster_model.get_ids()

    for ind in tqdm(all_grids, desc="pairs in clusters: "):

        network_obj = models.Network(gene.name, gene.chr, cluster_model.network_id)
        # function will return error string if the individual is already in a network
        err = cluster_model.construct_network(ind, network_obj)

        if len(network_obj.haplotypes) == 0:
            logger.info(
                f"There were no haplotypes found for network {network_obj.network_id} when individual, {ind}, were used as a seed"
            )
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
