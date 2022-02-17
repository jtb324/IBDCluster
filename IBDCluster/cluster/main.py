from typing import Dict, Generator, List, Tuple
from collections import namedtuple
import os
import models
import log

logger = log.get_logger(__name__)

def load_gene_info(filepath: str) -> Generator:
    """Function that will load in the information for each gene. This function will return a generator

    filepath : str
        filepath to a file that has the information for the genes of interest

    Returns

    Generator
        returns a generator of namedtuples that has the gene information
    """

    Genes = namedtuple("Genes", ["name", "chr", "start", "end"])

    with open(filepath, "r", encoding="utf-8") as gene_input:
        logger.info(f"Loaded in the gene information from the file, {filepath}, into a generator")
        for line in gene_input:
            split_line: List[str] = line.split()

            gene_tuple: Genes = Genes(*split_line)

            yield gene_tuple


def get_ibd_program(ibd_program: str):
    """Function that will take whatever ibd program the user pass and will get the correct info object for it

    Parameters

    ibd_program : str
        string of either 'hapibd' or 'ilash'

    Returns

    Object
        returns either a clusters.Hapibd_Info object or clusters.Ilash_Info object
    """
    logger.info(f"Building the model that has the necessary indices for the {ibd_program} output")

    if ibd_program == "hapibd":
        return models.Hapibd_Info(file_dir=os.environ.get("hapibd_files"))
    
    return models.Ilash_Info(file_dir=os.environ.get("ilash_files"))


def find_clusters(
    ibd_program: str, gene_info_filepath: str, cM_threshold: int
) -> Dict[Tuple[str, int], Dict]:
    """Main function that will handle the clustering into networks"""

    # create a dictionary that will have the gene name and chromosome as keys and the network information as values
    return_dict: Dict[Tuple[str, int], Dict] = {}

    # we will need the information for the correct ibd_program
    indices: models.File_Info = get_ibd_program(ibd_program)

    # gather all the ibd files into an attribute of the indice class called self.ibd_files
    logger.debug(f"gathering all the necessary files for the program: {ibd_program}")
    indices.gather_files()

    # creating a generator that returns the Genes namedtuple from the load_gene_info function
    gene_generator: Generator = load_gene_info(gene_info_filepath)

    # This for loop will iterate through each gene tuple in the generator. It then uses the
    # chromosome number to find the correct file. A Cluster object is then created which loads the pairs
    # that surround a certain location into a dataframe
    for gene_tuple in gene_generator:
        logger.debug(f"finding clusters for the gene: {gene_tuple.name}")
        hapibd_file: str = indices.find_file("".join(["chr", gene_tuple.chr]))

        cluster_model: models.Cluster = models.Cluster(gene_tuple.name, hapibd_file)

        cluster_model.load_file(
            gene_tuple.start, gene_tuple.end, indices.str_indx, indices.end_indx
        )
        # filtering the dataframe to >= specific centimorgan threshold
        cluster_model.filter_cM_threshold(cM_threshold, indices.cM_indx)
        # This line will filter the dataframe for only those pairs that have the same phase
        # cluster_model.filter_for_haplotype(indices.id1_phase_indx, indices.id2_phase_indx)

        network_info: Dict = cluster_model.find_networks(indices)

        return_dict[(gene_tuple.name, gene_tuple.chr)] = network_info

    return return_dict
