#!/usr/bin/env python
"""
This module is the main script for the IBDCluster program.
It contains the main cli and records inputs and creates
the typer app.
"""
import os
import pathlib
import shutil
import sys
from datetime import datetime
from enum import Enum
from typing import Dict, Optional

import analysis
import callbacks
import cluster
import log
import pandas as pd
import typer
from dotenv import load_dotenv
from models import DataHolder


class IbdProgram(str, Enum):
    """Enum used to define the options for the ibd_program flag in the cli"""

    HAPIBD = "hapibd"
    ILASH = "ilash"


class LogLevel(str, Enum):
    """Enum used to define the options for the log level in the cli"""

    WARNING = "warning"
    VERBOSE = "verbose"
    DEBUG = "debug"


app = typer.Typer(
    add_completion=False,
)


def load_phecode_descriptions(phecode_desc_file: str) -> Dict[str, Dict[str, str]]:
    """Function that will load the phecode_description file
    and then turn that into a dictionary

    Parameters
    ----------
    phecode_desc_file : str
        descriptions of each phecode

    Returns
    -------
    Dict[str, Dict[str, str]]
        returns a dictionary where the first key is the
        phecode and value is a dictionary where the inner key
        is 'phenotype' and the value is the descriptions
    """

    desc_df = pd.read_csv(phecode_desc_file, sep="\t", usecols=["phecode", "phenotype"])

    # making sure that the phecode keys are a string
    desc_df.phecode = desc_df.phecode.astype(str)

    # converting the dataframe into a dictionar
    return desc_df.set_index("phecode").T.to_dict()


@app.command()
def main(
    ibd_program: IbdProgram = typer.Option(
        IbdProgram.HAPIBD.value,
        "--ibd",
        "-i",
        help="IBD detection software that the output came from. The program expects these values to be either hapibd or ilash. The program also expects these values to be lowercase",
        case_sensitive=True,
    ),
    output: str = typer.Option(
        "./", "--output", "-o", help="directory to write the output files into."
    ),
    ibd_file: str = typer.Option(
        ...,
        "--ibd-file",
        "-f",
        help="path to either the hap-IBD or iLASH file that have the pairwise IBD sharing for each chromosome. This file should correspond to the chromosomes that are in the gene_info_file",
    ),
    env_path: str = typer.Option(
        "",
        "--env",
        "-e",
        help="path to a .env file that has variables for the hapibd files directory and the ilash file directory. These variables are called HAPIBD_PATH and ILASH_PATH respectively.",
        callback=callbacks.check_env_path,
    ),
    json_path: str = typer.Option(
        "",
        "--json-config",
        "-j",
        help="path to the json config file",
        callback=callbacks.check_json_path,
    ),
    gene_position: str = typer.Option(
        ...,
        "--gene-position",
        help="This will be the chromosome position in basepairs for the region. Example would be 10:12341234-12341234. The chromsome number comes first and then the start and end position of the region of interest. The chromosome position does not need to be prefixed with chr.",
        callback=callbacks.check_gene_pos_str,
    ),
    gene_name: str = typer.Option(
        "test", "--gene-name", "-n", help="name of the gene or region of interest"
    ),
    carriers: str = typer.Option(
        ...,
        "--carriers",
        "-c",
        help="Filepath to a text file that has the carrier status of all the individuals in the ibd dataset. The first column of these file should be a list of GRID ids and is expected to be called grids. If an individual has the phenotype they should be listed as a 1 otherwise they should be listed as 0.",
    ),
    cm_threshold: int = typer.Option(
        3,
        "--cM",
        help="Centimorgan threshold to filter the ibd segments",
    ),
    connection_threshold: int = typer.Option(
        0,
        "--connections",
        help="Threshold to filter out individuals who have fewer connections than the threshold. This parameter can be used if you are noticing that there are big networks of individuals (like 1000s of people connected)",
    ),
    steps: int = typer.Option(
        3, "-s", "--step", help="Number of steps used in the random walk"
    ),
    random_walk: bool = typer.Option(
        False,
        "--random-walk",
        help="Flag indicating that the user wishes to use a random walk instead of the grid search algorithm.",
        is_flag=True,
    ),
    max_network_size: int = typer.Option(
        30,
        "--max-network-size",
        help="Maximum number of individuals allowed in network before reclustering is performed. This argument should only be used when the user selects a random walk.",
    ),
    min_connectiveness: float = typer.Option(
        0.5,
        "--min-connectiveness",
        help="minimum threshold for connectiveness within the networks. This argument should only be used when the user selects a random walk.",
    ),
    phecode_descriptions: Optional[str] = typer.Option(
        None,
        "-d",
        "--phecode-desc",
        help="File that has the descriptions for each phecode. Expects two columns: 'phecode' and 'phenotype', that are tab separated.",
    ),
    sliding_window: bool = typer.Option(
        False,
        "--sliding-window",
        help="Optional flag that allows the user to run a 1MB sliding window along the region of interest if they want to. this method is recommended for large loci of interest",
        is_flag=True,
    ),
    loglevel: LogLevel = typer.Option(
        LogLevel.WARNING.value,
        "--loglevel",
        "-l",
        help="This argument sets the logging level for the program. Accepts values 'debug', 'warning', and 'verbose'.",
        case_sensitive=True,
    ),
    log_to_console: bool = typer.Option(
        False,
        "--log-to-console",
        help="Optional flag to log to only the console or also a file",
        is_flag=True,
    ),
    log_filename: str = typer.Option(
        "IBDCluster.log", "--log-filename", help="Name for the log output file."
    ),
    debug_iterations: int = typer.Option(
        3,
        "--debug-iterations",
        help="This argument will specify how many iterations the program should go through during the clustering step before it moves on. This argument should only be used if the loglevel is set to debug. If you wish to run in debug mode for a whole data set then set this argument to a high number. This practice is not recommended because the log file will get extremely large (Potentially TB's), so use with caution.",
    ),
    version: bool = typer.Option(  # pylint: disable=unused-argument
        False,
        "--version",
        help="version number of the IBDCluster program",
        callback=callbacks.display_version,
        is_eager=True,
        is_flag=True,
    ),
) -> None:
    """C.L.I. tool to identify networks of individuals at a biobank scale who share IBD segments overlapping a locus of interest and identify enrichment of phenotypes within each network"""
    # getting the programs start time
    start_time = datetime.now()

    # we are going to have the program append the $PLUGIN path to the interpreters search path
    sys.path.append(os.getenv("IBDCLUSTER_MAIN_PLUGINS"))

    if os.getenv("IBDCLUSTER_CUSTOM_PLUGINS"):
        sys.path.append(os.getenv("IBDCLUSTER_CUSTOM_PLUGINS"))

    # Now we can recreate the directory that the IBDCluster.log will be in
    pathlib.Path(output).mkdir(parents=True, exist_ok=True)

    # setting a path to the json file
    os.environ.setdefault("json_path", json_path)

    # loading the .env file
    load_dotenv(env_path)

    # creating the logger and then configuring it
    logger = log.create_logger()

    log.configure(
        logger,
        output,
        filename=log_filename,
        loglevel=loglevel,
        to_console=log_to_console,
    )

    # recording all the user inputs
    log.record_inputs(
        logger,
        ibd_program_used=ibd_program.value,
        ibd_filepath=ibd_file,
        output_path=output,
        environment_file=env_path,
        json_file=json_path,
        region_of_interest=gene_position,
        gene_name=gene_name,
        carrier_matrix=carriers,
        centimorgan_threshold=cm_threshold,
        connections_threshold=connection_threshold,
        sliding_window_enabled=sliding_window,
        loglevel=loglevel,
        log_filename=log_filename,
        random_walk_step_size=steps,
        random_walk_enabled=random_walk,
    )

    # adding the loglevel to the environment so that we can access it
    os.environ.setdefault("program_loglevel", str(log.get_loglevel(loglevel)))

    # adding the debug_iterations to the environment so that we can access it
    os.environ.setdefault("debug_iterations", str(debug_iterations))

    # need to first determine list of carriers for each phenotype
    carriers_df: pd.DataFrame = pd.read_csv(carriers, sep="\t")
    # forming a dictionary where the keys are phecodes and the
    # values are a list of indices in the carriers df that carry
    # the phecode
    carriers_dict = cluster.generate_carrier_dict(carriers_df)

    # loading the genes information into a generator object
    # pylint: disable-next="assignment-from-no-return"
    regions_of_interest = cluster.load_gene_info(
        gene_position, gene_name, sliding_window
    )

    # This section will handle preparing the phenocode
    # descriptions and the phenotype prevalances

    # loading in the phecode_descriptions
    if phecode_descriptions:
        phecode_desc = load_phecode_descriptions(phecode_descriptions)
    else:
        phecode_desc = {}

    # Now we will find the phecode percentages which will be used later
    phenotype_prevalances = cluster.get_phenotype_prevalances(
        carriers_dict, carriers_df.shape[0]
    )

    # We can then determine the different clusters for each gene. The genes_generator will always be an iterable so we can ignore that error
    for gene in regions_of_interest:

        network_generator = cluster.find_clusters(
            ibd_program.value, gene, cm_threshold, ibd_file, connection_threshold
        )
        # creating a specific output path that has the gene name
        gene_output = os.path.join(output, gene.name)

        # deleting the output directory incase there was already a
        # output there
        shutil.rmtree(gene_output, ignore_errors=True)
        # writing log messages for the networks and the allpairs.txt files
        logger.debug(
            "Information written to a networks.txt at: %s.",
            os.path.join(
                gene_output, "".join([ibd_program, "_", gene.name, "_networks.txt"])
            ),
        )

        logger.info(
            "Writing the allpairs.txt file to: %s",
            os.path.join(gene_output, "".join(["IBD_", gene.name, "_allpairs.txt"])),
        )

        # creating an object that holds useful information
        data_container = DataHolder(
            gene.name,
            gene.chr,
            carriers_dict,
            phenotype_prevalances,
            list(carriers_dict.keys()),
            ibd_program,
            phecode_desc,
        )
        # adding the networks, the carriers_df, the carriers_dict, and
        # the phenotype columns to a object that will be used in the analysis
        for network in network_generator:

            # This is the main function that will run the analysis of the networks
            analysis.analyze(data_container, network, gene_output)

    logger.info("analysis_finished")
    logger.info("Program Duration: %s", datetime.now() - start_time)


if __name__ == "__main__":
    app()
