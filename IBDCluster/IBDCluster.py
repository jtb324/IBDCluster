#!/usr/bin/env python
"""
This module is the main script for the IBDCluster program. It contains the main cli and records inputs and creates the typer app.
"""
import os
from enum import Enum
from typing import Dict, Optional
import pathlib
from dotenv import load_dotenv
import pandas as pd
import typer
import analysis
import callbacks
import shutil
import cluster
import log
from datetime import datetime
from models import DataHolder


class IbdProgram(str, Enum):
    hapibd = "hapibd"
    ilash = "ilash"


class LogLevel(str, Enum):
    warning = "warning"
    verbose = "verbose"
    debug = "debug"


app = typer.Typer(
    add_completion=False,
)


def load_phecode_descriptions(phecode_desc_file: str) -> Dict[str, Dict[str, str]]:
    """Function that will load the phecode_description file and then turn that into a dictionary

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
        IbdProgram.hapibd.value,
        "--ibd",
        "-i",
        help="IBD detection software that the output came from. The program expects these values to be hapibd or ilash. The program also expects these values to be lowercase",
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
    gene_info_file: str = typer.Option(
        ...,
        "--gene-file",
        "-g",
        help="Filepath to a text file that has information about the genes it should have four columns: Gene name, chromosome, gene start, and gene end (In this order). The file is expected to not have a header. The gene position should also correspond to the build used for the ibd data (GrCH37, etc...)",
        callback=callbacks.check_gene_file,
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
    phecode_descriptions: Optional[str] = typer.Option(
        None,
        "-d",
        "--phecode-desc",
        help="File that has the descriptions for each phecode. Expects two columns: 'phecode' and 'phenotype', that are tab separated.",
    ),
    loglevel: LogLevel = typer.Option(
        LogLevel.warning.value,
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
    debug_iterations: int = typer.Option(
        3,
        "--debug-iterations",
        help="This argument will specify how many iterations the program should go through durign the clustering step before it moves on. This argument should only be used if the loglevel is set to debug. If you wish to run in debug mode for a whole data set then set this argument to a high number. This practice is not recommended because the log file will get extremely large (Potentially TB's).",
    ),
    version: bool = typer.Option(
        False,
        "--version",
        help="version number of the IBDCluster program",
        callback=callbacks.display_version,
        is_eager=True,
        is_flag=True,
    ),
) -> None:
    """C.L.I. tool to identify networks of individuals who share IBD segments overlapping a locus of interest and identify enrichment of phenotypes within biobanks"""
    # getting the programs start time
    start_time = datetime.now()

    # Now we can recreate the directory that the IBDCluster.log will be in
    pathlib.Path(output).mkdir(parents=True, exist_ok=True)

    # setting a path to the json file
    os.environ.setdefault("json_path", json_path)

    # loading the .env file
    load_dotenv(env_path)

    # creating the logger and then configuring it
    logger = log.create_logger()

    log.configure(logger, output, loglevel=loglevel, to_console=log_to_console)

    # recording all the user inputs
    log.record_inputs(
        logger,
        ibd_program_used=ibd_program.value,
        ibd_filepath=ibd_file,
        output_path=output,
        environment_file=env_path,
        json_file=json_path,
        gene_info_file=gene_info_file,
        carrier_matrix=carriers,
        centimorgan_threshold=cm_threshold,
        loglevel=loglevel,
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
    genes_generator = cluster.load_gene_info(gene_info_file)

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

    # We can then determine the different clusters for each gene
    for gene in genes_generator:

        network_generator = cluster.find_clusters(
            ibd_program.value, gene, cm_threshold, ibd_file
        )
        # creating a specific output path that has the gene name
        gene_output = os.path.join(output, gene.name)

        # deleting the output directory incase there was already a
        # output there
        shutil.rmtree(gene_output, ignore_errors=True)
        # writing log messages for the networks and the allpairs.txt files
        logger.debug(
            f"Information written to a networks.txt at: {os.path.join( gene_output, ''.join([ibd_program, '_', gene.name, '_networks.txt']))}"
        )

        logger.info(
            f"Writing the allpairs.txt file to: {os.path.join(gene_output, ''.join(['IBD_', gene.name, '_allpairs.txt']))}"
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
    logger.info(f"Program Duration: {datetime.now() - start_time}")


if __name__ == "__main__":
    app()
