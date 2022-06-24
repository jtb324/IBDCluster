#!/usr/bin/env python
"""
This module is the main script for the IBDCluster program. It contains the main cli and records inputs and creates the typer app.
"""
import os
from typing import Dict, List, Optional
import pathlib
from dotenv import load_dotenv
import pandas as pd
import typer
import analysis
import callbacks
import cluster
import log
from datetime import datetime
from models import DataHolder, Network


app = typer.Typer(
    add_completion=False,
    help="Tool that identifies ibd sharing for specific loci for individuals within biobanks",
)


def record_inputs(logger, **kwargs) -> None:
    """function to record the user arguments that were passed to the
    program. Takes a logger and then a dictionary of the user
    arguments"""

    logger.setLevel(20)

    for parameter, value in kwargs.items():
        logger.info(f"{parameter}: {value}")

    # getting the correct log level to reset the logger
    logger.setLevel(log.get_loglevel(kwargs["loglevel"]))


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
    ibd_program: str = typer.Option(
        "hapibd",
        "--ibd",
        "-i",
        help="IBD detection software that the output came from. The program expects these values to be hapibd or ilash",
        callback=callbacks.check_ibd_program,
    ),
    output: str = typer.Option(
        "./", "--output", "-o", help="directory to write the output files into."
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
        help="Filepath to a text file that has information about the genes it should have four columns: Gene name, chromosome, gene start, and gene end. The file is expected to not have a header",
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
    loglevel: str = typer.Option(
        "warning",
        "--loglevel",
        "-l",
        help="This argument sets the logging level for the program. Accepts values 'debug', 'warning', and 'verbose'.",
        callback=callbacks.check_loglevel,
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
        help="This argument will specify how many iterations the program should go through durign the clustering step before it moves on. This argument should only be used if the loglevel is set to debug. If you wish to run in debug for all of the program then set this argument to a high number. This practice is not recommend because the log file will get quite large. The default value is 3",
    ),
    version: bool = typer.Option(
        False,
        "--version",
        callback=callbacks.display_version,
        is_eager=True,
        is_flag=True,
    ),
) -> None:
    """Main function for the program that has all the parameters that the user can use with type"""

    # getting the programs start time
    start_time = datetime.now()

    # create the directory that the IBDCluster.log will be in
    pathlib.Path(output).mkdir(parents=True, exist_ok=True)

    # loading the
    os.environ.setdefault("json_path", json_path)

    # loading the .env file
    load_dotenv(env_path)

    # adding the loglevel to the environment so that we can access it
    os.environ.setdefault("program_loglevel", str(log.get_loglevel(loglevel)))

    # adding the debug_iterations to the environment so that we can access it
    os.environ.setdefault("debug_iterations", str(debug_iterations))

    # creating the logger and then configuring it
    logger = log.create_logger()

    log.configure(logger, output, loglevel=loglevel, to_console=log_to_console)

    # recording all the user inputs
    record_inputs(
        logger,
        ibd_program_used=ibd_program,
        output_path=output,
        environment_file=env_path,
        json_file=json_path,
        gene_info_file=gene_info_file,
        carrier_matrix=carriers,
        centimorgan_threshold=cm_threshold,
        loglevel=loglevel,
    )

    # need to first determine list of carriers for each phenotype
    carriers_df: pd.DataFrame = pd.read_csv(carriers, sep="\t")

    carriers_dict = cluster.generate_carrier_list(carriers_df)

    genes_generator = cluster.load_gene_info(gene_info_file)

    # loading in the phecode_descriptions
    if phecode_descriptions:
        phecode_desc = load_phecode_descriptions(phecode_descriptions)
    else:
        phecode_desc = {}

    # We can then determine the different clusters for each gene
    for gene in genes_generator:

        networks_list: List[Network] = cluster.find_clusters(
            ibd_program, gene, cm_threshold, carriers_dict
        )

        # adding the networks, the carriers_df, the carriers_dict, and the
        # phenotype columns to a object that will be used in the analysis
        data_container = DataHolder(
            gene.name,
            gene.chr,
            networks_list,
            carriers_dict,
            carriers_df,
            carriers_df.columns[1:],
            ibd_program,
            phecode_desc,
        )

        # This is the main function that will run the analysis of the networks
        analysis.analyze(data_container, output)

    logger.info("analysis_finished")
    logger.info(f"Program Duration: {datetime.now() - start_time}")


if __name__ == "__main__":
    app()
