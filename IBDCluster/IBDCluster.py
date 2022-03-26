#!/usr/bin/env python

import callbacks
import typer
import cluster
import log
import os
import analysis
from dotenv import load_dotenv
import pandas as pd
from typing import Dict, Tuple, List
from models import Writer
import pathlib


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
        "./.env",
        "--env",
        "-e",
        help="path to a .env file that has variables for the hapibd files directory and the ilash file directory. These variables are called HAPIBD_PATH and ILASH_PATH respectively.",
    ),
    gene_info_file: str = typer.Option(
        ...,
        "--gene_file",
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
    loglevel: str = typer.Option(
        "warning",
        "--loglevel",
        "-l",
        help="This argument sets the logging level for the program. Accepts values 'debug', 'warning', and 'verbose'.",
        callback=callbacks.check_loglevel,
    ),
    log_to_console: bool = typer.Option(
        False,
        "--log_to_console",
        help="Optional flag to log to only the console or also a file",
        is_flag=True,
    ),
) -> None:
    """Main function for the program that has all the parameters that the user can use with type"""

    # adding the loglevel to the environment so that we can access it
    os.environ.setdefault("program_loglevel", str(log.get_loglevel(loglevel)))

    # loading in the environmental variables from the .env file
    load_dotenv(env_path)

    # creating the logger and then configuring it
    logger = log.create_logger()

    log.configure(logger, output, loglevel=loglevel, to_console=log_to_console)

    # recording all the user inputs
    record_inputs(
        logger,
        ibd_program_used=ibd_program,
        output_path=output,
        environment_file=env_path,
        gene_info_file=gene_info_file,
        carrier_matrix=carriers,
        centimorgan_threshold=cm_threshold,
        loglevel=loglevel,
    )

    # need to first determine list of carriers for each phenotype
    carriers_df: pd.DataFrame = pd.read_csv(carriers, sep="\t")

    phecode_list = carriers_df.columns[1:]

    carriers_dict = cluster.generate_carrier_list(carriers_df)

    # We can then determine the different clusters for each gene
    networks: Dict[Tuple[str, int], List] = cluster.find_clusters(
        ibd_program, gene_info_file, cm_threshold, carriers_dict, phecode_list
    )

    # iterate over each object
    for gene, networks_list in networks.items():

        gene_output = os.path.join(output, gene[0])

        pathlib.Path(gene_output).mkdir(parents=True, exist_ok=True)

        # create an object that will be used to write to an
        # appropriate file
        write_obj = Writer(gene_output, ibd_program)

        # This is the main function that will run the analysis of the networks
        analysis.analyze(gene, networks_list, carriers_df, write_obj, carriers_dict)

    logger.info("analysis_finished")


if __name__ == "__main__":
    app()
