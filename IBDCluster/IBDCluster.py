#!/usr/bin/env python

from dotenv import load_dotenv
import callbacks
import typer
import cluster
import log
import analysis
import pandas as pd
from typing import Dict, Tuple
from models import Writer


app = typer.Typer(
    add_completion=False,
    help="Tool that identifies ibd sharing for specific loci for individuals within biobanks",
)

def record_inputs(logger, **kwargs) -> None:
    """function to record the user arguments that were passed to the 
    program. Takes a logger and then a dictionary of the user 
    arguments"""

    for parameter, value in kwargs.items():
        logger.info(f"{parameter}: {value}")

@app.command()
def main(
    IBD_program: str = typer.Option(
        "hapibd",
        "--ibd",
        "-i",
        help="IBD detection software that the output came from. The program expects these values to be hapibd or ilash",
        callback=callbacks.check_ibd_program,
    ),
    output: str = typer.Option(
        "./", "--output", "-o", help="directory to write the output files into."
    ),
    env: str = typer.Option(
        ...,
        "--env",
        "-e",
        help="Filepath to an env file that has configuration setting. Program assumes that the default path is ../.env from the main file IBDCluster.py",
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
    cM_threshold: int = typer.Option(
        3,
        "--cM",
        help="Centimorgan threshold to filter the ibd segments",
    ),
    loglevel: str = typer.Option(
        "warning",
        "--loglevel",
        "-l",
        help="This argument sets the logging level for the program",
        callback=callbacks.check_loglevel,
    ),
    log_to_console: bool = typer.Option(
        False,
        "--log_to_console",
        help="Optional flag to log to only the console or also a file",
        is_flag=True,
    ),
) -> None:
    """Main function for the program that has all the parameters that the user can use with type
    """
    # loading in environmental variables from an env file
    load_dotenv(env)

    # creating the logger and then configuring it
    logger = log.create_logger()

    log.configure(logger, output, loglevel=loglevel, to_console=log_to_console)

    # recording all the user inputs
    record_inputs(
        logger,
        ibd_program_used=IBD_program,
        output_path=output,
        environment_file=env,
        gene_info_file=gene_info_file,
        carrier_matrix=carriers,
        centimorgan_threshold=cM_threshold,
        loglevel=loglevel
        )

    networks: Dict[Tuple[str, int], Dict] = cluster.find_clusters(
        IBD_program, gene_info_file, cM_threshold
    )
    
    # create an object that will be used to write to an
    # appropriate file
    write_obj = Writer(output, IBD_program)

    carriers_df: pd.DataFrame = pd.read_csv(carriers, sep="\t")

    for gene, networks_info in networks.items():
        # This is the main function that will run the analysis of the networks
        analysis.analyze(gene, networks_info, carriers_df, write_obj)

    logger.info("analysis_finished")

if __name__ == "__main__":
    app()
