#!/usr/bin/env python 

from dotenv import load_dotenv
import callbacks
import typer
import os
import cluster
from typing import Dict, Tuple
from models import Writer, Pair_Writer, Network_Writer


app = typer.Typer(
    add_completion=False,
    help="Tool that identifies ibd sharing for specific loci for individuals within biobanks",
)

def load_env_var(verbose: bool, loglevel: str) -> None:
    """Function that will add variables to the environment
    
    Parameters
    
    verbose : bool
        either true of False for if the user wants a more verbose setting
        
    loglevel : str
        level the user wishes to set for logging. It defaults to INFO
    """
    if verbose:
        print(f"adding verbosity settings and logging levels to environmental variables")

    os.environ["verbose"] = str(verbose)
    os.environ["loglevel"] = loglevel

    
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
        "./",
        "--output",
        "-o",
        help="directory to write the output files into."
    ),
    env: str = typer.Option(
        ...,
        "--env",
        "-e",
        help="Filepath to an env file that has configuration setting. Program assumes that the default path is ../.env from the main file IBDCluster.py"
    ),
    gene_info_file: str = typer.Option(
        ...,
        "--gene_file",
        "-g",
        help="Filepath to a text file that has information about the genes it should have four columns: Gene name, chromosome, gene start, and gene end. The file is expected to not have a header",
        callback=callbacks.check_gene_file,
    ),
    cM_threshold: int = typer.Option(
        3,
        "--cM",
        help="Centimorgan threshold to filter the ibd segments",
    ),
    verbose: bool = typer.Option(
        False, 
        "--verbose", 
        "-v", 
        help="Optional Flag to run the program in verbose mode", 
        is_flag=True
    ),
    loglevel: str = typer.Option(
        "INFO",
        "--loglevel",
        "-l",
        help="This argument sets the logging level for the program"
        )
) -> None:
    #loading in environmental variables from an env file
    load_dotenv(env)

    load_env_var(verbose, loglevel)
    
    networks: Dict[Tuple[str, int], Dict] = cluster.find_clusters(IBD_program, gene_info_file, cM_threshold)

    write_obj = Writer(IBD_program)
    
    for gene, networks_info in networks.items():
        write_obj.set_writer(Network_Writer(gene[0], gene[1]))

        write_obj.write_to_file(output, networks_info)

if __name__ == "__main__":
    app()

