#!/usr/bin/env python 

from dotenv import load_dotenv
import callbacks
import typer
import os
import cluster


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
    print(verbose)
    os.environ["verbose"] = str(verbose)
    os.environ["loglevel"] = loglevel

    
@app.command()
def main(
    IBD_programs: str = typer.Option(
        "hapibd",
        help="IBD detection software that the output came from. The program expects these values to be hapibd or ilash",
        callback=callbacks.check_ibd_program,
    ),
    env: str = typer.Option(
        ...,
        help="Filepath to an env file that has configuration setting. Program assumes that the default path is ../.env from the main file IBDCluster.py"
    ),
    gene_info_file: str = typer.Option(
        ...,
        "--gene_file",
        "-g",
        help="Filepath to a text file that has information about the genes it should have four columns: Gene name, chromosome, gene start, and gene end. The file is expected to not have a header",
        callback=callbacks.check_gene_file,
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
    cluster.find_clusters(IBD_programs, gene_info_file)

if __name__ == "__main__":
    app()

