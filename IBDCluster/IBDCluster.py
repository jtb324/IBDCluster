from dotenv import load_dotenv
import callbacks
import typer
import pandas as pd
from collections import namedtuple
from typing import List, Optional, Dict, Set
from glob import glob
import os
from dataclasses import dataclass, field
# import cluster
from numpy import where


app = typer.Typer(
    add_completion=False,
    help="Tool that identifies ibd sharing for specific loci for individuals within biobanks",
)

@app.command
def main(
    IBD_programs: str = typer.Option(
        "hapibd",
        help="IBD detection software that the output came from. The program expects these values to be hapibd or ilash",
        callback=callbacks.check_ibd_program,
    ),
    env: str = typer.Option(
        "../.env",
        help="Filepath to an env file that has configuration setting. Program assumes that the default path is ../.env from the main file IBDCluster.py"
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



if __name__ == "__main__":
    app()

