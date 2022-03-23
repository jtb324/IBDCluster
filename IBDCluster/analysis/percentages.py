# script to determine percentages of each phenotype in the subset
from typing import Dict
import os
import pandas as pd
from tqdm import tqdm
import log

logger = log.get_logger(__name__)


def _find_carrier_percentages(dataframe: pd.DataFrame) -> Dict[str, float]:
    """Function that will determine the percentages of carriers in each network

    Parameters

    col_series : pd.Series
        pandas series that has 0s and 1s for the carrier status of each
        individual for the specific phenotype

    """
    # get the carrier count for each column using sum and then dividing it by the total size of the
    # column to normalize
    normalized_carrier_counts: pd.Series = dataframe.sum(axis=0) / dataframe.count()

    return normalized_carrier_counts.to_dict()


def write_to_file(percentage_dict: Dict[str, float], output: str) -> None:
    """Writing the dictionary to a file

    Parameters

    percentage_dict : Dict[str, float]
        dictionary that has the phenotype name as the key and the percentage of
        individual carriers in the population as the value

    output : str
        filepath to write the output to
    """
    output_file_name = os.path.join(output, "percent_carriers_in_population.txt")

    logger.info(f"Writing file with population prevalences to {output_file_name}")

    with open(
        output_file_name,
        "w",
        encoding="utf-8",
    ) as output_file:
        output_file.write("phenotype\tpercentage_in_population\n")
        for phenotype, percent in percentage_dict.items():
            output_file.write(f"{phenotype}\t{percent}\n")


def check_phenotype_prevalence(percentage_dict: Dict[str, float]) -> None:
    """Function that will make sure that all the percentages are not 0. If they are then that will get logged to the output

    Parameters

    percentage_dict : Dict[str, float]
        dictionary where the key is the phenotype and the values are the prevalence in the population
    """

    if not any(percentage_dict.values()):
        logger.warning("All phenotypes have a population prevalence of 0%")


def get_percentages(carrier_matrix: pd.DataFrame) -> Dict[str, float]:
    """Function that will take the carrier matrix and determine the percentage of carriers out of the whole dataset

    Parameters

    carrier_matrix : pd.DataFrame
        dataframe that has 1s or 0s for whether each individual is a carrier for the
        phenotype or not. The first column needs to be named grids

    output : str
        filepath to the directory to write the output to

    Returns

    Dict[str, float]
        returns a dictionary that has the phenotype as a string and the percentage of carriers as the value
    """
    logger.info("Determining the dataset prevalance of each phenotype")

    # setting up a progress bar for this operation
    percent_carriers = _find_carrier_percentages(carrier_matrix.iloc[:, 1:])

    return percent_carriers
