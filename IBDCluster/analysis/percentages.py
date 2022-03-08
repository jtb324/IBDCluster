# script to determine percentages of each phenotype in the subset
from typing import Dict
import os
import pandas as pd
from tqdm import tqdm
import log

logger = log.get_logger(__name__)



def _find_carrier_percentages(
    col_series: pd.Series, percentage_dict: Dict[str, float]
) -> None:
    """Function that will determine the percentages of carriers in each network

    Parameters

    col_series : pd.Series
        pandas series that has 0s and 1s for the carrier status of each
        individual for the specific phenotype

    percentage_dict : Dict[str, float]
        dictionary that has the phenotype name as the key and the
        percentage of carriers out of the population as values
    """
    try:
        percentage_dict[col_series.name] = col_series.value_counts(normalize=True)[1]
        logger.debug(
            f"Phenotype prevalence for {col_series.name} is {col_series.value_counts(normalize=True)[1]}"
        )
    except KeyError:
        logger.warning(
            f"There were no individuals in the population affected by phenotype {col_series.name}"
        )
        percentage_dict[col_series.name] = 0


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

    percent_carriers: Dict[str, float] = {}

    # setting up a progress bar for this operation
    tqdm.pandas(desc="phenotypes population percentages calculated: ")
    carrier_matrix.iloc[:, 1:].progress_apply(
        lambda x: _find_carrier_percentages(x, percent_carriers)
    )

    return percent_carriers
