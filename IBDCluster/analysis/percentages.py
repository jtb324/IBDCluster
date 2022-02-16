# script to determine percentages of each phenotype in the subset
from typing import Dict
from black import out
import pandas as pd
import os


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
    percentage_dict[col_series.name] = col_series.value_counts(normalize=True)[1]


def write_to_file(percentage_dict: Dict[str, float], output: str) -> None:
    """Writing the dictionary to a file

    Parameters

    percentage_dict : Dict[str, float]
        dictionary that has the phenotype name as the key and the percentage of
        individual carriers in the population as the value

    output : str
        filepath to write the output to
    """
    with open(
        os.path.join(output, "percent_carriers_in_population.txt"),
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

    percent_carriers: Dict[str, float] = {}

    carrier_matrix.iloc[:, 1:].apply(
        lambda x: _find_carrier_percentages(x, percent_carriers)
    )

    return percent_carriers
