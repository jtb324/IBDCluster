from dataclasses import dataclass
from typing import Dict, Any, Tuple
import pandas as pd
import log
import os
from plugins import factory_register
import pathlib

logger = log.get_logger(__name__)


@dataclass
class PhecodePercentages:

    name: str = "Population Phecode Percentages"

    def analyze(self, **kwargs) -> Dict[int, Any]:
        """
        Function that will determine the percentages of each phenotype and then add it to the dataholder object for other analyses.
        """
        data_container = kwargs["data"]

        logger.info("Determining the dataset prevalance of each phenotype")

        # creating a dictionary that has the phecode value as the key and
        # the phecode percentage as the value
        prevalence_dict = self._find_carrier_percentages(data_container.phenotype_table)
        # adding the phenotype prevalence to the datacontainer
        # since this plugin calculates that. This will allow
        # the phenotype percentages to be used by other plugins
        data_container.phenotype_percentages = prevalence_dict

        self.check_phenotype_prevalence(prevalence_dict)

        return {
            "output": data_container.phenotype_percentages,
            "path": kwargs["output"],
        }

    @staticmethod
    def _find_carrier_percentages(dataframe: pd.DataFrame) -> Dict[str, float]:
        """Function that will determine the percentages of carriers in each network

        Parameters

        col_series : pd.Series
            pandas series that has 0s and 1s for the carrier status of each
            individual for the specific phenotype

        """
        # get the carrier count for each column using sum and then dividing it by the total size of the
        # column to normalize

        normalized_carrier_counts: pd.Series = (
            dataframe.iloc[:, 1:].sum(axis=0) / dataframe.iloc[:, 1:].count()
        )

        return normalized_carrier_counts.to_dict()

    def write(self, **kwargs) -> None:
        """Writing the dictionary to a file

        Parameters (Expected to be keys in kwargs)

        percentage_dict : Dict[str, float]
            dictionary that has the phenotype name as the key and the percentage of
            individual carriers in the population as the value

        output : str
            filepath to write the output to
        """
        percentage_dict = kwargs["input_data"]["output"]
        output_path = kwargs["input_data"]["path"]

        pathlib.Path(output_path).mkdir(parents=True, exist_ok=True)

        output_file_name = os.path.join(
            output_path, "percent_carriers_in_population.txt"
        )

        logger.info(f"Writing file with population prevalences to {output_file_name}")

        with open(
            output_file_name,
            "w",
            encoding="utf-8",
        ) as output_file:
            output_file.write("phenotype\tpercentage_in_population\n")
            # iterate over each item in the percentage dict and write the items to a file
            for phenotype, percent in percentage_dict.items():
                output_file.write(f"{phenotype}\t{percent}\n")

    def check_phenotype_prevalence(self, percentage_dict: Dict[str, float]) -> None:
        """Function that will make sure that all the percentages are not 0. If they are then that will get logged to the output

        Parameters

        percentage_dict : Dict[str, float]
            dictionary where the key is the phenotype and the values are the prevalence in the population
        """

        if not any(percentage_dict.values()):
            logger.warning("All phenotypes have a population prevalence of 0%")


def initialize() -> None:
    factory_register("phecode_percentages", PhecodePercentages)
