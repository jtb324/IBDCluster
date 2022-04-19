from dataclasses import dataclass
from typing import Dict, Any
import pandas as pd
import log
import os

logger = log.get_logger(__name__)


@dataclass
class PhecodePercentages:

    name: str = "Population Phecode Percentages"

    def analyze(self, **kwargs) -> Dict[str, Any]:
        """
        Function that will determine the percentages of each phenotype and then write them to a file
        """
        pheno_matrix: pd.DataFrame = kwargs["pheno_matrix"]

        logger.info("Determining the dataset prevalance of each phenotype")

        # creating a dictionary that has the phecode value as the key and
        # the phecode percentage as the value
        prevalence_dict = self.find_carrier_percentages(pheno_matrix)

        self.check_phenotype_prevalence(prevalence_dict)

        return prevalence_dict

    @staticmethod
    def find_carrier_percentages(dataframe: pd.DataFrame) -> Dict[str, float]:
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

    @staticmethod
    def write(**kwargs) -> None:
        """Writing the dictionary to a file

        Parameters (Expected to be keys in kwargs)

        percentage_dict : Dict[str, float]
            dictionary that has the phenotype name as the key and the percentage of
            individual carriers in the population as the value

        output : str
            filepath to write the output to
        """
        output = kwargs["output"]
        percentage_dict = kwargs["input_dict"]

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

    def check_phenotype_prevalence(self, percentage_dict: Dict[str, float]) -> None:
        """Function that will make sure that all the percentages are not 0. If they are then that will get logged to the output

        Parameters

        percentage_dict : Dict[str, float]
            dictionary where the key is the phenotype and the values are the prevalence in the population
        """

        if not any(percentage_dict.values()):
            logger.warning("All phenotypes have a population prevalence of 0%")
