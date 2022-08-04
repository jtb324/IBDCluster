"""
File with functions to determine the prevalence of the phenotype 
within the provided population
"""
import log

logger = log.get_logger(__name__)


def get_phenotype_prevalances(
    phenotype_carriers: dict[float, list[str]], individual_count: int
) -> dict[str, float]:
    """
    Function that will determine the percentages of each phenotype
    and then add it to the dataholder object for other analyses.

    Parameters
    ----------
    phenotype_carriers : dict[str, list[str]]
        dictionary where the keys are strings of phecodes and the values are the
        list of IIDs that are affected with the phenotype

    individual_count : int
        integer that is the total number of individuals in the phenotype dataset

    Returns
    -------
    dict[str, float]
        returns a dictionary
    """

    logger.info("Determining the dataset prevalance of each phenotype")

    # creating a dictionary that has the phecode value as the key and
    # the phecode percentage as the value
    prevalence_dict = _find_carrier_percentages(phenotype_carriers, individual_count)

    # checking to make sure that not all of the phenotypes have a 0% prevalance
    check_phenotype_prevalence(prevalence_dict)

    return prevalence_dict


@staticmethod
def _find_carrier_percentages(
    carriers_dict: dict[float, list[str]], all_individual_count: int
) -> dict[str, float]:
    """Function that will determine the percentages of carriers in each network

    Parameters
    ----------
    col_series : pd.Series
        pandas series that has 0s and 1s for the carrier status of each
        individual for the specific phenotype

    """
    return_dict = {}
    # iterating over each key, value pair and then returning a list of tuples
    # where the first value is the phecode and the second value is the phenotype
    # frequency
    phecode_percents = list(
        map(lambda x: (x[0], (len(x[1]) / all_individual_count)), carriers_dict.items())
    )
    # converting the list of tuples into a dictionary
    for percent_tuple in phecode_percents:
        return_dict[percent_tuple[0]] = percent_tuple[1]

    return return_dict


def check_phenotype_prevalence(percentage_dict: dict[str, float]) -> None:
    """Function that will make sure that all the percentages are not 0. If they are then that will get logged to the output

    Parameters
    ----------
    percentage_dict : Dict[str, float]
        dictionary where the key is the phenotype and the values are the prevalence in the population
    """

    if not any(percentage_dict.values()):
        logger.warning("All phenotypes have a population prevalence of 0%")
