from typing import Dict, List

from pandas import DataFrame, Series


def _identify_carriers_indx(
    carriers_series: Series, return_dict: Dict[str, List[str]]
) -> None:
    """Function that will insert the key, value pair of the phenotype and the carriers
    into the dictionary
    Parameters
    ----------
    carriers_series : pd.Series
        row of dataframe that is used in the apply function in
        the generate_carrier_list function

    return_dict : dict[str, list[str]]
        dictionary where the phecodes will be the keys and the values will be a list of
        indices indicating which grids are carriers
    """

    return_dict[carriers_series.name] = carriers_series[
        carriers_series == 1
    ].index.tolist()


def _identify_carriers_indx(
    carriers_series: Series, return_dict: Dict[str, List[str]]
) -> None:
    """Function that will insert the key, value pair of the phenotype and the carriers
    into the dictionary

    Parameters
    ----------
    carriers_series : pd.Series
        row of dataframe that is used in the apply function in
        the generate_carrier_list function

    return_dict : dict[str, list[str]]
        dictionary where the phecodes will be the keys and the values will be a list of
        indices indicating which grids are carriers
    """

    return_dict[carriers_series.name] = carriers_series[
        carriers_series == "NA"
    ].index.tolist()


def generate_carrier_dict(carriers_matrix: DataFrame) -> Dict[str, List[str]]:
    """Function that will take the carriers_pheno_matrix and generate a dictionary that
    has the list of indices for each carrier

    Parameters
    ----------
    carriers_matrix : pd.DataFrame
        dataframe where the columns are the phecodes and have 0's or 1's for whether or
        not they have the phecodes

    Returns
    -------
    Tuple[dict[str, list[str]], dict[str, List[str]]
        dictionary where the keys are phecodes and the values are list of integers
    """
    return_dict = {}

    new_indx_df = carriers_matrix.set_index("grids")

    # iterating over each phenotype which starts with the
    # second column
    new_indx_df.apply(lambda x: _identify_carriers_indx(x, return_dict), axis=0)

    return return_dict
