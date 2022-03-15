import numpy as np
import pandas as pd
from typing import List, Dict
from models import Network
from tqdm import tqdm
from collections import namedtuple


class CarriersInfo(
    namedtuple(
        "Carrier_Comp",
        ["ind_in_network", "percentage", "IIDs", "pvalue", "network_len"],
    )
):
    """This is an extension of the class tuple so that I can overwrite the __str__ method"""

    def __str__(self):
        return f"Carrier Info Object - Individuals in Network: {self.ind_in_network}, Percentages: {self.percentage}, IID List: {self.IIDs}, pvalue: {self.pvalue}"


def _gini(input_vector: List[str], weights=None) -> float:
    """Function to calculate the gini coeffficient"""
    # Array indexing requires reset indexes.
    input_vector = pd.Series(input_vector).reset_index(drop=True)
    if weights is None:
        weights = np.ones_like(input_vector)
    weights = pd.Series(weights).reset_index(drop=True)
    n = input_vector.size
    wxsum = sum(weights * input_vector)
    wsum = sum(weights)
    sxw = np.argsort(input_vector)
    sx = input_vector[sxw] * weights[sxw]
    sw = weights[sxw]
    pxi = np.cumsum(sx) / wxsum
    pci = np.cumsum(sw) / wsum
    g = 0.0
    for i in np.arange(1, n):
        g = g + pxi.iloc[i] * pci.iloc[i - 1] - pci.iloc[i] * pxi.iloc[i - 1]
    return g


def _process(phenotype_dict: Dict[str, CarriersInfo], carrier_dict: Dict) -> None:
    """Function that will create a dictionary where the keys are the phenotype and the keys are the the phecode and then the number of carriers for each network is stored in a list in an inner key"""

    for phenotype, carrier_obj in tqdm(
        phenotype_dict.items(), desc="gini coefficients calculated: "
    ):

        carrier_dict.setdefault(phenotype, {"num_of_carriers": []})

        carrier_dict[phenotype]["num_of_carriers"].append(carrier_obj.ind_in_network)


def _determine_gini_coefficient(
    phecode_carriers: Dict[str, List[int]]
) -> Dict[str, float]:
    """Function that will determine the gini coefficient for each phecode"""

    gini_coef: Dict[str, float] = {}

    for phecode, carrier_list in phecode_carriers.items():
        print(carrier_list)
        gini_coef[phecode] = _gini(carrier_list["num_of_carriers"])

    return gini_coef


def determine_phecode_genie_coefficient(
    networks_list: List[Network],
) -> Dict[str, float]:
    """Function that will determine the genie coefficient for each phenotype in each network and add that to an attribute called gini_coef"""

    carrier_dict: Dict[str, Dict[str, List[int]]] = {}

    for network in networks_list:

        _process(network.phenotype, carrier_dict)

    return _determine_gini_coefficient(carrier_dict)
