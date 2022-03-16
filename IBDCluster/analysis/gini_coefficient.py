from typing import List, Dict
from models import Network
from tqdm import tqdm
from collections import namedtuple
import log

logger = log.get_logger(__name__)


class CarriersInfo(
    namedtuple(
        "Carrier_Comp",
        ["ind_in_network", "percentage", "IIDs", "pvalue", "network_len"],
    )
):
    """This is an extension of the class tuple so that I can overwrite the __str__ method"""

    def __str__(self):
        return f"Carrier Info Object - Individuals in Network: {self.ind_in_network}, Percentages: {self.percentage}, IID List: {self.IIDs}, pvalue: {self.pvalue}"


def _gini(input_vector: List[str], phecode: str) -> float:
    """Function to calculate the gini coeffficient. This uses the biased calculation. Found at https://www.had2know.org/academics/gini-coefficient-calculator.html"""

    input_vector.sort()

    upper_val = 0

    for index, value in enumerate(input_vector):

        upper_val += (len(input_vector) + 1 - (index + 1)) * value

    try:
        return ((len(input_vector) + 1) / len(input_vector)) - (
            (2 * upper_val) / (len(input_vector) * sum(input_vector))
        )
    except ZeroDivisionError:
        logger.warning(
            f"Zero division error produced when determining the gini coefficient for {phecode}. The sum of the input vector was {sum(input_vector)}"
        )


def _process(phenotype_dict: Dict[str, CarriersInfo], carrier_dict: Dict) -> None:
    """Function that will create a dictionary where the keys are the phenotype and the keys are the the phecode and then the number of carriers for each network is stored in a list in an inner key"""

    for phenotype, carrier_obj in tqdm(
        phenotype_dict.items(), desc="gini coefficients calculated: "
    ):

        _ = carrier_dict.setdefault(phenotype, {"num_of_carriers": []})

        carrier_dict[phenotype]["num_of_carriers"].append(carrier_obj.ind_in_network)


def _determine_gini_coefficient(
    phecode_carriers: Dict[str, List[int]]
) -> Dict[str, float]:
    """Function that will determine the gini coefficient for each phecode"""

    gini_coef: Dict[str, float] = {}

    for phecode, carrier_list in phecode_carriers.items():
        print(carrier_list)
        gini_coef[phecode] = _gini(carrier_list["num_of_carriers"], phecode)

    return gini_coef


def determine_phecode_genie_coefficient(
    networks_list: List[Network],
) -> Dict[str, float]:
    """Function that will determine the genie coefficient for each phenotype in each network and add that to an attribute called gini_coef"""

    carrier_dict: Dict[str, Dict[str, List[int]]] = {}

    for network in networks_list:

        _process(network.phenotype, carrier_dict)

    return _determine_gini_coefficient(carrier_dict)
