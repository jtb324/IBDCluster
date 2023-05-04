from pathlib import Path
from typing import Dict, Union

from pandas import read_csv


def load_phenotype_descriptions(
    phecode_desc_file: Union[Path, str]
) -> Dict[str, Dict[str, str]]:
    """Function that will load the phecode_description file
    and then turn that into a dictionary

    Parameters
    ----------
    phecode_desc_file : str
        descriptions of each phecode

    Returns
    -------
    Dict[str, Dict[str, str]]
        returns a dictionary where the first key is the
        phecode and value is a dictionary where the inner key
        is 'phenotype' and the value is the descriptions
    """

    desc_df = read_csv(phecode_desc_file, sep="\t", usecols=["phecode", "phenotype"])

    # making sure that the phecode keys are a string
    desc_df.phecode = desc_df.phecode.astype(str)

    # converting the dataframe into a dictionar
    return desc_df.set_index("phecode").T.to_dict()
