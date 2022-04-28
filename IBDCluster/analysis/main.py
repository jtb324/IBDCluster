from typing import Dict, Any
from models import DataHolder
import log
import plugins
import os
import json


logger = log.get_logger(__name__)


def analyze(data_container: DataHolder, output: str) -> None:
    """Main function from the analyze module that will determine the
    pvalues and the number of haplotypes, and individuals and will
    write this all to file.

    Parameters

    data_container : DataHoldexr

    output : str
    """

    # This section will load in the analysis plugins
    with open(os.environ.get("json_path")) as json_config:

        config = json.load(json_config)

        plugins.load_plugins(config["plugins"])

        analysis_plugins = [plugins.factory_create(item) for item in config["modules"]]

        logger.info(
            f"Using plugins: {', '.join([obj.name for obj in analysis_plugins])}"
        )

        # iterating over every plugin and then running the analyze and write method
        for analysis_obj in analysis_plugins:

            logger.debug(analysis_obj.name)

            return_data: Dict[str, Any] = analysis_obj.analyze(
                data=data_container, output=output
            )

            # writing the output to a file
            analysis_obj.write(
                input_data=return_data,
                ibd_program=data_container.ibd_program,
                phenotype_list=data_container.phenotype_cols,
            )
