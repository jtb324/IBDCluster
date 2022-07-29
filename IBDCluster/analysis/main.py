from typing import Dict, Any
from models import DataHolder, Network
import log
import plugins
import os
import json


logger = log.get_logger(__name__)


def analyze(data_container: DataHolder, network_obj: Network, output: str) -> None:
    """Main function from the analyze module that will determine the
    pvalues and the number of haplotypes, and individuals and will
    write this all to file.

    Parameters

    data_container : DataHolder
        object that hass all the attributes that are used by
        different plugins

    network_obj : Network
        network object that has attributes about haplotypes and
        who is in the network

    output : str
        filepath to write the output to
    """

    # making sure that the output directory is created
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
            logger.debug(f"output being used in analysis: {output}")
            analysis_obj.analyze(
                data=data_container, network=network_obj, output=output
            )
