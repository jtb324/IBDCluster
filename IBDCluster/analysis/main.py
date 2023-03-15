"""
Main module that is used to call the different analysis plugins. This module takes
advantage of the factory to dynamically load the modules while also passing the 
appropriate information to the modules
"""
import json
import os

import factory
import log
from models import DataHolder, Network

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
    with open(os.environ.get("json_path"), encoding="utf-8") as json_config:

        config = json.load(json_config)

        factory.load_plugins(config["plugins"])

        analysis_plugins = [factory.factory_create(item) for item in config["modules"]]

        logger.debug(
            f"Using plugins: {', '.join([obj.name for obj in analysis_plugins])}"
        )

        # iterating over every plugin and then running the analyze and write method
        for analysis_obj in analysis_plugins:

            logger.debug(f"output path being used in analysis: {output}")
            analysis_obj.analyze(
                data=data_container, network=network_obj, output=output
            )
