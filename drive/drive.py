import json
import re
from datetime import datetime
from pathlib import Path
from typing import Optional

import typer

import drive.factory as factory
from drive.cluster import ClusterHandler, cluster
from drive.filters import IbdFilter
from drive.log import CustomLogger
from drive.models import Data, FormatTypes, Genes, OverlapOptions, create_indices
from drive.utilities.callbacks import check_input_exists, check_json_path
from drive.utilities.parser import PhenotypeFileParser, load_phenotype_descriptions

app = typer.Typer(add_completion=False)


def split_target_string(chromo_pos_str: str) -> Genes:
    """Function that will split the target string provided by the user.

    Parameters
    ----------
    chromo_pos_str : str
        String that has the region of interest in base pairs.
        This string will look like 10:1234-1234 where the
        first number is the chromosome number, then the start
        position, and then the end position of the region of
        interest.

    Returns
    -------
    Genes
        returns a namedtuple that has the chromosome number,
        the start position, and the end position

    Raises
    ------
    ValueError
        raises a value error if the string was formatted any
        other way than chromosome:start_position-end_position.
        Also raises a value error if the start position is
        larger than the end position
    """
    split_str = re.split(":|-", chromo_pos_str)

    if len(split_str) != 3:
        error_msg = f"Expected the gene position string to be formatted like chromosome:start_position-end_position. Instead it was formatted as {chromo_pos_str}"  # noqa: E501

        raise ValueError(error_msg)

    integer_split_str = [int(value) for value in split_str]

    if integer_split_str[1] > integer_split_str[2]:
        raise ValueError(
            f"expected the start position of the target string to be <= the end position. Instead the start position was {integer_split_str[1]} and the end position was {integer_split_str[2]}"  # noqa: E501
        )

    return Genes(*integer_split_str)


@app.command()
def main(
    input_file: Path = typer.Option(
        ..., "-i", "--input", help="IBD input file", callback=check_input_exists
    ),
    ibd_format: FormatTypes = typer.Option(
        FormatTypes.HAPIBD.value,
        "-f",
        "--format",
        help="IBD file format. Allowed values are hapibd, ilash, germline, rapid",
    ),
    target: str = typer.Option(
        ...,
        "-t",
        "--target",
        help="Target region or position, chr:start-end or chr:pos",
    ),
    output: Path = typer.Option(..., "-o", "--output", help="output file prefix"),
    min_cm: int = typer.Option(
        3, "-m", "--min-cm", help="minimum centimorgan threshold."
    ),
    step: int = typer.Option(3, "-k", "--step", help="steps for random walk"),
    max_check: int = typer.Option(
        5,
        "--max-recheck",
        help="Maximum number of times to re-perform the clustering. This value will not be used if the flag --no-recluster is used.",  # noqa: E501
    ),
    case_file: Optional[Path] = typer.Option(
        None,
        "-c",
        "--cases",
        help="A file containing individuals who are cases. This file expects for there to be two columns. The first column will have individual ids and the second has status where cases are indicated by a 1 and control are indicated by a 0.",  # noqa: E501
    ),
    segment_overlap: OverlapOptions = typer.Option(
        OverlapOptions.CONTAINS.value,
        "--segment-overlap",
        help="Indicates if the user wants the gene to contain the whole target region or if it just needs to overlap the segment.",  # noqa: E501
    ),
    phenotype_description_file: Optional[Path] = typer.Option(
        None,
        "-d",
        "--descriptions",
        help="tab delimited text file that has descriptions for each phecode. this file should have two columns called phecode and phenotype",  # noqa: E501
    ),
    max_network_size: int = typer.Option(
        30, "--max-network-size", help="maximum network size allowed"
    ),
    minimum_connected_thres: float = typer.Option(
        0.5,
        "--min-connected-threshold",
        help="minimum connectedness ratio required for the network",
    ),
    min_network_size: int = typer.Option(
        2,
        "--min-network-size",
        help="This argument sets the minimun network size that we allow. All networks smaller than this size will be filtered out. If the user wishes to keep all networks they can set this to 0",  # noqa: E501
    ),
    segment_dist_threshold: float = typer.Option(
        0.2,
        "--segment-distribution-threshold",
        help="Threshold to filter the network length to remove hub individuals",
    ),
    hub_threshold: float = typer.Option(
        0.01,
        "--hub-threshold",
        help="Threshold to determine what percentage of hubs to keep",
    ),
    json_path: Path = typer.Option(
        None,
        "--json-config",
        "-j",
        help="path to the json config file",
        callback=check_json_path,
    ),
    recluster: bool = typer.Option(
        True,
        help="whether or not the user wishes the program to automically recluster based on things lik hub threshold, max network size and how connected the graph is. ",  # noqa: E501
    ),
    verbose: int = typer.Option(
        0,
        "--verbose",
        "-v",
        help="verbose flag indicating if the user wants more information",
        count=True,
    ),
    log_to_console: bool = typer.Option(
        False,
        "--log-to-console",
        help="Optional flag to log to only the console or also a file",
        is_flag=True,
    ),
    log_filename: str = typer.Option(
        "drive.log", "--log-filename", help="Name for the log output file."
    ),
) -> None:
    # getting the programs start time
    start_time = datetime.now()

    # creating and configuring the logger and then recording user inputs
    logger = CustomLogger.create_logger()

    logger.configure(output.parent, log_filename, verbose, log_to_console)
    # record the input parameters
    logger.record_inputs(
        ibd_file=input_file,
        ibd_program_used=ibd_format,
        gene_target_region=target,
        output_prefix=output,
        phenotype_description_file=phenotype_description_file,
        phenotype_file=case_file,
        minimum_centimorgan_threshold=min_cm,
        random_walk_step_size=step,
        max_recheck_times=max_check,
        max_network_size=max_network_size,
        minimum_connection_threshold=minimum_connected_thres,
        min_network_size=min_network_size,
        log_to_console=log_to_console,
        log_filename=log_filename,
        recluster=recluster,
    )

    logger.debug(f"Parent directory for log files and output: {output.parent}")

    logger.info(f"Analysis start time: {start_time}")
    # we need to load in the phenotype descriptions file to get
    # descriptions of each phenotype
    if phenotype_description_file:
        logger.verbose(
            f"Using the phenotype descriptions file at: {phenotype_description_file}"
        )
        desc_dict = load_phenotype_descriptions(phenotype_description_file)
    else:
        logger.verbose("No phenotype descriptions provided")
        desc_dict = {}

    # if the user has provided a phenotype file then we will determine case/control/
    # exclusion counts. Otherwise we return an empty dictionary
    if case_file:
        with PhenotypeFileParser(case_file) as phenotype_file:
            phenotype_counts, cohort_ids = phenotype_file.parse_cases_and_controls()

            logger.info(
                f"identified {len(phenotype_counts.keys())} phenotypes within the file {case_file}"  # noqa: E501
            )
    else:
        logger.info(
            "No phenotype information provided. Only the clustering step of the analysis will be performed"  # noqa: E501
        )

        phenotype_counts = {}
        cohort_ids = []

    indices = create_indices(ibd_format.lower())

    logger.debug(f"created indices object: {indices}")

    ##target gene region or variant position
    target_gene = split_target_string(target)

    logger.debug(f"Identified a target region: {target_gene}")

    filter_obj: IbdFilter = IbdFilter.load_file(input_file, indices, target_gene)

    # choosing the proper way to filter the ibd files
    filter_obj.set_filter(segment_overlap)

    filter_obj.preprocess(min_cm, cohort_ids)

    # We need to invert the hapid_map dictionary so that the
    # integer mappings are keys and the values are the
    # haplotype string
    hapid_inverted = {value: key for key, value in filter_obj.hapid_map.items()}

    # creating the object that will handle clustering within the networks
    cluster_handler = ClusterHandler(
        minimum_connected_thres,
        max_network_size,
        max_check,
        step,
        min_network_size,
        segment_dist_threshold,
        hub_threshold,
        hapid_inverted,
        recluster,
    )

    networks = cluster(filter_obj, cluster_handler, indices.cM_indx)

    # creating the data container that all the plugins can interact with
    plugin_api = Data(networks, output, phenotype_counts, desc_dict)

    logger.debug(f"Data container: {plugin_api}")

    # making sure that the output directory is created
    # This section will load in the analysis plugins and run each plugin
    with open(json_path, encoding="utf-8") as json_config:
        config = json.load(json_config)

        factory.load_plugins(config["plugins"])

        analysis_plugins = [factory.factory_create(item) for item in config["modules"]]

        logger.debug(
            f"Using plugins: {', '.join([obj.name for obj in analysis_plugins])}"
        )

        # iterating over every plugin and then running the analyze and write method
        for analysis_obj in analysis_plugins:
            analysis_obj.analyze(data=plugin_api)

    end_time = datetime.now()

    logger.info(
        f"Analysis finished at {end_time}. Total runtime: {end_time - start_time}"
    )


if __name__ == "__main__":
    app()
