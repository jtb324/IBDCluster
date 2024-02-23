import json
import re
from datetime import datetime
from pathlib import Path
import argparse
from rich_argparse import RichHelpFormatter


import drive.factory as factory
from drive.cluster import ClusterHandler, cluster
from drive.filters import IbdFilter
from drive.log import CustomLogger
from drive.models import Data, Genes, create_indices
from drive.utilities.callbacks import CheckInputExist
from drive.utilities.parser import PhenotypeFileParser, load_phenotype_descriptions


def find_json_file() -> Path:
    """Method to find the default config file if the user does not provide one

    Returns
    -------
    Path
        returns the path to the json file

    Raises
    ------
    FileNotFoundError
        Raises a FileNotFoundError if the program can not locate a json file and the
        user does not provide the path to a file
    """

    src_dir = Path(__file__).parent

    config_path = src_dir / "config.json"

    if not config_path.exists():
        raise FileNotFoundError(
            f"Expected the user to either pass a configuration file path or for the config.json file to be present in the program root directory at {config_path}."  # noqa: E501
        )

    return config_path


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


def main() -> None:
    parser = argparse.ArgumentParser(
        description=" Distant Relatedness for Identification and Variant Evaluation (DRIVE) is a novel approach to IBD-based genotype inference used to identify shared chromosomal segments in dense genetic arrays",
        formatter_class=RichHelpFormatter,
    )

    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        help="IBD input file from ibd detection software",
        required=True,
        action=CheckInputExist,
    )

    parser.add_argument(
        "--format",
        "-f",
        default="hapibd",
        type=str,
        help="IBD program used to detect segments. Allowed values are hapibd, ilash, germline, rapid. Program expects for value to be lowercase. (default: %(default)s)",
        choices=["hapibd", "ilash", "germline", "rapid"],
    )

    parser.add_argument(
        "--target",
        "-t",
        type=str,
        help="Target region or position, chr:start-end or chr:pos",
        required=True,
    )

    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        help="output file prefix. The program will append .drive_networks.txt to the filename provided",
        required=True,
    )

    parser.add_argument(
        "--min-cm",
        "-m",
        default=3,
        type=int,
        help="minimum centimorgan threshold. The program expects this to be an integer value. (default: %(default)s)",
    )

    parser.add_argument(
        "--step",
        "-k",
        default=3,
        type=int,
        help="Minimum required number of steps for the community walktrap algorithm.(default: %(default)s)",
    )

    parser.add_argument(
        "--max-recheck",
        default=5,
        type=int,
        help="Maximum number of times to re-perform the clustering. This value will not be used if the flag --no-recluster is used.(default: %(default)s)",  # noqa: E501
    )

    parser.add_argument(
        "--cases",
        "-c",
        type=Path,
        help="A file containing individuals who are cases. This file expects for there to be two columns. The first column will have individual ids and the second has status where cases are indicated by a 1, control are indicated by a 0, and exclusions are indicated by NA.",  # noqa: E501
        action=CheckInputExist,
    )

    parser.add_argument(
        "--segment-overlap",
        default="contains",
        choices=["contains", "overlaps"],
        type=str,
        help="Indicates if the user wants the gene to contain the whole target region or if it just needs to overlap the segment. (default: %(default)s)",  # noqa: E501
    )

    parser.add_argument(
        "--descriptions",
        "-d",
        type=Path,
        help="tab delimited text file that has descriptions for each phecode. this file should have two columns called phecode and phenotype",  # noqa: E501
    )

    parser.add_argument(
        "--max-network-size",
        default=30,
        type=int,
        help="maximum network size allowed if the user has allowed the recluster option. (default: %(default)s)",
    )

    parser.add_argument(
        "--min-connected-threshold",
        default=0.5,
        type=float,
        help="minimum connectedness ratio required for the network. (default: %(default)s)",
    )

    parser.add_argument(
        "--min-network-size",
        default=2,
        type=int,
        help="This argument sets the minimun network size that we allow. All networks smaller than this size will be filtered out. If the user wishes to keep all networks they can set this to 0. (default: %(default)s)",  # noqa: E501
    )

    parser.add_argument(
        "--segment-distribution-threshold",
        default=0.2,
        type=float,
        help="Threshold to filter the network length to remove hub individuals. (default: %(default)s)",
    )

    parser.add_argument(
        "--phenotype-name",
        default=None,
        type=str,
        help="If the user wishes to only run 1 phenotype from a file with multiple phenotypes they can prioritize a column that has the phenotype name. The name must match with the column.",
    )

    parser.add_argument(
        "--hub-threshold",
        default=0.01,
        type=float,
        help="Threshold to determine what percentage of hubs to keep. (default: %(default)s)",
    )

    parser.add_argument(
        "--json-config",
        "-j",
        default=None,
        type=Path,
        help="path to the json config file",
    )

    parser.add_argument(
        "--recluster",
        type=bool,
        default=True,
        action=argparse.BooleanOptionalAction,
        help="whether or not the user wishes the program to automically recluster based on things like hub threshold, max network size and how connected the graph is. ",  # noqa: E501
    )

    parser.add_argument(
        "--verbose",
        "-v",
        default=0,
        help="verbose flag indicating if the user wants more information",
        action="count",
    )

    parser.add_argument(
        "--log-to-console",
        default=False,
        help="Optional flag to log to only the console or also a file",
        action="store_true",
    )

    parser.add_argument(
        "--log-filename",
        default="drive.log",
        type=str,
        help="Name for the log output file. (default: %(default)s)",
    )

    args = parser.parse_args()

    # getting the programs start time
    start_time = datetime.now()

    # We need to make sure that there is a configuration file
    json_config = args.json_config if args.json_config else find_json_file()

    # creating and configuring the logger and then recording user inputs
    logger = CustomLogger.create_logger()

    logger.configure(
        args.output.parent, args.log_filename, args.verbose, args.log_to_console
    )

    # record the input parameters using a method from the logger object that
    # takes the parser as an argument
    logger.record_namespace(args)

    logger.debug(f"Parent directory for log files and output: {args.output.parent}")

    logger.info(f"Analysis start time: {start_time}")
    # we need to load in the phenotype descriptions file to get
    # descriptions of each phenotype
    if args.descriptions:
        logger.verbose(f"Using the phenotype descriptions file at: {args.descriptions}")
        desc_dict = load_phenotype_descriptions(args.descriptions)
    else:
        logger.verbose("No phenotype descriptions provided")
        desc_dict = {}

    # if the user has provided a phenotype file then we will determine case/control/
    # exclusion counts. Otherwise we return an empty dictionary
    if args.cases:
        with PhenotypeFileParser(args.cases, args.phenotype_name) as phenotype_file:
            phenotype_counts, cohort_ids = phenotype_file.parse_cases_and_controls()

            logger.info(
                f"identified {len(phenotype_counts.keys())} phenotypes within the file {args.cases}"  # noqa: E501
            )
    else:
        logger.info(
            "No phenotype information provided. Only the clustering step of the analysis will be performed"  # noqa: E501
        )

        phenotype_counts = {}
        cohort_ids = []

    indices = create_indices(args.format.lower())

    logger.debug(f"created indices object: {indices}")

    ##target gene region or variant position
    target_gene = split_target_string(args.target)

    logger.debug(f"Identified a target region: {target_gene}")

    filter_obj: IbdFilter = IbdFilter.load_file(args.input, indices, target_gene)

    # choosing the proper way to filter the ibd files
    filter_obj.set_filter(args.segment_overlap)

    filter_obj.preprocess(args.min_cm, cohort_ids)

    # We need to invert the hapid_map dictionary so that the
    # integer mappings are keys and the values are the
    # haplotype string
    hapid_inverted = {value: key for key, value in filter_obj.hapid_map.items()}

    # creating the object that will handle clustering within the networks
    cluster_handler = ClusterHandler(
        args.min_connected_threshold,
        args.max_network_size,
        args.max_recheck,
        args.step,
        args.min_network_size,
        args.segment_distribution_threshold,
        args.hub_threshold,
        hapid_inverted,
        args.recluster,
    )

    networks = cluster(filter_obj, cluster_handler, indices.cM_indx)

    # creating the data container that all the plugins can interact with
    plugin_api = Data(networks, args.output, phenotype_counts, desc_dict)

    logger.debug(f"Data container: {plugin_api}")

    # making sure that the output directory is created
    # This section will load in the analysis plugins and run each plugin
    with open(json_config, encoding="utf-8") as json_config:
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
    main()
