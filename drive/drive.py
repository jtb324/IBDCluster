###

import re
from pathlib import Path

import log
import typer
from callbacks import check_input_exists
from cluster import cluster, ClusterHandler
from filters import IbdFilter
from models import FormatTypes, Genes, create_indices

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
        raise ValueError(
            f"Expected the gene position string to be formatted like chromosome:start_position-end_position. Instead it was formatted as {chromo_pos_str}"
        )

    integer_split_str = [int(value) for value in split_str]

    if integer_split_str[1] > integer_split_str[2]:
        raise ValueError(
            f"expected the start position of the target string to be <= the end position. Instead the start position was {integer_split_str[1]} and the end position was {integer_split_str[2]}"
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
        help="Maximum number of times to re-perform the clustering. If the user wishes to not re-perform the clustering, this value can just be set to 0.",
    ),
    case_file: Path = typer.Option(
        "",
        "-c",
        "--cases",
        help="A file containing individuals who are cases. This file expects for there to be two columns. The first column will have individual ids and the second has status where cases are indicated by a 1 and control are indicated by a 0.",
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
        help="This argument sets the minimun network size that we allow. All networks smaller than this size will be filtered out. If the user wishes to keep all networks they can set this to 0",
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
    logger = log.create_logger()

    log.configure(logger, output.parent, log_filename, verbose, log_to_console)

    indices = create_indices(ibd_format.lower())

    logger.debug(f"created indices object: {indices}")

    ##target gene region or variant position
    target_gene = split_target_string(target)

    logger.debug(f"Identified a target region: {target_gene}")

    # sys.exit()
    filter_obj = IbdFilter.load_file(input_file, indices, target_gene)

    filter_obj.preprocess(min_cm)

    # creating the object that will handle clustering within the networks
    cluster_handler = ClusterHandler(
        minimum_connected_thres, max_network_size, max_check, step, min_network_size
    )

    cluster(filter_obj, cluster_handler, indices.cM_indx)

    # with open("{}.DRIVE.txt".format(output), "w") as output:
    #     output.write(
    #         "## {0} IBD segments from {1} haplotypes\n".format(
    #             ibd_g.ecount(), ibd_g.vcount()
    #         )
    #     )
    #     output.write("## Identified {} IBD clusters\n".format(len(outclst)))
    #     output.write(
    #         "clstID\tn.total\tn.haplotype\ttrue.postive.n\ttrue.postive\tfalst.postive\tIDs\tID.haplotype\n"
    #     )
    #     for clst in outclst:
    #         clsthapid = list(
    #             ibdvs.loc[ibdvs["idnum"].isin(clst_info[clst]["memberID"])]["hapID"]
    #         )
    #         clstIID = set(
    #             ibdvs.loc[ibdvs["idnum"].isin(clst_info[clst]["memberID"])]["IID"]
    #         )
    #         n = len(clstIID)
    #         nhap = len(clsthapid)
    #         output.write(
    #             "clst{0}\t{1}\t{2}\t{3}\t{4:.4f}\t{5}\t{6}\t{7}\n".format(
    #                 clst,
    #                 n,
    #                 nhap,
    #                 clst_info[clst]["true_positive_n"],
    #                 clst_info[clst]["true_positive"],
    #                 clst_info[clst]["false_negative_n"],
    #                 ",".join(clstIID),
    #                 ",".join(clsthapid),
    #             )
    #         )


if __name__ == "__main__":
    app()
