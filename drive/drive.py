###

import itertools
import re
import sys
from pathlib import Path

import log
import pandas as pd
import typer
from callbacks import check_input_exists
from filter import IbdFilter
from models import FormatTypes, Genes, LogLevel, create_indices
from cluster import cluster, Networks

app = typer.Typer(add_completion=False)


def split_target_string(chromo_pos_str: str) -> Genes:
    """Function that will split the target string provided by the user.

    Parameters
    ----------
    chromo_pos_str : str
        string that has the region of interest in bbase pairs.
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
        5, "--max-recheck", help="Maximum number of times to check the clustering"
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
        3,
        "--min-network-size",
        help="This argument sets the minimun network size allowed for the clustering",
    ),
    loglevel: LogLevel = typer.Option(
        LogLevel.WARNING.value,
        "--loglevel",
        "-l",
        help="This argument sets the logging level for the program. Accepts values 'debug', 'warning', and 'verbose'.",
        case_sensitive=True,
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

    log.configure(logger, output.parent, log_filename, loglevel, log_to_console)

    indices = create_indices(ibd_format.lower())

    logger.debug(f"created indices object: {indices}")

    ##target gene region or variant position
    target_gene = split_target_string(target)

    logger.debug
    (f"Identified a target region: {target_gene}")

    # sys.exit()
    filter_obj = IbdFilter.load_file(input_file, indices, target_gene)

    filter_obj.preprocess(min_cm)

    # ibd_walktrap = ig.Graph.community_walktrap(ibd_g, weights="cm", steps=step)
    networks = Networks(minimum_connected_thres, max_network_size, max_check, step)
    # sys.exit()
    cluster(filter_obj, networks, min_network_size, indices.cM_indx)

    sys.exit()
    # allclst = [
    #     i for i, v in enumerate(ibd_walktrap_clusters.sizes()) if v > min_network_size
    # ]

    # def get_clst_info(name, clst, target_clsts, target_ig):
    #     clst_info[name] = {}
    #     clst_info[name]["memberID"] = [
    #         target_ig.vs()[c]["name"]
    #         for c, v in enumerate(target_clsts.membership)
    #         if v == clst
    #     ]
    #     clst_info[name]["member"] = [
    #         c for c, v in enumerate(target_clsts.membership) if v == clst
    #     ]
    #     clst_edge_n = (
    #         len(clst_info[name]["member"]) * (len(clst_info[name]["member"]) - 1) / 2
    #     )
    #     clst_info[name]["true_positive_edge"] = list(
    #         filter(
    #             lambda a: a != -1,
    #             target_ig.get_eids(
    #                 pairs=list(itertools.combinations(clst_info[name]["member"], 2)),
    #                 directed=False,
    #                 error=False,
    #             ),
    #         )
    #     )
    #     clst_info[name]["true_positive_n"] = len(clst_info[name]["true_positive_edge"])
    #     clst_info[name]["true_positive"] = (
    #         clst_info[name]["true_positive_n"] / clst_edge_n
    #     )
    #     all_edge = set([])
    #     for mem in clst_info[name]["member"]:
    #         all_edge = set(all_edge.union(set(target_ig.incident(mem))))
    #     clst_info[name]["false_negative_edge"] = list(
    #         all_edge.difference(
    #             list(
    #                 target_ig.get_eids(
    #                     pairs=list(
    #                         itertools.combinations(clst_info[name]["member"], 2)
    #                     ),
    #                     directed=False,
    #                     error=False,
    #                 )
    #             )
    #         )
    #     )
    #     clst_info[name]["false_negative_n"] = len(
    #         clst_info[name]["false_negative_edge"]
    #     )
    #     if (
    #         check_times < max_check
    #         and clst_info[name]["true_positive"] < minimum_connected_thres
    #         and len(clst_info[name]["member"]) > max_network_size
    #     ):
    #         recheck[check_times].append(name)
    #     else:
    #         outclst.append(name)

    # clst_info = {}
    # recheck = {}
    # check_times = 0
    # recheck[check_times] = []
    # outclst = []

    # for clst in allclst:
    #     get_clst_info(str(clst), clst, ibd_walktrap_clusters, ibd_g)

    def redo_clst(i):
        redopd = ibdpd.loc[
            (ibdpd["idnum1"].isin(clst_info[i]["memberID"]))
            & (ibdpd["idnum2"].isin(clst_info[i]["memberID"]))
        ]
        redo_g = ig.Graph.DataFrame(redopd, directed=False)
        redo_walktrap = ig.Graph.community_walktrap(redo_g, weights="cm", steps=step)
        redo_walktrap_clusters = redo_walktrap.as_clustering()
        if len(redo_walktrap_clusters.sizes()) == 1:
            clst_conn = pd.DataFrame(columns=["idnum", "conn", "conn.N", "TP"])
            for idnum in clst_info[i]["member"]:
                conn = sum(
                    list(
                        map(
                            lambda x: 1 / x,
                            redopd.loc[
                                (redopd["idnum1"] == idnum)
                                | (redopd["idnum2"] == idnum)
                            ]["cm"],
                        )
                    )
                )
                conn_idnum = list(
                    redopd.loc[(redopd["idnum1"] == idnum)]["idnum2"]
                ) + list(redopd.loc[(redopd["idnum2"] == idnum)]["idnum1"])
                conn_tp = len(
                    redopd.loc[
                        redopd["idnum1"].isin(conn_idnum)
                        & redopd["idnum2"].isin(conn_idnum)
                    ].index
                )
                if len(conn_idnum) == 1:
                    connTP = 1
                else:
                    connTP = conn_tp / (len(conn_idnum) * (len(conn_idnum) - 1) / 2)
                clst_conn.loc[idnum] = [idnum, conn, len(conn_idnum), connTP]
            rmID = list(
                clst_conn.loc[
                    (clst_conn["conn.N"] > (0.2 * len(clst_info[i]["member"])))
                    & (clst_conn["TP"] < minimum_connected_thres)
                    & (
                        clst_conn["conn"]
                        > sorted(clst_conn["conn"], reverse=True)[
                            int(0.01 * len(clst_info[i]["member"]))
                        ]
                    )
                ]["idnum"]
            )
            redopd = redopd.loc[
                (~redopd["idnum1"].isin(rmID)) & (~redopd["idnum2"].isin(rmID))
            ]
            redo_g = ig.Graph.DataFrame(redopd, directed=False)
            redo_walktrap = ig.Graph.community_walktrap(
                redo_g, weights="cm", steps=step
            )
            redo_walktrap_clusters = redo_walktrap.as_clustering()
        print(redo_walktrap_clusters.summary())
        allclst = [c for c, v in enumerate(redo_walktrap_clusters.sizes()) if v > 2]
        for clst in allclst:
            get_clst_info(
                "{0}.{1}".format(i, clst), clst, redo_walktrap_clusters, redo_g
            )

    while check_times < max_check and len(recheck[check_times]) > 0:
        check_times += 1
        print("recheck: {}".format(check_times))
        recheck[check_times] = []
        for redoclst in recheck[check_times - 1]:
            redo_clst(redoclst)
    #    print(recheck[check_times])

    with open("{}.DRIVE.txt".format(output), "w") as output:
        output.write(
            "## {0} IBD segments from {1} haplotypes\n".format(
                ibd_g.ecount(), ibd_g.vcount()
            )
        )
        output.write("## Identified {} IBD clusters\n".format(len(outclst)))
        output.write(
            "clstID\tn.total\tn.haplotype\ttrue.postive.n\ttrue.postive\tfalst.postive\tIDs\tID.haplotype\n"
        )
        for clst in outclst:
            clsthapid = list(
                ibdvs.loc[ibdvs["idnum"].isin(clst_info[clst]["memberID"])]["hapID"]
            )
            clstIID = set(
                ibdvs.loc[ibdvs["idnum"].isin(clst_info[clst]["memberID"])]["IID"]
            )
            n = len(clstIID)
            nhap = len(clsthapid)
            output.write(
                "clst{0}\t{1}\t{2}\t{3}\t{4:.4f}\t{5}\t{6}\t{7}\n".format(
                    clst,
                    n,
                    nhap,
                    clst_info[clst]["true_positive_n"],
                    clst_info[clst]["true_positive"],
                    clst_info[clst]["false_negative_n"],
                    ",".join(clstIID),
                    ",".join(clsthapid),
                )
            )


if __name__ == "__main__":
    app()
