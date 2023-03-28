###

from collections import namedtuple
from pathlib import Path
from enum import Enum
import re
import gzip
import igraph as ig
import pandas as pd
import itertools
import typer
from callbacks import check_input_exists
from generate_indices import create_indices

app = typer.Typer(add_completion=False)


Genes = namedtuple("Genes", ["chr", "start", "end"])


def split_target_string(chromo_pos_str: str) -> Genes:
    """Function that will split the target string provided by the user.

    Parameters
    ----------
    chromo_pos_str : str
        string that has the region of interest in bbase pairs. This string
        will look like 10:1234-1234 where the first number is the chromosome
        number, then the start position, and then the end position of the
        region of interest.

    Returns
    -------
    Genes
        returns a namedtuple that has the chromosome number, the start position, and the end position

    Raises
    ------
    ValueError
        raises a value error if the string was formatted any other way than chromosome:start_position-end_position
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
    input: Path = typer.Option(
        ..., "-i", "--input", help="IBD input file", callback=check_input_exists
    ),
    format: str = typer.Option(
        ..., "-f", "--format", help="IBD file format, e.g. hapIBD, iLASH"
    ),
    target: str = typer.Option(
        ...,
        "-t",
        "--target",
        help="Target region or position, chr:start-end or chr:pos",
    ),
    output: Path = typer.Option(..., "-o", "--output", help="output file prefix"),
    minCM: int = typer.Option(
        3, "-m", "--min-cm", help="minimum centimorgan threshold."
    ),
    step: int = typer.Option(3, "-k", "--step", help="steps for random walk"),
):
    # parser = argparse.ArgumentParser()
    # parser.add_argument('-i', '--input', type=str, required=True, help='IBD input file')
    # parser.add_argument('-f', '--format', type=str, required=True, help='IBD file format, e.g. hapIBD, iLASH, RaPID')
    # parser.add_argument('-t', '--target', type=str, required=True, help='Target region or position, chr:start-end or chr:pos')
    # parser.add_argument('-o', '--output', type=str, required=True, help='output file perfix')
    # parser.add_argument('-m', '--mincm', type=float, required=False, help='minimum cM, default=3', default=3)
    # parser.add_argument('-k', '--step', type=int, required=False, help='steps for random walk, default=3', default=3)

    # args = parser.parse_args()
    indices = create_indices(format.lower())
    ## set input format
    # id1_indx = 0
    # hap1_indx = 1
    # id2_indx = 2
    # hap2_indx = 3
    # chr_indx = 4
    # str_indx = 5
    # end_indx = 6

    # if format.lower() == "germline":
    #     cM_indx = 10
    #     unit = 11

    #     def getHAPID(IID, hapID):
    #         return hapID

    # elif format.lower() == "ilash":
    #     cM_indx = 9

    #     def getHAPID(IID, hapID):
    #         return hapID

    # elif format.lower() in ["hap-ibd", "hapibd"]:
    #     cM_indx = 7

    #     def getHAPID(IID, hapID):
    #         return "{0}.{1}".format(IID, hapID)

    # elif format.lower() == "rapid":
    #     id1_indx = 1
    #     hap1_indx = 3
    #     id2_indx = 2
    #     hap2_indx = 4
    #     chr_indx = 0
    #     cM_indx = 7

    #     def getHAPID(IID, hapID):
    #         return "{0}.{1}".format(IID, hapID)

    ##target gene region or variant position
    target_gene = split_target_string(target)
    # genechr = target.split(":")[0]
    # if len(target.split(":")[1].split("-")) == 2:
    #     genestr = int(target.split(":")[1].split("-")[0])
    #     geneend = int(target.split(":")[1].split("-")[1])
    # else:
    #     genestr, geneend = int(target.split(":")[1])

    ##other setting
    mincM = minCM
    kstep = step
    TP = 0.5
    maxN = 30
    maxcheck = 5

    ibdpd = pd.DataFrame(columns=["idnum1", "idnum2", "cm"])
    ibdvs = pd.DataFrame(columns=["idnum", "hapID", "IID"])
    hapid_to_int = {}
    allhapid = []
    idnum = int(0)

    with gzip.open(input, "rt") as ibdfile:
        for line in ibdfile:
            line = line.strip().split()
            CHR = line[indices.chr_indx]
            STR = min(int(line[indices.str_indx]), int(line[indices.end_indx]))
            END = max(int(line[indices.str_indx]), int(line[indices.end_indx]))
            iid1 = line[indices.id1_indx]
            iid2 = line[indices.id2_indx]
            hapid1 = str(
                indices.getHAPID(line[indices.id1_indx], line[indices.hap1_indx])
            )
            hapid2 = str(
                indices.getHAPID(line[indices.id2_indx], line[indices.hap2_indx])
            )
            cM = float(line[indices.cM_indx])
            if (
                CHR == target_gene.chr
                and STR <= target_gene.start
                and END >= target_gene.end
                and cM >= mincM
            ):
                if hapid1 != hapid2:
                    if hapid1 not in hapid_to_int:
                        hapid_to_int[hapid1] = int(idnum)
                        allhapid.append(hapid1)
                        ibdvs.loc[int(idnum)] = [int(idnum), hapid1, iid1]
                        idnum += 1
                    if hapid2 not in hapid_to_int:
                        hapid_to_int[hapid2] = int(idnum)
                        allhapid.append(hapid2)
                        ibdvs.loc[int(idnum)] = [int(idnum), hapid2, iid2]
                        idnum += 1
                    ibdpd = ibdpd.append(
                        {
                            "idnum1": int(hapid_to_int[hapid1]),
                            "idnum2": int(hapid_to_int[hapid2]),
                            "cm": cM,
                        },
                        ignore_index=True,
                    )
    ibdpd = ibdpd.astype({"idnum1": "int", "idnum2": "int"})
    # print(ibdpd)
    ibd_g = ig.Graph.DataFrame(ibdpd, directed=False, vertices=ibdvs, use_vids=True)
    ibd_walktrap = ig.Graph.community_walktrap(ibd_g, weights="cm", steps=kstep)
    ibd_walktrap_clusters = ibd_walktrap.as_clustering()
    print(ibd_walktrap_clusters.summary())
    allclst = [i for i, v in enumerate(ibd_walktrap_clusters.sizes()) if v > 2]

    def get_clst_info(name, clst, target_clsts, target_ig):
        clst_info[name] = {}
        clst_info[name]["memberID"] = [
            target_ig.vs()[c]["name"]
            for c, v in enumerate(target_clsts.membership)
            if v == clst
        ]
        clst_info[name]["member"] = [
            c for c, v in enumerate(target_clsts.membership) if v == clst
        ]
        clst_edge_n = (
            len(clst_info[name]["member"]) * (len(clst_info[name]["member"]) - 1) / 2
        )
        clst_info[name]["true_positive_edge"] = list(
            filter(
                lambda a: a != -1,
                target_ig.get_eids(
                    pairs=list(itertools.combinations(clst_info[name]["member"], 2)),
                    directed=False,
                    error=False,
                ),
            )
        )
        clst_info[name]["true_positive_n"] = len(clst_info[name]["true_positive_edge"])
        clst_info[name]["true_positive"] = (
            clst_info[name]["true_positive_n"] / clst_edge_n
        )
        all_edge = set([])
        for mem in clst_info[name]["member"]:
            all_edge = set(all_edge.union(set(target_ig.incident(mem))))
        clst_info[name]["false_negative_edge"] = list(
            all_edge.difference(
                list(
                    target_ig.get_eids(
                        pairs=list(
                            itertools.combinations(clst_info[name]["member"], 2)
                        ),
                        directed=False,
                        error=False,
                    )
                )
            )
        )
        clst_info[name]["false_negative_n"] = len(
            clst_info[name]["false_negative_edge"]
        )
        if (
            check_times < maxcheck
            and clst_info[name]["true_positive"] < TP
            and len(clst_info[name]["member"]) > maxN
        ):
            recheck[check_times].append(name)
        else:
            outclst.append(name)

    clst_info = {}
    recheck = {}
    check_times = 0
    recheck[check_times] = []
    outclst = []

    for clst in allclst:
        get_clst_info(str(clst), clst, ibd_walktrap_clusters, ibd_g)

    def redo_clst(i):
        redopd = ibdpd.loc[
            (ibdpd["idnum1"].isin(clst_info[i]["memberID"]))
            & (ibdpd["idnum2"].isin(clst_info[i]["memberID"]))
        ]
        redo_g = ig.Graph.DataFrame(redopd, directed=False)
        redo_walktrap = ig.Graph.community_walktrap(redo_g, weights="cm", steps=kstep)
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
                    & (clst_conn["TP"] < TP)
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
                redo_g, weights="cm", steps=kstep
            )
            redo_walktrap_clusters = redo_walktrap.as_clustering()
        print(redo_walktrap_clusters.summary())
        allclst = [c for c, v in enumerate(redo_walktrap_clusters.sizes()) if v > 2]
        for clst in allclst:
            get_clst_info(
                "{0}.{1}".format(i, clst), clst, redo_walktrap_clusters, redo_g
            )

    while check_times < maxcheck and len(recheck[check_times]) > 0:
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
