import itertools
import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Type, TypeVar

import igraph as ig
import log
from models import Filter
from pandas import DataFrame
from rich.progress import Progress

# creating a logger
logger: logging.Logger = log.get_logger(__name__)

# Create a generic variable that can represents the class from the
T = TypeVar("T", bound="Networks")


@dataclass
class Networks:
    graph: ig.Graph
    clusters: Dict[int, Any] = field(default_factory=dict)
    final_clusters: List[int] = field(default_factory=list)
    parent_cluster: Optional[int] = None

    @classmethod
    def generate_graph(
        cls: Type[T],
        ibd_edges: DataFrame,
        ibd_vertices: Optional[DataFrame] = None,
    ) -> T:
        """Method that will be responsible for creating the graph
        used in the network analysis

        Parameters
        ----------
        ibd_edges : DataFrame
            DataFrame that has the edges of the graph with the length
            of the edges.

        ibd_vertices : Optional[DataFrame]
            DataFrame that has information for each vertice in the
            graph. this value will be none when we are redoing the clustering

        redo : bool
            boolean flag indicating if the user is redoing the
            clustering to make the networks smaller
        """
        if ibd_vertices is not None:
            logger.debug("Generating graph with vertex labels.")
            graph = ig.Graph.DataFrame(
                ibd_edges, directed=False, vertices=ibd_vertices, use_vids=True
            )

        else:
            logger.debug(
                "No vertex metadata provided. Vertex ids will be nonnegative integers"
            )
            graph = ig.Graph.DataFrame(ibd_edges, directed=False)

        return cls(graph)

    def random_walk(self) -> ig.VertexClustering:
        self.logger.debug("performing the random walk")

        ibd_walktrap = ig.Graph.community_walktrap(
            self.graph, weights="cm", steps=self.random_walk_step_size
        )

        random_walk_clusters = ibd_walktrap.as_clustering()

        self.logger.info(random_walk_clusters.summary())

        return random_walk_clusters


@dataclass
class ClusterHandler:
    """Class responsible for performing the cluster on the network objects"""

    minimum_connected_thres: float
    max_network_size: int
    max_rechecks: int
    random_walk_step_size: int
    min_cluster_size: int
    check_times: int = 0
    recheck_clsts: Dict[int, List[int]] = field(default_factory=dict)

    def random_walk(self, network: Networks) -> ig.VertexClustering:
        """Method used to perform the random walk from igraph.community_walktrap

        Parameters
        ----------
        network : Networks
            object that has the ig.Graph

        Returns
        -------
        ig.VertexClustering
            result of the random walk cluster. This object has
            information about clusters and membership
        """
        logger.info(
            f"Using a random walk with a step size of {self.random_walk_step_size} to identify clusters within the provided graph"
        )

        ibd_walktrap = ig.Graph.community_walktrap(
            network.graph, weights="cm", steps=self.random_walk_step_size
        )

        random_walk_clusters = ibd_walktrap.as_clustering()

        logger.info(random_walk_clusters.summary())

        return random_walk_clusters

    def filter_cluster_size(
        self, random_walk_clusters: ig.VertexClustering
    ) -> List[int]:
        """Method to filter networks that are smaller than the min_cluster_size from the analysis

        Parameters
        ----------
        random_walk_clusters : ig.VertexClustering
            result from performing the random walk.

        Returns
        -------
        List[int]
            returns a list of integers where each integer represents
            a cluster id. Each cluster in the network will be >= to
            the min_cluster_size attribute.
        """
        return [
            i
            for i, v in enumerate(random_walk_clusters.sizes())
            if v > self.min_cluster_size
        ]

    def gather_cluster_info(
        self,
        network: Networks,
        cluster_ids: List[int],
        random_walk_clusters: ig.VertexClustering,
    ):
        """Method for getting the information about membership,
        true.positive, false.positives, etc... from the random
        walk

        Parameters
        ----------
        network : Networks
            object that has the graph and the clusters identified based on that graph.

        cluster_ids : List[int]
            list of integers for each cluster id

        random_walk_clusters : ig.VertexClustering
            result from performing the random walk.

        """

        with Progress(transient=True) as progress_bar:
            for clst_id in progress_bar.track(
                cluster_ids, description="clusters processed:", total=len(cluster_ids)
            ):
                network.clusters[clst_id] = {}

                for c, v in enumerate(random_walk_clusters.membership):
                    if v == clst_id:
                        # print(network.graph.vs()[c].attributes())

                        # network.clusters[clst_id].setdefault("memberID", []).append(network.graph.vs[c]["idnum"])

                        network.clusters[clst_id].setdefault("member", []).append(c)

                # getting the total number of edges possible
                theoretical_edge_count = len(
                    list(itertools.combinations(network.clusters[clst_id]["member"], 2))
                )

                cluster_edge_count = len(
                    random_walk_clusters.subgraph(clst_id).get_edgelist()
                )

                network.clusters[clst_id]["true_positive_n"] = cluster_edge_count

                network.clusters[clst_id]["true_positive"] = (
                    cluster_edge_count / theoretical_edge_count
                )
                all_edge = set([])

                for mem in network.clusters[clst_id]["member"]:
                    all_edge = all_edge.union(set(network.graph.incident(mem)))

                network.clusters[clst_id]["false_negative_edge"] = list(
                    all_edge.difference(
                        list(
                            network.graph.get_eids(
                                pairs=list(
                                    itertools.combinations(
                                        network.clusters[clst_id]["member"], 2
                                    )
                                ),
                                directed=False,
                                error=False,
                            )
                        )
                    )
                )

                network.clusters[clst_id]["false_negative_n"] = len(
                    network.clusters[clst_id]["false_negative_edge"]
                )

                if (
                    self.check_times < self.max_rechecks
                    and network.clusters[clst_id]["true_positive"]
                    < self.minimum_connected_thres
                    and len(network.clusters[clst_id]["member"]) > self.max_network_size
                ):
                    self.recheck_clsts.setdefault(self.check_times, []).append(clst_id)
                    # self.recheck_clsts[self.check_times].append(clst_id)
                else:
                    network.final_clusters.append(clst_id)

    def redo_clustering(
        self,
        orig_networks: Networks,
        ibd_pd: DataFrame,
        ibd_vs: DataFrame,
        clst_id: int,
        step: int,
    ) -> None:
        """Method that will redo the clustering, if the
        networks were too large or did not show a high degree
        of connectedness

        Parameters
        ----------
        orig_networks : Networks

        ibd_pd : pd.DataFrame

        clst_id : int

        step : int

        """
        # print(ibd_pd)
        # print(len(orig_networks.clusters[clst_id]["memberID"]))
        # filters for the specific cluster
        redopd = ibd_pd.loc[
            (ibd_pd["idnum1"].isin(orig_networks.clusters[clst_id]["member"]))
            & (ibd_pd["idnum2"].isin(orig_networks.clusters[clst_id]["member"]))
        ]
        redopd.to_csv("redopd_new.txt", sep="\t", index=None)
        print(redopd)

        # filters the vertex for specific members
        # redo_vs = ibd_vs[ibd_vs.idnum.isin(orig_networks.clusters[clst_id]["member"])].reset_index(drop=True)
        # We are going to generate a new Networks object using the redo graph
        # redo_networks = Networks.generate_graph(redopd, redo_vs)
        redo_networks = Networks.generate_graph(redopd.reset_index(drop=True))
        graph_take_2 = ig.Graph.DataFrame(redopd, directed=False)
        print(graph_take_2)
        # print(redo_networks.graph)
        redo_networks.parent_cluster = clst_id

        redo_walktrap_clusters = self.random_walk(redo_networks)

        # If only one cluster is found
        if len(redo_walktrap_clusters.sizes()) == 1:
            # creates an empty dataframe with these columns
            clst_conn = DataFrame(columns=["idnum", "conn", "conn.N", "TP"])
            # iterate over each member id
            for idnum in orig_networks.clusters[clst_id]["member"]:
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
                    (
                        clst_conn["conn.N"]
                        > (0.2 * len(self.clusters[clst_id]["member"]))
                    )
                    & (clst_conn["TP"] < self.minimum_connected_thres)
                    & (
                        clst_conn["conn"]
                        > sorted(clst_conn["conn"], reverse=True)[
                            int(0.01 * len(self.clusters[clst_id]["member"]))
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

        logger.info(redo_walktrap_clusters.summary())

        # if there are more
        allclst = self.filter_cluster_size(redo_walktrap_clusters)

        self.gather_cluster_info(redo_networks, allclst, redo_walktrap_clusters)


def cluster(
    filter_obj: Filter,
    cluster_obj: ClusterHandler,
    min_network_size: int,
    centimorgan_indx: int,
):
    """Main function that will perform the clustering using igraph

    Parameters
    ----------
    filter_obj : Filter
        Filter object that has two attributes for information about
        the edges and information about the vertices.

    cluster_obj : ClusterHandler
        Object that contains information about how the random walk
        needs to be performed. It will use the Networks.graph to
        construct clusters.

    min_network_size : int
        threshold so we can filter networks that are >= the threshold
    """
    filter_obj.ibd_pd = filter_obj.ibd_pd.rename(columns={centimorgan_indx: "cm"})

    ibd_pd = filter_obj.ibd_pd.loc[:, ["idnum1", "idnum2", "cm"]]
    print(ibd_pd)
    ibd_vs = filter_obj.ibd_vs.reset_index(drop=True)

    # Generate the first pass networks
    network = Networks.generate_graph(
        ibd_pd,
        ibd_vs,
    )

    random_walk_results = cluster_obj.random_walk(network)

    allclst = [
        i for i, v in enumerate(random_walk_results.sizes()) if v > min_network_size
    ]

    cluster_obj.gather_cluster_info(network, allclst, random_walk_results)

    while (
        cluster_obj.check_times < cluster_obj.max_rechecks
        and len(cluster_obj.recheck_clsts[cluster_obj.check_times]) > 0
    ):
        cluster_obj.check_times += 1

        logger.info(f"recheck: {cluster_obj.check_times}")

        _ = cluster_obj.recheck_clsts.setdefault(cluster_obj.check_times, [])

        for clst in cluster_obj.recheck_clsts.get(cluster_obj.check_times - 1):
            cluster_obj.redo_clustering(
                network,
                ibd_pd,
                ibd_vs,
                clst,
                cluster_obj.random_walk_step_size,
            )
