import itertools
import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import igraph as ig
import log
from models import Filter
from pandas import DataFrame
from rich.progress import Progress

# creating a logger
logger: logging.Logger = log.get_logger(__name__)

# Create a generic variable that can represents the class from the


@dataclass
class Network:
    clst_id: int
    true_positive_count: int
    true_positive_percent: float
    false_negative_edges: List[int]
    false_negative_count: int
    members: List[int] = field(default_factory=list)


@dataclass
class ClusterHandler:
    """Class responsible for performing the cluster on the network objects"""

    minimum_connected_thres: float
    max_network_size: int
    max_rechecks: int
    random_walk_step_size: int
    min_cluster_size: int
    check_times: int = 0
    recheck_clsts: Dict[int, List[Network]] = field(default_factory=dict)
    final_clusters: List[int] = field(default_factory=list)

    @staticmethod
    def generate_graph(
        ibd_edges: DataFrame,
        ibd_vertices: Optional[DataFrame] = None,
    ) -> ig.Graph:
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
        """
        if ibd_vertices is not None:
            logger.debug("Generating graph with vertex labels.")
            return ig.Graph.DataFrame(
                ibd_edges, directed=False, vertices=ibd_vertices, use_vids=True
            )

        else:
            logger.debug(
                "No vertex metadata provided. Vertex ids will be nonnegative integers"
            )
            return ig.Graph.DataFrame(ibd_edges, directed=False, use_vids=False)

    def random_walk(self, graph: ig.Graph) -> ig.VertexClustering:
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
            graph, weights="cm", steps=self.random_walk_step_size
        )

        random_walk_clusters = ibd_walktrap.as_clustering()

        logger.info(random_walk_clusters.summary())

        return random_walk_clusters

    def filter_cluster_size(
        self, random_walk_clusters_sizes: ig.VertexClustering
    ) -> List[int]:
        """Method to filter networks that are smaller than the min_cluster_size from the analysis

        Parameters
        ----------
        random_walk_clusters_sizes : List[int]
            size of each cluster from the random walk results

        Returns
        -------
        List[int]
            returns a list of integers where each integer
            represents a cluster id. Each cluster in the
            network will be >= to the min_cluster_size
            attribute.
        """
        return [
            i
            for i, v in enumerate(random_walk_clusters_sizes)
            if v > self.min_cluster_size
        ]

    @staticmethod
    def _gather_members(random_walk_members: List[int], clst_id: int) -> List[int]:
        """Generate a list of individuals in the network

        Parameters
        ----------
        random_walk_members : List[int]
            list of all members from the random walk results

        clst_id : int
            id for the original cluster

        Returns
        -------
        List[int]
            returns a list of ids of individuals in the network
        """
        member_list = []

        for c, v in enumerate(random_walk_members):
            if v == clst_id:
                member_list.append(c)

        return member_list

    @staticmethod
    def _determine_true_positive_edges(
        member_list: List[int], clst_id: int, random_walk_results: ig.VertexClustering
    ) -> Tuple[int, float]:
        """determining the number of true positive edges

        Parameters
        ----------
        member_list : List[int]
            list of ids within the specific network

        clst_id : int
            id for the original cluster

        random_walk_results : ig.VertexClustering
            vertexClustering object returned after the random
            walk that has the different clusters

        Returns
        -------
        Tuple[int, float]
            returns a tuple where the first element is the
            number of edges in the graph and the second
            element is the ratio of actually edges in the
            graph compared to the theoretical maximum number
            of edges in the graph.
        """
        # getting the total number of edges possible
        theoretical_edge_count = len(list(itertools.combinations(member_list, 2)))
        # Getting the number of edges within the graph and saving it
        # as a dictionary key, 'true_positive_n'
        cluster_edge_count = len(random_walk_results.subgraph(clst_id).get_edgelist())

        return cluster_edge_count, cluster_edge_count / theoretical_edge_count

    def _determine_false_positive_edges(
        graph: ig.Graph, member_list: List[int]
    ) -> Tuple[int, List[int]]:
        """determine the number of false positive edges

        Parameters
        ----------
        graph : ig.Graph
            graph object returned from ig.Graph.DataFrame

        member_list : List[int]
            list of ids within the specific network

        Returns
        -------
        Tuple[int, List[int]]
            returns a tuple where the first element is the
            number of edges in the graph and the second
            element is a list of false positive edges.
        """
        all_edge = set([])

        for mem in member_list:
            all_edge = all_edge.union(set(graph.incident(mem)))

        false_negative_edges = list(
            all_edge.difference(
                list(
                    graph.get_eids(
                        pairs=list(itertools.combinations(member_list, 2)),
                        directed=False,
                        error=False,
                    )
                )
            )
        )
        return len(false_negative_edges), false_negative_edges

    def gather_cluster_info(
        self,
        graph: ig.Graph,
        cluster_ids: List[int],
        random_walk_clusters: ig.VertexClustering,
        parent_cluster_id: Optional[int] = None,
    ) -> Network:
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

        parent_cluster_id : int
            id of the original cluster that is now being broken up. Child cluster ids will take the form parent_id.child_id
        """
        print(cluster_ids)
        for clst_id in cluster_ids:
            if parent_cluster_id:
                clst_name = f"{parent_cluster_id}.{clst_id}"
            else:
                clst_name = clst_id
            print(clst_name)

            network.clusters[clst_name] = {}

            member_list = ClusterHandler._gather_members(
                random_walk_clusters.membership, clst_id
            )

            (
                true_pos_count,
                true_pos_ratio,
            ) = ClusterHandler._determine_true_positive_edges(
                member_list, clst_id, random_walk_clusters
            )

            (
                false_neg_count,
                false_neg_list,
            ) = ClusterHandler._determine_false_positive_edges(graph, member_list)

            network = Network(
                clst_id,
                true_pos_count,
                true_pos_ratio,
                false_neg_list,
                false_neg_count,
                member_list,
            )

            if (
                self.check_times < self.max_rechecks
                and network.true_positive_percent < self.minimum_connected_thres
                and len(network.members) > self.max_network_size
            ):
                self.recheck_clsts.setdefault(self.check_times, []).append(network)

            else:
                network.final_clusters.append(clst_name)
            print(network.clusters[clst_name])

    def redo_clustering(
        self,
        ibd_pd: DataFrame,
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
            Id  from the original clustering that is being
            broken up

        step : int

        """
        # filters for the specific cluster
        redopd = ibd_pd[
            (ibd_pd["idnum1"].isin(orig_networks.clusters[clst_id]["member"]))
            & (ibd_pd["idnum2"].isin(orig_networks.clusters[clst_id]["member"]))
        ]

        # filters the vertex for specific members
        # redo_vs = ibd_vs[ibd_vs.idnum.isin(orig_networks.clusters[clst_id]["member"])].reset_index(drop=True)
        # We are going to generate a new Networks object using the redo graph
        # redo_networks = Networks.generate_graph(redopd, redo_vs)
        redo_networks = Networks.generate_graph(redopd)
        # print(redo_networks.graph)
        # redo_networks.parent_cluster = clst_id

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

        # if there are more
        allclst = self.filter_cluster_size(redo_walktrap_clusters)

        self.gather_cluster_info(
            redo_networks, allclst, redo_walktrap_clusters, clst_id
        )


def cluster(
    filter_obj: Filter,
    cluster_obj: ClusterHandler,
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
    # filtering the edges dataframe to the correct columns
    ibd_pd = filter_obj.ibd_pd.loc[:, ["idnum1", "idnum2", "cm"]]

    ibd_vs = filter_obj.ibd_vs.reset_index(drop=True)

    # Generate the first pass networks
    network_graph = cluster_obj.generate_graph(
        ibd_pd,
        ibd_vs,
    )

    random_walk_results = cluster_obj.random_walk(network_graph)

    allclst = cluster_obj.filter_cluster_size(random_walk_results)

    cluster_obj.gather_cluster_info(network_graph, allclst, random_walk_results)

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
                clst,
                cluster_obj.random_walk_step_size,
            )
