import itertools
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

import igraph as ig
from pandas import DataFrame

from drive.log import CustomLogger
from drive.models import Filter, Network, Network_Interface

# creating a logger
logger: logging.Logger = CustomLogger.get_logger(__name__)

# Create a generic variable that can represents the class from the


@dataclass
class ClusterHandler:
    """Class responsible for performing the cluster on the network objects"""

    minimum_connected_thres: float
    max_network_size: int
    max_rechecks: int
    random_walk_step_size: int
    min_cluster_size: int
    segment_dist_threshold: int
    hub_threshold: float
    haplotype_mappings: Dict[str, int]
    recluster: bool
    check_times: int = 0
    recheck_clsts: Dict[int, List[Network_Interface]] = field(default_factory=dict)
    final_clusters: List[Network_Interface] = field(default_factory=list)

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
            graph. This value will be none when we are redoing the clustering
        """
        if ibd_vertices is not None:
            logger.debug("Generating graph with vertex labels.")
            return ig.Graph.DataFrame(
                ibd_edges, directed=False, vertices=ibd_vertices, use_vids=False
            )

        else:
            logger.debug(
                "No vertex metadata provided. Vertex ids will be nonnegative integers"
            )
            # return ig.Graph.DataFrame(ibd_edges, directed=False, use_vids=False)
            return ig.Graph.DataFrame(ibd_edges, directed=False)

    def random_walk(self, graph: ig.Graph) -> ig.VertexClustering:
        """Method used to perform the random walk from igraph.community_walktrap

        Parameters
        ----------
        graph : ig.Graph
            graph object created by ig.Graph.DataFrame

        Returns
        -------
        ig.VertexClustering
            result of the random walk cluster. This object has
            information about clusters and membership
        """
        logger.debug(
            f"Performing the community walktrap algorithm with a random walk step size of: {self.random_walk_step_size}"
        )

        ibd_walktrap = ig.Graph.community_walktrap(
            graph, weights="cm", steps=self.random_walk_step_size
        )

        random_walk_clusters = ibd_walktrap.as_clustering()

        logger.verbose(random_walk_clusters.summary())

        return random_walk_clusters

    def filter_cluster_size(self, random_walk_clusters_sizes: List[int]) -> List[int]:
        """Method to filter networks that are smaller than the min_cluster_size from
        the analysis

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
    def _gather_members(
        random_walk_members: List[int], clst_id: int, graph: ig.Graph
    ) -> Tuple[List[int], List[int]]:
        """Generate a list of individuals ids in the network

        Parameters
        ----------
        random_walk_members : List[int]
            list of all members from the random walk results

        clst_id : int
            id for the cluster. This value is used to pull out members
            that belong to the cluster

        graph : ig.Graph
            Graph object formed by ig.Graph.DataFrame

        Returns
        -------
        Tuple[List[int], List[int]]
            returns a list of ids of individuals in the network and then
            a list of vertex ids. The individual ids are just the index
            of the element in the membership list and the vertex ids are
            the list of ids provided by nme label in the vs() property.
        """
        member_list = []
        # this list has the ids. It is sometimes the same as the
        # member_list but it will not be the same in the redo_networks
        # graph
        vertex_ids = []

        for member_id, assigned_clst_id in enumerate(random_walk_members):
            if assigned_clst_id == clst_id:
                member_list.append(graph.vs()[member_id]["name"])
                vertex_ids.append(member_id)

        return member_list, vertex_ids

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
        graph: ig.Graph, vertex_list: List[int]
    ) -> Tuple[int, List[int]]:
        """determine the number of false positive edges

        Parameters
        ----------
        graph : ig.Graph
            graph object returned from ig.Graph.DataFrame

        vertex_list : List[int]
            list of vertex ids within the specific network

        Returns
        -------
        Tuple[int, List[int]]
            returns a tuple where the first element is the
            number of edges in the graph and the second
            element is a list of false positive edges.
        """
        all_edge = set([])

        for mem in vertex_list:
            all_edge = all_edge.union(set(graph.incident(mem)))

        false_negative_edges = list(
            all_edge.difference(
                list(
                    graph.get_eids(
                        pairs=list(itertools.combinations(vertex_list, 2)),
                        directed=False,
                        error=False,
                    )
                )
            )
        )
        return len(false_negative_edges), false_negative_edges

    def _map_ids_back_to_haplotypes(
        self, members: List[int]
    ) -> Tuple[List[str], Set[str]]:
        """remap the haplotype integer ids back to the haplotype
        strings for the output file

        Parameters
        ----------
        members : List[int]
            list of integers that represent the name of each vertex in the network
            corresponding the the graph

        Returns
        -------
        Tuple[List[str], List[str]]
            returns a list of haplotype id strings in the network. These
            strings include the phase number. Also returns a list of ids
            within the networks. This will be the same as the haplotype
            id without the phase number
        """
        haplotypes = [self.haplotype_mappings[value] for value in members]

        member_ids = {value[:-2] for value in haplotypes}

        return haplotypes, member_ids

    def gather_cluster_info(
        self,
        graph: ig.Graph,
        cluster_ids: List[int],
        random_walk_clusters: ig.VertexClustering,
        parent_cluster_id: Optional[str] = None,
    ) -> None:
        """Method for getting the information about membership,
        true.positive, false.positives, etc... from the random
        walk

        Parameters
        ----------
        graph : ig.Graph
            Graph object generated by ig.Graph.DataFrame

        cluster_ids : List[int]
            list of integers for each cluster id

        random_walk_clusters : ig.VertexClustering
            result from performing the random walk.

        parent_cluster_id : Optional[str]
            id of the original cluster that is now being broken up. Child
            cluster ids will take the form parent_id.child_id
        """

        for clst_id in cluster_ids:
            # We need to form the appropriate id if the cluster has a
            # parent otherwise they get the value of the clst_id argument
            if parent_cluster_id:
                clst_name = f"{parent_cluster_id}.{clst_id}"
            else:
                clst_name = f"{clst_id}"

            # We are going to get the vertex id and member id of each
            # graph
            member_list, vertex_ids = ClusterHandler._gather_members(
                random_walk_clusters.membership, clst_id, graph
            )

            # Next we get the number of edges/ ratio of actual edges to
            # the potential edges
            (
                true_pos_count,
                true_pos_ratio,
            ) = ClusterHandler._determine_true_positive_edges(
                member_list, clst_id, random_walk_clusters
            )
            # next we determine the number of false positive edges
            (
                false_neg_count,
                false_neg_list,
            ) = ClusterHandler._determine_false_positive_edges(graph, vertex_ids)

            # If the graph is too sparse and it is too large and the max
            # number of rechecks has not been reached then we will put
            # the network into a recluster dictionary. Otherwise it is
            # added to the final_clst list
            if (
                self.check_times < self.max_rechecks
                and true_pos_ratio < self.minimum_connected_thres
                and len(member_list) > self.max_network_size
                and self.recluster
            ):
                # We can put all of this information into a network class. Here the
                # member list will still be in integers
                network = Network(
                    clst_name,
                    true_pos_count,
                    true_pos_ratio,
                    false_neg_list,
                    false_neg_count,
                    member_list,
                    vertex_ids,
                )
                # debug statement if we want to see the members and the haplotypes
                logger.debug(f"members: {member_list}\nhaplotypes: {vertex_ids}")

                # we are going to append the network to the list of networks that needs
                # to be rechecked
                self.recheck_clsts.setdefault(self.check_times, []).append(network)

            else:
                # we need to convert the integer ids back to strings
                haplotype_ids, member_ids = self._map_ids_back_to_haplotypes(
                    member_list
                )

                network = Network(
                    clst_name,
                    true_pos_count,
                    true_pos_ratio,
                    false_neg_list,
                    false_neg_count,
                    member_ids,
                    haplotype_ids,
                )
                # logger.info(f"members: {member_ids}\nhaplotypes: {haplotype_ids}")

                self.final_clusters.append(network)

    def redo_clustering(
        self, network: Network_Interface, ibd_pd: DataFrame, ibd_vs: DataFrame
    ) -> None:
        """Method that will redo the clustering, if the
        networks were too large or did not show a high degree
        of connectedness

        Parameters
        ----------
        network : Network_InterFace
            object that represents each cluster. These objects have information
            about the cluster id, number and ratio of edges, true_positive_percent,
            false_negative_edges, false_negative_count

        ibd_pd : pd.DataFrame
            DataFrame that has information about the edges that a pair shares
        """
        # pulling the id from the original cluster
        original_id = network.clst_id

        # filters for the specific cluster
        redopd = ibd_pd[
            (ibd_pd["idnum1"].isin(network.haplotypes))
            & (ibd_pd["idnum2"].isin(network.haplotypes))
        ]

        redo_vs = ibd_vs[ibd_vs.idnum.isin(network.haplotypes)]

        # If the redopd or redo_vs is empty it causes strange behavior and the code will
        # usually fail. The desired behavior is for the program to tell teh user that
        # the graph could not be constructed and then for it to move on.
        if not redopd.empty and not redo_vs.empty:
            # We are going to generate a new Networks object using the redo graph
            redo_networks = ClusterHandler.generate_graph(redopd, redo_vs)
            # redo_networks = ClusterHandler.generate_graph(redopd)
            # performing the random walk
            redo_walktrap_clusters = self.random_walk(redo_networks)

            # If only one cluster is found
            if len(redo_walktrap_clusters.sizes()) == 1:
                # creates an empty dataframe with these columns
                clst_conn = DataFrame(columns=["idnum", "conn", "conn.N", "TP"])
                # iterate over each member id
                # for idnum in network.haplotypes:
                for idnum in network.members:
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
                    # assert 1 == 0
                    if len(conn_idnum) == 1:
                        connTP = 1
                    else:
                        try:
                            connTP = conn_tp / (
                                len(conn_idnum) * (len(conn_idnum) - 1) / 2
                            )
                        except ZeroDivisionError:
                            raise ZeroDivisionError(
                                f"There was a zero division error encountered when looking at the network with the id {idnum}"  # noqa: E501
                            )  # noqa: E501

                    clst_conn.loc[idnum] = [idnum, conn, len(conn_idnum), connTP]
                rmID = list(
                    clst_conn.loc[
                        (
                            clst_conn["conn.N"]
                            > (self.segment_dist_threshold * len(network.members))
                        )
                        & (clst_conn["TP"] < self.minimum_connected_thres)
                        & (
                            clst_conn["conn"]
                            > sorted(clst_conn["conn"], reverse=True)[
                                int(self.hub_threshold * len(network.members))
                            ]
                        )
                    ]["idnum"]
                )

                redopd = redopd.loc[
                    (~redopd["idnum1"].isin(rmID)) & (~redopd["idnum2"].isin(rmID))
                ]
                redo_graph = self.generate_graph(
                    redopd,
                )
                # redo_g = ig.Graph.DataFrame(redopd, directed=False)
                redo_walktrap_clusters = self.random_walk(redo_graph)
                # redo_walktrap = ig.Graph.community_walktrap(
                #     redo_g, weights="cm", steps=self.random_walk_step_size
                # )
                # redo_walktrap_clusters = redo_walktrap.as_clustering()

            # Filter to the clusters that are llarger than the minimum size
            allclst = self.filter_cluster_size(redo_walktrap_clusters.sizes())

            self.gather_cluster_info(
                redo_networks, allclst, redo_walktrap_clusters, original_id
            )
        else:
            logger.debug(
                f"A graph was not able to be generated when we attempted to recluster the network: {original_id}. This error probably indicates that there were There were none of the {len(network.haplotypes)} individuals in that specific network that shared ibd segments with one another."
            )


def cluster(
    filter_obj: Filter,
    cluster_obj: ClusterHandler,
    centimorgan_indx: int,
) -> List[Network_Interface]:
    """Main function that will perform the clustering using igraph

    Parameters
    ----------
    filter_obj : Filter
        Filter object that has two attributes: ibd_pd and ibd_vs. These
        attributes are two dataframes that have information about the
        edges and information about the vertices.

    cluster_obj : ClusterHandler
        Object that contains information about how the random walk
        needs to be performed. It will use the construct networks and return those
        values in a list.

    min_network_size : int
        threshold so we can filter networks that are >= the threshold

        Returns
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

    allclst = cluster_obj.filter_cluster_size(random_walk_results.sizes())

    cluster_obj.gather_cluster_info(network_graph, allclst, random_walk_results)

    while (
        cluster_obj.check_times < cluster_obj.max_rechecks
        and len(cluster_obj.recheck_clsts.get(cluster_obj.check_times, [])) > 0
    ):
        cluster_obj.check_times += 1
        logger.verbose(f"recheck: {cluster_obj.check_times}")

        _ = cluster_obj.recheck_clsts.setdefault(cluster_obj.check_times, [])

        for network in cluster_obj.recheck_clsts.get(cluster_obj.check_times - 1):
            cluster_obj.redo_clustering(network, ibd_pd, ibd_vs)
    # logginng the number of segments, haplotypes, and clusters
    # identified in the analysis
    logger.info(
        f"Identified {network_graph.ecount()} IBD segments from {network_graph.vcount()} haplotypes"  # noqa: E501
    )

    logger.info(f"Identified {len(cluster_obj.final_clusters)} IBD clusters")

    return cluster_obj.final_clusters
