from dataclasses import dataclass, field
import itertools
import logging
from typing import Any, Dict, List, Optional
from pandas import DataFrame
from rich.progress import Progress

import igraph as ig
from models import Filter
import log


@dataclass
class Networks:
    minimum_connected_thres: float
    max_network_size: int
    max_rechecks: int
    random_walk_step_size: int
    graph: Optional[Any] = None
    redo_graph: Optional[Any] = None
    clusters: Dict[int, Any] = field(default_factory=dict)
    recheck_clsts: dict[int, List[int]] = field(default_factory=dict)
    outside_of_clusters: List[int] = field(default_factory=list)
    check_times: int = 0
    logger: logging.Logger = log.get_logger(__name__)

    def generate_graph(
        self,
        ibd_edges: DataFrame,
        ibd_vertices: Optional[DataFrame] = None,
        redo: bool = False,
    ) -> None:
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
        if redo:
            self.logger.debug("creating the dataframe for redo-ing the cluster")
            self.redo_graph = ig.Graph.DataFrame(ibd_edges, directed=False)
        else:
            self.logger.debug("Generating the original graph")
            self.graph = ig.Graph.DataFrame(
                ibd_edges, directed=False, vertices=ibd_vertices, use_vids=True
            )

    def random_walk(self) -> ig.VertexClustering:
        self.logger.debug("performing the random walk")

        ibd_walktrap = ig.Graph.community_walktrap(
            self.graph, weights="cm", steps=self.random_walk_step_size
        )

        random_walk_clusters = ibd_walktrap.as_clustering()

        self.logger.info(random_walk_clusters.summary())

        return random_walk_clusters

    def gather_cluster_info(
        self, cluster_ids: List[int], random_walk_clusters: ig.VertexClustering
    ):
        """Method for getting the information about membership,
        true.positive, false.positives, etc... from the random
        walk

        Parameters
        ----------
        cluster_ids : List[int]
            list of integers for each cluster id

        random_walk_clusters : ig.VertexClustering

        """
        with Progress(transient=True) as progress_bar:
            for clst_id in progress_bar.track(
                cluster_ids, description="clusters_processed:", total=len(cluster_ids)
            ):
                # print(self.graph)
                # print(self.graph.vs)
                # print(self.graph.vs.attribute_names())
                # print(self.graph.vs)
                # print(self.graph.vs()[1])
                self.clusters[clst_id] = {}
                self.clusters[clst_id]["memberID"] = [
                    self.graph.vs[c]["idnum"]
                    for c, v in enumerate(random_walk_clusters.membership)
                    if v == clst_id
                ]

                self.clusters[clst_id]["member"] = [
                    c
                    for c, v in enumerate(random_walk_clusters.membership)
                    if v == clst_id
                ]
                clst_edge_n = (
                    len(self.clusters[clst_id]["member"])
                    * (len(self.clusters[clst_id]["member"]) - 1)
                    / 2
                )
                self.clusters[clst_id]["true_positive_edge"] = list(
                    filter(
                        lambda a: a != -1,
                        self.graph.get_eids(
                            pairs=list(
                                itertools.combinations(
                                    self.clusters[clst_id]["member"], 2
                                )
                            ),
                            directed=False,
                            error=False,
                        ),
                    )
                )
                self.clusters[clst_id]["true_positive_n"] = len(
                    self.clusters[clst_id]["true_positive_edge"]
                )
                self.clusters[clst_id]["true_positive"] = (
                    self.clusters[clst_id]["true_positive_n"] / clst_edge_n
                )
                all_edge = set([])
                for mem in self.clusters[clst_id]["member"]:
                    all_edge = all_edge.union(set(self.graph.incident(mem)))

                self.clusters[clst_id]["false_negative_edge"] = list(
                    all_edge.difference(
                        list(
                            self.graph.get_eids(
                                pairs=list(
                                    itertools.combinations(
                                        self.clusters[clst_id]["member"], 2
                                    )
                                ),
                                directed=False,
                                error=False,
                            )
                        )
                    )
                )

                self.clusters[clst_id]["false_negative_n"] = len(
                    self.clusters[clst_id]["false_negative_edge"]
                )

                if (
                    self.check_times < self.max_rechecks
                    and self.clusters[clst_id]["true_positive"]
                    < self.minimum_connected_thres
                    and len(self.clusters[clst_id]["member"]) > self.max_network_size
                ):
                    self.recheck_clsts.setdefault(self.check_times, []).append(clst_id)
                    # self.recheck_clsts[self.check_times].append(clst_id)
                else:
                    self.outside_of_clusters.append(clst_id)


def cluster(
    filter_obj: Filter, network: Networks, min_network_size: int, centimorgan_indx: int
):
    """Main function that will perform the clustering using igraph

    Parameters
    ----------
    filter_obj : Filter
        Filter object that has two attributes for information about the edges and information about the vertices.

    network : Networks
        class that has information about the networks identified in the random walk

    min_network_size : int
        threshold so we can filter networks that are >= the threshold
    """
    filter_obj.ibd_pd = filter_obj.ibd_pd.rename(columns={centimorgan_indx: "cm"})

    network.generate_graph(
        filter_obj.ibd_pd.loc[:, ["idnum1", "idnum2", "cm"]],
        filter_obj.ibd_vs.reset_index(drop=True),
    )

    random_walk_results = network.random_walk()

    allclst = [
        i for i, v in enumerate(random_walk_results.sizes()) if v > min_network_size
    ]

    network.gather_cluster_info(allclst, random_walk_results)
