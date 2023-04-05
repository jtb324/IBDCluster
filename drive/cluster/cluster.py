from dataclasses import dataclass
from typing import Any, Dict, List

import igraph as ig
from filter import Filter


@dataclass
class Networks:
    graph: Any

    def random_walk(self, step_size: int) -> ig.VertexClustering:
        ibd_walktrap = ig.Graph.community_walktrap(
            self.graph, weights="cm", steps=step_size
        )

        random_walk_clusters = ibd_walktrap.as_clustering()

        print(random_walk_clusters.summary())

        return random_walk_clusters


def cluster(
    filter_obj: Filter,
):
    """Main function that will perform the clustering using igraph"""

    ibd_graph = ig.Graph.DataFrame(
        filter_obj.ibd_pd, directed=False, vertices=filter_obj.ibd_vs, use_vids=True
    )
