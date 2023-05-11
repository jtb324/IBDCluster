DRIVE Output
============

DRIVE currently outputs two files: a file providing a network level description of the clustering results with the suffix ".DRIVE.txt" and a log file that has the suffix ".log". This  file has the results from the clustering analysis with information such as number of members, members ids, how connected the graph is internally, and the binomial test statistics. There are eight columns in the output file like the table below. 

-------------

Column descriptions:
--------------------

* **clustID**: ID given to each network identified. This value will have the form "clst#".

---
* **n.total**: Total number of individuals in the network.

---
* **n.haplotype**: The number of haplotypes in the network. This value may be different than n.total due to inbreeding.

---
* **true.positive.n**: Number of shared IBD segments that are identified in the network. 

---
* **true.postive**: Proportion of identified IBD segments in networks vs the total number of possible IBD segments that could exist between all individuals in the network.

---
* **false.postive**: Proportion of individuals within the cluster that share an IBD segment with another individual outside of the cluster.

---
* **IDs**: List of ids that are in the network. 

---
* **ID.haplotype**: List of haplotypes that are in the network. These will be equivalent to the ids in the "IDs" column except each id will have a phase value attached to it.