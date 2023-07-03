.. raw:: html

    <style> .yellow {color:yellow; font-weight:bold;} </style>

.. role:: yellow

Output File Structure
=====================

DRIVE outputs a file that ends in .DRIVE.txt. This file has the results of the clustering analysis. The file begins with two lines that tell how many IBD segments and haplotypes were identified as well as how many clusters were identified. There are eight columns in the output file are described below

Column descriptions:
--------------------
* :yellow:`clustID` ID given to each network identified. This value will have the form "clst#".

----

* :yellow:`n.total` Total number of individuals in the network.

----

* :yellow:`n.haplotype` The number of haplotypes in the network. This value may be different than n.total due to inbreeding.

----

* :yellow:`true.positive.n` Number of shared IBD segments that are identified in the network. 

----

* :yellow:`true.postive` Proportion of identified IBD segments in networks vs the total number of possible IBD segments that could exist between all individuals in the network.

----

* :yellow:`false.postive` Proportion of individuals within the cluster that share an IBD segment with another individual outside of the cluster.

----

* :yellow:`IDs` List of ids that are in the network. 

----

* :yellow:`ID.haplotype` List of haplotypes that are in the network. These will be equivalent to the ids in the "IDs" column except each id will have a phase value attached to it.