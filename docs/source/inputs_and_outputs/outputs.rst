DRIVE Output
============

DRIVE currently outputs two files. These files are described below.

Networks File
-------------
DRIVE creates a file with the suffix ".drive_networks.txt". This file has the results from the clustering analysis with information such as number of members, members ids, haplotype ids how connected the graph is internally, and the binomial test statistics. This file has at a minimum of 11 columns depending on where the user provides a phenotype file or not. These columns, plus the possible additional columns are described below.

Column descriptions:
^^^^^^^^^^^^^^^^^^^^

* **clustID**: ID given to each network identified. This value will have the form "clst#".

----

* **n.total**: Total number of individuals in the network.

----

* **n.haplotype**: The number of haplotypes in the network. This value may be different than n.total due to inbreeding.

----

* **true.positive.n**: Number of shared IBD segments that are identified in the network. 

----

* **true.postive**: Proportion of identified IBD segments in networks vs the total number of possible IBD segments that could exist between all individuals in the network.

----

* **false.postive**: Proportion of individuals within the cluster that share an IBD segment with another individual outside of the cluster.

----

* **IDs**: List of ids that are in the network. 

----

* **ID.haplotype**: List of haplotypes that are in the network. These will be equivalent to the ids in the "IDs" column except each id will have a phase value attached to it.

----

* **min_pvalue**: Value of the smallest pvalue calculated for the network from the binomial test. If a phenotype file is not provided then this value will be N/A.

----

* **min_phenotype**: Name of the phenotype that corresponds to the smallest pvalue. This value will also be N/A if a phenotype file is not provided.

----

* **min_phenotype_description**: Description of what the phenotype is. This value will be N/A if a descriptions file is not provided, if the phenotype doesn't have a description, or if a phenotype file is not provided.

----

* **\*_cases_in_network**: Number of individuals in the network that are affected by the phenotype. 

----

* **\*_excluded_in_network**: Number of individuals in the network that are excluded from the statistics analysis. 

----

* **\*_pvalue**: Pvalue determined for the specific phenotype

.. note::

    The final three columns, "\*_cases_in_network, \*_excluded_in_network, \*_pvalue" are only created if the user provides a case file, otherwise the output file will only have the first 11 columns. If the user provides a case file then these three columns will be created for each phenotype so if you provided 3 phenotypes then 9 columns would be added to the output file.


Log File
--------

DRIVE creates a log file with whatever name the user provides. This file has the suffix ".log". This file has information about the arguments the user passed and then runtime information from the program such as how many networks were identified and how many haplotypes were identified. The amount of information written to this file will vary depending on what level of verbosity the user chooses.