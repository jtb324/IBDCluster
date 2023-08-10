Pre-Installed Plugins 
=====================


"pvalues" plugin
----------------
This plugin performs the binomial test to determine enrichment of phenotypes of interest within the networks identified by DRIVE's clustering algorithm. This plugin updates two attributes of the network class: "min_pvalue_str" and "pvalues". The plugin will construct a string that has the most statistically enriched phenotype id, the corresponding pvalue, and the phenotype description if it is provided. The plugin also generates a dictionary where the keys are the phenotype ids and the values are strings that describe how the number of carriers in the network, the ids of the carriers in the network, how many individuals are excluded from the network, and then the pvalue calculated from the binomial test.


"network writer" plugin
-----------------------
This plugin is responsible for writing all the output about what clusters were detected and if those networks are enriched for a phenotype of interest to a file. This plugin will create a file with the extension "\*drive.networks.txt". This plugin has to be run after the pvalues plugin because it uses some of the results of the pvalues plugin.

