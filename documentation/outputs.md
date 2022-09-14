---
layout: default 
title: Outputs
parent: Documentation
nav_order: 3
---
# Output Files:
---
IBDCluster produces three files using the stock plugins. The Networks file and the Allpairs file will be in a subdirectory that is named for whatever Gene you are interested in. This subdirectory will be in the user specified output directory. The phenotype prevalence file will be in the main output directory.

**Networks.txt File:**

 The first file is the \*networks.txt file. This file gives information on a per network level. Each row is a different network identified in the analysis. There are 11 columns that will be in each file. These first 11 columns are shown below in the table and then they are described below.

| program |    gene     | network_id | chromosome | IIDs_count | haplotypes_count |       IIDs       |     haplotypes      | min_pvalue | min_pvalue_phecode |    min_phecode_desc    |
|:--------|:------------|:-----------|:-----------|:-----------|:-----------------|:-----------------|:--------------------|:-----------|:-------------------|:-----------------------|
| hapibd  |  gene name  |      #     |      #     |      #     |         #        | grids in network |haplotypes in network|      #     |    phecode label   | description of phecode |


* <span style="color: #F0FF00">**program:**</span> This will have the name of the IBD program that was used to identify the IBD segments. This value currently will be hapibd because it is the only supported program.

* <span style="color: #F0FF00">**gene:**</span> This value will be the gene that you are interested in for the analysis.

* <span style="color: #F0FF00">**network_id:**</span> This is the number of the network.

* <span style="color: #F0FF00">**chromosome:**</span> This is the chromosome number that the gene is on.

* <span style="color: #F0FF00">**IIDs_count:**</span> This is the number of different individuals in the network. This number does not take into account haplotypes.

* <span style="color: #F0FF00">**haplotypes_count:**</span> This is the number of haplotypes in the network. This value will not be higher than the IIDs_count but could be lower if there is inbreeding in the population.

* <span style="color: #F0FF00">**IIDs:**</span> This value is a list of the individual IDs in the network. 

* <span style="color: #F0FF00">**haplotypes:**</span> This value is a list of the haplotypes in the network. Each haplotype consider of the individual ID concatenated with the phase number.

* <span style="color: #F0FF00">**min_pvalue:**</span> This value is the lowest pvalue found during the binomial test.

* <span style="color: #F0FF00">**min_pvalue_phecode:**</span> This value is the PheCode label for that corresponds to the lowest pvalue.

* <span style="color: #F0FF00">**min_phecode_desc:**</span> This value is the description of the PheCode that corresponds to the lowest pvalue.

Every other column in the file depends on what phenotypes the user supplies. For each phenotype provided, two columns will be formed: a pvalue columns and an ind_in_network columns. The pvalue column will have the pvalue calculated during the binomial test and the ind_in_network columns will have the number of individuals in the network affected by the phenotype.

**Allpairs.txt File:**

The second file produced will an *allpairs.txt file. This file gives information about the segment shared between each pair in a network. There are eleven columns shown in the table below and then each column is described in more detail after the table.

| program |  network_id  | pair_1 | pair_2 |    phase_1    |      phase_2     |  chromosome  | gene_name |     start     |      end     |       length        |
|:--------|:-------------|:-------|:-------|:--------------|:-----------------|:-------------|:----------|:--------------|:-------------|:--------------------|
| hapibd  |      1       | grid 1 | grid 2 | id1 with phase|  id2 with phase  |       1      |    name   | start position| end position |  length of segment  |

* <span style="color: #F0FF00">**program:**</span> This will have the name of the IBD program that was used to identify the IBD segments. This value currently will be hapibd because it is the only supported program.

* <span style="color: #F0FF00">**network_id:**</span> This value has the network number that the pair is in. This value starts at one and goes to how ever many networks were identified

* <span style="color: #F0FF00">**pair_1:**</span> This value will have the id of the first individual in the pair.

* <span style="color: #F0FF00">**pair_2:**</span> This value will have the id of the second individual in the pair

* <span style="color: #F0FF00">**phase_1:**</span> This value will have the haplotype phase for the first individual in the pair. This value will be the pair_1 id plus the haplotype number attached to it (This value is 1 or 2 for hapibd). An example of this is R12341234.1

* <span style="color: #F0FF00">**phase_2:**</span> This value will have the haplotype phase for the second individual. It follows the same format as the phase_1 value

* <span style="color: #F0FF00">**chromosome:**</span> This value has the chromosome number that the gene of interest is on

* <span style="color: #F0FF00">**gene_name:**</span> This value has the gene name of interest

* <span style="color: #F0FF00">**start:**</span> This value is an integer that has the start base position of the gene. You can use gnomAD to find this, just make sure that the build matches what was providied to hap-IBD to identify segments.

* <span style="color: #F0FF00">**end:**</span> This value is an integer that has the end base position of the gene. This value can also be found using gnomAD, just again be aware of the genome build.

* <span style="color: #F0FF00">**length:**</span> This value will be an integer that has the total length of the segment in cM.

**Log File:**

The program also produces a log file in whatever directory the specified as output. This file is called "IBDCluster.log" by default. It will contain all of the parameters passed to the program and it will have all of the logging messages for whatever log level the user supplied.