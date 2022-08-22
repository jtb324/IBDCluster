---
layout: default 
title: Inputs
parent: Documentation
nav_order: 2
---
# Inputs:
---
IBDCluster has a number of inputs, some of which are required and some that are optional. These are listed below. The required input files will be colored yellow while the optional ones will be grey:

---

* <span style="color: #F0FF00">**ibd file**:</span> This input file is the result of running hap-IBD to determine pairwise IBD sharing between individuals in a BioBank. This file has to be formed before because IBDCluster does not determine pairwise ibd sharing. The table below shows the format of the file (The file doesn't have a header line but I have added one to the table for ease of explaining the file).

| Pair 1 ID | Pair 1 Phase | Pair 2 ID | Pair 2 Phase | Chromosome | Segment Start Pos | Segment End Pos | Segment Length (cM) |
|:----------|:-------------|:----------|:-------------|:-----------|:------------------|:----------------|:---------------|
| ID1   | Phase 1 | ID2 | Phase 2 | Chromosome # | Segment start base position | Segment end base position | Segment length |

There should be a total of eight columns that are all tab separated. The file should have information about each individual in the pair and the haplotype phase for each variant, the chromosome that the segment is on, the location of the segment, and the length of the segment. All of this information will be used by the program so it is mandatory. 

This file will be supplied through the --ibd-file or -f flag. It is a required file.

---

* <span style="color: #A0A0A0">**.env file**:</span> This input file is just a environment file that has two variables "HAPIBD_PATH" and "JSON_PATH". The HAPIBD_PATH variable has the directory where the ibd files are. The suffix of the file has to be .env. This file will be supplied through the --env or the -e argument. This file is required by the program but the the program comes with a default file within the install directory  so you do not have to supply an argument unless you are using a custom .env file.

---

* <span style="color: #A0A0A0">**config.json**:</span> This input file is a json file that has information about the plugins that are going to be used in the analysis. The default config.json for the stock plugins are shown below. 

```json
{
    "plugins":["plugins.pvalues", "plugins.network_writer", "plugins.allpair_writer"],
    "modules": [
        {
            "name":"pvalues"
        },
        {
            "name": "network_writer"
        },
        {
            "name": "allpair_writer"
        }
    ]
}
```
The json files should have the keys "plugins" and "modules". 

The value for the "plugins" key is a list of strings that have the module path of the plugin. To break this down each plugin is located in the plugins module in the IBDCluster plugin, therefore the first part of the module path in the default config.json is "plugins". The name of the plugin .py file is the second part of the module path. For example there is a default plugin pvalues.py within the plugin directory, therefore the module path is "plugins.pvalues".  This structure of the plugins directory is explained in more detail within the section called "Stock Plugins". 

The value for the "modules" key is a list of key, value pairs where the key is "name" and the value is the name of the plugin which can be found in the initialize function in the plugin .py file (This is also described in more detail in the "Stock Plugins" section). 

This file will be supplied through the optional --json-config or the -j flag. There is a default config file that IBDCluster reads from but you can provide your own config file if you wish to turn on or off a plugin or if you wish to use a custom plugin (This process will be described in the "Custom Plugin Design" section).

---
* <span style="color: #F0FF00">**Gene Information File**:</span> This input file is a tab separated text file that has information about the gene region of interest such as the gene name, the chromosome the gene is on, the start base position of the gene, and the end base position of the gene. This information can usually be found at the [GnomAD genome browser](https://gnomad.broadinstitute.org/). An example of this file for the gene rbm20 is shown below (Once again I've added a header line for readability. The file should not have a header line or else IBDCluster will give an error that will not be straightforward to debug). 

| gene name | chromosome | gene start position | gene end position |
|:----------|:-----------|:--------------------|:------------------|
|   RBM20   |     10     |      112404115      |      112599227    |

***Note on what information can be proivded in the file:*** <br>
So far we only discussed how to provide information about a gene through the Gene Information File. There is no reason though that IDBCluster can only be used to cluster around a gene. If you were interested in a variant instead of a gene then you could just replace the information in the file as shown below for the variant chr10:110812304 (using build Hg38) 

|     gene name    | chromosome | gene start position | gene end position |
|:-----------------|:-----------|:--------------------|:------------------|
|   10:110812304   |     10     |      110812284      |      110812284    | 

The program does not necessarily use the gene name for anything but record keeper so you could use the variant name or another id, and then in this example we just made the start and end position the variant location. The program will see this as valid input.

IBDCluster can also be used on just any section of interest in the gene name. You just need to provide a "name" for the gene name, and then you need to provide the start and end position of the segment of interest.

This file will be supplied through the --gene-file or -g flag. It is required by IBDCluster. 

{: .warning }
***Warning about the genome position:***<br>
Make sure that the base position in the file corresponds to the same build of the human genome as what you used in the IBD detection software. If the builds are different then you will get inaccurate results.

---

* <span style="color: #F0FF00">**carrier file**:</span> This input file is a tab separated text file that indicates which individuals are affected by the phenotypes of interest. The first column is expected to be titled "grids". Every other column after this should be either a 0 or 1 for each phenotype of interest. Currently this program only supports binary phenotypes. An example of this format is shown below.

| grids | Phenotype A | Phenotype B |
|:------|:------------|:------------|
| grid 1|      1      |      0      |
| grid 2|      0      |      1      |
| grid 3|      0      |      0      |

This required input file will be supplied to the --carriers or -c flag. 

{: .optional}
Some people may be familiar with a PheCode matrix, and thats all this file really is. You're not restricted to only using PheCodes though. Any phenotype will work for this program as long as you can use binary phenotyping to determine cases and controls.