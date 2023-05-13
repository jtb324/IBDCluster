DRIVE Inputs
============

The DRIVE program has several command line arguments shown in the image below. 

.. ![image](https://belowlab.github.io/drive/assets/images/DRIVE_cli_options.png)

----------

Required Inputs:
----------------
* **input**: This input file describes the pairwise shared IBD segments within the cohort. The file is formed as the result of running hap-IBD, iLASH, GERMLINE, or RapID.

----

* **format**: This argument indicates which ibd program that was used to identify IBD segments. Currently the program supports values of hapibd, hap-IBD, iLASH, RapID, and GERMLINE.

----

* **target**: This argument is the region of interest that you wish to cluster around. The program will filter out segments that don't contain this entire region. This argument should be of the format chromosome_number:start-end (eg: 7:1234-2345).


.. warning::

    Make sure that the base position in the file corresponds to the same build of the human genome as what you used in the IBD detection software. If the builds are different then you will get inaccurate results.


----

* **output**: This argument will indicate a path to write the output file to. The user should provide a file path without an extension and the program will add the extension *.drive_networks.txt*.

Optional Arguments:
-------------------

* **min-cm**: DRIVE will use this argument to filter out all pairwise IBD segments that are shorter than the provided threshold. This value defaults to 3cM.

----

* **step**: This argument indicates the number of minimum steps that the random walk will use to generate a network. By default this value is 3.

----

* **max-check**: This value indicates the maximum number of times the program will redo the clustering in an effort to perform tree pruning. This argument defaults to 5 if the user doesn't provide any value. This argument is ignored if the flag --no-recluster is provided. 

----

* **cases**: A tab separated text file containing individuals who are either cases, controls, or exclusions. This file expects for there to be at least columns. The first column will have individual ids. All other columns in the file are for each phenotype being analyzed. Each column is expected to have the individual's status where cases are indicated by a 1, controls are indicated by a 0, and excluded individuals are indiciated by -1, -1.0, N/A, or a blank space. Excluded individuals will not be included in the binomial test but will be included in the clustering analysis. The file is expected to have a header where the first column is grid or grids (case insensitive) and the remaining columns are the phenotype names. If this argument is not supplied than the program will not perform the binomial test that determines enrichment of phenotypes within the identified networks. 

----

* **descriptions**: This argument provides the path to a tab separated text file that contains a text description of each phenotype. This file is typically supplied when the user is running the program phenomewide with PheCodes but can be provided for any phenotype. This file is expected to have two columns called "phecode" and "phenotype"

----

* **max-network-size**: This argument provides one of the thresholds that determines if a network needs to be re-clustered. max-network-size indicates the largest size of a network that is permitted. If a network is bigger than the max-network-size argument and it is smaller than the min-connected-threshold argument then the network will be re-clustered. 

-----

* **min-connected-threshold**: This argument provides one of the thresholds that determines if a network needs to be re-clustered. min-connected-threshold indicates how connected we want the network to be. If the network is sparser than the min-connected-threshold and is larger than the max-network-size then the network will be re-clustered.

----

* **min-network-size**: This argument sets a threshold for the minimum size a network has to be to be included in the analysis. Users can filter out pairs or trios but changing this value. By default this argument is set to 2

----

* **segment-distribution-threshold**: Threshold to filter the network length to remove hub individuals

----

* **hub-threshold**: Threshold to determine what percentage of hubs to keep

----

* **json_path**: This argument provides a path to the config.json file that DRIVE uses to identify what plugins the user wishes to use. DRIVE uses the plugin architecture to allow users to add new functionality to the program. The program comes with default plugins to determine enrichment pvalues and to write networks to a file. These plugins are listed in the config.json file found within the programs install directory. If the user wishes to use DRIVE in its standard form then they should just use the default value. If the user wishes to extend DRIVE then they can read more about the plugin system and how to add new plugins HERE.

----

* **recluster**: Flag indicating if the user wishes to recluster networks or if they wish to only return the initial networks without redoing the clustering algorithm. If the user wishes to use reclustering then DRIVE will use the values for hub-threshold, segment-distribution-threshold, min-connected-threshold, and max-network-size. If the user wishes to not perform reclustering, then they should pass the flag "*--no-recluster*".

----

* **verbose**: Flag indicating how verbose the user wants the drive program to be. The flag can be combined with itself to indicating more verbosity (-v = verbose while -vv = debug mode). By default the program will provide minimum information. If the user passes -v the program will run in verbose mode. If the user passes -vv then it will run in debug mode. Debug mode will generate a lot of logging output so use with caution.

----

* **log_to_console**: Flag indicating if the user wishes to write output to the console. If this flag is passed then output will be written to the console. Otherwise the output will only be written to the log file.

----

* **log_filename**: Argument providing a name for the log file generated by drive. This log file will be written to the parent directory from the "output" argument