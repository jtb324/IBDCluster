.. raw:: html

    <style> .yellow {color:yellow; font-weight:bold;} </style>

.. role:: yellow

DRIVE program inputs
=====================

The DRIVE program has several command line arguments shown in the image below. 

.. ![image](https://belowlab.github.io/drive/assets/images/DRIVE_cli_options.png)

Required Inputs:
----------------

----

* :yellow:`INPUT`: This input file is formed as the result of running hap-IBD, iLASH, GERMLINE, or RapID. These program output a file where pairwise shared IBD segments have been identified.

----

* :yellow:`FORMAT`: This argument indicates which ibd program that was used to identify IBD segments. Currently the program supports values of hapibd, hap-IBD, iLASH, RapID, and GERMLINE.

----

* :yellow:`TARGET`: This argument is the region of interest that you wish to cluster around. The program will filter out segments that don't contain this entire region. This argument should be of the format chromosome_number:start-end (Ex: 7:1234-2345).


.. caution:: 

    Make sure that the base position in the file corresponds to the same build of the human genome as what you used in the IBD detection software. If the builds are different then you will get inaccurate results.


----

* :yellow:`OUTPUT`: This argument will indicate a path to write output to. The user should provide a file path without an extension and the program will add the extension.

Optional Arguments:
-------------------

* :yellow:`MINCM`: This argument indicates the minimum centimorgan threshold that the program will use to filter out IBD segments smaller than the threshold. By default this value is 3.

----

* :yellow:`STEP`: This argument indicates the number of minimum steps that the random walk will use to generate a network. By default this value is 3.