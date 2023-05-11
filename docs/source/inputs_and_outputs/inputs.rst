DRIVE Inputs
============

The DRIVE program has several command line arguments shown in the image below. 

.. ![image](https://belowlab.github.io/drive/assets/images/DRIVE_cli_options.png)

----------

Required Inputs:
----------------
* <span style="color: #F0FF00">**INPUT**:</span> This input file is formed as the result of running hap-IBD, iLASH, GERMLINE, or RapID. These program output a file where pairwise shared IBD segments have been identified.

---

* <span style="color: #F0FF00">**FORMAT**:</span> This argument indicates which ibd program that was used to identify IBD segments. Currently the program supports values of hapibd, hap-IBD, iLASH, RapID, and GERMLINE.

---

* <span style="color: #F0FF00">**TARGET**:</span> This argument is the region of interest that you wish to cluster around. The program will filter out segments that don't contain this entire region. This argument should be of the format chromosome_number:start-end (Ex: 7:1234-2345).


.. warning::

    Make sure that the base position in the file corresponds to the same build of the human genome as what you used in the IBD detection software. If the builds are different then you will get inaccurate results.


---

* <span style="color: #F0FF00">**OUTPUT**:</span> This argument will indicate a path to write output to. The user should provide a file path without an extension and the program will add the extension.

Optional Arguments:
-------------------

* <span style="color: #A0A0A0">**MINCM**:</span> This argument indicates the minimum centimorgan threshold that the program will use to filter out IBD segments smaller than the threshold. By default this value is 3.

---
* <span style="color: #A0A0A0">**STEP**:</span> This argument indicates the number of minimum steps that the random walk will use to generate a network. By default this value is 3.