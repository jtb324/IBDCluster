.. .. raw:: html

..     <style> .yellow {color:yellow; font-weight:bold;} </style>

.. .. role:: yellow

Example command format
============================

The following command assumes that you have either installed DRIVE using pip or that you have installed it from github. The following command will show you how to call it if you install DRIVE using pip. This example command has only the required arguments.

.. code::

    drive -i {input ibd filepath} -t {chromosome position to cluster around} -o {output filepath}

.. note::

    If you installed DRIVE from github then you can replace the 'drive' portion with "python /path_to_drive.py" or you can add drive to your path. The rest of the command will be the same.

Explanation of command:
-----------------------

* **input ibd filepath**: filepath to the output from a program like hap-IBD, iLASH, or GERMLINE. This file should correspond to the chromosome that the locus of interest is on.


* **chromosome position to cluster around**: string indicating the target region of interest should be of the form chromosome:start position-end position (An example is chrX:XXXX-XXXX).


* **output filepath**: filepath to write an output file to. This value should not include a suffix. DRIVE will automatically append the suffix ".DRIVE.txt".