.. DRIVE documentation master file, created by

   sphinx-quickstart on Sat May  6 10:53:12 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to DRIVE's documentation!
=================================

Distant Relatedness for Identification and Variant Evaluation (DRIVE) is a novel approach to IBD-based genotype inference used to identify shared chromosomal segments in dense genetic arrays. DRIVE implemented a random walk algorithm that identifies clusters of individuals who pairwise share an IBD segment overlapping a locus of interest. This tool was developed in python by the Below Lab at Vanderbilt University.

Quick Installation:
-------------------
The easiest way to install DRIVE is from the PYPI register. Users can use the following command to install the program. This method ensures that all the necessary dependencies are installed.

.. code:: bash

   pip install drive-ibd

.. note::

   To read more about this install go to the Pip Installation section

Citation:
---------
The paper that originally discusses DRIVE can be found here: `<https://www.medrxiv.org/content/10.1101/2023.04.19.23288831v1>`_


Contact:
--------
If you have any questions about DRIVE, you can either post an issue on the Github issues page or you can contact us at the email address, insert email here.

.. toctree::
   :maxdepth: 2
   :caption: User's Guide
   :hidden:

   /installation/installation
   /inputs_and_outputs/inputs
   /inputs_and_outputs/outputs
   
   

.. toctree::
   :maxdepth: 2
   :caption: Tutorial
   :hidden:

   /examples/tutorial

.. toctree::
   :maxdepth: 2
   :caption: Developer's Guide
   :hidden:

   /contributing/contributing
   /plugin_descriptions/extending_drive
   /modules/modules



