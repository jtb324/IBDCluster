# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
# import drive


project = 'DRIVE'
copyright = '2023, James Baker, Hung-Hsin Chen, David Samuels, Jennifer Piper-Below'
author = 'James Baker, Hung-Hsin Chen, David Samuels, Jennifer Piper-Below'
release = '2.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

sys.path.insert(0, os.path.abspath('../../drive/'))
sys.path.insert(0, os.path.abspath('../../'))

github_url = "https://github.com/belowlab/drive"

extensions = [
    'sphinx.ext.napoleon',
    'sphinx_copybutton',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
    # 'sphinxcontrib_autodocgen'
]

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
master_doc = 'index'
extensions = []
