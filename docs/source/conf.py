# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'DRIVE'
copyright = '2023, James Baker, Hung-Hsin Chen, David Samuels, Jennifer Piper-Below'
author = 'James Baker, Hung-Hsin Chen, David Samuels, Jennifer Piper-Below'
release = '2.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

github_url = "https://github.com/belowlab/drive"

extensions = [
    'sphinx.ext.napoleon',
    'sphinx_copybutton',
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "piccolo_theme"

html_theme_options = {
    "source_url": 'https://github.com/belowlab/drive',
    "source_icon": 'github'
}
html_static_path = ['_static']
