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

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = "piccolo_theme"

html_theme = 'sphinx_book_theme'
html_title = "DRIVE documentation"
html_theme_options = {
    "repository_url": 'https://github.com/belowlab/drive',
    "repository_provider": 'github',
    "use_edit_page_button": True,
    "use_source_button": True,
    # "use_issues_button": True,
    "repository_branch": "refactor",
    "path_to_docs": "/docs/source/",
    "use_repository_button": True,
    # 'navbar_pagenav': True
}
html_static_path = ['_static']

# Napolean Settings
napoleon_numpy_docstring = True

# autodocgen_config = [{
#         # 'modules':[drive],
#         # 'generated_source_dir': DOC_SOURCE_DIR+'/autodocgen-output/',

#         # if module matches this then it and any of its submodules will be skipped
#         # 'skip_module_regex': '(.*[.]__|myskippedmodule)',

#         # produce a text file containing a list of everything documented. you can use this in a test to notice when you've
#         # intentionally added/removed/changed a documented API
#         'write_documented_items_output_file': 'autodocgen_documented_items.txt',

#         # customize autodoc on a per-module basis
#         # 'autodoc_options_decider': {
#         #         'mymodule.FooBar':    { 'inherited-members':True },
#         # },

#         # # choose a different title for specific modules, e.g. the toplevel one
#         # 'module_title_decider': lambda modulename: 'API Reference' if modulename=='mymodule' else modulename,
# }]