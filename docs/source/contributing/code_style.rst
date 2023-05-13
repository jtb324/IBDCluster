Code Styling Tips
=================

When contributing to this project we ask that you use a few development tools to maintain consistent code style across the code base. Each of these tools will be described in the sections below. These tools can be install with poetry when you run the following command:

.. note:: 
    :class: sidebar

    If you are only contributing to the documentation of DRIVE then you can jump to the section titled "Commitizen"

.. code:: bash

    poetry install --with dev

All of these tools can also be installed separately with pip if the user did not originally clone the repository and does not have the pyproject.toml file.

Black
-----
For this project we use the famous "Uncompromising Code Formatter", Black. Black makes most of the style decisions for us so that we don't have to. Individuals who contribute to DRIVE are expected to run black before each push so that the code style is maintained. For more information about Black, you can read the documentation here `Black's documentation <https://black.readthedocs.io/en/stable/>`_.

Isort
-----
Isort is used to alphabetically sort imports. Isort also arranges imports by type. Users are expected to run Isort before pushing any proposed changes. For more information about Isort and how to run the program, the documentation is here: `Isort's documentation <https://pycqa.github.io/isort/index.html>`_. For compatibility with Black, ensure that the following lines are within the pyproject.toml file.

.. code::

    [tool.isort]
    profile = "black"

Pylint
------
Pylint is used as a static code analyzer for the project. This tool will detect common code smell, catch common errors, indicate if code blocks are too complex, and suggest ways to refactor code blocks. For more information about Pylint, the documentation can be found here: `Pylint documentation <https://pylint.readthedocs.io/en/latest/>`_. For compatibility with Black, ensure that the following lines are within the pyproject.toml file.

.. code::

    [tool.pylint.format]
    max-line-length = "88"

Commitizen
----------
Commitizen is used to enforce consistency across commits and help with the automatic generation of changelogs. Commitizen relies on the concept of `conventional commits <https://www.conventionalcommits.org/en/v1.0.0/>`_ to make commit messages more readable and meaning for both humans and machines. For more information on how to use commitizen, the documentation can be found here `Commitizen documentation <https://commitizen-tools.github.io/commitizen/>`_. If the repository wasn't cloned from github then commitizen will require the user to create a pyproject.toml file for configuration options. If this file needs to be made then the user should copy and paste the following code into the pyproject.toml file to maintain consistency:

.. code::
    
    [tool.commitizen]
    name = "cz_conventional_commits"
    version = "1.1.0"
    tag_format = "$version"