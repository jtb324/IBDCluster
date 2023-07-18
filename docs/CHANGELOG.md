## Unreleased

### Feat

- **logging**: added a logging submodule

### Fix

- fixed bug where the loglevel was not being reset to the original one in the record_inputs method
- fixed cases when there are no controls
- Fixed a bug where the header line was appending a new line incorrectly
- Fixed typos that were messing with the pvalue calculation
- Removed duplicate values from the ibd_vs attribute of the IbdFilter

### Refactor

- log situation where ibd_pd is empty and early termination of program
- Add a check to see if dataframe is empty after filtering for individuals. If so then it continues to next loop
- **logging**: Added more informative debug logging to the case_file_parser.py
- switched list in phenotype dictionary to sets
- adding more informative logging statements
- make sure that the ibd file matches the chromosome of the chromosome target
- updated the exclusion criteria so it includes blank spaces
- **logging**: Fixed the record_inputs method of the logger so that now it records the inputs and writes it to a file
- **logging**: Changed the code to adapt to the new OOP logger style
- Finished refactoring drive
- refactor filter to remove unnecessary function and to determine all individuals in cohort
- **plugins**: Added a config file for the plugin systemm to specify what plugins exist
- **plugins**: Setup the plugin architecture
- Updated poetry lock file for dependencies
- **clustering**: Finished refactoring the clustering module
- Created models for the Data class and the network class
- removed the accidental file
- moved callbacks into utilities modules
- Added the load_phenotype_descriptions to the __init__.py
- Added a function called load_phenotype_descriptions that will read in the phecode descriptions file if the user passes it
- Added attribute to keep a list of all individuals in cohort
- **phenotyping**: adjusted the phenotype parser to support multiple phenotypes
- **phenotypeing**: refactored how phenotype files are read in to determine case, control, or exclusion individuals
- fixed merge conflict in .gitmodules
- merging the proper changes
- **clustering**: finished first refactor of the clustering algorithm
- Added jupyter notebook as dependency group to test.
- **igraph-clustering**: Fixed clustering bug
- **logging**: changed log levels to be 1 or 2 to reflect typers use of -v or -vv
- started adding in the redo clustering steps to the code
- Adjusted the imports
- Made sure to reset the index so that the dataframe is the same as the original code
- **PhenotypeFileParser**: added functionality to the parser so it can determine cases/controls/exclusions
- **gitignore**: Add the /tests/test_input directory to version controls
- **logging**: Changed the logging in the filter.py file
- **name-change**: Changed the filter module to filters to avoid name collision with filter function
- Added a directory for utility functions.
- **clustering**: Broke the Networks class into two classes
- **clustering**: Finished the initial cluster step and fixed the missing attribute in igraph 0.10.4
- Created the clustering class and wrote functionality for the initial clustering step
- Starting refactoring the clustering and reorganizing the models
- **type-hints**: using type hints from Typing module
- **Error-Handling**: updated the Filter._remove_dups method to raise a KeyError if columns aren't present
- **logging**: Added a submodule for logging
- **filtering-ibd-file**: Restructured hung-hsin's filtering into a module
- **logging**: added logging to the program
- **log**: removed the messed up log folder
- removed the log from the .gitmodules
- **debugging**: added a __str__ method to indices classes

### Perf

- updated set operations

## v1.0.0 (2023-03-31)

### Refactor

- **DRIVE.py**: moved DRIVE.py to drive.py and added a pre-commit hook
- **global-variables**: removed more global variables
- **global-variable-removal**: removed some global variables that were constants and converted them to user options in the commandline
- **drive/generate_indices.py**: remove branching in create_indices
- **file-indices**: Switched to a strategy pattern
- **split_target_string**: Added a check to the split_target_string function
- **drive/DRIVE.py**: refactor how the target string was split and how to indices are created
- **drive/callbacks.py**: created callback function to make sure that the input ibd file exist
- **DRIVE.py**: Refactored CLI to use typer-cli
- **profiles/python_3_8_drive_dcm_hh.prof**: addeed a profile for python 3.8
- **pyproject.toml**: Added static type checkers and linters

### Perf

- **profiles**: added directory to keep track of profiles
