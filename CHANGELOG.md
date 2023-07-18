## Unreleased

### Feat

- **logging**: added a logging submodule

### Refactor

- **logging**: Added a submodule for logging
- **filtering-ibd-file**: Restructured hung-hsin's filtering into a module
- **logging**: added logging to the program
- **log**: removed the messed up log folder
- removed the log from the .gitmodules
- **debugging**: added a __str__ method to indices classes
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

## 1.0.1 (2023-03-26)

### Perf

- **profiles**: added directory to keep track of profiles
