from pathlib import Path
import argparse
import sys


class CheckInputExist(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs) -> None:
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(CheckInputExist, self).__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: Path,
        option_string: str = None,
    ) -> None:
        if values.exists():
            setattr(namespace, self.dest, values)
        else:
            print(
                f"ERROR: The file, {values}, was not found. Please make sure that there is not a typo in the file name."
            )
            sys.exit(1)


# class CheckJsonPath(argparse.Action):

#     def __init__(self, option_strings, dest, nargs=None, **kwargs) -> None:
#         if nargs is not None:
#             raise ValueError("nargs not allowed")
#         super(CheckJsonPath, self).__init__(option_strings, dest, **kwargs)


#     def __call__(self, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: Path, option_string: str =None) -> None:
#         print(f"This is the content in values: {values}")
#         if values:
#             setattr(namespace, self.dest, values)
#         else:
#             src_dir = Path(__file__).parent.parent.parent

#             config_path = src_dir / "config.json"

#             if not config_path.exists():
#                 raise FileNotFoundError(
#                     f"Expected the user to either pass a configuration file path or for the config.json file to be present in the program root directory at {config_path}."  # noqa: E501
#                 )

#             setattr(namespace, self.dest, config_path)
