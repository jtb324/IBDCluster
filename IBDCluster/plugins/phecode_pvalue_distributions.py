from dataclasses import dataclass
import plugins


@dataclass
class Phecode_Pvalue_Distributions:

    name: str = "Phecode_Pvalue_Distributions"

    def analyze(self, **kwargs) -> None:
        """execute method"""
        print(f"networks_file: {kwargs['network_filepath']}")


def initialize() -> None:
    plugins.factory_register(
        "phecode_pvalue_distributions", Phecode_Pvalue_Distributions
    )
