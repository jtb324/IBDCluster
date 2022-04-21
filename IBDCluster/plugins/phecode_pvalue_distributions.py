from dataclasses import dataclass
from typing import Dict, Any
from plugins import factory_register


@dataclass
class Phecode_Pvalue_Distributions:

    name: str = "Phecode_Pvalue_Distributions"

    def analyze(self, **kwargs) -> Dict[str, Any]:
        """execute method"""
        print("Not implemented")

    def write(self, **kwargs) -> None:
        print("not implemented")


def initialize() -> None:
    factory_register("phecode_pvalue_distributions", Phecode_Pvalue_Distributions)
