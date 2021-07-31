import os
from dataclasses import dataclass, field
from typing import Dict, List

from biomass.exec_model import ModelObject


@dataclass
class WeightingFactors(object):
    """
    Prepare for adding information about gene expression data to model.

    Attributes
    ----------
    model : biomass.exec_model.ModelObject
        BioMASS model object.
    gene_expression : dict
        Pairs of proteins and their related genes.
    prefix : str (default: "w_")
        Prefix of weighting factors on gene expression levels.
    indentation : str (default: 4 spaces)
        How many spaces as indentation.
    """
    model: ModelObject
    gene_expression: Dict[str, List[str]]
    prefix: str = field(default="w_", init=False)
    indentation: str = field(default=" " * 4, init=False)

    def add(self) -> None:
        """
        Add weighting factors to model parameters.
        """
        weighting_factors: List[str] = []
        for genes in self.gene_expression.values():
            for gene in genes:
                if self.prefix + gene not in self.model.parameters:
                    weighting_factors.append(self.prefix + gene)
        if weighting_factors:
            with open(
                os.path.join(self.model.path, "name2idx", "parameters.py"),
                mode="r",
                encoding="utf-8",
            ) as f:
                lines = f.readlines()
            for line_num, line in enumerate(lines):
                if line.startswith("NUM: int"):
                    lines[line_num] = "NAMES.extend(\n"
                    lines[line_num] += (
                        f"{self.indentation}[\n"
                        + f'{2 * self.indentation}"'
                        + f'",\n{2 * self.indentation}"'.join(weighting_factors)
                        + '",\n'
                        f"{self.indentation}]\n"
                    )
                    lines[line_num] += ")\n\nNUM: int = len(NAMES)\n"
            with open(
                os.path.join(self.model.path, "name2idx", "parameters.py"),
                mode="w",
                encoding="utf-8",
            ) as f:
                f.writelines(lines)


    def set_search_bounds(self, lb: float = 0.01, ub: float = 100.0) -> None:
        """
        Set search bounds for weighting factors.

        Parameters
        ----------
        lb : float
            Lower bound.
        ub : float
            Upper bound.
        """
        weighting_factors: List[str] = []
        for param in self.model.parameters:
            if param.startswith(self.prefix):
                weighting_factors.append(param)
        if weighting_factors:
            search_bounds = [f"search_rgn[:, C.{wf}] = [{lb}, {ub}]" for wf in weighting_factors]
            with open(
                os.path.join(self.model.path, "set_search_param.py"),
                mode="r",
                encoding="utf-8",
            ) as f:
                lines = f.readlines()
            for line_num, line in enumerate(lines):
                if line.startswith(f"{2 * self.indentation}search_rgn = convert_scale("):
                    lines[line_num] = (
                        2 * self.indentation
                        + f"\n{2 * self.indentation}".join(search_bounds)
                    )
                    lines[line_num] += f"\n\n{2 * self.indentation}search_rgn = convert_scale(\n"
            with open(
                os.path.join(self.model.path, "set_search_param.py"),
                mode="w",
                encoding="utf-8",
            ) as f:
                f.writelines(lines)
