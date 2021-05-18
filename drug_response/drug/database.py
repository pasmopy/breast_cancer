import os
from dataclasses import dataclass, field
from typing import List, NamedTuple, NoReturn, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import brunnermunzel


class DrugResponse(NamedTuple):
    ccle_cell_line_name: str
    primary_cell_line_name: str
    compound: str
    target: str
    doses: np.ndarray
    activity_data: np.ndarray
    activity_sd: np.ndarray
    num_data: int
    fit_type: str
    ec50: float
    ic50: float
    amax: float
    act_area: float


@dataclass
class CancerCellLineEncyclopedia(object):
    """
    Cancer Cell Line Encyclopedia (CCLE)
    https://portals.broadinstitute.org/ccle

    Attributes
    ----------
    _drug_response_data : pandas DataFrame
        Pharmacologic profiles for 24 anticancer drugs across 504 cell lines.
    """
    _drug_response_data: pd.DataFrame = field(
        default=pd.read_csv(
            "https://data.broadinstitute.org/ccle_legacy_data/"
            "pharmacological_profiling/CCLE_NP24.2009_Drug_data_2015.02.24.csv"
        ),
        init=False,
    )

    @property
    def drug_response_data(self) -> pd.DataFrame:
        return self._drug_response_data

    @staticmethod
    def _convert_drug_name(name: str) -> str:
        if name == "AZD6244":
            return "Selumetinib"
        elif name == "ZD-6474":
            return "Vandetanib"
        else:
            return name
    
    def _drug2target(self, drug: str) -> str:
        target = list(
            self.drug_response_data[self.drug_response_data["Compound"] == drug]["Target"]
        )
        return target[0]

    def _get_drug_response(
        self,
        cell_line: Optional[List[str]] = None,
        compound: Optional[List[str]] = None,
    ) -> List[DrugResponse]:
        """
        Extract drug-response data.

        Parameters
        ----------
        cell_line : list of strings, optional
            List of CCLE Cell Line Names.
        compound : list of strings, optional
            List of drug names.

        Returns
        -------
        drug_response_info : list of DrugResponse
            Drug response information.
        """

        if cell_line is None:
            cell_line = list(self.drug_response_data.loc[:, "CCLE Cell Line Name"])
        if compound is None:
            compound = list(self.drug_response_data.loc[:, "Compound"])
        df = self.drug_response_data[
            (self.drug_response_data["CCLE Cell Line Name"].isin(cell_line))
            & (self.drug_response_data["Compound"].isin(compound))
        ]
        drug_response_info = []
        for i in df.index:
            drug_response_info.append(
                DrugResponse(
                    df.at[i, "CCLE Cell Line Name"],
                    df.at[i, "Primary Cell Line Name"],
                    df.at[i, "Compound"],
                    df.at[i, "Target"],
                    np.array([float(val) for val in df.at[i, "Doses (uM)"].split(",")]),
                    np.array(
                        [float(val) for val in df.at[i, "Activity Data (median)"].split(",")]
                    ),
                    np.array([float(val) for val in df.at[i, "Activity SD"].split(",")]),
                    df.at[i, "Num Data"],
                    df.at[i, "FitType"],
                    df.at[i, "EC50 (uM)"],
                    df.at[i, "IC50 (uM)"],
                    df.at[i, "Amax"],
                    df.at[i, "ActArea"],
                )
            )
        return drug_response_info
    
    def _check_args(self, drug: str, save_format: str) -> Optional[NoReturn]:
        if drug not in set(self.drug_response_data["Compound"]):
            raise ValueError(
                f"{drug} doesn't exist in the CCLE drug response data.\n"
                f"Should be one of {', '.join(list(set(self.drug_response_data['Compound'])))}"
            )
        if save_format not in ["pdf", "png"]:
            raise ValueError("save_format must be either 'pdf' of 'png'.")

    def save_dose_response_curve(
        self,
        erbb_expression_ratio: pd.DataFrame,
        drug: str,
        *,
        show_individual: bool = False,
        error: str = "std",
        save_format: str = "pdf",
    ) -> None:
        self._check_args(drug, save_format)
        if error not in ["std", "sem"]:
            raise ValueError("error must be either 'std' or 'sem'.")

        os.makedirs(
            os.path.join(
                "dose_response",
                f"{self._drug2target(drug)}",
            ),
            exist_ok=True
        )

        top30 = []
        bottom30 = []
        for i, level in enumerate(list(erbb_expression_ratio.loc[:, "value"])):
            if level == "high":
                top30.append(erbb_expression_ratio.index[i])
            elif level == "low":
                bottom30.append(erbb_expression_ratio.index[i])

        egfr_high = self._get_drug_response(
            top30,
            [drug] * len(top30),
        )
        egfr_low = self._get_drug_response(
            bottom30,
            [drug] * len(bottom30),
        )
        
        plt.figure(figsize=(7, 5))
        plt.rcParams['font.size'] = 24
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['axes.linewidth'] = 2.4
        plt.rcParams['xtick.major.width'] = 2.4
        plt.rcParams['ytick.major.width'] = 2.4
        plt.rcParams['lines.linewidth'] = 3
        plt.rcParams["lines.markersize"] = 1
        plt.gca().spines["right"].set_visible(False)
        plt.gca().spines["top"].set_visible(False)

        DOSE = [2.50e-03, 8.00e-03, 2.50e-02, 8.00e-02, 2.50e-01, 8.00e-01, 2.53e+00, 8.00e+00]

        if show_individual:
            for i, _ in enumerate(egfr_high):
                plt.plot(
                    DOSE,
                    np.interp(DOSE, egfr_high[i].doses, egfr_high[i].activity_data) + 100,
                    "o",
                    color="darkmagenta",
                )
        y_mean = np.mean(
            [
                (np.interp(DOSE, egfr_high[i].doses, egfr_high[i].activity_data) + 100)
                for i, _ in enumerate(egfr_high)
            ],
            axis=0
        )
        y_err = np.std(
            [
                (np.interp(DOSE, egfr_high[i].doses, egfr_high[i].activity_data) + 100)
                for i, _ in enumerate(egfr_high)
            ],
            axis=0,
            ddof=1
        ) / (1 if error == "std" else np.sqrt(len(egfr_high)))
        plt.plot(
            DOSE,
            y_mean,
            "-",
            color="darkmagenta",
            label="EGFR high",
        )
        plt.fill_between(
            DOSE,
            y_mean - y_err,
            y_mean + y_err,
            lw=0,
            color="darkmagenta",
            alpha=0.1,
        )

        if show_individual:
            for i, _ in enumerate(egfr_low):
                plt.plot(
                    DOSE,
                    np.interp(DOSE, egfr_low[i].doses, egfr_low[i].activity_data) + 100,
                    "o",
                    color="goldenrod",
                )
        y_mean = np.mean(
            [
                (np.interp(DOSE, egfr_low[i].doses, egfr_low[i].activity_data) + 100)
                for i, _ in enumerate(egfr_low)
            ],
            axis=0
        )
        y_err = np.std(
            [
                (np.interp(DOSE, egfr_low[i].doses, egfr_low[i].activity_data) + 100)
                for i, _ in enumerate(egfr_low)
            ],
            axis=0,
            ddof=1
        ) / (1 if error == "std" else np.sqrt(len(egfr_high)))
        plt.plot(
            DOSE,
            y_mean,
            "-",
            color="goldenrod",
            label="EGFR low",
        )
        plt.fill_between(
            DOSE,
            y_mean - y_err,
            y_mean + y_err,
            lw=0,
            color="goldenrod",
            alpha=0.1,
        )

        plt.xscale("log")
        plt.xlabel("Concentration (μM)", fontsize=28)
        plt.title(f"{self._convert_drug_name(drug)}", fontsize=36)

        plt.ylabel("Relative viability (%)", fontsize=28)
        plt.yticks([0, 25, 50, 75, 100])

        plt.legend(loc="lower left", frameon=False, labelspacing=1)
        #plt.show()
        plt.savefig(
            os.path.join(
                "dose_response",
                f"{self._drug2target(drug)}",
                f"{self._convert_drug_name(drug)}.{save_format}",
            ),
            dpi=1200 if save_format == "png" else None,
            bbox_inches="tight",
        )
        plt.close()

    
    def save_activity_area(
        self,
        erbb_expression_ratio: pd.DataFrame,
        drug: str,
        *,
        save_format: str = "pdf",
    ) -> None:
        self._check_args(drug, save_format)

        os.makedirs(
            os.path.join(
                "activity_area",
                f"{self._drug2target(drug)}",
            ),
            exist_ok=True
        )
        
        top30 = []
        bottom30 = []
        for i, level in enumerate(list(erbb_expression_ratio.loc[:, "value"])):
            if level == "high":
                top30.append(erbb_expression_ratio.index[i])
            elif level == "low":
                bottom30.append(erbb_expression_ratio.index[i])
        
        egfr_high = self._get_drug_response(
            top30, 
            [drug]*len(top30),
        )
        egfr_low = self._get_drug_response(
            bottom30,
            [drug]*len(bottom30),
        )

        activity_area = np.array(
            [
                [egfr_high[i].act_area for i, _ in enumerate(egfr_high)],
                [egfr_low[i].act_area for i, _ in enumerate(egfr_low)],
            ],
        )

        plt.figure(figsize=(3, 5))
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['mathtext.it'] = 'Arial:italic'
        plt.rcParams['font.size'] = 24
        plt.rcParams['axes.linewidth'] = 2.4
        plt.rcParams['xtick.major.width'] = 2.4
        plt.rcParams['ytick.major.width'] = 2.4
        plt.rcParams['lines.linewidth'] = 1.8
        plt.rcParams['lines.markersize'] = 11
        plt.gca().spines["right"].set_visible(False)
        plt.gca().spines["top"].set_visible(False)

        p_value = brunnermunzel(activity_area[0], activity_area[1]).pvalue

        sns.boxplot(
            data=(activity_area[0], activity_area[1]),
            palette=sns.color_palette(['darkmagenta', 'goldenrod']),
        )
        #sns.stripplot(data=(activity_area[0], activity_area[1]),color=".26")
        sns.swarmplot(data=(activity_area[0], activity_area[1]),color=".48")
        #sns.boxenplot(data=(activity_area[0], activity_area[1]), palette=sns.color_palette(['#76106cff', '#f7a426ff']))
        #sns.violinplot(data=(activity_area[0], activity_area[1]), palette=sns.color_palette(['mediumslateblue', 'tomato']))
        plt.xticks([0, 1], ["EGFR\nhigh", "EGFR\nlow"])
        plt.ylabel("Activity area", fontsize=28)
        plt.suptitle(
            (r"$\it{p}$" + " = {:.2e}".format(p_value)) if p_value < 0.05 else "n.s.",
            fontsize=22,
        )
        #plt.show()
        plt.savefig(
            os.path.join(
                "activity_area",
                f"{self._drug2target(drug)}",
                f"{self._convert_drug_name(drug)}.{save_format}",
            ),
            dpi=1200 if save_format == "png" else None,
            bbox_inches="tight",
        )
        plt.close()
    
    def save_all(
        self,
        erbb_expression_ratio: pd.DataFrame,
        drug: str,
        *,
        save_format: str = "pdf",
    ) -> None:
        self.save_dose_response_curve(erbb_expression_ratio, drug, save_format=save_format)
        self.save_activity_area(erbb_expression_ratio, drug, save_format=save_format)