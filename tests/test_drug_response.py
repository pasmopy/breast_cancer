import os
import shutil

import pandas as pd
from drug_response.drug.database import CancerCellLineEncyclopedia


def test_create_figs():
    for dir in ["dose_response", "activity_area"]:
        if os.path.isdir(dir):
            shutil.rmtree(dir)
    ccle = CancerCellLineEncyclopedia()
    erbb_expression_ratio = pd.read_csv(
        os.path.join("drug_response", "data", "ErbB_expression_ratio.csv"),
        index_col=0,
    )
    compounds = list(set(ccle.drug_response_data["Compound"]))
    for compound in compounds:
        ccle.save_all(erbb_expression_ratio, compound)
        for dir in ["dose_response", "activity_area"]:
            assert os.path.isfile(
                os.path.join(
                    f"{dir}",
                    f"{ccle._drug2target(compound)}",
                    f"{ccle._convert_drug_name(compound)}.pdf",
                )
            )


def test_cleanup():
    for dir in ["dose_response", "activity_area"]:
        shutil.rmtree(dir)
