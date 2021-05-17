import os
import shutil

os.chdir("drug_response")

import pandas as pd
from drug.database import CancerCellLineEncyclopedia


def test_create_figs():
    ccle = CancerCellLineEncyclopedia()
    erbb_expression_ratio = pd.read_csv(
        os.path.join("data", "CCLE_receptor_ratio.csv"),
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
            shutil.rmtree(dir)


os.chdir("..")
