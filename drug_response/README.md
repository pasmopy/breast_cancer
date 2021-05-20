# CCLE drug response data analysis and visualization

## Drug response data

#### Description:

> Pharmacologic profiles for 24 anticancer drugs across 504 CCLE lines.

#### Source:

- https://data.broadinstitute.org/ccle_legacy_data/pharmacological_profiling/CCLE_NP24.2009_Drug_data_2015.02.24.csv

#### Reference:

- Barretina, J., Caponigro, G., Stransky, N. _et al_. The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity. _Nature_ **483**, 603–607 (2012). https://doi.org/10.1038/nature11003

## Usage

1. Calculation of the ErbB receptor expression ratio

   ```bash
   $ Rscript data/calc_erbb_ratio.R
   ```

1. Drug response analysis and visualization

   ```python
   import pandas as pd
   from drug.database import CancerCellLineEncyclopedia

   ccle = CancerCellLineEncyclopedia()

   erbb_expression_ratio = pd.read_csv("./data/ErbB_expression_ratio.csv", index_col=0)

   compounds = ["Erlotinib", "Lapatinib", "AZD6244", "PD-0325901"]

   for compound in compounds:
       ccle.save_all(erbb_expression_ratio, compound)
   ```
