---
layout: default
title: Investigation of patient specific pathway activities
parent: Getting started
nav_order: 5
---

# Investigation of patient specific pathway activities

## Sensitivity analysis

- Calculate sensitivity coefficients by varying the amount of each nonzero species

  ```python
  import os
  from pathlib import Path

  from pasmopy import PatientModelAnalyses

  import models.breast


  TNBC_ID = [
      l.strip() for l in Path("models", "breast", "selected_tnbc.txt").read_text("utf-8").splitlines()
  ]
  analyses = PatientModelAnalyses(
      models.breast.__package__,
      TNBC_ID,
      biomass_kws={
          "metric": "maximum", "style": "heatmap", "options": {"excluded_initials": ["PIP2"]}
      },
  )
  analyses.run()
  ```

## Drug response data analysis

1. Calculation of the ErbB receptor expression ratio

   ```bash
   $ cd drug_response
   $ Rscript data/calc_erbb_ratio.R
   ```

1. Drug response analysis and visualization

   ```python
   import os

   import pandas as pd
   from drug.database import CancerCellLineEncyclopedia


   ccle = CancerCellLineEncyclopedia()

   erbb_expression_ratio = pd.read_csv(
       os.path.join("data", "ErbB_expression_ratio.csv"),
       index_col=0
   )
   compounds = ["Erlotinib", "Lapatinib", "AZD6244", "PD-0325901"]
   for compound in compounds:
       ccle.save_all(erbb_expression_ratio, compound)
   ```