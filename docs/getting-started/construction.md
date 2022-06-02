---
layout: default
title: Construction of a comprehensive model of the ErbB signaling network
parent: Getting started
nav_order: 2
---

# Construction of a comprehensive model of the ErbB signaling network

## From text into executable models

1. Use [**pasmopy.Text2Model**](https://pasmopy.readthedocs.io/en/latest/model_development.html) to build a mechanistic model

   ```python
   import os

   from pasmopy import Text2Model


   Text2Model(os.path.join("models", "erbb_network.txt")).convert()
   ```

1. Rename `erbb_network/` to CCLE_name or TCGA_ID, e.g., `MCF7_BREAST` or `TCGA_3C_AALK_01A`

   ```python
   import shutil

   shutil.move(
       os.path.join("models", "erbb_network"),
       os.path.join("models", "breast", "TCGA_3C_AALK_01A")
   )
   ```

## Other tasks for incorporating gene expression levels

1. Add weighting factors for each gene (prefix: `"w_"`) to `name2idx/parameters.py`

   ```python
   from pasmopy import Model
   from pasmopy.preprocessing import WeightingFactors

   from models import erbb_network


   model = Model(erbb_network.__package__).create()

   gene_expression = {
       "ErbB1": ["EGFR"],
       "ErbB2": ["ERBB2"],
       "ErbB3": ["ERBB3"],
       "ErbB4": ["ERBB4"],
       "Grb2": ["GRB2"],
       "Shc": ["SHC1", "SHC2", "SHC3", "SHC4"],
       "RasGAP": ["RASA1", "RASA2", "RASA3"],
       "PI3K": ["PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"],
       "PTEN": ["PTEN"],
       "SOS": ["SOS1", "SOS2"],
       "Gab1": ["GAB1"],
       "RasGDP": ["HRAS", "KRAS", "NRAS"],
       "Raf": ["ARAF", "BRAF", "RAF1"],
       "MEK": ["MAP2K1", "MAP2K2"],
       "ERK": ["MAPK1", "MAPK3"],
       "Akt": ["AKT1", "AKT2"],
       "PTP1B": ["PTPN1"],
       "GSK3b": ["GSK3B"],
       "DUSP": ["DUSP5", "DUSP6", "DUSP7"],
       "cMyc": ["MYC"],
   }

   weighting_factors = WeightingFactors(model, gene_expression)
   weighting_factors.add_to_params()
   weighting_factors.set_search_bounds()
   ```

1. Edit `set_search_param.py`

   ```python
   import os

   import numpy as np

   from pasmopy import Individualization

   from . import __path__
   from .name2idx import C, V
   from .set_model import initial_values, param_values


   incorporating_gene_expression_levels = Individualization(
       parameters=C.NAMES,
       species=V.NAMES,
       transcriptomic_data=os.path.join("transcriptomic_data", "TPM_RLE_postComBat_BRCA_BREAST.csv"),
       gene_expression={
           "ErbB1": ["EGFR"],
           "ErbB2": ["ERBB2"],
           "ErbB3": ["ERBB3"],
           "ErbB4": ["ERBB4"],
           "Grb2": ["GRB2"],
           "Shc": ["SHC1", "SHC2", "SHC3", "SHC4"],
           "RasGAP": ["RASA1", "RASA2", "RASA3"],
           "PI3K": ["PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"],
           "PTEN": ["PTEN"],
           "SOS": ["SOS1", "SOS2"],
           "Gab1": ["GAB1"],
           "RasGDP": ["HRAS", "KRAS", "NRAS"],
           "Raf": ["ARAF", "BRAF", "RAF1"],
           "MEK": ["MAP2K1", "MAP2K2"],
           "ERK": ["MAPK1", "MAPK3"],
           "Akt": ["AKT1", "AKT2"],
           "PTP1B": ["PTPN1"],
           "GSK3b": ["GSK3B"],
           "DUSP": ["DUSP5", "DUSP6", "DUSP7"],
           "cMyc": ["MYC"],
       },
       read_csv_kws={"index_col": "Description"}
   )

   ...

   def update(self, indiv):
       x = param_values()
       y0 = initial_values()

       for i, j in enumerate(self.idx_params):
           x[j] = indiv[i]
       for i, j in enumerate(self.idx_initials):
           y0[j] = indiv[i + len(self.idx_params)]

       # As maximal transcription rate
       x[C.V291] = incorporating_gene_expression_levels.as_reaction_rate(
           __path__[0].split(os.sep)[-1], x, "V291", "DUSP"
       )
       x[C.V310] = incorporating_gene_expression_levels.as_reaction_rate(
           __path__[0].split(os.sep)[-1], x, "V310", "cMyc"
       )
       # As initial conditions
       y0 = incorporating_gene_expression_levels.as_initial_conditions(
           __path__[0].split(os.sep)[-1], x, y0
       )

       ...
   ```