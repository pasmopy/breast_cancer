# Breast cancer [![Actions Status](https://github.com/pasmopy/breast_cancer/workflows/Tests/badge.svg)](https://github.com/pasmopy/breast_cancer/actions)

This repository contains analysis code for the following paper:



## Requirements

| Language      | Dependent packages                                             |
| ------------- | -------------------------------------------------------------- |
| Python >= 3.7 | See [`requirements.txt`](requirements.txt)                     |
| Julia >= 1.5  | [`BioMASS.jl==0.5.0`](https://github.com/biomass-dev/BioMASS.jl)        |
| R >= 4.0      | TCGAbiolinks, sva, biomaRt, edgeR, ComplexHeatmap, viridisLite |

## Table of contents

- [Integration of TCGA and CCLE data](#integration-of-tcga-and-ccle-data)

- [Construction of a comprehensive model of the ErbB signaling network](#construction-of-a-comprehensive-model-of-the-ErbB-signaling-network)

- [Individualization of the mechanistic model](#individualization-of-the-mechanistic-model)

- [Subtype classification based on the ErbB signaling dynamics](#subtype-classification-based-on-the-ErbB-signaling-dynamics)

- [Investigation of patient-specific pathway activities](#investigation-of-patient-specific-pathway-activities)

## Integration of TCGA and CCLE data

- Run [`transcriptomic_data_integration.R`](transcriptomic_data/transcriptomic_data_integration.R)

  ```bash
  $ cd transcriptomic_data
  $ Rscript transcriptomic_data_integration.R
  ```

  Output: `TPM_RLE_postComBat.csv`

## Construction of a comprehensive model of the ErbB signaling network

1. Use `pasmopy.Text2Model` to build a mechanistic model

   ```python
   import os

   from pasmopy import Text2Model


   Text2Model(os.path.join("models", "erbb_network.txt")).convert()
   ```

1. Rename `erbb_network/` to CCLE_name or TCGA_ID, e.g., `MCF7_BREAST` or `TCGA_3C_AALK_01A`

1. Add weighting factors for each gene (prefix: `"w_"`) to [`name2idx/parameters.py`](models/breast/TCGA_3C_AALK_01A/name2idx/parameters.py)

1. Edit [`set_search_param.py`](models/breast/TCGA_3C_AALK_01A/set_search_param.py)

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
       transcriptomic_data=os.path.join("transcriptomic_data", "TPM_RLE_postComBat.csv"),
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

## Individualization of the mechanistic model

### Use time-course datasets to train kinetic constants and weighting factors

1. Build a mechanistic model to identify model parameters

   ```python
   import os

   from pasmopy import Text2Model


   Text2Model(os.path.join("models", "erbb_network.txt"), lang="julia").convert()
   ```

1. Add time-series data to [`experimental_data.jl`](training/erbb_network_jl/experimental_data.jl)

1. Set an objective function to be minimized in [`fitness.jl`](training/erbb_network_jl/fitness.jl)

1. Run [`optimize_parallel.sh`](training/optimize_parallel.sh)

   ```bash
   $ mv erbb_network_jl training
   $ cd training
   $ mkdir errout
   $ sh optimize_parallel.sh  # It will take more than a few days to optimize parameters.
   $ cd ..
   ```

   When finished, run:

   ```bash
   $ julia
   ```

   ```julia
   using BioMASS

   param2biomass("training")
   ```

   And you will get [`dat2npy/out/`](https://github.com/pasmopy/breast_cancer/tree/master/training/erbb_network_jl/dat2npy/out).
   This is the estimated parameter sets that [`biomass`](https://github.com/biomass-dev/biomass) can recognize and read.
   Copy `out/` to each patient-specific model folder via:

   ```python
   import os
   import shutil


   breast_cancer_models = []
   path_to_models = os.path.join("models", "breast")
   for model in os.listdir(path_to_models):
        if os.path.isdir(os.path.join(path_to_models, model)) and (
            model.startswith("TCGA_") or model.endswith("_BREAST")
        ):
            breast_cancer_models.append(model)
   # Set optimized parameters
   for model in breast_cancer_models:
       shutil.copytree(
           os.path.join("training", "erbb_network_jl", "dat2npy", "out"),
           os.path.join(path_to_models, f"{model}", "out"),
       )
   ```

### Execute patient-specific models

- Use `pasmopy.PatientModelSimulations`

  ```python
  import os
  import shutil

  from pasmopy import PatientModelSimulations

  import models.breast


  with open (os.path.join("models", "breast", "sample_names.txt"), mode="r") as f:
      TCGA_ID = f.read().splitlines()
  # Create patient-specific models
  for patient in TCGA_ID:
      if patient != "TCGA_3C_AALK_01A":
          shutil.copytree(
              os.path.join("models", "breast", "TCGA_3C_AALK_01A"),
              os.path.join("models", "breast", f"{patient}"),
          )
  # Execute patient-specific models
  simulations = PatientModelSimulations(models.breast.__package__, TCGA_ID)
  simulations.run()
  ```

## Subtype classification based on the ErbB signaling dynamics

1. Extract response characteristics from patient-specific simulations

   ```python
   simulations.subtyping(
       None,
       {
           "Phosphorylated_Akt": {"EGF": ["max"], "HRG": ["max"]},
           "Phosphorylated_ERK": {"EGF": ["max"], "HRG": ["max"]},
           "Phosphorylated_c-Myc": {"EGF": ["max"], "HRG": ["max"]},
       }
   )
   ```

1. Visualize patient classification by executing [`brca_heatmap.R`](classification/brca_heatmap.R)

   ```bash
   $ cd classification
   # $ Rscript brca_heatmap.R [n_cluster: int] [figsize: tuple]
   $ Rscript brca_heatmap.R 6 8,5
   ```

## Investigation of patient-specific pathway activities

### Sensitivity analysis

- Calculate sensitivity coefficients by varying the amount of each nonzero species

  ```python
  import os

  from pasmopy import PatientModelAnalyses

  import models.breast
  

  with open (os.path.join("models", "breast", "selected_tnbc.txt"), mode="r") as f:
      TNBC_ID = f.read().splitlines()
  analyses = PatientModelAnalyses(
      models.breast.__package__,
      TNBC_ID,
      biomass_kws={
          "metric": "maximum", "style": "heatmap", "options": {"excluded_initials": ["PIP2"]}
      },
  )
  ```

### Drug response data analysis

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

## Author

- Hiroaki Imoto
- Sawa Yamashiro

## License

[Apache License 2.0](LICENSE)
