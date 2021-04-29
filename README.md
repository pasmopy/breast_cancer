# breast_cancer

[![Actions Status](https://github.com/dyaus-dev/breast_cancer/workflows/Tests/badge.svg)](https://github.com/dyaus-dev/breast_cancer/actions)
[![License](https://img.shields.io/badge/License-Apache%202.0-green.svg)](https://github.com/dyaus-dev/breast_cancer/blob/master/LICENSE)

Workflow for classifying breast cancer subtypes based on intracellular signaling dynamics.

## Requirements

| Language      | Dependent packages                                 |
| ------------- | -------------------------------------------------- |
| Python >= 3.7 | [dyaus](https://github.com/dyaus-dev/dyaus)        |
| Julia >= 1.5  | [BioMASS.jl](https://github.com/himoto/BioMASS.jl) |
| R >= 4.0      | TCGAbiolinks, sva, biomaRt, edgeR                  |

## Table of contents

- [Integration](#integration-of-tcga-and-ccle-data)

- [Construction](#construction-of-a-comprehensive-model-of-the-ErbB-signaling-network)

- [Individualization](#individualization-of-the-mechanistic-model)

- [Classification](#subtype-classification-based-on-the-ErbB-signaling-dynamics)

## Integration of TCGA and CCLE data

- Run [`transcriptomic_data_integration.R`](transcriptomic_data/transcriptomic_data_integration.R)

  ```bash
  $ cd transcriptomic_data
  $ Rscript transcriptomic_data_integration.R
  ```

## Construction of a comprehensive model of the ErbB signaling network

1. Use `dyaus.Text2Model` to build a mechanistic model

   ```python
   from dyaus import Text2Model

   Text2Model("models/erbb_network.txt").to_biomass_model()
   ```

1. Rename `erbb_network/` to CCLE_name or TCGA_ID, e.g., `MCF7_BREAST` or `TCGA_3C_AALK_01A`

1. Add weighting factors for each gene (prefix: `"w_"`) to [`name2idx/parameters.py`](models/breast/TCGA_3C_AALK_01A/name2idx/parameters.py)

1. Edit [`set_search_param.py`](models/breast/TCGA_3C_AALK_01A/set_search_param.py)

   ```python
   import os

   import numpy as np

   from dyaus import Individualization

   from . import __path__
   from .name2idx import C, V
   from .set_model import initial_values, param_values

   incorporating_gene_expression_levels = Individualization(
       parameters=C.NAMES,
       species=V.NAMES,
       tpm_values="transcriptomic_data/TPM_RLE_postComBat.csv",
       relationship={
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

1. Build a mechanistic model for parameter estimation with BioMASS.jl

   ```python
   from dyaus import Text2Model

   Text2Model("models/erbb_network.txt", lang="julia").to_biomass_model()
   ```

1. Add time-series data to [`experimental_data.jl`](training/erbb_network_jl/experimental_data.jl)

1. Set an objective function to be minimized in [`fitness.jl`](training/erbb_network_jl/fitness.jl)

1. Run [`optimize_parallel.sh`](training/optimize_parallel.sh)

   ```bash
   $ mv erbb_network_jl training
   $ cd training
   $ mkdir errout
   $ sh optimize_parallel.sh
   $ cd ..
   ```

   When finished, run:

   ```julia
   using BioMASS

   param2biomass("training")
   ```

   And you will get **dat2npy/out/**. This is the optimized parameter sets that biomass can load.
   Copy **out/** to each biomass model folder via:

   ```python
   import os
   import shutil

   models = []
   path_to_models = os.path.join("models", "breast")
   for f in os.listdir(path_to_models):
        if os.path.isdir(os.path.join(path_to_models, f)) and (
            f.startswith("TCGA_") or f.endswith("_BREAST")
        ):
            models.append(f)
   # Set optimized parameters
   for model in models:
       shutil.copytree(
           os.path.join("training", "dat2npy", "out"),
           os.path.join(path_to_models, f"{model}", "out"),
       )
   ```

### Execute patient-specific models

- Use `dyaus.PatientModelSimulations`

  ```python
  import os
  import shutil

  from dyaus import PatientModelSimulations


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
  simulations = PatientModelSimulations("models.breast", TCGA_ID)
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

1. Visualize patient classification

## License

[Apache License 2.0](https://github.com/dyaus-dev/dyaus/blob/master/LICENSE)
