# dyaus_breast

[![License](https://img.shields.io/badge/License-Apache%202.0-brightgreen.svg)](https://opensource.org/licenses/Apache-2.0)

Workflow for classifying breast cancer subtypes based on intracellular signaling dynamics.

## Requirements

| Language      | Dependent packages                                                                             |
| ------------- | ---------------------------------------------------------------------------------------------- |
| Python >= 3.7 | [biomass](https://github.com/okadalabipr/biomass), [dyaus](https://github.com/dyaus-dev/dyaus) |
| Julia >= 1.5  | [BioMASS.jl](https://github.com/himoto/BioMASS.jl)                                             |
| R             | [TODO] Write dependent packages here.                                                          |

## Individualization of the mechanistic model

### Integrate TCGA and CCLE data

[TODO] Write analysis procedure here.

### Build an executable model of the ErbB signaling network

1. Use `biomass.Text2Model` to build a mechanistic model

   ```python
   from biomass import Text2Model

   Text2Model("models/erbb_network.txt").to_biomass()
   ```

1. Rename **erbb_network/** to CCLE_name or TCGA_ID

1. Edit **name2idx/parameters.py**
   - Add weighting factors (`'w_XXX'`)
1. Edit **set_search_param.py**

   - Import `dyaus.Individualization`

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
       structure={
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
   )

   ...

   def update(self, indiv):
       x = param_values()
       y0 = initial_values()

       for i, j in enumerate(self.idx_params):
           x[j] = indiv[i]
       for i, j in enumerate(self.idx_initials):
           y0[j] = indiv[i + len(self.idx_params)]

       x[C.V291] = incorporating_gene_expression_levels.as_maximal_transcription_rate(
           __path__[0].split(os.sep)[-1], x, "V291", "DUSP"
       )
       x[C.V310] = incorporating_gene_expression_levels.as_maximal_transcription_rate(
           __path__[0].split(os.sep)[-1], x, "V310", "cMyc"
       )
       y0 = incorporating_gene_expression_levels.as_initial_condition(
           __path__[0].split(os.sep)[-1], x, y0
       )

       ...
   ```

### Train model parameters against time-course datasets obtained from breast cancer cell lines

1. Use `biomass.Text2Model` to build a mechanistic model for parameter estimation with BioMASS.jl

   ```python
   from biomass import Text2Model

   Text2Model("models/erbb_network.txt", lang="julia").to_biomass()
   ```

1. Add time-series data to **experimental_data.jl**

1. Set an objective function to be minimized in **fitness.jl**

1. Run **optimize_parallel.sh**

   ```bash
   $ mv erbb_network_jl training
   $ cd training
   $ mkdir errout
   $ sh optimize_parallel.sh
   ```

   When finished, run:

   ```julia
   using BioMASS

   param2biomass(".")
   ```

   And you will get **dat2npy/out/**. This is the optimized parameter sets that biomass can load. Copy **out/** to each biomass model folder.

### Execute patient-specific models

1. Use `dyaus.PatientModelSimulations`

   ```python
   import os
   import shutil

   from dyaus import PatientModelSimulations


   with open (
       os.path.join(
           "models",
           "breast",
           "sample_names.txt",
       ), mode="r"
   ) as f:
       TCGA_ID = f.read().splitlines()

   # Create patient-specific models
   for patient in TCGA_ID:
       if patient != "TCGA_3C_AALK_01A":
           shutil.copytree(
               os.path.join(
                   "models",
                   "breast",
                   "TCGA_3C_AALK_01A",
               ),
               os.path.join(
                   "models",
                   "breast",
                   f"{patient}",
               )
           )

   simulations = PatientModelSimulations("models.breast", TCGA_ID)

   simulations.run()
   ```

## Subtype classification based on the ErbB signaling dynamics

[TODO] Write analysis procedure here.

## License

[Apache-2.0 License](https://opensource.org/licenses/Apache-2.0)
