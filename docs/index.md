# Requirements

## Manual installation of package requirements

The following packages are required for transcriptomic data integration, parameter estimation, patient-specific simulations, and result visualization.

### General:

- Python >= 3.7
- Julia >= 1.5
- R >= 4.0

### Python:

- [pasmopy==0.1.0](https://github.com/pasmopy/pasmopy)
- [biomass>=0.5.2,<0.6](https://github.com/biomass-dev/biomass)
- matplotlib==3.3.4
- numpy==1.19.2
- pandas==1.2.4
- seaborn==0.11.2
- scipy==1.6.0

### Julia:

- [BioMASS.jl==0.5.0](https://github.com/biomass-dev/BioMASS.jl)

### R:

- TCGAbiolinks
- sva
- biomaRt
- ComplexHeatmap
- viridisLite
- dplyr
- edgeR
- sva
- tibble
- data.table
- stringr
- biomaRt

# Getting started
## Integration of TCGA and CCLE data
### Download TCGA clinical/subtype information

- Move to `transcriptomic_data/` and start `R`

  ```bash
  $ cd transcriptomic_data
  $ R
  ```

- Read `integration.R`

  ```R
  source("integration.R")
  ```

- Run `outputClinical()` or `outputSubtype()`

  ```R
  outputClinical("BRCA")
  outputSubtype("BRCA")
  ```

  Output: `<TCGA Study Abbreviation>_clinic.csv` or `<TCGA Study Abbreviation>_subtype.csv`

### Select samples in reference to clinical or subtype data

- You can select the patient's state based on the clinical or subtype data obtained above.

  ```R
  patientSelection(type = subtype,
                   ID = "patient",
                   pathologic_stage %in% c("Stage_I", "Stage_II"),
                   age_at_initial_pathologic_diagnosis < 60)
  ```

### Download TCGA gene expression data (HTSeq-Counts)

- Download the gene expression data of the specified sample types ([Sample Type Codes](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes)) in the cancer type specified by `outputClinical()` or `outputSubtype()`. By running this code, you can get data of only the patients selected by `sampleSelection()`.

  ```R
  downloadTCGA(cancertype = "BRCA",
               sampletype = c("01", "06"),
               outputresult = FALSE)
  ```

  Output: Number of selected samples

### Download CCLE transcriptomic data

- Download CCLE transcriptomic data. You can select cell lines derived from [one specific cancer type](https://github.com/pasmopy/breast_cancer/blob/master/transcriptomic_data/CCLE_cancertype.txt).

  ```R
  downloadCCLE(cancertype = "BREAST",
               outputresult = FALSE)
  ```

  Output: Number of selected samples

### Merge TCGA and CCLE data

1.  Merge TCGA data download with `downloadTCGA()` and CCLE data download with `downloadCCLE()`.
1.  Run ComBat-seq program to remove batch effects between TCGA and CCLE datasets.
1.  Output total read counts of all samples in order to decide the cutoff value of total read counts for `normalization()`.

    ```R
    mergeTCGAandCCLE(outputesult = FALSE)
    ```

    Output : `totalreadcounts.csv`

### Normalize RNA-seq counts data

- Conduct normalization of RNA-seq.
- You can specify min and max value for truncation of total read counts.
- If you do not want to specify values for truncation, please set `min=F` or `max=F`.

  ```R
  normalization(min=40000000, max=140000000)
  ```

  Output : `TPM_RLE_postComBat_<TCGA>_<CCLE>.csv`

## Construction of a comprehensive model of the ErbB signaling network

### From text into executable models

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

### Other tasks for incorporating gene expression levels

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

       ## As maximal transcription rate
       x[C.V291] = incorporating_gene_expression_levels.as_reaction_rate(
           __path__[0].split(os.sep)[-1], x, "V291", "DUSP"
       )
       x[C.V310] = incorporating_gene_expression_levels.as_reaction_rate(
           __path__[0].split(os.sep)[-1], x, "V310", "cMyc"
       )
       ## As initial conditions
       y0 = incorporating_gene_expression_levels.as_initial_conditions(
           __path__[0].split(os.sep)[-1], x, y0
       )

       ...
   ```

## Individualization of the mechanistic model
### Parameter estimation

Here, we use phospho-protein time-course datasets to train kinetic constants and weighting factors.

📒 **Note:** In this study, we used `BioMASS.jl`, but in most cases you can use `pasmopy.optimize()` function for parameter estimation.

1. Build a mechanistic model to identify model parameters

   ```python
   import os

   from pasmopy import Text2Model


   Text2Model(os.path.join("models", "erbb_network.txt"), lang="julia").convert()
   ```

1. Add time-series data to [`experimental_data.jl`](https://github.com/pasmopy/breast_cancer/blob/master/training/erbb_network_jl/experimental_data.jl)

1. Set an objective function to be minimized in [`fitness.jl`](https://github.com/pasmopy/breast_cancer/blob/master/training/erbb_network_jl/fitness.jl)

1. Run [`optimize_parallel.sh`](https://github.com/pasmopy/breast_cancer/blob/master/training/optimize_parallel.sh)

   ```bash
   $ mv erbb_network_jl training
   $ cd training
   $ mkdir errout
   $ sh optimize_parallel.sh  ## It will take more than a few days to optimize parameters.
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
   ## Set optimized parameters
   for model in breast_cancer_models:
       shutil.copytree(
           os.path.join("training", "erbb_network_jl", "dat2npy", "out"),
           os.path.join(path_to_models, f"{model}", "out"),
       )
   ```

### Patient-specific simulations

- Use `pasmopy.PatientModelSimulations`

  ```python
  import os
  import shutil
  from pathlib import Path

  from pasmopy import PatientModelSimulations

  import models.breast


  TCGA_ID = [
      l.strip() for l in Path("models", "breast", "sample_names.txt").read_text("utf-8").splitlines()
  ]
  ## Create patient-specific models
  for patient in TCGA_ID:
      if patient != "TCGA_3C_AALK_01A":
          shutil.copytree(
              os.path.join("models", "breast", "TCGA_3C_AALK_01A"),
              os.path.join("models", "breast", f"{patient}"),
          )
  ## Execute patient-specific models
  simulations = PatientModelSimulations(models.breast.__package__, TCGA_ID)
  simulations.run()
  ```

## Subtype classification based on the ErbB signaling dynamics

### Extraction of response characteristics from patient-specific simulations

- Execute `subtyping()`

  ```python
  simulations.subtyping(
      fname=None,
      dynamical_features={
          "Phosphorylated_Akt": {"EGF": ["max"], "HRG": ["max"]},
          "Phosphorylated_ERK": {"EGF": ["max"], "HRG": ["max"]},
          "Phosphorylated_c-Myc": {"EGF": ["max"], "HRG": ["max"]},
      }
  )
  ```

### Model-based patient stratification

- Run [`brca_heatmap.R`](https://github.com/pasmopy/breast_cancer/tree/master/classification/brca_heatmap.R)

  ```bash
  $ cd classification
  ## $ Rscript brca_heatmap.R [n_cluster: int] [figsize: tuple]
  $ Rscript brca_heatmap.R 6 8,5
  ```

## Investigation of patient specific pathway activities

### Sensitivity analysis

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