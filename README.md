# Breast cancer [![Actions Status](https://github.com/pasmopy/breast_cancer/workflows/Tests/badge.svg)](https://github.com/pasmopy/breast_cancer/actions)

Workflow for classifying breast cancer subtypes based on intracellular signaling dynamics.

## Requirements

| Language      | Dependent packages                                             |
| ------------- | -------------------------------------------------------------- |
| Python >= 3.7 | See [`requirements.txt`](requirements.txt)                     |
| Julia >= 1.5  | [BioMASS.jl](https://github.com/biomass-dev/BioMASS.jl)        |
| R >= 4.0      | TCGAbiolinks, sva, biomaRt, ComplexHeatmap, viridisLite, dplyr, edgeR, sva, tibble, data.table, stringr, biomaRt |

## Table of contents

- [Integration of TCGA and CCLE data](#integration-of-tcga-and-ccle-data)

- [Construction of a comprehensive model of the ErbB signaling network](#construction-of-a-comprehensive-model-of-the-ErbB-signaling-network)

- [Individualization of the mechanistic model](#individualization-of-the-mechanistic-model)

- [Subtype classification based on the ErbB signaling dynamics](#subtype-classification-based-on-the-ErbB-signaling-dynamics)

- [Investigation of patient-specific pathway activities](#investigation-of-patient-specific-pathway-activities)

## Integration of TCGA and CCLE data


### Download package for run following R script
- Run code.R

  ```R
  source("code.R")
  ```

### Download TCGA clinical/subtype information

- Run `outputclinical()` or `outputsubtype()`  
**outputclinical()** :  You can select all cancer types in [TCGA Study Abbreviations](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations).  
**outputsubtype()** :  You can select "LGG", "LUAD", "STAD", "BRCA", COAD", "READ"  

  ```R
  outputclinical("BRCA")
  outputsubtype("BRCA")
  ```

  Output: `"TCGA Study Abbreviation"_clinic.csv` or `"TCGA Study Abbreviation"_subtype.csv`


### Select samples in reference to clinical or subtype data

- You can select the patient's state based on the clinical or subtype data obrained above.   

  ```R
  sampleselection(type = subtype, 
                  ID = "patient",
                  pathologic_stage %in% c("Stage_I", "Stage_II"),
                  age_at_initial_pathologic_diagnosis < 60)
  ```

    **type** :   
    You can choose `clinic` or `subtype`. If you specify `clinic`, refer to `"TCGA Study Abbreviation"_clinic.csv`, and if you specify `subtype`, refer to `"TCGA Study Abbreviation"_subtype.csv` to select the patient. In order to select each one, you need to run outputclinical() or outputsubtype() before running this code.  

    **ID** :   
    Column name that contains the patient's ID (ex. TCGA-E2-A14U, TCGA-E9-A1RC, ...) in the .csv file referenced by "type". 

    **After line 3** :  
    You can set multiple conditions for selecting samples. 

    | Terms      | How to write                                             |
    | ------------- | -------------------------------------------------------------- |
    | all patients meet "A" in column x | x == "A" |
    | all patients meet "A" or "B" or "C" in column x | x %in% c("A", "B", "C") |
    | all patients have a value greater than A in column x | x > A  |
    | all patients have a value less than A in column x | x < A  |




### Download TCGA gene expression data (HTSeq-Counts)

 - Download the gene expression data of the specified sample types [(Sample Type Codes)](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes) in the cancer type specified by `outputclinical()` or `outputsubtype()`. By running this code, you can get data of only the patients selected by `sampleselection()`.

   ```R
   downloadTCGA(cancertype = "BRCA", 
                sampletype = c("01", "06"))
   ```  
   Output: Number of selected samples


### Download CCLE transcriptomic data


- Download CCLE transcriptomic data. You can select cell lines derived from [`one specific cancer type`](CCLE_cancertype.txt).

  ```R
  downloadCCLE(cancertype = "BREAST")
  ```  
  Output: Number of selected samples
 

### Merge TCGA and CCLE data
 1. Merge TCGA data download with `downloadTCGA()` and CCLE data download with `downloadCCLE()`
 2. Run ComBat-seq program to remove batch effects between TCGA and CCLE datasets  
 3. Output total read counts of all samples in order to decide the cutoff value of total read counts for `Normalization()`

    ```R
     MergeTCGAandCCLE()
     ```  

    Output : totalreadcounts.csv 

### Normalization of RNA-seq counts data
 - Conduct noramlization of RNA-seq .
 - You can specify min and max value for truncation of total read counts.
 - If you do not want to specify values for truncation, please set "min = F" or "max = F"

   ```R
   Normalization(min = 40000000, max = 140000000)
   ```  
   Output : TPM_RLE_postComBat.csv

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
   for f in os.listdir(path_to_models):
        if os.path.isdir(os.path.join(path_to_models, f)) and (
            f.startswith("TCGA_") or f.endswith("_BREAST")
        ):
            breast_cancer_models.append(f)
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
  from pasmopy import PatientModelAnalyses

  import models.breast

  with open ("selected_tnbc.txt" mode="r") as f:
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
