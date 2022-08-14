## CCLE drug response data analysis and visualization

### Drug response data

#### Description:

> Pharmacologic profiles for 24 anticancer drugs across 504 human cancer lines.

#### Source:

- https://data.broadinstitute.org/ccle_legacy_data/pharmacological_profiling/CCLE_NP24.2009_Drug_data_2015.02.24.csv

#### Reference:

- Barretina, J., Caponigro, G., Stransky, N. _et al_. The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity. _Nature_ **483**, 603–607 (2012). https://doi.org/10.1038/nature11003

### Usage

1. Prepare [`sample.txt`](data/sample.txt) and [`gene.txt`](data/gene.txt), which contain the names of the samples and genes, respectively.

1. Open R

   ```bash
   $ R
   ```

1. Load `CCLE_normalization.R`

   ```R
   source("calc_erbb_ratio.R")
   ```

1. Read [`sample.txt`](data/sample.txt) and [`gene.txt`](data/gene.txt)

   If you want data for a specific samples or genes, you will need to create a list of those samples as [`sample.txt`](data/sample.txt) or [`gene.txt`](data/gene.txt) and read it here.

   ```R
   gene <- scan("gene.txt", what="character")
   #sample <- scan("sample.txt", what="character")
   ```

1. Create TPM/RLE normalized RNA-seq data matrix of selected samples and genes

   ```R
   CCLEnormalization(gene, sample = NULL)
   ```

   Output: [`CCLE_normalized.csv`](data/CCLE_normalized.csv) (TPM/RLE normalized RNA-seq data for specific samples or genes)

1. Calculate EGFR/(ERBB2+ERBB3+ERBB4) using [`CCLE_normalized.csv`](data/CCLE_normalized.csv)

   If you want to run this code, you need to prepare `gene.txt` containing the names of these 4 genes.

   ```R
   CCLE_normalized <- read.csv("CCLE_normalized.csv", row.names = 1)
   receptor_ratio(data = CCLE_normalized, num = 30)
   ```

   Output: [`ErbB_expression_ratio.csv`](data/ErbB_expression_ratio.csv)

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
