---
layout: default
title: Integration of TCGA and CCLE data
parent: Getting started
nav_order: 1
---

# Integration of TCGA and CCLE data

## Download TCGA clinical/subtype information

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

## Select samples in reference to clinical or subtype data

- You can select the patient's state based on the clinical or subtype data obtained above.

  ```R
  patientSelection(type = subtype,
                   ID = "patient",
                   pathologic_stage %in% c("Stage_I", "Stage_II"),
                   age_at_initial_pathologic_diagnosis < 60)
  ```

## Download TCGA gene expression data (HTSeq-Counts)

- Download the gene expression data of the specified sample types ([Sample Type Codes](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes)) in the cancer type specified by `outputClinical()` or `outputSubtype()`. By running this code, you can get data of only the patients selected by `sampleSelection()`.

  ```R
  downloadTCGA(cancertype = "BRCA",
               sampletype = c("01", "06"),
               outputresult = FALSE)
  ```

  Output: Number of selected samples

## Download CCLE transcriptomic data

- Download CCLE transcriptomic data. You can select cell lines derived from [one specific cancer type](https://github.com/pasmopy/breast_cancer/blob/master/transcriptomic_data/CCLE_cancertype.txt).

  ```R
  downloadCCLE(cancertype = "BREAST",
               outputresult = FALSE)
  ```

  Output: Number of selected samples

## Merge TCGA and CCLE data

1.  Merge TCGA data download with `downloadTCGA()` and CCLE data download with `downloadCCLE()`.
1.  Run ComBat-seq program to remove batch effects between TCGA and CCLE datasets.
1.  Output total read counts of all samples in order to decide the cutoff value of total read counts for `normalization()`.

    ```R
    mergeTCGAandCCLE(outputesult = FALSE)
    ```

    Output : `totalreadcounts.csv`

## Normalize RNA-seq counts data

- Conduct normalization of RNA-seq.
- You can specify min and max value for truncation of total read counts.
- If you do not want to specify values for truncation, please set `min=F` or `max=F`.

  ```R
  normalization(min=40000000, max=140000000)
  ```

  Output : `TPM_RLE_postComBat_<TCGA>_<CCLE>.csv`