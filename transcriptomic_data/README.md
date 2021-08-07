## Transcriptomic data processing

### Requirements

| Language      | Dependent packages                                             |
| ------------- | -------------------------------------------------------------- |
| R | dplyr, edgeR, sva, tibble, data.table, stringr, TCGAbiolinks ,  biomaRt |
<br>

## Transcriptomic data integration
Integrating TCGA and CCLE data for parameterization and individualization of the mechanistic model.  
<br>
### Download TCGA clinical/subtype information

- Read `integration.R`
  ```bash
  $ R
  ```  

  ```R
  source("integration.R")
  ```

- Run `outputclinical()` or `outputsubtype()`

  | Function | Description |
  | ---      | ---         |
  | `outputClinical()` | You can select all cancer types in [TCGA Study Abbreviations](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations). |
  | `outputSubtype()` |  You can select "ACC", "BRCA", "BLCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "PAAD", "PCPG", "PRAD", "READ", "SKCM", "SARC", "STAD", THCA", "UCEC", "UCS", "UVM".

  ```R
  outputClinical("COAD")
  outputSubtype("COAD")
  ```

  Output: `{TCGA Study Abbreviation}_clinic.csv` or `{TCGA Study Abbreviation}_subtype.csv`


### Select samples in reference to clinical or subtype data

- You can select the patient's state based on the clinical or subtype data obtained above.   

  ```R
  patientSelection(type = clinical, 
                   ID = "submitter_id",
                   pathologic_stage %in% c("Stage_Ia", "Stage_II"),
                   age_at_initial_pathologic_diagnosis < 80)
  ```

    **type** :   
    You can choose `clinical` or `subtype`. If you specify `clinical`, refer to `{TCGA Study Abbreviation}_clinical.csv`, and if you specify `subtype`, refer to `{TCGA Study Abbreviation}_subtype.csv` to select the patient. In order to select each one, you need to run outputClinical() or outputSubtype() before running this code.  

    **ID** :   
    Column name that contains the patient's ID (ex. TCGA-E2-A14U, TCGA-E9-A1RC, ...) in the .csv file referenced by "type". 

    **After line 3** :  
    You can set multiple conditions for selecting samples. 

    | Expression | Meaning |
    | ---------- | ------- |
    | x == "A" | All patients meeting "A" in column x |
    | x %in% c("A", "B", "C") | All patients meeting "A" or "B" or "C" in column x |
    | x > A | All patients having a value greater than A in column x |
    | x < A | All patients having a value less than A in column x |




### Download TCGA gene expression data (HTSeq-Counts)

 - Download the gene expression data of the specified sample types [(Sample Type Codes)](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes) in the cancer type specified by `outputClinical()` or `outputSubtype()`. By running this code, you can get data of only the patients selected by `sampleSelection()`.

   ```R
   downloadTCGA(cancertype = "COAD", 
                sampletype = "01",
                outputresult = FALSE)
   ```  
   Output: Number of selected samples


### Download CCLE transcriptomic data


- Download CCLE transcriptomic data. You can select cell lines derived from [`one specific cancer type`](CCLE_cancertype.txt).

  ```R
  downloadCCLE(cancertype = "BREAST",
               outputresult = FALSE)
  ```  
  Output: Number of selected samples
 

### Merge TCGA and CCLE data
 1. Merge TCGA data download with `downloadTCGA()` and CCLE data download with `downloadCCLE()`.
 1. Run ComBat-seq program to remove batch effects between TCGA and CCLE datasets.
 1. Output total read counts of all samples in order to decide the cutoff value of total read counts for `normalization()`.

    ```R
    mergeTCGAandCCLE(outputesult = FALSE)
    ```  

    Output : `totalreadcounts.csv `

### Normalization of RNA-seq counts data

- Conduct noramlization of RNA-seq.
- You can specify min and max value for truncation of total read counts.
- If you do not want to specify values for truncation, please set `min = F` or `max = F`.

  ```R
  normalization(min = 40000000, max = 140000000)
  ```  
  Output : `TPM_RLE_postComBat_{TCGA}_{CCLE}.csv`  

<br>
<br>  

## Data transformation for reusing parameter sets  
Predicting the ErbB signaling dynamics of the new samples using parameter sets optimized in the previous

  ```R
  source(value_transformation.R)

  BRCA_TPM <- fread("TPM_RLE_postComBat_BRCA_BREAST.csv")
  COAD_TPM <- fread("TPM_RLE_postComBat_COAD_BREAST.csv")
  gene <- scan(gene.txt, what="character")
  cellline <- c("MCF7_BREAST", "BT474_BREAST", "SKBR3_BREAST", "MDAMB231_BREAST")

  value_transformation(upstream = BRCA_TPM, 
                       origin = COAD_TPM,
                       gene = gene,
                       cellline = cellline,
                       outputfile = "TPM_RLE_postComBat_COAD_BREAST.csv")
  ```  

  **upstream** :   
    Gene expression data used for parameter estimation

  **origin** :   
    Gene expression data for the samples for which you want to predict dynamics using the parameters estimated by `optimizeddata`.

  **gene** :  
    List of genes used in the ODE model
    
  **cellline** :  
    Cell lines for which time-series experimental data were collected　　

  **filename** :   
    Name of output file

   Output : `"outputfile"`