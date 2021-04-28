library(dplyr)
library(stringr)
library(tibble)

# TCGA ---------------------------------------------------------------------------------------

# download TCGA data from TCGA biolinks -----------------------------------
# https://bioconductor.org/packages/release/bioc/manuals/TCGAbiolinks/man/TCGAbiolinks.pdf

library(TCGAbiolinks)

# download counts file
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
results <- getResults(query)

GDCdownload(query)
df <- GDCprepare(query, save = T, save.filename = "TCGABRCA_counts.rda",
                 summarizedExperiment = F)

# download clinical data
BRCA_path_subtypes <- TCGAquery_subtype(tumor = "brca")
BRCA_path_subtypes_age_stage <- data.frame(patient = BRCA_path_subtypes$patient,
                                           stage = BRCA_path_subtypes$pathologic_stage,
                                           subtype = BRCA_path_subtypes$BRCA_Subtype_PAM50,
                                           age = BRCA_path_subtypes$age_at_initial_pathologic_diagnosis)



# exclude columns other than counts data ----------------------------------
df_counts <- df[grep("ENSG", df$X1),]


# exclude some samples ----------------------------------------------------
# exclude normal(non-cancerous) and duplicated samples
# https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
pattern <- c("-11A-", "-11B-", "-01B-", "-01C-")
TCGA_counts <- df_counts[,-grep(paste(pattern, collapse = "|"), colnames(df_counts))]
colnames(TCGA_counts) <- str_sub(colnames(TCGA_counts), start = 1, end = 16)
TCGA_counts <- TCGA_counts[,!duplicated(colnames(TCGA_counts))] #exclude duplicated samples
TCGA_counts$X1 <- str_sub(TCGA_counts$X1, start = 1, end = 15) 


# screen samples with patient's age and pathologic stage ------------------
# patients's age <= 59 & pathological stage <- stage I and II
BRCA_path_subtypes_age_stage_filtered <- BRCA_path_subtypes_age_stage[BRCA_path_subtypes_age_stage$age <= 59 &
                                                                        BRCA_path_subtypes_age_stage$stage %in% c("Stage_I", "Stage_II"),]
BRCA_path_subtypes_age_stage_filtered$patient <- gsub("-", ".", BRCA_path_subtypes_age_stage_filtered$patient)

# make counts file with filtered patients ---------------------------------
TCGA_counts_filtered <- data.frame(Name = TCGA_counts$X1,
                                   TCGA_counts[,grep(paste(BRCA_path_subtypes_age_stage_filtered$patient, collapse = "|"),
                                                     colnames(TCGA_counts))])


# CCLE ---------------------------------------------------------------------------------------

# download counts file
# https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_counts_20180929.gct.gz
url <- "https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_counts_20180929.gct.gz"
tmp <- tempfile()
download.file(url, tmp)
CCLE_counts <- read.table(gzfile(tmp), skip = 2, header = T)

# get only cell lines of breast cancer from CCLE datasets
CCLE_BRCA <- data.frame(Name = CCLE_counts$Name,
                        Description = CCLE_counts$Description,
                        CCLE_counts[,grep("BREAST", colnames(CCLE_counts))])
CCLE_BRCA[,1] <- str_sub(CCLE_BRCA[,1], start = 1, end = 15)




# Merge TCGA and CCLE ------------------------------------------------------------------------

TCGA_CCLE_counts <- inner_join(CCLE_BRCA, TCGA_counts_filtered, 
                               by = "Name")
TCGA_CCLE_counts <- TCGA_CCLE_counts %>%
  dplyr::distinct(Description, .keep_all=T) %>% column_to_rownames(var = "Name")
TCGA_CCLE_counts_numeric <- sapply(TCGA_CCLE_counts[,-1], as.numeric)
rownames(TCGA_CCLE_counts_numeric) <- rownames(TCGA_CCLE_counts)




# ComBat-seq ---------------------------------------------------------------------------------
# https://academic.oup.com/nargab/article/2/3/lqaa078/5909519
# https://rdrr.io/bioc/sva/man/ComBat_seq.html
# make batch file for ComBat-seq

library(sva)

batch <- data.frame(sample = colnames(TCGA_CCLE_counts_numeric),
                    batch = c(rep(1,length(grep("BREAST", colnames(TCGA_CCLE_counts_numeric)))), 
                              rep(2,length(grep("TCGA", colnames(TCGA_CCLE_counts_numeric))))))
TCGA_CCLE_counts_matrix <- data.matrix(TCGA_CCLE_counts_numeric)
postComBat <- ComBat_seq(TCGA_CCLE_counts_matrix, batch = batch$batch, group = NULL)
postComBat_ensembl <- data.frame(Name = rownames(postComBat),
                                 Description = TCGA_CCLE_counts$Description,
                                 postComBat)



# Normalization ------------------------------------------------------------------------------

# get gene length for calculating TPM
library(biomaRt)

mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
listAttributes(mart)

gene.annotations <- biomaRt::getBM(mart = mart, attributes=c("ensembl_gene_id", "start_position", "end_position"))
gene.annotations <- dplyr::transmute(gene.annotations, ensembl_gene_id, length = end_position - start_position)   

# excluding samples with total read counts below 40000000 and above 140000000
allread <- apply (postComBat_ensembl[,-c(1:2)],2,sum)
allread_sub <- allread[allread >= 40000000 & allread <= 140000000]


postComBat_sub_ensembl <- data.frame(ensembl_gene_id = postComBat_ensembl$Name,
                                     Description = postComBat_ensembl$Description,
                                     postComBat_ensembl[,colnames(postComBat_ensembl) %in% names(allread_sub)])

countfile_length <- inner_join(gene.annotations,
                               postComBat_sub_ensembl,
                               by = "ensembl_gene_id")



# Count TPM ---------------------------------------------------------------

library(edgeR)

counts_tpm <- (countfile_length[,4:ncol(countfile_length)]/countfile_length$length) *1000
counts_tpm <- cbind(ensembl_gene_id = countfile_length$ensembl_gene_id, counts_tpm) %>% as_tibble()

normfactor <- DGEList(counts = counts_tpm[,2:ncol(counts_tpm)],
                      group = colnames(counts_tpm[,2:ncol(counts_tpm)]))
normfactor <- calcNormFactors(normfactor, method="RLE")
normfactor_samples <- normfactor$samples
normfactor_samples$normlib <- normfactor_samples$lib.size*normfactor_samples$norm.factors


for(i in 1:(dim(counts_tpm)[2]-1)){
  counts_tpm[,i+1] <- (counts_tpm[,i+1]/normfactor_samples$normlib[i])*1000000
}
ensembl_symbol <- countfile_length[,c(1,3)]
counts_tpm_symbol <- inner_join(ensembl_symbol,
                                counts_tpm,
                                by = "ensembl_gene_id")

write.csv(counts_tpm_symbol, "TPM_RLE_postComBat.csv")


# counts_tpm_symbol_log2 <- data.frame(ensembl_gene_id = counts_tpm_symbol$ensembl_gene_id,
#                                      Description = counts_tpm_symbol$Description,
#                                      log2(counts_tpm_symbol[,-c(1,2)] + 1))

# write.csv(counts_tpm_symbol_log2, "TPM_RLE_log2_postComBat.csv")
