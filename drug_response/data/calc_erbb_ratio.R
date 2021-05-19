suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(stringr))




# download counts file
# https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_counts_20180929.gct.gz
url <- "https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_counts_20180929.gct.gz"
tmp <- tempfile()
download.file(url, tmp)
CCLE_counts <- read.table(gzfile(tmp), skip = 2, header = T)
CCLE_counts[,1] <- str_sub(CCLE_counts[,1], start = 1, end = 15)



# get RNA length 
dataset="hsapiens_gene_ensembl"
mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = dataset)
gene.annotations <- biomaRt::getBM(mart = mart, attributes=c("ensembl_gene_id", "start_position", "end_position"))
gene.annotations <- dplyr::transmute(gene.annotations, ensembl_gene_id, length = end_position - start_position)                        


# merge count file with length by ensembl_gene_id_version
colnames(CCLE_counts)[1] <- "ensembl_gene_id"
countfile_length <- data.frame(merge(gene.annotations,CCLE_counts))



# data normalization
counts_tpm <- (countfile_length[,4:ncol(countfile_length)]/countfile_length$length) *1000
counts_tpm <- cbind(Geneid = countfile_length$ensembl_gene_id, counts_tpm) %>% as_tibble()
normfactor <- DGEList(counts = counts_tpm[,2:ncol(counts_tpm)],
                      group = colnames(counts_tpm[,2:ncol(counts_tpm)]))
normfactor <- calcNormFactors(normfactor, method="RLE")
normfactor_samples <- normfactor$samples
normfactor_samples$normlib <- normfactor_samples$lib.size*normfactor_samples$norm.factors


for(i in 1:(dim(counts_tpm)[2]-1)){
  counts_tpm[,i+1] <- (counts_tpm[,i+1]/normfactor_samples$normlib[i])*1000000
}
counts_tpm[,1] <- countfile_length$Description


# log2(TPM+1) calculation
CCLE_symbol <- counts_tpm %>% distinct(Geneid, .keep_all=T) %>% column_to_rownames("Geneid") %>% t() %>% as.data.frame()
CCLE_symbol <- log2(CCLE_symbol + 1)


# receptor gene expression
CCLE_receptor <- data.frame(EGFR = CCLE_symbol$EGFR,
                            ERBB234 = CCLE_symbol$ERBB2 + CCLE_symbol$ERBB3 + CCLE_symbol$ERBB4,
                            ratio = CCLE_symbol$EGFR/(CCLE_symbol$ERBB2 + CCLE_symbol$ERBB3 + CCLE_symbol$ERBB4),
                            row.names = rownames(CCLE_symbol))




# download drug sensitivity file
# https://data.broadinstitute.org/ccle_legacy_data/pharmacological_profiling/CCLE_NP24.2009_Drug_data_2015.02.24.csv
url_drug <- "https://data.broadinstitute.org/ccle_legacy_data/pharmacological_profiling/CCLE_NP24.2009_Drug_data_2015.02.24.csv"
tmp_drug <- tempfile()
download.file(url_drug, tmp_drug)
CCLE_drug <- read.table(gzfile(tmp_drug), header = T, sep = ",")

drug_cellline <- unique(CCLE_drug$CCLE.Cell.Line.Name)

CCLE_receptor_drug <- CCLE_receptor[rownames(CCLE_receptor) %in% drug_cellline,]
CCLE_receptor_drug_median <- CCLE_receptor_drug[CCLE_receptor_drug$EGFR > median(CCLE_receptor_drug$EGFR),]
CCLE_receptor_order <- CCLE_receptor_drug_median[order(CCLE_receptor_drug_median$ratio, decreasing = T),]


CCLE_receptor_order$value <- c(rep("high", 30),
                               rep("middle",nrow(CCLE_receptor_order) - 2*30),
                               rep("low", 30))



write.csv(CCLE_receptor_order, file.path("data", "ErbB_expression_ratio.csv"))
