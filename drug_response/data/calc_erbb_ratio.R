
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(stringr))



# CCLE normalization ------------------------------------------------------

CCLEnormalization <- function(gene, sample = NULL){
  # download CCLE RNA-seq counts data
  url <- "https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_counts_20180929.gct.gz"
  tmp <- tempfile()
  download.file(url, tmp)
  CCLE_counts <- read.table(gzfile(tmp), skip = 2, header = T)
  CCLE_counts[,1] <- str_sub(CCLE_counts[,1], start=1, end=15)
  colnames(CCLE_counts)[1] <- "ensembl_gene_id"
  
  # get gene length
  mart <- useEnsembl(biomart = "ensembl", 
                     dataset = "hsapiens_gene_ensembl", 
                     mirror = "asia")
  gene.annotations <- biomaRt::getBM(mart=mart,
                                     attributes=c(
                                       "ensembl_gene_id", "start_position", "end_position"))
  gene.annotations <- dplyr::transmute(gene.annotations, ensembl_gene_id, 
                                       length = end_position - start_position)    
  
  # merge count file with length by ensembl_gene_id_version
  countsfile_length <- inner_join(
    gene.annotations,
    CCLE_counts,
    by = "ensembl_gene_id"
  )
  
  # data normalization
  counts_tpm <- (countsfile_length[,4:ncol(countsfile_length)]/countsfile_length$length)*1000
  counts_tpm <- cbind(Description=countsfile_length$Description, counts_tpm) %>% as_tibble()
  normfactor <- DGEList(counts=counts_tpm[,2:ncol(counts_tpm)],
                        group=colnames(counts_tpm[,2:ncol(counts_tpm)]))
  normfactor <- calcNormFactors(normfactor, method="RLE")
  normfactor_samples <- normfactor$samples
  normfactor_samples$normlib <- normfactor_samples$lib.size*normfactor_samples$norm.factors
  
  for(i in 1:(dim(counts_tpm)[2]-1)){
    counts_tpm[,i+1] <- (counts_tpm[,i+1]/normfactor_samples$normlib[i])*1000000
  } 
  
  if(is.null(sample)){
    counts_tpm_selected <- data.frame(Description = counts_tpm[counts_tpm$Description %in% gene,]$Description,
                                      counts_tpm[counts_tpm$Description %in% gene,-1]) %>% 
      column_to_rownames(var = "Description") %>% t() 
  } else {
    counts_tpm_selected <- data.frame(Description = counts_tpm[counts_tpm$Description %in% gene,]$Description,
                                      counts_tpm[counts_tpm$Description %in% gene ,colnames(counts_tpm) %in% sample]) %>% 
      column_to_rownames(var = "Description") %>% t() 
  }
  write.csv(counts_tpm_selected, "CCLE_normalized.csv")
}

CCLEnormalization(gene = gene)



# calculate receptor ratio ------------------------------------------------


receptor_ratio <- function(data, num){
  data <- log2(data + 1)
  CCLE_receptor <- data.frame(EGFR = data$EGFR,
                              ERBB234 = data$ERBB2 + data$ERBB3 + data$ERBB4,
                              ratio = data$EGFR/(data$ERBB2 + data$ERBB3 + data$ERBB4),
                              row.names = rownames(data))
  
  #filtering cell-lines (only cell-lines in CCLE drug response data)
  url_drug <- "https://data.broadinstitute.org/ccle_legacy_data/pharmacological_profiling/CCLE_NP24.2009_Drug_data_2015.02.24.csv"
  tmp_drug <- tempfile()
  download.file(url_drug, tmp_drug)
  CCLE_drug <- read.table(gzfile(tmp_drug), header = T, sep = ",")
  
  drug_cellline <- unique(CCLE_drug$CCLE.Cell.Line.Name)
  
  #filtering cell-lines (only cell-lines whose EGFR expression level is higher than the median expression level of EGFR in all cell-lines)
  CCLE_receptor_drug <- CCLE_receptor[rownames(CCLE_receptor) %in% drug_cellline,]
  CCLE_receptor_drug_median <- CCLE_receptor_drug[CCLE_receptor_drug$EGFR > median(CCLE_receptor_drug$EGFR),]
  CCLE_receptor_order <- CCLE_receptor_drug_median[order(CCLE_receptor_drug_median$ratio, decreasing = T),]
  
  #evaluation of cell-lines
  CCLE_receptor_order$value <- c(rep("high", num),
                                 rep("middle",nrow(CCLE_receptor_order) - 2*num),
                                 rep("low", num))
  
  write.csv(CCLE_receptor_order, "ErbB_expression_ratio.csv")
}


