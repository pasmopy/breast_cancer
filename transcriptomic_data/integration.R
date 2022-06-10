
library(dplyr)
library(edgeR)
library(sva)
library(tibble)
library(data.table)
library(stringr)
library(TCGAbiolinks)
library(biomaRt)
library(cluster)




# Download TCGA clinical/subtype information ------------------------------
#outputclinical()
outputClinical <- function(type){
  clinical <<- GDCquery_clinic(paste("TCGA", type, sep="-"), type = "clinical")
  write.csv(clinical, paste0(type, "_clinical.csv"))
}

#outputsubtype()
outputSubtype <- function(type){
  subtype <<- TCGAquery_subtype(paste(tolower(type)))
  write.csv(subtype, paste0(type, "_subtype.csv"))
}



# Select samples in reference to clinical or subtype data -----------------
patientSelection <- function(type,ID,
                            ...){
  object <- as.list(substitute(list(...)))
  condition <- rep(TRUE, nrow(type))
  for(obj in object[-c(1)]){
    condition <- condition & eval(obj, type)
  }
  selected_samples <<- type[condition, ID][!is.na(type[condition, ID])]
  message("********************************************************\n",
          paste(" patientSelection(): Number of selected samples is ", length(selected_samples),"\n"),
          "********************************************************")
}



# Download TCGA gene expression data (HTSeq-Counts) -----------------------
downloadTCGA <- function(cancertype, sampletype,
                         outputresult = FALSE){
  query <- GDCquery(project = paste("TCGA", cancertype, sep = "-"),
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts")
  cancertype_TCGA <<- cancertype
  GDCdownload(query)
  df<- GDCprepare(query, save = F, summarizedExperiment = F)
  
  df_TCGA <- df[grep("ENSG", df$gene_id),]
  df_TCGA<-as.data.frame(df_TCGA,check.names=F)
  pattern <- paste0("-",
                    paste0(sampletype, "A"),
                    "-")
  TCGA_counts <- df_TCGA[,grepl(paste(pattern, collapse = "|"), colnames(df_TCGA))]
  TCGA_counts <- TCGA_counts[,grepl("unstranded",colnames(TCGA_counts))]
  colnames(TCGA_counts) <- str_sub(colnames(TCGA_counts), start = 1, end = 27)
  TCGA_counts <- TCGA_counts[,!duplicated(colnames(TCGA_counts))]
  TCGA_counts$gene_id <- str_sub(df_TCGA$gene_id, start = 1, end = 15)
  TCGA_counts_filtered <<- data.frame(Name = TCGA_counts$gene_id,
                                      TCGA_counts[,grep(paste(selected_samples, collapse = "|"),
                                                        colnames(TCGA_counts))],check.names=F)
  colnames(TCGA_counts_filtered)<-gsub("-","_",colnames(TCGA_counts_filtered))
  if (outputresult == TRUE){
    fwrite(TCGA_counts_filtered, "TCGA_counts.csv")
  }
  message("********************************************************\n",
          paste("  downloadTCGA(): Number of selected samples is ", ncol(TCGA_counts_filtered)-1,"\n"),
          "********************************************************")
}




# Download CCLE transcriptomic data ---------------------------------------
downloadCCLE <- function(cancertype,
                         outputresult = FALSE){
  url <- "https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_counts_20180929.gct.gz"
  tmp <- tempfile()
  download.file(url, tmp)
  cancertype_CCLE <<- cancertype
  CCLE_counts <- read.table(gzfile(tmp), skip = 2, header = T)
  CCLE_cancer <<- data.frame(Name = str_sub(CCLE_counts$Name, start = 1, end = 15),
                             Description = CCLE_counts$Description,
                             CCLE_counts[,grep(cancertype, colnames(CCLE_counts))])
  message("********************************************************\n",
          paste("  downloadCCLE(): Number of selected samples is ", ncol(CCLE_cancer)-2,"\n"),
          "********************************************************")
  if (outputresult == TRUE){
    write.csv(CCLE_cancer, "CCLE_counts.csv")
  }
}


# Combine TCGA and CCLE data ----------------------------------------------
mergeTCGAandCCLE <- function(outputresult = FALSE){
  TCGA_CCLE_counts <- inner_join(CCLE_cancer, TCGA_counts_filtered, 
                                 by = "Name")
  TCGA_CCLE_counts <- TCGA_CCLE_counts %>%
    dplyr::distinct(Description, .keep_all=T) %>% column_to_rownames(var = "Name")
  TCGA_CCLE_counts_numeric <- sapply(TCGA_CCLE_counts[,-1], as.numeric)
  rownames(TCGA_CCLE_counts_numeric) <- rownames(TCGA_CCLE_counts)
  batch <- data.frame(sample = colnames(TCGA_CCLE_counts_numeric),
                      batch = c(rep(1,length(grep("BREAST", colnames(TCGA_CCLE_counts_numeric)))), 
                                rep(2,length(grep("TCGA", colnames(TCGA_CCLE_counts_numeric))))),check.names=F)
  TCGA_CCLE_counts_matrix <- data.matrix(TCGA_CCLE_counts_numeric)
  postComBat <- ComBat_seq(TCGA_CCLE_counts_matrix, batch = batch$batch, group = NULL)
  postComBat_ensembl <<- data.frame(Name = rownames(postComBat),
                                    Description = TCGA_CCLE_counts$Description,
                                    postComBat,check.names=F)
  allread <- apply(postComBat_ensembl[,-c(1:2)],2,sum)
  write.csv(data.frame(Sample = names(allread),
                       total_read_counts = allread) %>% remove_rownames(),
            "totalreadcounts.csv",
            row.names = F)
  colnames(postComBat_ensembl)<-gsub("-","_",colnames(postComBat_ensembl))
  if (outputresult == TRUE){
    write.csv(postComBat_ensembl, "merged_TCGA_CCLE.csv")
  }
}



# Normalization -----------------------------------------------------------
normalization <- function(min, max){
  mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  gene.annotations <- biomaRt::getBM(mart = mart, attributes=c("ensembl_gene_id", "start_position", "end_position"))
  gene.annotations <- dplyr::transmute(gene.annotations, ensembl_gene_id, length = end_position - start_position)
  allread <- apply (postComBat_ensembl[,-c(1:2)],2,sum)
  allread_sub <- allread[allread >= min & allread <= max]
  postComBat_sub_ensembl <- data.frame(ensembl_gene_id = postComBat_ensembl$Name,
                                       Description = postComBat_ensembl$Description,
                                       postComBat_ensembl[,colnames(postComBat_ensembl) %in% names(allread_sub)],check.names=F)
  countfile_length <<- inner_join(gene.annotations,
                                  postComBat_sub_ensembl,
                                  by = "ensembl_gene_id")
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
  counts_tpm_symbol <<- inner_join(ensembl_symbol,
                                  counts_tpm,
                                  by = "ensembl_gene_id")
  colnames(counts_tpm_symbol)<-gsub("-","_",colnames(counts_tpm_symbol))
  fwrite(counts_tpm_symbol, paste0("TPM_RLE_postComBat_", cancertype_TCGA, "_", cancertype_CCLE,".csv"))
  message("********************************************************\n",
          paste("  normalization(): Number of TCGA samples is ", length(grep("TCGA", colnames(counts_tpm_symbol))),"\n"),
          paste("  normalization(): Number of CCLE samples is ", ncol(counts_tpm_symbol)-length(grep("TCGA", colnames(counts_tpm_symbol))),"\n"),
          "********************************************************")
}

