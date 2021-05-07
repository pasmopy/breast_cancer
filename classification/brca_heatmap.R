
#use following packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(grDevices))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(TCGAbiolinks))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(cluster))


#read arguments
args <- commandArgs(trailingOnly = TRUE)


#read csv file of features
path <- getwd()

for(i in 1:length(list.files(path, pattern = "csv"))){
  if (i == 1){
    temp <- list.files(path, pattern = "csv")[i]
    read_temp <- data.frame(lapply(temp, read.csv))
    colnames(read_temp)[-1] <- paste(str_sub(temp, end = -5), colnames(read_temp)[-1], sep = "_")
    dynamics_info <- read_temp
  } else {
    temp <- list.files(path, pattern = "csv")[i]
    read_temp <- data.frame(lapply(temp, read.csv))
    colnames(read_temp)[-1] <- paste(str_sub(temp, end = -5), colnames(read_temp)[-1], sep = "_")
    dynamics_info <- cbind(dynamics_info, read_temp[,-1])
  }
}



#read clinical information of TCGA-BRCA
BRCA_path_subtypes <- TCGAquery_subtype(tumor = "brca")
BRCA_path_subtypes_annotation <- data.frame(patient = gsub("-", "_", BRCA_path_subtypes$patient),
                                            days_to_death = BRCA_path_subtypes$days_to_death,
                                            subtype = BRCA_path_subtypes$BRCA_Subtype_PAM50)


#fix data frame of annotation for heatmap
annotation <- c()
for(i in 1:nrow(dynamics_info)){
  temp <- BRCA_path_subtypes_annotation[BRCA_path_subtypes_annotation$patient == 
                                          str_sub(dynamics_info$Sample, end = -5)[i],]
  annotation <- rbind(annotation, temp)
}


suppressWarnings(annotation$days_to_death <- as.numeric(annotation$days_to_death))
annotation[is.na(annotation)] <- 0


#add prognosis information
for (i in 1:nrow(annotation)){
  if (annotation$days_to_death[i] == 0){
    annotation$prognosis[i] <- 20
  } else if (annotation$days_to_death[i] < 365 & annotation$days_to_death[i] > 0) {
    annotation$prognosis[i] <- 1
  } else if (annotation$days_to_death[i] > 366 & annotation$days_to_death[i] < 730){
    annotation$prognosis[i] <- 2
  } else if (annotation$days_to_death[i] > 731 & annotation$days_to_death[i] < 1095){
    annotation$prognosis[i] <- 3
  } else if (annotation$days_to_death[i] > 1096 & annotation$days_to_death[i] < 1825){
    annotation$prognosis[i] <- 4
  } else if (annotation$days_to_death[i] > 1826 & annotation$days_to_death[i] < 2190){
    annotation$prognosis[i] <- 5
  } else if (annotation$days_to_death[i] > 2191 & annotation$days_to_death[i] < 2555){
    annotation$prognosis[i] <- 6
  } else if (annotation$days_to_death[i] > 2556 & annotation$days_to_death[i] < 2920){
    annotation$prognosis[i] <- 7
  } else if (annotation$days_to_death[i] > 2921 & annotation$days_to_death[i] < 3650){
    annotation$prognosis[i] <- 8
  } else if (annotation$days_to_death[i] > 3651 & annotation$days_to_death[i] < 4015){
    annotation$prognosis[i] <- 9
  } else if (annotation$days_to_death[i] > 4016 & annotation$days_to_death[i] < 4380){
    annotation$prognosis[i] <- 10
  } else if (annotation$days_to_death[i] > 4381){
    annotation$prognosis[i] <- 11
  }
}


#heatmap
dynamics_info_zscore <- dynamics_info  %>% column_to_rownames(var = "Sample") %>% scale()
km <- pam(as.matrix(dynamics_info_zscore), metric = "euclidean", k = as.numeric(args[1]))

gradPal <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(7)
set.seed(1234)
hm <- Heatmap(as.matrix(dynamics_info_zscore), show_row_names = F, name = "Z-score",
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "ward.D2",
              row_split = km$clustering,
              na_col = "black",
              col = colorRamp2(c(-3:3), gradPal),
              row_dend_reorder = TRUE,
              cluster_row_slices = F) +
  Heatmap(annotation$prognosis, name = "Prognosis",
          show_row_names = F,
          width = unit(0.5,"cm"),
          col = viridis(100)) +
  Heatmap(annotation$subtype, name = "Subtype",
          show_row_names = F,
          width = unit(0.5,"cm"),
          col = c("LumA" = "#304698", "LumB" = "#68c6ea", "Her2" = "#f489b9", "Basal" = "#e02828", "Normal" = "#00882c")) 


#save heatmap as pdf file
hw <- as.numeric(strsplit(args[2], ",")[[1]])
pdf("Heatmap.pdf", height = hw[1], width = hw[2])
hm
invisible(dev.off())

