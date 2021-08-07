
#value_transformation

value_transformation <- function(upstream, origin, gene, cellline, outputfile){
  optimizeddata_genelist <- upstream[upstream$Description %in% genelist,]
  newdata_genelist <- origin[origin$Description %in% genelist,]
  
  newdata_converted <- c()
  for(i in 1:length(gene)){
    temp_gene <- optimizeddata_genelist$Description[i]
    max_optimized <- max(optimizeddata_genelist[i, colnames(optimizeddata_genelist) %in% cellline])
    
    column_new <- which(newdata_genelist$Description == temp_gene)
    max_new <- max(newdata_genelist[column_new, colnames(newdata_genelist) %in% cellline])
    
    newdata_converted_row <- c(newdata_genelist[column_new,c(1,2)],newdata_genelist[column_new,-c(1,2)]*max_optimized/max_new)
    newdata_converted <- data.frame(rbind(newdata_converted, newdata_converted_row))
  }
  rownames(newdata_converted) <- c(1:length(gene))
  fwrite(newdata_converted, paste(outputfile))
}


