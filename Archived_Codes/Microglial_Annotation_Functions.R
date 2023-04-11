library(dplyr)


combineSeuratWOBatch <- function(seurat_list) {
  # Leaving only the intersected genes
  seurat_list <- lapply(seurat_list, function(x) rownames(x)) %>% 
    {Reduce(intersect, .)} %>%
    {lapply(seurat_list, function(x) x[.,])}
  
  # Roughly merged the genes together
  seurat_combi <- merge(seurat_list[[1]], seurat_list[-1], project = "CombinedSeurat") 
  return(seurat_combi)
}