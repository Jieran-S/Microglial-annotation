# -------------------Internship - Max von Kolczynski -------------------
# ------------sc-RNA analysis of microglia in Alzheimer context--------

#install packaged 
install.packages("tidyverse")
install.packaged("seurat")
install.packages("xlsx")
install.packages('BiocManager') # makes finding the markers more efficient 
BiocManager::install('limma')

#load the libraries
library(hdf5r) 
library(tidyverse) 
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(gridExtra)
library(openxlsx)
library(xlsx)
library(xlsxjars)

#############################################--- read the data and set up seurat objects ---############################

# three different samples: 100-WT, 101 - GpnmbKO, 102-GpnmbKO)
# Read10X() reads in the output of the cellranger pipeline from 10X, returning a unique molecular identified (UMI) count matrix
# The values in this matrix represent the number of molecules for each feature that are detected in each cell 

# INPUT: filtered feature matrix .h5
seurat_list <- lapply(list.files("Data/GSE98969_RAW", full.name =TRUE), 
  function(file){
    file_ind <- which(list.files("Data/GSE98969_RAW", full.name =TRUE) == file)
    seurat_object <- read.delim2(file) %>% 
      CreateSeuratObject(project = c("WT", "GpnmbKO_1", "GpnmbKO_2")[file_ind], min.cells = 3, min.features = 200)
    seurat_object$stim <- sprintf("cond%s", file_ind-1)
    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
    seurat_object[["percent.rb"]] <- PercentageFeatureSet(seurat_object, pattern = "^Rp[sl]")
    seurat_object %>% subset( subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 15 & percent.rb < 25 ) %>%
      NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 
    seurat_object
})

# visualization by violin plot
lapply(seurat_list, function(seurat_object){
  Vplot <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
  print(Vplot)
})

if (FALSE){
  cond0.data <- Read10X_h5(object = , use.names = TRUE, unique.features = TRUE)
  cond1.data <- Read10X_h5(object = , use.names = TRUE, unique.features = TRUE)
  cond2.data <- Read10X_h5(object = , use.names = TRUE, unique.features = TRUE)
  
  # set up Seurat objects. The object serves as a container that contains both data (like the count matrix) and analysis (like PCA, or clustering results) for a single-cell dataset
  # the count matrix is stored in pbmc[["RNA"]]@counts.
  cond0 <- CreateSeuratObject(counts = cond0.data, project = "condition0", min.cells = 3, min.features = 200) # you could rename project to WT    
  cond1 <- CreateSeuratObject(counts = cond1.data, project = "condition1", min.cells = 3, min.features = 200) # you could rename project to GpnmbKO_1  
  cond2 <- CreateSeuratObject(counts = cond2.data, project = "condition2", min.cells = 3, min.features = 200) # you could rename project to GpnmbKO_2 
  
  cond0$stim <- "cond0"  #WT                                     
  cond1$stim <- "cond1"  #Gnpmb KO                               
  cond2$stim <- "cond2"  #Gnpmn KO                             
  # what is this for again?? bilge knows 
  
  # View(cond0@meta.data)
  # View(cond1@meta.data)
  # View(cond2@meta.data)
  
  # calculate mitochondrial and ribsosmmal RNA percentage for each condition
  cond0[["percent.mt"]] <- PercentageFeatureSet(cond0, pattern = "^mt-")
  cond0[["percent.rb"]] <- PercentageFeatureSet(cond0, pattern = "^Rp[sl]")
  cond1[["percent.mt"]] <- PercentageFeatureSet(cond1, pattern = "^mt-")
  cond1[["percent.rb"]] <- PercentageFeatureSet(cond1, pattern = "^Rp[sl]")
  cond2[["percent.mt"]] <- PercentageFeatureSet(cond2, pattern = "^mt-")
  cond2[["percent.rb"]] <- PercentageFeatureSet(cond2, pattern = "^Rp[sl]")


############################################# --- QC and Filtering --- ##################################################### 

# filter the data seperately 7   additional values for n_Count? filter before normalize? 
cond0_filt  <- subset(cond0, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 15 & percent.rb < 25) 
cond1_filt  <- subset(cond1, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 15 & percent.rb < 25) 
cond2_filt  <- subset(cond2, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 15 & percent.rb < 25) 

# View(cond0@meta.data)
# View(cond1@meta.data)
# View(cond2@meta.data)

####################################################--- Integration ---#####################################################

# normalize and identify variable features for each dataset independently 
# Normalize: ensuring expression values across cells are on a comparable scale.
cond.list  = list(cond0_filt, cond1_filt, cond2_filt) 
cond.list <- lapply(X = cond.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000) 
})
}

# Perform Integration: select integration features, find anchor genes, integrate
seurat_integrated <- seurat_list %>% SelectIntegrationFeatures() %>% 
  {FindIntegrationAnchors(object.list = seurat_list, anchor.features = .)} %>% 
  IntegrateData()



features_cond.list <- SelectIntegrationFeatures(object.list = seurat_list)
anchors_cond.list  <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features_cond.list)
mg.combined        <- IntegrateData(anchorset = anchors_cond.list) 

#  warining message: unique cell names

# specify that we will perform downstream analysis on the corrected data. note that the original unmodified data still resides in the 'RNA' assay
DefaultAssay(mg.combined) <- "integrated"  

#View(mg.combined@meta.data)

######################################################### --- Analysis --- ###################################################

# (after filtering and integrating)

#PCA (Linear dimensional reduction) and Clustering 
all.genes  <- rownames(mg.combined)  
mg.combined <- ScaleData(mg.combined, features = all.genes)  
mg.combined <- RunPCA(mg.combined, features = VariableFeatures(object = mg.combined), npcs = 25) # 25 is a default value, dicussed with Bilge  
mg.combined <- FindNeighbors(mg.combined, dims = 1:25)
mg.combined <- FindClusters(mg.combined, resolution = 0.35) # changed the resoultion from 0.5 -> 0.35 (Paper)  important?   

# determine dimensionality if you like with ElbowPlot or JackStrawPlot

# Assign subpopulation cell type identities to clusters (here we don't know the subpopulations yet, so just microglia (mg) 1-13)
levels(mg.combined)
new.cluster.ids <- c("mg1", "mg2", "mg3", "mg4", "mg5", "mg6","mg7","mg8","mg9","mg10","mg11","mg12", "mg13")
names(new.cluster.ids) <- levels(mg.combined)
mg.combined <- RenameIdents(mg.combined, new.cluster.ids)

# run tSNE and UMAP  (non linear dimensional reduction)
mg.combined <- RunTSNE(object = mg.combined, dims = 1:25)
mg.combined <- RunUMAP(object = mg.combined, dims = 1:25)

# find differentially expressed features // only with the merged? split.by condition?  (s. Unten)
allMarkers <- FindAllMarkers(mg.combined, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
# only.pos = FALSE, otherwise mg2 is empty
allMarkers %>%
  group_by(cluster) %>%
  top_n(n = 25) # keep it 25 - like Bilge

# export to Exel 
write.xlsx(allMarkers, file= "AllMarkers.xlsx", row.names = FALSE)  #sep = "," ??? 

# Vector of well studied Genes: read about them (refer to papers!!!!)
heatmapMarkers <- c("Mrc1","Ms4a7","Pf4","Trem2","Spp1","Itgax","Cd83","Ifit1","Ifit2","Isg15","Ifitm3","Ccl3","Ccl4","Cst7","Il1b","Rtp4","Cd68","Lamp1","Usp18","Mki67","Birc5","Cd36","Cd74","H2-Aa")

# important: mrc1, trem2, cd68, Mki67, lamp1

# some features were omitted as they were not found in the scale.data slot for the integrated assay: 
# Cdk1, Ctsd, Cd9, Hexb, Csf1r, Tmem119, P2ry13, Cx3cr1, P2ry12 (from the well studied genes. I earsed them here)

#export to exel
write.table(heatmapMarkers, file = "heatmapMarkers.csv", row.names = FALSE, sep = ",") 
write.table(top25_VarFeat, file = "top25_VarFeat.csv", row.names = FALSE, sep = ",") 
# is this already ready? -> just the markers with no further entries

# SubsetGenes
heatmapMarkers1 <- c("Mrc1","Ms4a7","Pf4","Trem2")
heatmapMarkers2 <- c("Spp1", "Itgax","Cd83","Ifit1")
heatmapMarkers3 <- c("Ifit2","Isg15","Ifitm3","Ccl3")
heatmapMarkers4 <- c("Ccl4","Cst7","Il1b","Rtp4")
heatmapMarkers5 <- c("Cd68","Lamp1","Usp18","Mki67")
heatmapMarkers6 <- c("Birc5","Cd36","Cd74","H2-Aa")

heatmapMarkers1_3 <- c("Mrc1","Ms4a7","Pf4")
heatmapMarkers2_3 <- c("Trem2", "Spp1", "Itgax")
heatmapMarkers3_3 <- c("Cd83","Ifit1", "Ifit2")
heatmapMarkers4_3 <- c("Isg15","Ifitm3","Ccl3")
heatmapMarkers5_3 <- c("Ccl4","Cst7","Il1b")
heatmapMarkers6_3 <- c("Rtp4","Cd68","Lamp1")
heatmapMarkers7_3 <- c("Usp18","Mki67","Birc5")
heatmapMarkers8_3 <- c("Cd36", "Cd74","H2-Aa")

top15_VarFeat <- head(VariableFeatures(mg.combined), 15)
top20_VarFeat <- head(VariableFeatures(mg.combined), 20)
top25_VarFeat <- head(VariableFeatures(mg.combined), 25)

top15_VarFeat_subset1 <- c("Spp1", "H2-Eb1", "F13a1", "Mrc1")
top15_VarFeat_subset2 <- c("H2-Aa", "H2-Ab1", "Cd209f", "Cd74")
top15_VarFeat_subset3 <- c("Lyve1", "Cd209a", "Pf4", "Ly6c1") 
top15_VarFeat_subset4 <- c("Cldn5", "S100a6", "Aoah")

top15_VarFeat_subset1_3 <- c("Spp1", "H2-Eb1", "F13a1")
top15_VarFeat_subset2_3 <- c("Mrc1", "H2-Aa", "H2-Ab1")
top15_VarFeat_subset3_3 <- c("Cd209f", "Cd74", "Lyve1") 
top15_VarFeat_subset4_3 <- c("Cd209a", "Pf4", "Ly6c1") 
top15_VarFeat_subset5_3 <- c("Cldn5", "S100a6", "Aoah")

# #other Features:  Gm42418  Cst3  Slc2a5  Slco2b1  Gpr34  Sall1  Tmem119  P2ry13  P2ry12  Hpgds  Olfml3  Fcrls Elmo1 Apoe Etl4 Tmem119 Fth1 Cd74 Mt1 Cldn5


# percentatges of suspopulations by different condition 

#View(mg.combined@meta.data)

n_cellsM <- FetchData(mg.combined, vars = c("ident", "stim")) %>%
  dplyr::count(ident, stim) %>% 
  spread(ident, n) 

View(n_cellsM)   # cell number detected per cluster per condition

write.xlsx(n_cellsM, "CellNumberByCluster.xlsx") # basically the same as write.table with .csv? 

ncells  <- read.xlsx("CellNumberByCluster.xlsx", sheetIndex = 2)

View(ncells)

# Percentage of Cluster 
ggplot(ncells, aes(x = Cluster , y = Percentage, fill = cond)) +
  geom_col(position = "dodge") + labs(title = "Percentage of detcted cells per cluster by condition")

# ------  Total number of cells per condition ------ 
# cond0_total =  14648    cond1_total =  11427    cond2_total =  16323

# Notes
# once with seperate conditions, once with cond1 + cond2 (see in the exel sheets)
# Normalize -> normalize it based on how many cells I have in that certain condition 
# Formular:   %_mg(i) in cond j =  (cell number in cluster i) / total number of cells in cond_j 
# put percentages in a table
# plot the percentages between wild and alzheimer   


################################################# -- Plots -- ##################################################  

# NOTE: merge has nothing to do with "integrating". It just concatenates the data sets -> makes plotting easier


# unfiltered (cond0, cond1, cond2)

merged <- merge(cond0, y = c(cond1, cond2), 
                add.cell.ids = c("cond0", "cond1", "cond2"),
                project = 'mg_merged') 

# QC metrics Violin Plots 

VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, group.by = 'stim')

# if not using merged 
# QC_unfilt_c0 <- VlnPlot(cond0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
# QC_unfilt_c1 <- VlnPlot(cond1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
# QC_unfilt_c2 <- VlnPlot(cond2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)


# Feature Scatter Plots

plot_featureScatter_mt  <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot_featureScatter_rb  <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.rb")
plot_featureScatter_RNA <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_featureScatter_RNA + plot_featureScatter_mt +plot_featureScatter_rb 


# filtered (cond0_filt, cond1_filt, cond2_filt)

merged_filt <- merge(cond0_filt, y = c(cond1_filt, cond2_filt), 
                     add.cell.ids = c("cond0_filt", "cond1_filt", "cond2_filt"),
                     project = 'mg_merged_filt') 

# QC metrics Violin plots

VlnPlot(merged_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, group.by = 'stim')

# VlnPlot(cond0_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
# VlnPlot(cond1_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
# VlnPlot(cond2_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)


# Feature Scatter Plots

plot_featureScatter_mt_f  <- FeatureScatter(merged_filt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot_featureScatter_rb_f  <- FeatureScatter(merged_filt, feature1 = "nCount_RNA", feature2 = "percent.rb")
plot_featureScatter_RNA_f <- FeatureScatter(merged_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_featureScatter_RNA_f + plot_featureScatter_mt_f + plot_featureScatter_rb_f 


# filtered and integrated data 

# PCA visualization 
# Examine and Visualize PCA results a few different ways  -> Seurat Vignette Guided Clustering
print(mg.combined[["pca"]], dims = 1:5, nfeatures = 5) 
VizDimLoadings(mg.combined, dims = 1:4, reduction = "pca") 
DimPlot(mg.combined, reduction = "pca") 
DimHeatmap(mg.combined, dims = 1:9, cells = 500, balanced = TRUE) 

# t-SNE (with integrated data) 
tSNE_all        <- DimPlot(mg.combined, reduction = "tsne", label = TRUE, pt.size = 0.2) + labs(title = "microglial cell subpopulations clusterd on a two-dimensional tSNE") 
tSNE_by_cond    <- DimPlot(mg.combined, reduction = "tsne", split.by = "orig.ident", label = TRUE, pt.size = 0.5) + labs(title = "microglial cell subpopulations clusterd on a two-dimensional tSNE \n- grouped by condition ")# tSNE of microglial cell clusterstSNE #by stim // from three different mice
tSNE_all_byCond <- DimPlot(mg.combined, reduction = "tsne", group.by = 'orig.ident', label = FALSE) + labs(title = "tSNE - grouped by condition") # "three different conditiosn "by stim? //   from three different mice

# UMAP (with integrated data) 
UMAP_all        <- DimPlot(mg.combined, reduction = "umap", label = TRUE, pt.size = 0.2) + labs(title = "microglial cell subpopulations clusterd on a two-dimensional UMAP") #UMAP of microglia clusters 1-13
UMAP_by_cond    <- DimPlot(mg.combined, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.5) + labs(title = "microglial cell subpopulations clusterd on a two-dimensional UMAP \n- grouped by sample from three different mice")# tSNE of microglial cell clusterstSNE #by stim 
UMAP_all_byCond <- DimPlot(mg.combined, reduction = "umap", group.by = 'orig.ident', label = FALSE) + labs(title = "UMAP - grouped by sample from three different mice") # "three different conditiosn "by stim? 


# create Heat-Maps 

allMarkers %>%
  group_by(cluster) %>% 
  top_n(n = 40,  wt = avg_log2FC) -> top40
HeatMap1 <- DoHeatmap(mg.combined, features = top40$gene, angle = 90, ) + NoLegend()  # additional flags?  # dont show feature names
# I changed 60 to 40 for computational issues
HeatMap2 <- DoHeatmap(mg.combined, features = heatmapMarkers, angle = 90) + NoLegend() # additional flags? 
HeatMap3 <- DoHeatmao(mg.combined, features = top15_VarFeat, angle = 90) + NoLegend() 


# Feature Scatter Plots of marker-expression in clusters

# for well studied genes (you may wanna change the features: change the vectors)
FeaturePlot1 <- FeaturePlot(mg.combined, features = heatmapMarkers1)
FeaturePlot2 <- FeaturePlot(mg.combined, features = heatmapMarkers2)
FeaturePlot3 <- FeaturePlot(mg.combined, features = heatmapMarkers3)
FeaturePlot4 <- FeaturePlot(mg.combined, features = heatmapMarkers4)
FeaturePlot5 <- FeaturePlot(mg.combined, features = heatmapMarkers5)
FeaturePlot6 <- FeaturePlot(mg.combined, features = heatmapMarkers6)

# for high variable features (you may wanna change the features: change the vectors)
FeaturePlot_HVF1 <- FeaturePlot(mg.combined, features = top15_VarFeat_subset1)
FeaturePlot_HVF2 <- FeaturePlot(mg.combined, features = top15_VarFeat_subset2)
FeaturePlot_HVF3 <- FeaturePlot(mg.combined, features = top15_VarFeat_subset3)
FeaturePlot_HVF4 <- FeaturePlot(mg.combined, features = top15_VarFeat_subset4)


# Violoin Plots of marker-expression in clusters 

# well studied genes (you may wanna change the features: change the vectors)
VlnPlot1_ <- VlnPlot(mg.combined, features = heatmapMarkers1_3) 
VlnPlot2_ <- VlnPlot(mg.combined, features = heatmapMarkers2_3) 
VlnPlot3_ <- VlnPlot(mg.combined, features = heatmapMarkers3_3) 
VlnPlot4_ <- VlnPlot(mg.combined, features = heatmapMarkers4_3) 
VlnPlot5_ <- VlnPlot(mg.combined, features = heatmapMarkers5_3) 
VlnPlot6_ <- VlnPlot(mg.combined, features = heatmapMarkers6_3) 
VlnPlot7_ <- VlnPlot(mg.combined, features = heatmapMarkers7_3) 
VlnPlot8_ <- VlnPlot(mg.combined, features = heatmapMarkers8_3) 

# high variable features (you may wanna change the features: change the vectors)
VlnPlot_HFV1_ <- VlnPlot(mg.combined, features = top15_VarFeat_subset1_3) 
VlnPlot_HFV2_ <- VlnPlot(mg.combined, features = top15_VarFeat_subset2_3) 
VlnPlot_HFV3_ <- VlnPlot(mg.combined, features = top15_VarFeat_subset3_3) 
VlnPlot_HFV4_ <- VlnPlot(mg.combined, features = top15_VarFeat_subset4_3) 
VlnPlot_HFV5_ <- VlnPlot(mg.combined, features = top15_VarFeat_subset5_3)


# mean logcounts by cluster:
pb <- aggregateData(sce, "logcounts", by=c("cluster"), fun="mean")

assayNames(pb) <- "logcounts"
rowData(pb)$marker4 <- NA
rowData(pb)[unlist(km),"marker4"] <- rep(names(km),lengths(km))
sechm::sechm(pb, unlist(km), assayName = "logcounts", gaps_row = "marker4", show_colnames = TRUE, do.scale = TRUE, breaks=1, row_title_rot=0)

# build a heatmap of the mean logcounts of the known markers:
h1 <- pheatmap(assay(pb)[unlist(km),], annotation_row=data.frame(row.names=unlist(km), type=rep(names(km), lengths(km))), split=rep(names(km), lengths(km)), cluster_rows=FALSE, scale="row")

# heatmap for the de-novo markers:
h2 <- pheatmap(assay(pb)[markers,], scale="row")

# we will assign clusters to the cell type whose markers are the most expressed
# we get rid of the unspecific neuronal markers:
km <- km[names(km)!="neuronal"]
# we extract the pseudo-bulk counts of the markers for each cluster
mat <- assay(pb)[unlist(km),]
# we aggregate across markers of the same type
mat <- aggregate(t(scale(t(mat))), by=list(type=rep(names(km), lengths(km))), FUN=sum)
# for each column (cluster), we select the row (cell type) which has the maximum aggregated value
cl2 <- mat[,1][apply(mat[,-1], 2, FUN=which.max)]
# we convert the cells' cluster labels to cell type labels:
sce$cluster2 <- cl2[sce$cluster]

# we aggregate again to pseudo-bulk using the new clusters
pb <- aggregateData(sce, "logcounts", by=c("cluster2"), fun="mean")
# we plot again the expression of the markers as a sanity check
h1 <- pheatmap(assay(pb)[unlist(km),], annotation_row=data.frame(row.names=unlist(km), type=rep(names(km), lengths(km))), split=rep(names(km), lengths(km)), cluster_rows=FALSE, scale="row")

# we aggregate by cluster x sample to perform pseudo-bulk differential state analysis
sce <- muscat::prepSCE(sce, kid="cluster2")
pb <- aggregateData(sce)
assays(pb)

# top genes in a given cell type
pbHeatmap(sce, res, k="astrocytes")


######################## END ############################
