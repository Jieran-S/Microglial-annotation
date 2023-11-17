# Microglial-annotation
For annotating scRNA-seq data for microglial cells

## Workflow 

![Single_cell_pipeline - Microglial Annotation](https://user-images.githubusercontent.com/91852421/229467316-5538a78e-ff2a-4811-a184-5e90aa13b3db.png)

1. Import and preprocess the data

2. Quality control and visualization 

3. Normalization, Scaling and PCA/integration 

4. Clustering and visualization 

5. Annotation 

6. Visualization

## Annotation breakdown

1. Import all markers and plot the heatmap for overlapping heatmap population distribution (High resolution)

2. Pseudo-bulk analysis: mean log expression data for each cluster -> Heatmap again (low resolution)

3. Summarizing each cell type based on the mean expression value for its corresponding marker genes

4. Try hirearchical classification on labelling (Heatmap)

5. Label and visualization again

6. Subsequent visulization of barplot (in case of larger sample sizes), vlnplot (expression profile), and ridgeplot (similar comparision data)

<img width="594" alt="P11_UMAP_AllCases" src="https://github.com/Jieran-S/Microglial-annotation/assets/91852421/e31ea752-a835-401a-9a12-6fb94d6f0be5">

