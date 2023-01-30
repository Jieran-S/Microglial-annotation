# Microglial-annotation
For annotating scRNA-seq data for microglial cells

## Workflow 

1. Import and preprocess the data

2. Quality control and visualization 

3. Normalization, Scaling and PCA/integration 

> Note: Here should try different integration methods + Different sequenced between integration and PCA

4. Clustering and visualization 

-------- Annotation --------

5. Import all markers and plot the heatmap for overlapping heatmap population distribution (High resolution)

6. Pseudo-bulk analysis: mean log expression data for each cluster -> Heatmap again (low resolution)

7. Selecting labels for each cluster based on the most expressed gene and their corresponding cell types

> Discuss: is that a reliable way for labeling? Should we tune the clustering resolution to avoid inaccurate labelling?

8. Label and visualization again

9. Subsequent visulization of barplot (cluster percentage), vlnplot (expression profile), and ridgeplot (similar comparision data)
