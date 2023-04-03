# Microglial-annotation
For annotating scRNA-seq data for microglial cells

## Workflow 

<img width="828" alt="image" src="https://user-images.githubusercontent.com/91852421/229467110-89ac1f96-93b5-4b16-861a-9e8fea0d4476.png">

1. Import and preprocess the data

2. Quality control and visualization 

3. Normalization, Scaling and PCA/integration 

> Note: Here should try different integration methods + Different sequenced between integration and PCA

4. Clustering and visualization 

5. Annotation and visualization

## Annotation breakdown

1. Import all markers and plot the heatmap for overlapping heatmap population distribution (High resolution)

2. Pseudo-bulk analysis: mean log expression data for each cluster -> Heatmap again (low resolution)

3. Selecting labels for each cluster based on the most expressed gene and their corresponding cell types

4. Try hirearchical classification on labelling (Heatmap)

> Discuss: is that a reliable way for labeling? Should we tune the clustering resolution to avoid inaccurate labelling?

5. Label and visualization again

6. Subsequent visulization of barplot (cluster percentage), vlnplot (expression profile), and ridgeplot (similar comparision data)

![image](https://user-images.githubusercontent.com/91852421/215492226-2d0a6f63-bf0c-4420-b0d7-bf9a457cb2fc.png)
Visualization of example datasets
