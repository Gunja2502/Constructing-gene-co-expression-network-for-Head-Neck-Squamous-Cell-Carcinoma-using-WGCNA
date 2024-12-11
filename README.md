# **Constructing-gene-co-expression-network-for-Head-Neck-Squamous-Cell-Carcinoma-using-WGCNA**

This repository contains the analysis and scripts used for constructing a weighted gene co-expression network for Head-Neck Squamous Cell Carcinoma (HNSC) using data from The Cancer Genome Atlas (TCGA). The project employs the Weighted Gene Co-expression Network Analysis (WGCNA) framework in R to identify gene modules and investigate their potential biological significance.

## **Overview**

Head-Neck Squamous Cell Carcinoma (HNSC) is a significant global health concern, arising from epithelial cells in the head and neck region. This study aims to explore gene expression patterns and identify co-expressed gene modules that may provide insights into the biological processes underlying HNSC.

---

## **Objectives**
- Perform quality control, normalization, and filtering of RNA-seq data from TCGA.
- Construct a co-expression network using WGCNA.
- Identify gene modules and visualize their eigengenes.
- Generate and interpret heatmaps to uncover patterns in gene expression.
- Investigate potential biological associations between gene modules and clinical traits.

---

## **Key Features**

- Data Source: TCGA RNA-Seq gene expression data for HNSC.
- Tools: R packages such as WGCNA, DESeq2, TCGAbiolinks, tidyverse, org.Hs.eg.db, and others.
- Techniques:
   - Quality control (outlier detection using PCA and hierarchical clustering).
   - Normalization using variance stabilization (DESeq2).
   - Network construction with soft-thresholding to ensure scale-free topology.
   - Module eigengene analysis to summarize gene expression patterns.
 
 ---

## **Data source**

TCGA: Head and neck squamous cell carcinomas (HNSC)

Req_ColData (Data plotted on heatmap)
![Data HNSC](https://github.com/user-attachments/assets/5ae82375-1eb2-4847-bd1a-7db1e05b252a)

 ---

## **Project Workflow**

- Data Preparation:
   - RNA-Seq data fetched from TCGA using TCGAbiolinks.
   - Selected a subset of tumor primary (TP) and normal tissue (NT) samples for demonstration.
- Quality Control:
   - Outlier samples and genes detected and filtered using goodSamplesGenes() and PCA.
   - Hierarchical clustering to identify relationships among samples.
- Normalization:
   - Applied variance stabilization to RNA-Seq data.
   - Filtered genes with low expression for reliable results.
- Network Construction:
   - Soft-thresholding power selected to achieve scale-free topology.
   - Blockwise module detection to group co-expressed genes into modules.
- Heatmap Analysis:
   - Heatmaps of specific modules (e.g., turquoise, cyan, pink) visualized gene expression across samples.
   - Interpretation of patterns concerning tumor vs. normal tissue, gender, and alcohol history.
- Module Analysis:
   - Gene symbols and descriptions mapped to module genes.
   - Investigated biological significance of specific modules (e.g., pink module showing distinct patterns).
 
 --- 

## **Plots and parameters**

**Cluster Dendrogram**
![Cluster dendogram](https://github.com/user-attachments/assets/d4f21f9c-d22f-4391-b8c6-4400ff541870)

This is a cluster dendrogram, typically used in hierarchical clustering to display the relationships between samples or features based on a distance metric.

**PCA** 
![PCA](https://github.com/user-attachments/assets/8b40fc1c-f4d7-4613-b82e-9d3853e2f588) 

This graph helps identify outliers. The distant samples are outliers. 

**Soft thresholding power selection plot based on R²** 
![Soft threshold](https://github.com/user-attachments/assets/4ae907b1-4a84-40c0-bad2-705a0c86c51c)

In WGCNA, one crucial step is determining the soft-thresholding power to achieve a scale-free topology in the network (followed by most biological networks).  The goal is to identify a power where the network approximates a scale-free topology (fit index ≈ 0.85) while maintaining sufficient connectivity among genes.

Soft threshold I chose: 16

**List of colors for different modules** 
![Modules](https://github.com/user-attachments/assets/89c5117a-787b-4e85-b613-e48ee4862a5c)

For the analysis and heatmap plotting, I first chose turquoise but didn’t get a heatmap with a clear pattern. Then I randomly plotted for different colors. 

--- 

## **Results**

- Identified distinct gene modules representing co-expressed genes in HNSC.
- The pink module exhibited differential expression patterns between tumor and normal samples, potentially highlighting pathways affected in HNSC.
- Insights obtained are crucial for understanding HNSC at the molecular level and can guide future research.

**Heatmap: Pink module**
![Heatmap](https://github.com/user-attachments/assets/020983cd-8d0c-45f5-baf8-078857908118)

Sample Types:
- TP: Primary tumors.
- NT: Normal adjacent tissues.
- Gender and alcohol history are also noted for each sample, providing demographic and clinical stratification.

Gene Clustering: Genes are hierarchically clustered based on expression patterns. Some genes show distinct expression in certain clusters of samples, suggesting differential expression between tumor and normal samples or based on gender/alcohol consumption history.

The first 10 columns are for primary tumor type and the next 10 are for normal tumor type. 

Primary tumor gene expression pattern: There are three main blocks of expression pattern. There are sets of genes that are highly expressed (red/orange), and another set where another set of genes are less expressed comparatively (white/yellow/orange. 

Normal tumor gene expression pattern: In comparison to the TP columns, for the NT, gene expression is high overall. With a small block where gene expression is low. 

Overall, there are certain sets of genes that are highly expressed in normal tumors and those same genes are less expressed in cancer types. And this is a major trend. Biologically, this pattern could be because this cancer type could result in the downregulation of certain pathways that lead to the downregulation of genes involved in the particular pathway. 

--- 

##**Conclusion**

The gene expression heatmap analysis for Head-Neck Squamous Cell Carcinoma (HNSC) showed distinct patterns among different colors (modules). For certain modules, I could not draw a clear gene expression pattern. However, in the Pink module, there was a clear pattern in the primary tumor and normal tumor. There were sets of genes differentially expressed between the two. 
Further exploration using Weighted Gene Co-expression Network Analysis (WGCNA) highlights co-expressed gene modules, which could uncover meaningful biological associations. Integrating clinical data with gene expression profiles is crucial to understanding these discrepancies. Identifying patient-specific factors that contribute to this variability is essential for advancing personalized medicine approaches in HNSC.

---
