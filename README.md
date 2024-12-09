# **Transcriptomic-profiling-of-Hepatocellular-carcinoma-HCC-dataset-from-TCGA**

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

## **Results**

- Identified distinct gene modules representing co-expressed genes in HNSC.
- The pink module exhibited differential expression patterns between tumor and normal samples, potentially highlighting pathways affected in HNSC.
- Insights obtained are crucial for understanding HNSC at the molecular level and can guide future research.
