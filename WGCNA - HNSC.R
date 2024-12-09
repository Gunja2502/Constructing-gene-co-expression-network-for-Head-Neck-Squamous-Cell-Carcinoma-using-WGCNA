library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

allowWGCNAThreads()          # allow multi-threading (optional)

# 1. Fetch Data from TCGA ------------------------------------------------
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  BiocManager::install("TCGAbiolinks")
}
library(TCGAbiolinks)

query_TCGA_HNSC <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification")
#get a full list of barcodes
samplesDown <- getResults(query_TCGA_HNSC, cols=c("cases"))

#identify barcodes with Tumor Primary samples
dataSmTP <- TCGAquery_SampleTypes(
  barcode = samplesDown,
  typesample = "TP"
)

#identify barcodes with Normal Tissue
dataSmNT <- TCGAquery_SampleTypes(
  barcode = samplesDown,
  typesample = "NT"
)

#Choose a small subset of 10 from each.
#Omit this step for full analysis.
dataSmTP_short <- dataSmTP[1:10]
dataSmNT_short <- dataSmNT[1:10]

#Create a query for only the selected samples
query.selected.samples <- GDCquery(
  project = "TCGA-HNSC", 
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts", 
  barcode = c(dataSmTP_short, dataSmNT_short) #only for a subset of barcodes
)

#Download the data
GDCdownload(query = query.selected.samples)
##Confirm in the console if the download starts

#prepare a summary of all the data
# Prepare a summary of all the data using the temp directory
dataPrep <- GDCprepare(
  query = query.selected.samples,
  save = TRUE,
  save.filename = "TCGA-HNSCTranscriptome_ProfilingSun_Nov_24_22_01_44_2024.RData", # specified the filename to solve the error
  summarizedExperiment = TRUE  # this ensures proper data structure
)

data_counts <- TCGAanalyze_Preprocessing(
  object = dataPrep, 
  cor.cut = 0.6,  # filter out if correlation between samples is below this threshold
  datatype = "unstranded",
  filename = file.path( "preprocessing_plot.pdf"),  # save plot to temp directory
  width = 1000, 
  height = 1000
)  

# detect outlier genes

gsg <- goodSamplesGenes(t(data_counts))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detected as outliers
data_counts <- data_counts[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering
htree <- hclust(dist(t(data_counts)), method = "average")
plot(htree)

# PCA for detecting sample outliers
pca <- prcomp(t(data_counts))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

### NOTE: If there are batch effects observed, correct for them before moving ahead

##If no outlier samples to be excluded
data.subset <- data_counts

# 3. Normalization (using variance stabilization)
# create a deseq2 dataset
# get the required colData table for dds

reqBarcodes <- colnames(data.subset)
req_clin_info <- c("shortLetterCode", "tumor_grade", "gender", "alcohol_history")
req_colData <- as.data.frame(colData(dataPrep)[reqBarcodes, req_clin_info])

# create dds
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = req_colData,
                              design = ~ 1) # not specifying model

## remove all genes with counts < 15 in more than 75% of samples (31*0.75=23.25)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 15,]
nrow(dds75) # 13284 genes

# perform variance stabilization
dds_norm <- vst(dds75)

# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices

# visualization to pick power
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 16
temp_cor <- cor
cor <- WGCNA::cor

# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

cor <- temp_cor

# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs

# Print out a preview
head(module_eigengenes)

# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# grey module = all genes that doesn't fall into other modules were assigned to the grey module

##Get gene symbols and names from IDs
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)

# List of Ensembl IDs
col_oi <- "cyan"
ensembl_ids <- names(which(bwnet$colors==col_oi))
ensembl_ids <- sapply(ensembl_ids, function(x) strsplit(x, split = "\\.")[[1]][1])

# Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",   # Retrieve HGNC symbols
  keytype = "ENSEMBL", # Input type
  multiVals = "first"  # If multiple matches, return the first
)

# Map Ensembl IDs to gene names
gene_names <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "GENENAME", # Retrieve gene names/descriptions
  keytype = "ENSEMBL", # Input type
  multiVals = "first"  # If multiple matches, return the first
)

# Combine results into a data frame
gene_info <- data.frame(
  Ensembl_ID = ensembl_ids,
  Gene_Symbol = gene_symbols,
  Gene_Name = gene_names,
  stringsAsFactors = FALSE
)

# View the results
View(gene_info)

## expression patterns of module genes
expression_oi <- norm.counts[,names(which(bwnet$colors==col_oi))]
colLabels_oi <- apply(req_colData[rownames(req_colData) %in% rownames(norm.counts), ], 1, paste, collapse = "_")
heatmap(margins = c(15,15), x = t(expression_oi), labCol = colLabels_oi, labRow = gene_info$Gene_Name)

# Set row names of the expression matrix to gene symbols
rownames(expression_oi) <- gene_info$Gene_Symbol

# Create a color palette for the heatmap
heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# Plot the heatmap
pheatmap(
  expression_oi,
  color = heatmap_colors,
  scale = "row",                # Normalize expression by row (gene-wise z-score)
  annotation_col = req_colData, # Add sample metadata annotations
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = paste("Heatmap of Gene Expression -", col_oi, "Module")
)
