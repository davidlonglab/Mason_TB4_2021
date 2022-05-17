## scRNAseq analysis pipeline for Mason et al.
## Data derived from 2 x Ctrl (control) and 2 x ADR (doxorubicin) glomeruli scRNAseq datsets (C57Bl/6 background) See paper for experimental details: https://jasn.asnjournals.org/content/early/2020/07/10/ASN.2020020220.
## Daniyal Jafree & Gideon Pomeranz | 6th May 2022 | Version 5 | Acquisition of post-aligned files from GEO accession (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146912) and comparison of TB4 expression with statistics.
## Post-revision at Scientific reports with violin plots for other molecules within thymosin family

#----------------------------------------------------------------------------------------------------------------#

# Load up packages.
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(Matrix)
library(ggrepel)
library(patchwork)
library(tidyselect)

# Set working directory for things to save to. Please change as required.
setwd("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/glom_analysis/Will/New_analysis")

#----------------------------------------------------------------------------------------------------------------#

## PART 1: Read in count matrices and generate Seurat objects. Please download the .tar file from the GEO Accession. Here, we analyse control numbers 2 and 3 and doxorubicin 1 and 2.
# Load data from file (please change file location as required), deletes Ensembl IDs and makes gene names unique and create Seurat object for control number 2.
# Add 'project = "Ctrl"' to differentiate between control and ADR samples.
ctrl2.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/glom_analysis/Will/New_analysis/GSE146912_RAW/GSM4409508_control_2.txt.gz"),sep="\t", header = T, row.names = NULL)
ctrl2.data <- subset(ctrl2.data, select = -c(IGIS))
names <- make.unique(as.character(ctrl2.data$SYMBOL))
rownames(ctrl2.data) <- names
ctrl2.data <- ctrl2.data[,-1]
ctrl2 <- CreateSeuratObject(counts = ctrl2.data, min.cells = 2, min.features = 200, project = "Ctrl")
ctrl2

# Load data from file (please change file location as required), deletes Ensembl IDs and makes gene names unique and create Seurat object for control number 3.
# Add 'project = "Ctrl"' to differentiate between control and ADR samples.
ctrl3.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/glom_analysis/Will/New_analysis/GSE146912_RAW/GSM4409509_control_3.txt.gz"),sep="\t", header = T, row.names = NULL)
ctrl3.data <- subset(ctrl3.data, select = -c(IGIS))
name <- make.unique(as.character(ctrl3.data$SYMBOL))
rownames(ctrl3.data) <- names
ctrl3.data <- ctrl3.data[,-1]
ctrl3 <- CreateSeuratObject(counts = ctrl3.data, min.cells = 2, min.features = 200, project = "Ctrl")
ctrl3

# Merge control datasets and add cell IDs for each group.
Ctrl <- merge(ctrl2, y = ctrl3, add.cell.ids = c("Ctrl2", "Ctrl3"))
Ctrl

# Load data from file (please change file location as required), deletes Ensembl IDs and makes gene names unique and create Seurat object for ADR number 1.
# Add 'project = "Adr"' to differentiate between control and ADR samples.
adr1.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/glom_analysis/Will/New_analysis/GSE146912_RAW/GSM4409514_doxorubicin_1.txt.gz"),sep="\t", header = T, row.names = NULL)
adr1.data <- subset(adr1.data, select = -c(IGIS))
name <- make.unique(as.character(adr1.data$SYMBOL))
rownames(adr1.data) <- names
adr1.data <- adr1.data[,-1]
adr1 <- CreateSeuratObject(counts = adr1.data, min.cells = 2, min.features = 200, project = "Adr")
adr1

# Load data from file (please change file location as required), deletes Ensembl IDs and makes gene names unique and create Seurat object for ADR number 2.
# Add 'project = "Adr"' to differentiate between control and ADR samples.
adr2.data <- read.table(file=paste0("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Collaborations/glom_analysis/Will/New_analysis/GSE146912_RAW/GSM4409515_doxorubicin_2.txt.gz"),sep="\t", header = T, row.names = NULL)
adr2.data <- subset(adr2.data, select = -c(IGIS))
name <- make.unique(as.character(adr2.data$SYMBOL))
rownames(adr2.data) <- names
adr2.data <- adr2.data[,-1]
adr2 <- CreateSeuratObject(counts = adr2.data,  min.cells = 2, min.features = 200, project = "Adr")
adr2

# Merge ADR datasets and add cell IDs for each group.
Adr <- merge(adr1, y = adr2, add.cell.ids = c("Adr1", "Adr2"))
Adr

# Merge Seurat objects into one to create single, concatenated Seurat object.
Combined <- merge(Ctrl, y = Adr)
Combined

#----------------------------------------------------------------------------------------------------------------#

## PART 2: Quality control. Removes potential doublets and dead / low quality cells.
# Check mitochondrial content and features per barcode.
Combined[["percent.mt"]] <- PercentageFeatureSet(Combined, pattern = "^mt-")
VlnPlot(Combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1
plot2 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

# Check plots above and process data to exclude cells with less than 200 or more than 4000 transcripts and less than 10% mitochondrial transcripts.
Combined <- subset(Combined, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)

# Rerun mitochondrial content and features per barcode after QC.
Combined[["percent.mt"]] <- PercentageFeatureSet(Combined, pattern = "^mt-")
VlnPlot(Combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1
plot2 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

#----------------------------------------------------------------------------------------------------------------#

## PART 3: Data processing. Normalisation, finding variables genes, scaling and principle component analysis.
# Normalise, find variable features, scale data (based on all genes) and perform PCA.
Combined <- NormalizeData(Combined, normalization.method = "LogNormalize", scale.factor = 10000)
Combined <- FindVariableFeatures(Combined, selection.method = "vst", nfeatures = 2000)
top20 <- head(VariableFeatures(Combined), 20)
plot1 <- VariableFeaturePlot(Combined)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot2
all.genes <- rownames(Combined)
Combined <- ScaleData(Combined, features = all.genes)
Combined <- RunPCA(Combined, features = VariableFeatures(object = Combined))

# Draw elbow plot to assess how many PCs to use for analysis.
ElbowPlot(Combined)

#----------------------------------------------------------------------------------------------------------------#

## PART 4: Data integration using Harmony package. Alters PCs so that the same cell types in Ctrl vs Adr cluster together.
# Create plot for PCA, grouped by "orig.ident", and compare expression levels between  datatsets.
p1 <- DimPlot(object = Combined, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = Combined, features = "PC_1", group.by = "orig.ident", pt.size = .1)
p1
p2

# Run Harmony and perform batch correction, returns plot with number of iterations required for convergence.
Combined <- Combined %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)

# Access and show first five Harmony embeddings for each barcode.
Combined_harmony <- Embeddings(Combined, 'harmony')
Combined_harmony[1:5, 1:5]

# Replot PCA and expression comparison, this time after Harmony.
p1 <- DimPlot(object = Combined, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = Combined, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
p1
p2

#----------------------------------------------------------------------------------------------------------------#

## PART 5: Clustering, production of UMAP plots, differential gene expression between clusters and cell-type identification.
# Compute UMAP using Harmony embeddings by setting reduction technique to "harmony". Change number of dimensions and resolution of Louvain / Leiden clustering as desired.
Combined <- FindNeighbors(Combined, reduction = "harmony", dims = 1:9) 
Combined <- FindClusters(Combined, resolution = 0.4)
Combined <- RunUMAP(Combined, reduction = "harmony", dims = 1:9)

# Reorder identities so that Ctrl comes before Adr in all successive graphs.
my_levels <- c("Ctrl", "Adr")
Combined@meta.data$orig.ident <- factor(x = Combined@meta.data$orig.ident, levels = my_levels)

# Plot UMAP assigned by experimental group (Figure S1B).
p1 <- DimPlot(Combined, reduction = "umap", group.by = "orig.ident", pt.size = 0.5)
p1

# Plot UMAP assigned by cluster
p2 <- DimPlot(Combined, reduction = "umap", label = TRUE)
p2

# Calculate top 10 differentially expressed genes per cluster and save as CSV file in the current directory.
#Combined.markers <- FindAllMarkers(Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#Combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#top10 <- Combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#write.csv(Combined.markers,"celltypemarkers.csv", row.names = FALSE)

# Assign cell type identity to clusters after recognising them using DE genes and redo UMAP. Here, we group subtypes of major cell-types (e.g. glomerular endothelial cells or podocytes) together.
new.cluster.ids <- c("GEC", "GEC", "Mesangium", "Podocyte", "AEC", "PEC", "SMC", "Macrophage", "Monocyte", "TEC", "PEC", "T lymphocyte",  "Mesangium", "Podocyte")
names(new.cluster.ids) <- levels(Combined)
Combined <- RenameIdents(Combined, new.cluster.ids)

# Plot UMAP assigned by cell type (Figure 1A).
DimPlot(Combined, reduction = "umap", label = T, pt.size = 1, label.size = 2)

# Save / load RDS file for exchange with Gideon
#saveRDS(Combined, file = "/Users/daniyaljafree/Desktop/Collaborations/Will TB4 project/Will 2021 paper/Dataset/Combined.rds")
#Combined <- readRDS("/Users/daniyaljafree/Desktop/Collaborations/Will TB4 project/Will 2021 paper/Dataset/Combined.rds")

# Validate podocyte identity using Featureplots for Nphs1 and Nphs2 (Figure 1B).
FeaturePlot(Combined, features = "Nphs1", pt.size = 0.5)
FeaturePlot(Combined, features = "Nphs2", pt.size = 0.5)

# Draw dotplot for cell type-specific markers for each condition and featureplot for those markers across the whole dataset (Figure S1A).
glomerular.cell.markers <- c("Emcn", "Ehd3", "Ptn", "Pdgfrb", "Wt1", "Nphs2", "Fbln2", "Ptprc", "Cd52", "Pax8", "Acta2", "Myh11", "Cdh1", "Aqp2", "Cd3e")
DotPlot(Combined, features = glomerular.cell.markers)

# Retrieve number of cells per cluster by condition and saves as CSV file in current directory (use as input numbers for Figure S1C).
cell_numbers <- table(Combined@active.ident, Combined@meta.data$group)
cell_numbers
write.csv(cell_numbers, "cell_numbers.csv")


p1 <- VlnPlot(Combined, features = "Tmsb10", split.by = "orig.ident", pt.size = 0.5)
p1

#----------------------------------------------------------------------------------------------------------------#

## PART 6: Differential expression between Ctrl and Adr for Tmsb4x.
# Visualize Tmsb4x expression between all Ctrl and Adr podocytes using violin plot (Figure 1D)
p1 <- VlnPlot(Combined, features = "Tmsb4x", split.by = "orig.ident", idents = "Podocyte",  pt.size = 0.5)
p1

# Use FindAllMarkers function to compare cell type-specific TB4 between conditions. Write CSV file and save to current directory. Includes data used for Wilcoxon Rank Sum test used in Figure 1D.
compare_tb4_clusters <- FindAllMarkers(Combined, features = "Tmsb4x", logfc.threshold = 0)
compare_tb4_clusters
write.csv(compare_tb4_clusters, "compare_tb4_clusters.csv")

# Use FindAllMarkers function to compare cell type-specific TB10 between conditions. Write CSV file and save to current directory. Includes data used for Wilcoxon Rank Sum test used in Figure 1F.
compare_tb10_clusters <- FindAllMarkers(Combined, features = "Tmsb10", logfc.threshold = 0)
compare_tb10_clusters
write.csv(compare_tb4_clusters, "compare_tb10_clusters.csv")

# Subset all glomerular cells (podocyte, glomerular endothelial and mesangial cells). 
podo_endo_mes_subset <- subset(Combined, idents = c("GEC", "Mesangium", "Podocyte"))
podo_endo_mes_subset

# Violin plot of TB4 levels between ADR and control for these cell types (Figure 1C).
Idents(podo_endo_mes_subset) <- podo_endo_mes_subset@meta.data$orig.ident
DefaultAssay(podo_endo_mes_subset) <- "RNA"
p1 <- VlnPlot(podo_endo_mes_subset, features = "Tmsb4x", split.by = "orig.ident",  pt.size = 0.5)
p1

# Compare TB4 expression in grouped cells between conditions and save as CSV in current directory. Includes data used for Wilcoxon Rank Sum test used in Figure 1E.
compare_tb4_podo_endo_mes <- FindAllMarkers(podo_endo_mes_subset, features = "Tmsb4x", logfc.threshold = 0)
compare_tb4_podo_endo_mes
write.csv(compare_tb4_podo_endo_mes, "compare_tb4_podo_endo_mes.csv")


## END ##
