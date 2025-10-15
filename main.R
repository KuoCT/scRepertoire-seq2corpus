# Setup working environment
working_dir <- file.path(this.path::here(),"..")
setwd(working_dir)

source("scripts/setup.R")
source("scripts/functions.R")

set.seed(1234)

# Create out2 folder
exp_name <- tools::file_path_sans_ext(basename(this.path::this.path()))
out2 <- file.path(out, exp_name)
if (!dir.exists(out2)) dir.create(out2)

# scRNA-seq ----
## Import data ----
h5_S1_HT7 <- Read10X_h5("data/20240419_scRNAseq/demultiplexed_samples/outs/per_sample_outs/HT7/count/sample_filtered_feature_bc_matrix.h5")
h5_S1_HT8 <- Read10X_h5("data/20240419_scRNAseq/demultiplexed_samples/outs/per_sample_outs/HT8/count/sample_filtered_feature_bc_matrix.h5")
h5_S1_HT7 <- h5_S1_HT7[["Gene Expression"]]
h5_S1_HT8 <- h5_S1_HT8[["Gene Expression"]]

h5_S3_unsort <- Read10X_h5("data/20250326_scRNAseq/Unsort_22HNT5LT4/outs/filtered_feature_bc_matrix.h5")
h5_S3_sorted <- Read10X_h5("data/20250326_scRNAseq/Sort_22HNT5LT4/outs/filtered_feature_bc_matrix.h5")

## Modify barcode suffix ----
colnames(h5_S1_HT8) <- sub("-1$", "-2", colnames(h5_S1_HT8))
colnames(h5_S3_unsort) <- sub("-1$", "-3", colnames(h5_S3_unsort))
colnames(h5_S3_sorted) <- sub("-1$", "-4", colnames(h5_S3_sorted))

## Create Seurat Object ----
seu_S1_HT7 <- CreateSeuratObject(h5_S1_HT7, assay = "RNA", project = "S1_HT7")
seu_S1_HT8 <- CreateSeuratObject(h5_S1_HT8, assay = "RNA", project = "S1_HT8")

seu_S3_unsort <- CreateSeuratObject(h5_S3_unsort, assay = "RNA", project = "S3_unsort")
seu_S3_sorted <- CreateSeuratObject(h5_S3_sorted, assay = "RNA", project = "S3_sorted")

## Make Seurat object list ----
seu_list <- list(
  seu_S1_HT7,
  seu_S1_HT8,
  seu_S3_unsort,
  seu_S3_sorted
)

## Pre-processing ----
for (i in 1:length(seu_list)) {
  # Set different cutoff values for sample 1 (S1)
  cutoff_vals <- if (i %in% c(1, 2)) {
    c(200, 6000, 10)
  } else {
    c(200, 10000, 10)
  }
  
  # Apply the Seurat preprocessing workflow
  seu_list[[i]] <- seurat_workflow1(
    seu_list[[i]], 
    out_path = file.path(out2, "QC"),
    cutoff = cutoff_vals, 
    lim_nFeature_RNA = c(NA_real_, NA_real_),
    lim_nCount_RNA = c(NA_real_, NA_real_),
    lim_percent_mt = c(NA_real_, NA_real_)
  )
}

## S1: Merge Seurat Object ----
seu_S1_merged <- merge(seu_list[[1]], y = seu_list[[2]])
head(seu_S1_merged)
tail(seu_S1_merged)

### Normalize ----
# seu_S1_merged <- NormalizeData(seu_S1_merged) %>% FindVariableFeatures(nfeatures = 4000) %>%  ScaleData()
options(future.globals.maxSize = 4 * 1024^3) 
seu_S1_merged <- SCTransform(seu_S1_merged, variable.features.n = 4000)

### PCA/UMAP ----
seu_S1_merged <- RunPCA(seu_S1_merged)
ElbowPlot(seu_S1_merged, ndims = 50)

seu_S1_merged <- RunUMAP(seu_S1_merged, dims = 1:30, seed.use = 417, return.model = T) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.3)

### Rename clusters ----
clusters_ver1 <- c("C4", "C1", "C7", "C2", "C3", "C8", "C5", "C6")
clusters_ver1_order <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")
seu_S1_merged <- RenameSeuratClusters(seu_S1_merged, clusters_ver1, clusters_ver1_order, new_column_name = "clusters_ver1")

DimPlot(seu_S1_merged, group.by = "orig.ident", pt.size = 1, shuffle = T) + UMAP_preset("UMAP: RNA (Jurkat + MDA-MB-231)", nrow_legend = 1)
DimPlot(seu_S1_merged, group.by = "clusters_ver1", pt.size = 1, shuffle = T) + UMAP_preset("UMAP: RNA (Jurkat + MDA-MB-231)")

### Add module score ----
markers_MDA_MB_231 <- read.csv("data/MDA_MB_231_high_geneset.csv") %>% pull(Symbol)
markers_JURKAT <- read.csv("data/JURKAT_high_geneset.csv") %>% pull(Symbol)

seu_S1_merged <- AddModuleScore2(seu_S1_merged, markers_MDA_MB_231, ScoreName = "MDA-MB-231", assay = "RNA")
seu_S1_merged <- AddModuleScore2(seu_S1_merged, markers_JURKAT, ScoreName = "Jurkat", assay = "RNA")

FeaturePlot3(seu_S1_merged, features = "MDA-MB-231")
FeaturePlot3(seu_S1_merged, features = "Jurkat")
FeaturePlot3(seu_S1_merged, features = "EGFP")

### << Checkpoint (Save RDS) >> ----
# Setup working environment
working_dir <- file.path(this.path::here(),"..")
setwd(working_dir)

source("scripts/setup.R")
source("scripts/functions.R")

set.seed(1234)

# Create out2 folder
exp_name <- tools::file_path_sans_ext(basename(this.path::this.path()))
out2 <- file.path(out, exp_name)
if (!dir.exists(out2)) dir.create(out2)

overwrite_checkpoint <- F
if (overwrite_checkpoint) saveRDS(seu_S1_merged, file.path(rds, "20250604_seu_S1_merged.rds"))
if (overwrite_checkpoint) saveRDS(seu_list, file.path(rds, "20250604_seu_list.rds"))

if (!exists("seu_list")) seu_list <- readRDS(file.path(rds, "20250604_seu_list.rds"))
if (!exists("seu_S1_merged")) seu_S1_merged <- readRDS(file.path(rds, "20250604_seu_S1_merged.rds"))
if (!exists("markers_MDA_MB_231")) markers_MDA_MB_231 <- read.csv("data/MDA_MB_231_high_geneset.csv") %>% pull(Symbol)
if (!exists("markers_JURKAT")) markers_JURKAT <- read.csv("data/JURKAT_high_geneset.csv") %>% pull(Symbol)
if (!exists("geneset_list")) geneset_list <- readRDS("rds/geneset_list.rds")
if (!exists("geneset_1")) geneset_1 <- geneset_list[[1]]
if (!exists("geneset_2")) geneset_2 <- geneset_list[[2]]
if (!exists("geneset_3")) geneset_3 <- geneset_list[[3]]
if (!exists("geneset_4")) geneset_4 <- geneset_list[[4]]
if (!exists("geneset_5")) geneset_5 <- Reduce(intersect, list(geneset_1, geneset_2, geneset_3, geneset_4)) # intersect
if (!exists("geneset_6")) geneset_6 <- Reduce(union, list(geneset_1, geneset_2, geneset_3, geneset_4)) # union
if (!exists("geneset_ctl")) geneset_ctl <- c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1")

### Subset Jurkat cells ----
Jurkat_S1_merged <- subset(seu_S1_merged, idents = "C8", invert = T)
DimPlot(Jurkat_S1_merged, group.by = "orig.ident", pt.size = 1, shuffle = T) + UMAP_preset("UMAP: RNA (Jurkat)", nrow_legend = 1)
DimPlot(Jurkat_S1_merged, group.by = "clusters_ver1", pt.size = 1, shuffle = T) + UMAP_preset("UMAP: RNA (Jurkat)")

### Plots ----
if (!dir.exists(file.path(out2, "UMAP"))) dir.create(file.path(out2, "UMAP"), recursive = T)

svglite(file.path(out2, "UMAP", "S1_all_orig.svg"), width = 5, height = 4.15)
DimPlot(seu_S1_merged, group.by = "orig.ident", pt.size = 1, shuffle = T) + UMAP_preset("UMAP: S1 (Jurkat + MDA-MB-231)")
dev.off()

svglite(file.path(out2, "UMAP", "S1_all_clusters.svg"), width = 5, height = 4.15)
DimPlot(seu_S1_merged, group.by = "clusters_ver1", pt.size = 1, shuffle = T) + UMAP_preset("UMAP: S1 (Jurkat + MDA-MB-231)")
dev.off()

svglite(file.path(out2, "UMAP", "S1_Jurkat_orig.svg"), width = 4, height = 5.3)
DimPlot(Jurkat_S1_merged, group.by = "orig.ident", pt.size = 1, shuffle = T) + UMAP_preset("UMAP: S1 (Jurkat)")
dev.off()

svglite(file.path(out2, "UMAP", "S1_Jurkat_clusters.svg"), width = 4, height = 5.3)
DimPlot(Jurkat_S1_merged, group.by = "clusters_ver1", pt.size = 1, shuffle = T) + UMAP_preset("UMAP: S1 (Jurkat)")
dev.off()

# Build the ggplot object
p <- DimPlot(seu_S1_merged, group.by = "orig.ident", pt.size = 1, shuffle = T) + UMAP_preset("UMAP: S1 (Jurkat + MDA-MB-231)")
p <- ggplot_build(p)

# Extract x/y limits from panel ranges
x_limits <- p$layout$panel_params[[1]]$x$limits
x_breaks <- p$layout$panel_params[[1]]$x$breaks
y_limits <- p$layout$panel_params[[1]]$y$limits
y_breaks <- p$layout$panel_params[[1]]$y$breaks

FeaturePlot2(
  seu_S1_merged, 
  features = c(
    "MDA-MB-231", 
    "Jurkat"
  ),
  outfolder = file.path(out2, "FeaturePlot"),
  x_limits = x_limits,
  x_breaks = x_breaks,
  y_limits = y_limits,
  y_breaks = y_breaks,
  plot_width = 5,
  plot_height = 4
)

## S3 Unsort: reference mapping ----
### Normalize/PCA/UMAP ----
options(future.globals.maxSize = 4 * 1024^3) 
seu_list[[3]] <- SCTransform(seu_list[[3]], variable.features.n = 4000)

# seu_list[[3]] <- SCTransform(seu_list[[3]], variable.features.n = 4000) %>% RunPCA()
# ElbowPlot(seu_list[[3]], ndims = 50)
# 
# seu_list[[3]] <- RunUMAP(seu_list[[3]], dims = 1:30, seed.use = 417) %>%
#   FindNeighbors(dims = 1:30) %>%
#   FindClusters(resolution = 0.3)
# DimPlot(seu_list[[3]], pt.size = 1, shuffle = T) + UMAP_preset("UMAP: RNA (Jurkat + MDA-MB-231)")

### Finding anchors ----
seu_S3_unsort_anchors <- FindTransferAnchors(
  reference = seu_S1_merged,
  query = seu_list[[3]],
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:50
)

# predictions <- TransferData(anchorset = seu_S3_unsort_anchors, refdata = seu_S1_merged$clusters_ver1, dims = 1:50)
# query <- AddMetaData(seu_list[[3]], metadata = predictions)
# table(query$predicted.id)

### Unimodal UMAP Projection ----
seu_S3_unsort_mapped <- MapQuery(
  anchorset = seu_S3_unsort_anchors,
  query = seu_list[[3]],
  reference = seu_S1_merged,
  refdata = list(
    clusters_ver1 = "clusters_ver1"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

### Rename clusters ----
seu_S3_unsort_mapped$clusters_ver1 <- seu_S3_unsort_mapped$predicted.clusters_ver1
seu_S3_unsort_mapped <- SetIdent(seu_S3_unsort_mapped, value = "clusters_ver1")

Idents(seu_S3_unsort_mapped) <- factor(
  Idents(seu_S3_unsort_mapped),
  levels = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")
)
levels(seu_S3_unsort_mapped)

DimPlot(seu_S3_unsort_mapped, group.by = "clusters_ver1", pt.size = 1, shuffle = T) +
  UMAP_preset("Reference mapping: S3 Unsort", nrow_legend = 1)

### Add module score ----
seu_S3_unsort_mapped <- AddModuleScore2(seu_S3_unsort_mapped, markers_MDA_MB_231, ScoreName = "MDA-MB-231", assay = "RNA")
seu_S3_unsort_mapped <- AddModuleScore2(seu_S3_unsort_mapped, markers_JURKAT, ScoreName = "Jurkat", assay = "RNA")

FeaturePlot3(seu_S3_unsort_mapped, features = "MDA-MB-231")
FeaturePlot3(seu_S3_unsort_mapped, features = "Jurkat")
FeaturePlot3(seu_S3_unsort_mapped, features = "EGFP")

### Subset Jurkat cells ----
seu_S3_unsort_mapped[["UMAP_1"]] <- Embeddings(seu_S3_unsort_mapped, "ref.umap")[, 1]
seu_S3_unsort_mapped[["UMAP_2"]] <- Embeddings(seu_S3_unsort_mapped, "ref.umap")[, 2]
Jurkat_S3_unsort <- subset(seu_S3_unsort_mapped, cells = WhichCells(seu_S3_unsort_mapped, expression = UMAP_1 < 5))

DimPlot(seu_S3_unsort_mapped, group.by = "clusters_ver1", pt.size = 1, shuffle = T) + UMAP_preset("UMAP: RNA (Jurkat)")
DimPlot(Jurkat_S3_unsort, group.by = "clusters_ver1", pt.size = 1, shuffle = T) + UMAP_preset("UMAP: RNA (Jurkat)")

## S3 Sort: reference mapping ----
### Normalize/PCA/UMAP ----
options(future.globals.maxSize = 4 * 1024^3) 
seu_list[[4]] <- SCTransform(seu_list[[4]], variable.features.n = 4000)

### Finding anchors ----
seu_S3_sorted_anchors <- FindTransferAnchors(
  reference = seu_S1_merged,
  query = seu_list[[4]],
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:50
)

### Unimodal UMAP Projection ----
seu_S3_sorted_mapped <- MapQuery(
  anchorset = seu_S3_sorted_anchors,
  query = seu_list[[4]],
  reference = seu_S1_merged,
  refdata = list(
    clusters_ver1 = "clusters_ver1"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

### Rename clusters ----
seu_S3_sorted_mapped$clusters_ver1 <- seu_S3_sorted_mapped$predicted.clusters_ver1
seu_S3_sorted_mapped <- SetIdent(seu_S3_sorted_mapped, value = "clusters_ver1")

Idents(seu_S3_sorted_mapped) <- factor(
  Idents(seu_S3_sorted_mapped),
  levels = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")
)
levels(seu_S3_sorted_mapped)

DimPlot(seu_S3_sorted_mapped, group.by = "clusters_ver1", pt.size = 1, shuffle = T) +
  UMAP_preset("Reference mapping: S3 Unsort", nrow_legend = 1)

### Add module score ----
seu_S3_sorted_mapped <- AddModuleScore2(seu_S3_sorted_mapped, markers_MDA_MB_231, ScoreName = "MDA-MB-231", assay = "RNA")
seu_S3_sorted_mapped <- AddModuleScore2(seu_S3_sorted_mapped, markers_JURKAT, ScoreName = "Jurkat", assay = "RNA")

FeaturePlot3(seu_S3_sorted_mapped, features = "MDA-MB-231")
FeaturePlot3(seu_S3_sorted_mapped, features = "Jurkat")
FeaturePlot3(seu_S3_sorted_mapped, features = "EGFP")

### Subset Jurkat cells ----
seu_S3_sorted_mapped[["UMAP_1"]] <- Embeddings(seu_S3_sorted_mapped, "ref.umap")[, 1]
seu_S3_sorted_mapped[["UMAP_2"]] <- Embeddings(seu_S3_sorted_mapped, "ref.umap")[, 2]
Jurkat_S3_sorted <- subset(seu_S3_sorted_mapped, cells = WhichCells(seu_S3_sorted_mapped, expression = UMAP_1 < 5))
Jurkat_S3_sorted <- subset(Jurkat_S3_sorted, idents = "C8", invert = T)

DimPlot(seu_S3_sorted_mapped, group.by = "clusters_ver1", pt.size = 1, shuffle = T) + UMAP_preset("UMAP: RNA (Jurkat)")
DimPlot(Jurkat_S3_sorted, group.by = "clusters_ver1", pt.size = 1, shuffle = T) + UMAP_preset("UMAP: RNA (Jurkat)")

### Compare between Unsort and Sort ----
p1 <- DimPlot(Jurkat_S3_unsort, reduction = "ref.umap", group.by = "predicted.clusters_ver1", pt.size = 1, shuffle = T) +
  UMAP_preset("Reference mapping: S3 Unsort", nrow_legend = 1)

p2 <- DimPlot(Jurkat_S3_sorted, reduction = "ref.umap", group.by = "predicted.clusters_ver1", pt.size = 1, shuffle = T) +
  UMAP_preset("Reference mapping: S3 sort", nrow_legend = 1)

p1 + p2

### << Checkpoint (Save RDS) >> ----
# Setup working environment
working_dir <- file.path(this.path::here(),"..")
setwd(working_dir)

source("scripts/setup.R")
source("scripts/functions.R")

set.seed(1234)

# Create out2 folder
exp_name <- tools::file_path_sans_ext(basename(this.path::this.path()))
out2 <- file.path(out, exp_name)
if (!dir.exists(out2)) dir.create(out2)

overwrite_checkpoint <- F
if (overwrite_checkpoint) saveRDS(seu_S3_unsort_mapped, file.path(rds, "20250604_seu_S3_unsort_mapped.rds"))
if (overwrite_checkpoint) saveRDS(seu_S3_sorted_mapped, file.path(rds, "20250604_seu_S3_sorted_mapped.rds"))
if (overwrite_checkpoint) saveRDS(Jurkat_S1_merged, file.path(rds, "20250604_Jurkat_S1_merged.rds"))
if (overwrite_checkpoint) saveRDS(Jurkat_S3_unsort, file.path(rds, "20250604_Jurkat_S3_unsort.rds"))
if (overwrite_checkpoint) saveRDS(Jurkat_S3_sorted, file.path(rds, "20250604_Jurkat_S3_sorted.rds"))

if (!exists("Jurkat_S1_merged")) Jurkat_S1_merged <- readRDS(file.path(rds, "20250604_Jurkat_S1_merged.rds"))
if (!exists("Jurkat_S3_unsort")) Jurkat_S3_unsort <- readRDS(file.path(rds, "20250604_Jurkat_S3_unsort.rds"))
if (!exists("Jurkat_S3_sorted")) Jurkat_S3_sorted <- readRDS(file.path(rds, "20250604_Jurkat_S3_sorted.rds"))
if (!exists("geneset_list")) geneset_list <- readRDS("rds/geneset_list.rds")
if (!exists("geneset_1")) geneset_1 <- geneset_list[[1]]
if (!exists("geneset_2")) geneset_2 <- geneset_list[[2]]
if (!exists("geneset_3")) geneset_3 <- geneset_list[[3]]
if (!exists("geneset_4")) geneset_4 <- geneset_list[[4]]
if (!exists("geneset_5")) geneset_5 <- Reduce(intersect, list(geneset_1, geneset_2, geneset_3, geneset_4)) # intersect
if (!exists("geneset_6")) geneset_6 <- Reduce(union, list(geneset_1, geneset_2, geneset_3, geneset_4)) # union
if (!exists("geneset_ctl")) geneset_ctl <- c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1")

### Add module score ----
gene_sets <- list(geneset_ctl, geneset_1, geneset_2, geneset_3, geneset_4, geneset_5, geneset_6)
for (i in seq_along(gene_sets)) {
  if (i == 1) {
    Jurkat_S1_merged <- AddModuleScore2(Jurkat_S1_merged, gene_sets[[i]], ScoreName = "Gene set CTL", assay = "RNA")
    Jurkat_S3_unsort <- AddModuleScore2(Jurkat_S3_unsort, gene_sets[[i]], ScoreName = "Gene set CTL", assay = "RNA")
    Jurkat_S3_sorted <- AddModuleScore2(Jurkat_S3_sorted, gene_sets[[i]], ScoreName = "Gene set CTL", assay = "RNA")
  } else {
    Jurkat_S1_merged <- AddModuleScore2(Jurkat_S1_merged, gene_sets[[i]], ScoreName = paste0("Gene set ", (i - 1)), assay = "RNA")
    Jurkat_S3_unsort <- AddModuleScore2(Jurkat_S3_unsort, gene_sets[[i]], ScoreName = paste0("Gene set ", (i - 1)), assay = "RNA")
    Jurkat_S3_sorted <- AddModuleScore2(Jurkat_S3_sorted, gene_sets[[i]], ScoreName = paste0("Gene set ", (i - 1)), assay = "RNA")
  }
}

### Analysis ----
#### DEG + PseudoBulk ----
S1_markers <- FindAllMarkers(Jurkat_S1_merged, logfc.threshold = 0)
S1_pseud_expr_list <- AverageExpression(Jurkat_S1_merged)
S1_pseud_expr <- as.data.frame(S1_pseud_expr_list[["SCT"]])

write.csv(S1_markers, file.path(out2, "S1_markers.csv"), row.names = F)
write.csv(S1_pseud_expr, file.path(out2, "S1_pseud_bulk.csv"))

#### ORA ----
if (!dir.exists(file.path(out2, "ORA"))) dir.create(file.path(out2, "ORA"))
for (i in unique(S1_markers$cluster)) {
  WebGestaltR(
    enrichMethod = "ORA",
    organism = "hsapiens",
    enrichDatabase = "pathway_KEGG",
    interestGene = S1_markers %>% filter(cluster == i, p_val_adj < 0.05, avg_log2FC > 0.1) %>% pull(.),
    interestGeneType = "genesymbol",
    referenceSet = "genome_protein-coding",
    minNum = 5,
    maxNum = 2000,
    sigMethod = "top",
    fdrMethod = "BH",
    useWeightedSetCover = T,
    usekMedoid = F,
    fdrThr = 0.05,
    topThr = 40,
    setCoverNum = 10,
    kMedoid_k = 10,
    reportNum = 40,
    isOutput = T,
    outputDirectory = file.path(out2, "ORA"),
    projectName = i
  )
}

#### GSEA ----
if (!dir.exists(file.path(out2, "GSEA"))) dir.create(file.path(out2, "GSEA"))
for (i in unique(S1_markers$cluster)) {
  WebGestaltR(
    enrichMethod = "GSEA",
    organism = "hsapiens",
    enrichDatabase = "pathway_KEGG",
    interestGene = S1_markers %>% 
      filter(cluster == i) %>% 
      select(gene, avg_log2FC),
    interestGeneType = "genesymbol",
    referenceSet = "genome_protein-coding",
    minNum = 5,
    maxNum = 2000,
    sigMethod = "top",
    fdrMethod = "BH",
    useWeightedSetCover = T,
    usekMedoid = F,
    fdrThr = 0.05,
    topThr = 40,
    setCoverNum = 10,
    kMedoid_k = 10,
    reportNum = 40,
    isOutput = T,
    outputDirectory = file.path(out2, "GSEA"),
    projectName = i
  )
}

# Add GSEA color bar
svg_files <- list.files(
  path = "out/seurat_exp_20250603_reanalysis/GSEA/Project_C7/Project_C7_GSEA/",
  pattern = "\\.svg$",
  full.names = T,
  recursive = F
)
svg_files <- svg_files[!grepl("color\\.svg$", svg_files)] # Avoid reprocessing previously completed tasks

for (svg in svg_files) {
  add_GSEA_heatmap(
    svg_path = svg,
    score_path = "out/seurat_exp_20250603_reanalysis/GSEA/Project_C7/interestingID_mappingTable_C7.txt",
    output_folder = "out/seurat_exp_20250603_reanalysis/GSEA/Project_C7/Project_C7_GSEA/",
  )
}

### Plots ----
# Build the ggplot object
p <- DimPlot(Jurkat_S1_merged, group.by = "orig.ident", pt.size = 1, shuffle = T) + UMAP_preset("UMAP: S1 (Jurkat)")
p <- ggplot_build(p)

# Extract x/y limits from panel ranges
x_limits <- p$layout$panel_params[[1]]$x$limits
x_breaks <- p$layout$panel_params[[1]]$x$breaks
y_limits <- p$layout$panel_params[[1]]$y$limits
y_breaks <- p$layout$panel_params[[1]]$y$breaks

if (!dir.exists(file.path(out2, "UMAP"))) dir.create(file.path(out2, "UMAP"), recursive = T)
svglite(file.path(out2, "UMAP", "S3_Jurkat_unsort_clusters.svg"), width = 4, height = 5.3)
DimPlot(Jurkat_S3_unsort, group.by = "clusters_ver1", pt.size = 1, shuffle = T) + UMAP_preset("UMAP: S3 Unsort (Jurkat)")
dev.off()

svglite(file.path(out2, "UMAP", "S3_Jurkat_sorted_clusters.svg"), width = 4, height = 5.3)
DimPlot(Jurkat_S3_sorted, group.by = "clusters_ver1", pt.size = 1, shuffle = T) + UMAP_preset("UMAP: S3 Sorted (Jurkat)")
dev.off()

svglite(file.path(out2, "UMAP", "S3_Jurkat_unsort_density.svg"), width = 4, height = 5.3)
plot_umap_density(Jurkat_S3_unsort, Embedding = "ref.umap") + UMAP_preset("UMAP density: S3 Unsort (Jurkat)") +
  theme(plot.margin = margin(t = 10.5, r = 10.5, b = 10.5, l = 10.5))
dev.off()

svglite(file.path(out2, "UMAP", "S3_Jurkat_sorted_density.svg"), width = 4, height = 5.3)
plot_umap_density(Jurkat_S3_sorted, Embedding = "ref.umap") + UMAP_preset("UMAP density: S3 Sorted (Jurkat)") +
  theme(plot.margin = margin(t = 10.5, r = 10.5, b = 10.5, l = 10.5))
dev.off()

FeaturePlot2(
  Jurkat_S1_merged,
  features = c(
    geneset_5,
    "NFATC1",
    "Gene set 1",
    "Gene set 2",
    "Gene set 3",
    "Gene set 4",
    "Gene set 5",
    "Gene set 6",
    "Gene set CTL"
  ),
  outfolder = file.path(out2, "FeaturePlot", "S1"),
  x_limits = x_limits,
  x_breaks = x_breaks,
  y_limits = y_limits,
  y_breaks = y_breaks,
  plot_width = 4,
  plot_height = 5.17
)

FeaturePlot2(
  Jurkat_S3_unsort,
  features = c(
    geneset_5,
    "Gene set 1",
    "Gene set 2",
    "Gene set 3",
    "Gene set 4",
    "Gene set 5",
    "Gene set 6",
    "Gene set CTL"
  ),
  outfolder = file.path(out2, "FeaturePlot", "S2"),
  x_limits = x_limits,
  x_breaks = x_breaks,
  y_limits = y_limits,
  y_breaks = y_breaks,
  plot_width = 4,
  plot_height = 5.17
)

FeaturePlot2(
  Jurkat_S3_sorted,
  features = c(
    geneset_5,
    "Gene set 1",
    "Gene set 2",
    "Gene set 3",
    "Gene set 4",
    "Gene set 5",
    "Gene set 6",
    "Gene set CTL"
  ),
  outfolder = file.path(out2, "FeaturePlot", "S3"),
  x_limits = x_limits,
  x_breaks = x_breaks,
  y_limits = y_limits,
  y_breaks = y_breaks,
  plot_width = 4,
  plot_height = 5.17
)

# scRepertoire-seq ----
## S1: Add Repertoire Annotation ----
Jurkat_S1_merged <- AddRepertoireAnnotation(
  seu_obj = Jurkat_S1_merged,
  contig_path = "data/20240618_scREPERTOIREseq/scFv_denovo/outs/all_contig.fasta", 
  contig_anno_path = "data/20240618_scREPERTOIREseq/scFv_denovo/outs/all_contig_annotations.csv",
  output_path = file.path(out2, "Repertoire"),
  intermediate_path = file.path(out2, "Repertoire/intermediate"),
  prefix = "S1"
)

scFv_colors <- scFv_color_list(Jurkat_S1_merged)
DimPlot(Jurkat_S1_merged, group.by = "scFv_type", pt.size = 1, order = names(scFv_colors), cols = scFv_colors) + 
  UMAP_preset("UMAP: RNA (Jurkat)")

### Assemble Training data ----
S1_scFv_dict <- read.csv(file.path(out2, "Repertoire", "S1_scFv_dict.csv"))

#### Set 1-6 (loop)
gene_sets <- list(geneset_1, geneset_2, geneset_3, geneset_4, geneset_5, geneset_6)
for (i in seq_along(gene_sets)) {
  geneset <- c(paste0("Gene set ", i), gene_sets[[i]])
  
  cat(crayon::blue(geneset[1], ":\n"))
  print(gene_sets[[i]])
  cat("\n\n")
  
  training_data <- AssembleTrainingData(
    Jurkat_S1_merged,
    S1_scFv_dict,
    meatdata = c("orig.ident", "clusters_ver1"), 
    genes = geneset
  )
  
  write.csv(training_data, file = file.path(out2, "Repertoire", paste0("S1_TrainingData_set_", i, ".csv")), row.names = F)
}

## S3 Unsort: Add Repertoire Annotation ----
Jurkat_S3_unsort <- AddRepertoireAnnotation(
  seu_obj = Jurkat_S3_unsort,
  contig_path = list(
    "data/20250416_scREPERTOIREseq/Unsort-Target_22HNT5LT4/outs/all_contig.fasta",
    "data/20250416_scREPERTOIREseq/Unsort-Target_22LWTWLT4/outs/all_contig.fasta",
    "data/20250416_scREPERTOIREseq/Unsort-Target_22LWVNLT4/outs/all_contig.fasta"
  ), 
  contig_anno_path = list(
    "data/20250416_scREPERTOIREseq/Unsort-Target_22HNT5LT4/outs/all_contig_annotations.csv",
    "data/20250416_scREPERTOIREseq/Unsort-Target_22LWTWLT4/outs/all_contig_annotations.csv",
    "data/20250416_scREPERTOIREseq/Unsort-Target_22LWVNLT4/outs/all_contig_annotations.csv"
  ),
  output_path = file.path(out2, "Repertoire"),
  intermediate_path = file.path(out2, "Repertoire/intermediate"),
  prefix = "S3_unsort"
)

scFv_colors <- scFv_color_list(Jurkat_S3_unsort)
DimPlot(Jurkat_S3_unsort, group.by = "scFv_type", pt.size = 1, order = names(scFv_colors), cols = scFv_colors) + 
  UMAP_preset("UMAP: RNA (Jurkat)")

### Assemble Training data ----
S3_unsort_scFv_dict <- read.csv(file.path(out2, "Repertoire", "S3_unsort_scFv_dict.csv"))

#### Set 1-6 (loop)
gene_sets <- list(geneset_1, geneset_2, geneset_3, geneset_4, geneset_5, geneset_6)
for (i in seq_along(gene_sets)) {
  geneset <- c(paste0("Gene set ", i), gene_sets[[i]])
  
  cat(crayon::blue(geneset[1], ":\n"))
  print(gene_sets[[i]])
  cat("\n\n")
  
  training_data <- AssembleTrainingData(
    Jurkat_S3_unsort,
    S3_unsort_scFv_dict,
    meatdata = c("orig.ident", "clusters_ver1"), 
    genes = geneset
  )
  
  write.csv(training_data, file = file.path(out2, "Repertoire", paste0("S3_unsort_TrainingData_set_", i, ".csv")), row.names = F)
}

## S3 Sorted: Add Repertoire Annotation ----
Jurkat_S3_sorted <- AddRepertoireAnnotation(
  seu_obj = Jurkat_S3_sorted,
  contig_path = list(
    "data/20250326_scREPERTOIREseq/Sort-Target_22HNT5LT4/outs/all_contig.fasta",
    "data/20250502_scREPERTOIREseq/Sort-Target_22L2CTLT4/outs/all_contig.fasta",
    "data/20250602_scREPERTOIREseq/Sort-Target_22LWKVLT4/outs/all_contig.fasta"
  ), 
  contig_anno_path = list(
    "data/20250326_scREPERTOIREseq/Sort-Target_22HNT5LT4/outs/all_contig_annotations.csv",
    "data/20250502_scREPERTOIREseq/Sort-Target_22L2CTLT4/outs/all_contig_annotations.csv",
    "data/20250602_scREPERTOIREseq/Sort-Target_22LWKVLT4/outs/all_contig_annotations.csv"
  ),
  output_path = file.path(out2, "Repertoire"),
  intermediate_path = file.path(out2, "Repertoire/intermediate"),
  prefix = "S3_sorted"
)

scFv_colors <- scFv_color_list(Jurkat_S3_sorted)
DimPlot(Jurkat_S3_sorted, group.by = "scFv_type", pt.size = 1, order = names(scFv_colors), cols = scFv_colors) + 
  UMAP_preset("UMAP: RNA (Jurkat)")

### Assemble Training data ----
S3_sorted_scFv_dict <- read.csv(file.path(out2, "Repertoire", "S3_sorted_scFv_dict.csv"))

#### Set 1-6 (loop)
gene_sets <- list(geneset_1, geneset_2, geneset_3, geneset_4, geneset_5, geneset_6)
for (i in seq_along(gene_sets)) {
  geneset <- c(paste0("Gene set ", i), gene_sets[[i]])
  
  cat(crayon::blue(geneset[1], ":\n"))
  print(gene_sets[[i]])
  cat("\n\n")
  
  training_data <- AssembleTrainingData(
    Jurkat_S3_sorted,
    S3_sorted_scFv_dict,
    meatdata = c("orig.ident", "clusters_ver1"), 
    genes = geneset
  )
  
  write.csv(training_data, file = file.path(out2, "Repertoire", paste0("S3_sorted_TrainingData_set_", i, ".csv")), row.names = F)
}

# << Checkpoint (Save RDS) >> ----
# Setup working environment
working_dir <- file.path(this.path::here(),"..")
setwd(working_dir)

source("scripts/setup.R")
source("scripts/functions.R")

set.seed(1234)

# Create out2 folder
exp_name <- tools::file_path_sans_ext(basename(this.path::this.path()))
out2 <- file.path(out, exp_name)
if (!dir.exists(out2)) dir.create(out2)

overwrite_checkpoint <- F
if (overwrite_checkpoint) saveRDS(Jurkat_S1_merged, file.path(rds, "20250604_Jurkat_S1_merged.rds"))
if (overwrite_checkpoint) saveRDS(Jurkat_S3_unsort, file.path(rds, "20250604_Jurkat_S3_unsort.rds"))
if (overwrite_checkpoint) saveRDS(Jurkat_S3_sorted, file.path(rds, "20250604_Jurkat_S3_sorted.rds"))

if (!exists("Jurkat_S1_merged")) Jurkat_S1_merged <- readRDS(file.path(rds, "20250604_Jurkat_S1_merged.rds"))
if (!exists("Jurkat_S3_unsort")) Jurkat_S3_unsort <- readRDS(file.path(rds, "20250604_Jurkat_S3_unsort.rds"))
if (!exists("Jurkat_S3_sorted")) Jurkat_S3_sorted <- readRDS(file.path(rds, "20250604_Jurkat_S3_sorted.rds"))
if (!exists("geneset_list")) geneset_list <- readRDS("rds/geneset_list.rds")
if (!exists("geneset_1")) geneset_1 <- geneset_list[[1]]
if (!exists("geneset_2")) geneset_2 <- geneset_list[[2]]
if (!exists("geneset_3")) geneset_3 <- geneset_list[[3]]
if (!exists("geneset_4")) geneset_4 <- geneset_list[[4]]
if (!exists("geneset_5")) geneset_5 <- Reduce(intersect, list(geneset_1, geneset_2, geneset_3, geneset_4)) # intersect
if (!exists("geneset_6")) geneset_6 <- Reduce(union, list(geneset_1, geneset_2, geneset_3, geneset_4)) # union
if (!exists("geneset_ctl")) geneset_ctl <- c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1")

# Plots ----
if (!dir.exists(file.path(out2, "UMAP"))) dir.create(file.path(out2, "UMAP"), recursive = T)

svglite(file.path(out2, "UMAP", "S1_Jurkat_scFv.svg"), width = 4, height = 5.3)
scFv_colors <- scFv_color_list(Jurkat_S1_merged)
DimPlot(Jurkat_S1_merged, group.by = "scFv_type", pt.size = 1, order = names(scFv_colors), cols = scFv_colors) + 
  UMAP_preset("UMAP: S1 (Jurkat)")
dev.off()

# Clonal analysis ----
file_list <- list(
  S1 = "S1_scFv_dict.csv",
  S3_unsort = "S3_unsort_scFv_dict.csv",
  S3_sorted = "S3_sorted_scFv_dict.csv"
)

scFv_dict_list <- lapply(names(file_list), function(name) {
  df <- read.csv(file.path(out2, "Repertoire", file_list[[name]]))
  return(df)
})

scFv_dict <- bind_rows(scFv_dict_list) %>% mutate(clone_ID = as.integer(factor(Full_V_domain)))
write.csv(scFv_dict, file.path(out2, "Repertoire/scFv_dict_all.csv"), row.names = F)

# Add `clone_ID` back to seurat object
meta_clone_id <- scFv_dict %>% select(seu_barcode, clone_ID)
meta_clone_id_unique <- meta_clone_id %>%
  group_by(seu_barcode) %>%
  filter(n_distinct(clone_ID) == 1) %>%
  slice(1) %>% 
  ungroup()

Jurkat_S1_merged <- AppendSeuratMetadata(Jurkat_S1_merged, meta_clone_id_unique)
Jurkat_S3_unsort <- AppendSeuratMetadata(Jurkat_S3_unsort, meta_clone_id_unique)
Jurkat_S3_sorted <- AppendSeuratMetadata(Jurkat_S3_sorted, meta_clone_id_unique)

# Extract Seurat data
names(Jurkat_S1_merged[[]])
ext_meta <- c("clusters_ver1", "scFv_type", "clone_ID", "Gene set 2", "orig.ident")
df_Jurkat_S1_merged <- ExtractSeuratDataFrame(Jurkat_S1_merged, reduction = "umap", metadata = ext_meta)
df_Jurkat_S3_unsort <- ExtractSeuratDataFrame(Jurkat_S3_unsort, reduction = "ref.umap", metadata = ext_meta) %>%
  rename(umap_1 = refUMAP_1, umap_2 = refUMAP_2)
df_Jurkat_S3_sorted <- ExtractSeuratDataFrame(Jurkat_S3_sorted, reduction = "ref.umap", metadata = ext_meta) %>%
  rename(umap_1 = refUMAP_1, umap_2 = refUMAP_2)

# Build the ggplot object
p <- DimPlot(Jurkat_S1_merged, group.by = "orig.ident", pt.size = 1, shuffle = T) + UMAP_preset("UMAP: S1 (Jurkat + MDA-MB-231)")
p <- ggplot_build(p)

# Extract x/y limits from panel ranges
x_limits <- p$layout$panel_params[[1]]$x$limits
x_breaks <- p$layout$panel_params[[1]]$x$breaks
y_limits <- p$layout$panel_params[[1]]$y$limits
y_breaks <- p$layout$panel_params[[1]]$y$breaks

# Get max clone number
max_clone_num <- max(c(df_Jurkat_S1_merged$clone_ID,
                       df_Jurkat_S3_unsort$clone_ID,
                       df_Jurkat_S3_sorted$clone_ID), na.rm = TRUE)
max_score <- max(c(df_Jurkat_S1_merged$`Gene set 2`,
                   df_Jurkat_S3_unsort$`Gene set 2`,
                   df_Jurkat_S3_sorted$`Gene set 2`), na.rm = TRUE)
min_score <- min(c(df_Jurkat_S1_merged$`Gene set 2`,
                   df_Jurkat_S3_unsort$`Gene set 2`,
                   df_Jurkat_S3_sorted$`Gene set 2`), na.rm = TRUE)

# Clone S1 + S3
df_Jurkat <- rbind(df_Jurkat_S1_merged, df_Jurkat_S3_unsort, df_Jurkat_S3_sorted)
df_Jurkat_filtered <- df_Jurkat %>% filter(!is.na(clone_ID))
stats_orig <- table(df_Jurkat_filtered$orig.ident)
stats_all <- df_Jurkat_filtered %>%
  summarise(
    total_clone = n(),
    unique_clone = n_distinct(clone_ID)
  )
stats_act <- df_Jurkat_filtered %>%
  filter(`Gene set 2` > 0) %>%
  summarise(
    activated_clone = n(),
    unique_activated_clone = n_distinct(clone_ID)
  )
stats <- bind_cols(stats_all, stats_act)

cat(glue::glue(
  "Clone S1 + S3:\n",
  "{paste(glue::glue('{names(stats_orig)}: {as.integer(stats_orig)}'), collapse = '\n')}\n",
  "total_clone: {stats$total_clone}\n",
  "unique_clone: {stats$unique_clone}\n",
  "activated_clone: {stats$activated_clone}\n",
  "unique_activated_clone: {stats$unique_activated_clone}\n\n"
))

df_Jurkat_filtered <- df_Jurkat_filtered %>%
  mutate(scFv_type2 = case_when(scFv_type == "STM004_HL" ~ "aPDL1", 
                               scFv_type == "Atezolizumab_LH" ~ "aPDL1",
                               TRUE ~ "Novel scFvs"))

p <- ggplot(df_Jurkat_filtered, aes(x = umap_1, y = umap_2, color = scFv_type2)) +
  geom_point(size = 1) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) + 
  coord_cartesian(xlim = x_limits, ylim = y_limits) +
  labs(title = "UMAP colored by clone_ID",
       color = "scFv_type") +
  theme_prism(base_size = 10) +
  theme(plot.title = element_text(hjust = 0),
        legend.title = element_text(),
        legend.position = "bottom")
svglite(file.path(out2, "UMAP", "scFv_type_S1+S3.svg"), width = 4, height = 5.3)
print(p)
dev.off()

p <- ggplot(df_Jurkat_filtered, aes(x = umap_1, y = umap_2, color = clone_ID)) +
  geom_point(size = 1) +
  viridis::scale_color_viridis(option = "turbo", limits = c(1, max_clone_num), discrete = F) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) + 
  coord_cartesian(xlim = x_limits, ylim = y_limits) +
  labs(title = "UMAP colored by clone_ID",
       color = "clone_ID") +
  theme_prism(base_size = 10) +
  theme(plot.title = element_text(hjust = 0),
        legend.title = element_text(),
        legend.position = "bottom")
svglite(file.path(out2, "UMAP", "Clone_S1+S3.svg"), width = 4, height = 5.3)
print(p)
dev.off()

p <- ggplot(df_Jurkat_filtered, aes(x = umap_1, y = umap_2, color = `Gene set 2`)) +
  geom_point(size = 1) +
  viridis::scale_color_viridis(option = "turbo", limits = c(min_score, max_score)) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) + 
  coord_cartesian(xlim = x_limits, ylim = y_limits) +
  labs(title = "UMAP colored by gene set 2",
       color = "gene set 2") +
  theme_prism(base_size = 10) +
  theme(plot.title = element_text(hjust = 0),
        legend.title = element_text(),
        legend.position = "bottom")
svglite(file.path(out2, "UMAP", "Clone_S1+S3_score.svg"), width = 4, height = 5.3)
print(p)
dev.off()

# Clone S1 only
df_Jurkat <- rbind(df_Jurkat_S1_merged)
df_Jurkat_filtered <- na.omit(df_Jurkat)
stats_orig <- table(df_Jurkat_filtered$orig.ident)
stats_all <- df_Jurkat_filtered %>%
  summarise(
    total_clone = n(),
    unique_clone = n_distinct(clone_ID)
  )
stats_act <- df_Jurkat_filtered %>%
  filter(`Gene set 2` > 0) %>%
  summarise(
    activated_clone = n(),
    unique_activated_clone = n_distinct(clone_ID)
  )
stats <- bind_cols(stats_all, stats_act)

cat(glue::glue(
  "Clone S1:\n",
  "{paste(glue::glue('{names(stats_orig)}: {as.integer(stats_orig)}'), collapse = '\n')}\n",
  "total_clone: {stats$total_clone}\n",
  "unique_clone: {stats$unique_clone}\n",
  "activated_clone: {stats$activated_clone}\n",
  "unique_activated_clone: {stats$unique_activated_clone}\n\n"
))

p <- ggplot(df_Jurkat_filtered, aes(x = umap_1, y = umap_2, color = clone_ID)) +
  geom_point(size = 1) +
  viridis::scale_color_viridis(option = "turbo", limits = c(1, max_clone_num), discrete = F) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) + 
  coord_cartesian(xlim = x_limits, ylim = y_limits) +
  labs(title = "UMAP colored by clone_ID",
       color = "clone_ID") +
  theme_prism(base_size = 10) +
  theme(plot.title = element_text(hjust = 0),
        legend.title = element_text(),
        legend.position = "bottom")
svglite(file.path(out2, "UMAP", "Clone_S1.svg"), width = 4, height = 5.3)
print(p)
dev.off()

p <- ggplot(df_Jurkat_filtered, aes(x = umap_1, y = umap_2, color = `Gene set 2`)) +
  geom_point(size = 1) +
  viridis::scale_color_viridis(option = "turbo", limits = c(min_score, max_score)) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) + 
  coord_cartesian(xlim = x_limits, ylim = y_limits) +
  labs(title = "UMAP colored by gene set 2",
       color = "gene set 2") +
  theme_prism(base_size = 10) +
  theme(plot.title = element_text(hjust = 0),
        legend.title = element_text(),
        legend.position = "bottom")
svglite(file.path(out2, "UMAP", "Clone_S1_score.svg"), width = 4, height = 5.3)
print(p)
dev.off()

# Clone S3 only
df_Jurkat <- rbind(df_Jurkat_S3_unsort, df_Jurkat_S3_sorted)
df_Jurkat_filtered <- na.omit(df_Jurkat)
stats_orig <- table(df_Jurkat_filtered$orig.ident)
stats_all <- df_Jurkat_filtered %>%
  summarise(
    total_clone = n(),
    unique_clone = n_distinct(clone_ID)
  )
stats_act <- df_Jurkat_filtered %>%
  filter(`Gene set 2` > 0) %>%
  summarise(
    activated_clone = n(),
    unique_activated_clone = n_distinct(clone_ID)
  )
stats <- bind_cols(stats_all, stats_act)

cat(glue::glue(
  "Clone S3:\n",
  "{paste(glue::glue('{names(stats_orig)}: {as.integer(stats_orig)}'), collapse = '\n')}\n",
  "total_clone: {stats$total_clone}\n",
  "unique_clone: {stats$unique_clone}\n",
  "activated_clone: {stats$activated_clone}\n",
  "unique_activated_clone: {stats$unique_activated_clone}\n\n"
))

p <- ggplot(df_Jurkat_filtered, aes(x = umap_1, y = umap_2, color = clone_ID)) +
  geom_point(size = 1) +
  viridis::scale_color_viridis(option = "turbo", limits = c(1, max_clone_num), discrete = F) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) + 
  coord_cartesian(xlim = x_limits, ylim = y_limits) +
  labs(title = "UMAP colored by clone_ID",
       color = "clone_ID") +
  theme_prism(base_size = 10) +
  theme(plot.title = element_text(hjust = 0),
        legend.title = element_text(),
        legend.position = "bottom")
svglite(file.path(out2, "UMAP", "Clone_S3.svg"), width = 4, height = 5.3)
print(p)
dev.off()

p <- ggplot(df_Jurkat_filtered, aes(x = umap_1, y = umap_2, color = `Gene set 2`)) +
  geom_point(size = 1) +
  viridis::scale_color_viridis(option = "turbo", limits = c(min_score, max_score)) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) + 
  coord_cartesian(xlim = x_limits, ylim = y_limits) +
  labs(title = "UMAP colored by gene set 2",
       color = "gene set 2") +
  theme_prism(base_size = 10) +
  theme(plot.title = element_text(hjust = 0),
        legend.title = element_text(),
        legend.position = "bottom")
svglite(file.path(out2, "UMAP", "Clone_S3_score.svg"), width = 4, height = 5.3)
print(p)
dev.off()

# For RNA Velocity ----
if (!dir.exists(file.path(out2, "metadata"))) dir.create(file.path(out2, "metadata"))

head(Jurkat_S1_merged)
tail(Jurkat_S1_merged)

## Extracting metadata
write.csv(data.frame("x" = Cells(Jurkat_S1_merged)), file = file.path(out2, "metadata", "cell_id.csv"))
write.csv(Embeddings(Jurkat_S1_merged, reduction = "umap"), file = file.path(out2, "metadata", "umap.csv"))
write.csv(Jurkat_S1_merged$clusters_ver1, file = file.path(out2, "metadata", "clusters.csv"))

# GSEA bubble plot ----
## Load GSEA result
clusters <- paste0("C", 1:7)
GSEA_all_clusters <- map_dfr(clusters, read_GSEA, out_folder = out2) %>%
  # Standardize/rename key columns for clarity
  rename(
    description = description, 
    NES = normalizedEnrichmentScore,
    FDR = FDR,
    leading_edge_size = leadingEdgeNum
  ) %>%
  mutate(
    # label for y-axis
    label = paste(cluster, description, sep = " | ")
  ) %>%
  filter(FDR <= 0.05)

## Check GSEA result by cluster
GSEA_all_clusters %>%
  filter(cluster == "C1") %>%
  select(description, NES, FDR, leading_edge_size) %>%
  arrange(desc(NES)) %>%
  print()

## User-selected enriched gene set descriptions for each cluster
user_selection <- list(
  C1 = c(
    "ATP-dependent chromatin remodeling",
    "Cytokine-cytokine receptor interaction",
    "Neutrophil extracellular trap formation",
    "Motor proteins"
  ),
  C2 = c(
    "Motor proteins",
    "Proteasome",
    "DNA replication",
    "Ribosome",
    "Mismatch repair"
  ),
  C3 = c(
    "Ribosome",
    "PI3K-Akt signaling pathway",
    "Chemokine signaling pathway",
    "Cytokine-cytokine receptor interaction",
    "Neutrophil extracellular trap formation",
    "ATP-dependent chromatin remodeling"
  ),
  C4 = c(
    "Proteasome",
    "DNA replication",
    "ATP-dependent chromatin remodeling",
    "Cytokine-cytokine receptor interaction",
    "Necroptosis"
  ),
  C5 = c(
    "Graft-versus-host disease",
    "Cell adhesion molecules",
    "Proteasome",
    "Cell cycle",
    "DNA replication"
  ),
  C6 = c(
    "Cell adhesion molecules",
    "Rap1 signaling pathway",
    "Wnt signaling pathway",
    "Serotonergic synapse",
    "Ribosome",
    "Oxidative phosphorylation",
    "Systemic lupus erythematosus"
  ),
  C7 = c(
    "NF-kappa B signaling pathway",
    "T cell receptor signaling pathway",
    "Cytokine-cytokine receptor interaction",
    "Natural killer cell mediated cytotoxicity",
    "JAK-STAT signaling pathway",
    "ATP-dependent chromatin remodeling",
    "Neutrophil extracellular trap formation"
  )
)
GSEA_selected <- filter_GSEA_by_selection(GSEA_all_clusters, user_selection)

# Factor order on Y-axis
# Order labels within each cluster by NES (or by FDR). Here we use NES.
GSEA_selected <- GSEA_selected %>%
  group_by(cluster) %>%
  arrange(NES, .by_group = TRUE) %>%
  mutate(label = forcats::fct_inorder(label)) %>%
  ungroup()

# Add Y-breaks
GSEA_spaced <- add_GSEA_spacers(GSEA_selected)
spacer_pos <- GSEA_spaced %>%
  mutate(y_num = as.numeric(label)) %>%
  filter(str_starts(as.character(label), "<<SPACER_")) %>%
  arrange(cluster) %>%
  pull(y_num)

p <- ggplot(GSEA_spaced, aes(x = NES, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_point(
    aes(size = leading_edge_size, fill = FDR),
    shape = 21
  ) +
  scale_fill_viridis_c(direction = -1) +
  geom_hline(yintercept = spacer_pos, color = "black", linewidth = 0.4) +
  scale_y_discrete(
    labels = function(lab) ifelse(grepl("^<<SPACER_", lab), "", lab),
    expand = expansion(add = c(1, 1))
  ) +
  scale_x_continuous(limits = c(-5, 5)) +
  theme_prism(base_size = 10, border = T) +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.3, linetype = "dashed"),
    axis.ticks.y = element_blank(),
    legend.title = element_text(),
    legend.ticks = element_line(colour = "black"),
    legend.frame = element_rect(colour = "black"),
    panel.border = element_rect(size = 1)
  ) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Cluster | Gene set",
    size = "Leading Edge\nSize",
    fill = "FDR"
  )
print(p)

svglite(file.path(out2, "GSEA", "GSEA_bubble_plot.svg"), width = 5.5, height = 7)
print(p)
dev.off()

# GEX bubble plot ----
# used_features <- read.csv("data/Gene_list_0822_GSEA.csv") %>% pull(Gene)
# 
# # Get wide format for clustering
# df <- DotPlot(Jurkat_S1_merged, features = used_features)$data
# df_wide <- df %>%
#   select(id, features.plot, avg.exp.scaled) %>%
#   pivot_wider(names_from = id, values_from = avg.exp.scaled)
# mat <- as.matrix(df_wide[,-1])
# rownames(mat) <- df_wide$features.plot
# 
# # Hierarchical clustering
# hc_rows <- hclust(dist(mat))
# row_order <- rownames(mat)[hc_rows$order]
# 
# # Dendrogram plot (rows)
# dendro_data_rows <- ggdendro::dendro_data(hc_rows)
# dendro_plot_rows <- ggplot(dendro_data_rows$segments) +
#   geom_segment(aes(x = y, y = x, xend = yend, yend = xend)) +
#   scale_y_continuous(
#     breaks = seq_along(row_order),
#     labels = row_order,
#     expand = c(0.02, 0.02)
#   ) +
#   scale_x_continuous(trans = "reverse") +
#   theme_void()
# 
# df$features.plot <- factor(df$features.plot, levels = row_order)
# 
# hc_columns <- hclust(dist(t(mat)))
# col_order <- colnames(mat)[hc_columns$order]
# 
# # Dendrogram plot (columns)
# dendro_data_columns <- ggdendro::dendro_data(hc_columns)
# dendro_plot_columns <- ggplot(dendro_data_columns$segments) +
#   geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
#   scale_x_continuous(
#     breaks = seq_along(col_order),
#     labels = col_order,
#     expand = c(0.1, 0.1)
#   ) + 
#   theme_void()
#   
# df$id <- factor(df$id, levels = col_order)
# 
# # Plot
# dotplot <- ggplot(df, aes(x = id, y = features.plot)) +
#   geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21, color = "black", stroke = 0.4) +
#   scale_fill_viridis_c() +
#   scale_y_discrete(expand = c(0.02, 0.02), position = "right") +
#   scale_x_discrete(expand = c(0.1, 0.1)) +
#   theme_prism(base_size = 10, border = T) +
#   theme(
#     panel.grid.major = element_line(color = "grey80", size = 0.3, linetype = "dashed"),
#     axis.text.x = element_text(angle = 90, hjust = 1),
#     axis.ticks.y.right = element_blank(),
#     axis.title = element_blank(),
#     axis.ticks.length = unit(2, "pt"),
#     legend.title = element_text(),
#     panel.border = element_rect(size = 1), 
#     plot.margin = margin(c(0, 0, 0, 0))
# )
# dendro_plot_columns + dendro_plot_rows + dotplot + plot_layout(design = "\n#A\nBC\n", widths = c(1, 5), heights = c(1, 20))

# Use function
if (!dir.exists(file.path(out2, "GEX"))) dir.create(file.path(out2, "GEX"))
used_features <- read.csv("data/Gene_list_0822_GSEA.csv") %>% pull(1) %>% c("EGFP") %>% unique()

p <- BubblePlot(Jurkat_S1_merged, used_features)
svglite(file.path(out2, "GEX", "Gene_list_0822_GSEA_bubble.svg"), width = 4.5, height = 9)
print(p)
dev.off()

used_features <- read.csv("data/Gene_list_0825_GSEA.csv") %>% pull(1) %>% c("EGFP") %>% unique()
p <- BubblePlot(Jurkat_S1_merged, used_features, col_hclust = F)
svglite(file.path(out2, "GEX", "Gene_list_0825_GSEA_bubble.svg"), width = 4.5, height = 9)
print(p)
dev.off()

p <- BubblePlot(Jurkat_S1_merged, used_features)
svglite(file.path(out2, "GEX", "Gene_list_0825_GSEA_bubble_2.svg"), width = 4.5, height = 9)
print(p)
dev.off()

# GEX scatter plot ----
exp_xy_candidate <- c("EGFP", "NR4A1")
exp_data <- FetchData(Jurkat_S1_merged, vars = exp_xy_candidate)
ggplot(exp_data, aes(x = .data[[exp_xy_candidate[1]]], y = .data[[exp_xy_candidate[2]]])) +
  geom_point(size = 1) + 
  UMAP_preset()