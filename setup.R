# Set up environment
# Set working directory
working_dir <- file.path(this.path::here(),"..")
setwd(working_dir)

# Load libraries
# GEX
library(Seurat)
options(Seurat.object.assay.version = "v3")
# GSEA/ORA
library(WebGestaltR)
# Tool
library(dplyr)
library(ggplot2)
library(ggprism)
library(ComplexHeatmap)
library(patchwork)
library(scales)
library(tidyr)
library(stringr)
library(purrr)
library(tibble)
library(svglite)
library(fs)

# Create folders
out <- file.path(getwd(), "out")
data <- file.path(getwd(), "data")
rds <- file.path(getwd(), "rds")
if (!dir.exists(out)) dir.create(out)
if (!dir.exists(data)) dir.create(data)
if (!dir.exists(rds)) dir.create(rds)

# Copy .h5 file from "02_cellrangerData"
path_cellranger <- "../02_cellrangerData"
# S1_HT7 (GEX 1st) ----
h5 <- "20240419_scRNAseq/demultiplexed_samples/outs/per_sample_outs/HT7/count/sample_filtered_feature_bc_matrix.h5"
if (!file.exists(file.path(data, h5)) & file.exists(file.path(path_cellranger, h5))) {
  dir_create(path_dir(file.path(data, h5)))
  file_copy(file.path(path_cellranger, h5), file.path(data, h5), overwrite = T)
}

# S1_HT8 (GEX 1st) ----
h5 <- "20240419_scRNAseq/demultiplexed_samples/outs/per_sample_outs/HT8/count/sample_filtered_feature_bc_matrix.h5"
if (!file.exists(file.path(data, h5)) & file.exists(file.path(path_cellranger, h5))) {
  dir_create(path_dir(file.path(data, h5)))
  file_copy(file.path(path_cellranger, h5), file.path(data, h5), overwrite = T)
}

# S1_scFv (Repertoire 1st) ----
fasta <- "/20240618_scREPERTOIREseq/scFv_denovo/outs/all_contig.fasta"
csv <- "/20240618_scREPERTOIREseq/scFv_denovo/outs/all_contig_annotations.csv"
if (!file.exists(file.path(data, fasta)) & file.exists(file.path(path_cellranger, fasta))) {
  dir_create(path_dir(file.path(data, fasta)))
  file_copy(file.path(path_cellranger, fasta), file.path(data, fasta), overwrite = T)
}
if (!file.exists(file.path(data, csv)) & file.exists(file.path(path_cellranger, csv))) {
  dir_create(path_dir(file.path(data, csv)))
  file_copy(file.path(path_cellranger, csv), file.path(data, csv), overwrite = T)
}

# S2_WT (GEX 2nd) ----
h5 <- "/20241018_scRNAseq/WT_22CYF7LT4/outs/filtered_feature_bc_matrix.h5"
if (!file.exists(file.path(data, h5)) & file.exists(file.path(path_cellranger, h5))) {
  dir_create(path_dir(file.path(data, h5)))
  file_copy(file.path(path_cellranger, h5), file.path(data, h5), overwrite = T)
}

# S2_KO (GEX 2nd) ----
h5 <- "/20241018_scRNAseq/KO_22CYF7LT4/outs/filtered_feature_bc_matrix.h5"
if (!file.exists(file.path(data, h5)) & file.exists(file.path(path_cellranger, h5))) {
  dir_create(path_dir(file.path(data, h5)))
  file_copy(file.path(path_cellranger, h5), file.path(data, h5), overwrite = T)
}

# S2_scFv_WT (Repertoire 2nd) ----
fasta <- "/20241111_scREPERTOIREseq/scFv_denovo_WT/outs/all_contig.fasta"
csv <- "/20241111_scREPERTOIREseq/scFv_denovo_WT/outs/all_contig_annotations.csv"
if (!file.exists(file.path(data, fasta)) & file.exists(file.path(path_cellranger, fasta))) {
  dir_create(path_dir(file.path(data, fasta)))
  file_copy(file.path(path_cellranger, fasta), file.path(data, fasta), overwrite = T)
}
if (!file.exists(file.path(data, csv)) & file.exists(file.path(path_cellranger, csv))) {
  dir_create(path_dir(file.path(data, csv)))
  file_copy(file.path(path_cellranger, csv), file.path(data, csv), overwrite = T)
}

# S2_scFv_KO (Repertoire 2nd) ----
fasta <- "/20241111_scREPERTOIREseq/scFv_denovo_KO/outs/all_contig.fasta"
csv <- "/20241111_scREPERTOIREseq/scFv_denovo_KO/outs/all_contig_annotations.csv"
if (!file.exists(file.path(data, fasta)) & file.exists(file.path(path_cellranger, fasta))) {
  dir_create(path_dir(file.path(data, fasta)))
  file_copy(file.path(path_cellranger, fasta), file.path(data, fasta), overwrite = T)
}
if (!file.exists(file.path(data, csv)) & file.exists(file.path(path_cellranger, csv))) {
  dir_create(path_dir(file.path(data, csv)))
  file_copy(file.path(path_cellranger, csv), file.path(data, csv), overwrite = T)
}

# Unsort_22HNT5LT4 (GEX 3rd) ----
h5 <- "/20250326_scRNAseq/Unsort_22HNT5LT4/outs/filtered_feature_bc_matrix.h5"
if (!file.exists(file.path(data, h5)) & file.exists(file.path(path_cellranger, h5))) {
  dir_create(path_dir(file.path(data, h5)))
  file_copy(file.path(path_cellranger, h5), file.path(data, h5), overwrite = T)
}

# Sort_22HNT5LT4 (GEX 3rd) ----
h5 <- "/20250326_scRNAseq/Sort_22HNT5LT4/outs/filtered_feature_bc_matrix.h5"
if (!file.exists(file.path(data, h5)) & file.exists(file.path(path_cellranger, h5))) {
  dir_create(path_dir(file.path(data, h5)))
  file_copy(file.path(path_cellranger, h5), file.path(data, h5), overwrite = T)
}

# # Unsort-Target_CONCATENATE (Repertoire 3rd) ----
# fasta <- "/20250416_scREPERTOIREseq/Unsort-Target_CONCATENATE/outs/all_contig.fasta"
# csv <- "/20250416_scREPERTOIREseq/Unsort-Target_CONCATENATE/outs/all_contig_annotations.csv"
# if (!file.exists(file.path(data, fasta)) & file.exists(file.path(path_cellranger, fasta))) {
#   dir_create(path_dir(file.path(data, fasta)))
#   file_copy(file.path(path_cellranger, fasta), file.path(data, fasta), overwrite = T)
# }
# if (!file.exists(file.path(data, csv)) & file.exists(file.path(path_cellranger, csv))) {
#   dir_create(path_dir(file.path(data, csv)))
#   file_copy(file.path(path_cellranger, csv), file.path(data, csv), overwrite = T)
# }

# # Sort-Target_CONCATENATE (Repertoire 3rd) ----
# fasta <- "/20250502_scREPERTOIREseq/Sort-Target_CONCATENATE/outs/all_contig.fasta"
# csv <- "/20250502_scREPERTOIREseq/Sort-Target_CONCATENATE/outs/all_contig_annotations.csv"
# if (!file.exists(file.path(data, fasta)) & file.exists(file.path(path_cellranger, fasta))) {
#   dir_create(path_dir(file.path(data, fasta)))
#   file_copy(file.path(path_cellranger, fasta), file.path(data, fasta), overwrite = T)
# }
# if (!file.exists(file.path(data, csv)) & file.exists(file.path(path_cellranger, csv))) {
#   dir_create(path_dir(file.path(data, csv)))
#   file_copy(file.path(path_cellranger, csv), file.path(data, csv), overwrite = T)
# }

# Unsort-Target_22HNT5LT4 (Repertoire 3rd) ----
fasta <- "/20250416_scREPERTOIREseq/Unsort-Target_22HNT5LT4/outs/all_contig.fasta"
csv <- "/20250416_scREPERTOIREseq/Unsort-Target_22HNT5LT4/outs/all_contig_annotations.csv"
if (!file.exists(file.path(data, fasta)) & file.exists(file.path(path_cellranger, fasta))) {
  dir_create(path_dir(file.path(data, fasta)))
  file_copy(file.path(path_cellranger, fasta), file.path(data, fasta), overwrite = T)
}
if (!file.exists(file.path(data, csv)) & file.exists(file.path(path_cellranger, csv))) {
  dir_create(path_dir(file.path(data, csv)))
  file_copy(file.path(path_cellranger, csv), file.path(data, csv), overwrite = T)
}

# Unsort-Target_22LWTWLT4 (Repertoire 3rd) ----
fasta <- "/20250416_scREPERTOIREseq/Unsort-Target_22LWTWLT4/outs/all_contig.fasta"
csv <- "/20250416_scREPERTOIREseq/Unsort-Target_22LWTWLT4/outs/all_contig_annotations.csv"
if (!file.exists(file.path(data, fasta)) & file.exists(file.path(path_cellranger, fasta))) {
  dir_create(path_dir(file.path(data, fasta)))
  file_copy(file.path(path_cellranger, fasta), file.path(data, fasta), overwrite = T)
}
if (!file.exists(file.path(data, csv)) & file.exists(file.path(path_cellranger, csv))) {
  dir_create(path_dir(file.path(data, csv)))
  file_copy(file.path(path_cellranger, csv), file.path(data, csv), overwrite = T)
}

# Unsort-Target_22LWVNLT4 (Repertoire 3rd) ----
fasta <- "/20250416_scREPERTOIREseq/Unsort-Target_22LWVNLT4/outs/all_contig.fasta"
csv <- "/20250416_scREPERTOIREseq/Unsort-Target_22LWVNLT4/outs/all_contig_annotations.csv"
if (!file.exists(file.path(data, fasta)) & file.exists(file.path(path_cellranger, fasta))) {
  dir_create(path_dir(file.path(data, fasta)))
  file_copy(file.path(path_cellranger, fasta), file.path(data, fasta), overwrite = T)
}
if (!file.exists(file.path(data, csv)) & file.exists(file.path(path_cellranger, csv))) {
  dir_create(path_dir(file.path(data, csv)))
  file_copy(file.path(path_cellranger, csv), file.path(data, csv), overwrite = T)
}

# Sort-Target_22HNT5LT4 (Repertoire 3rd) ----
fasta <- "/20250326_scREPERTOIREseq/Sort-Target_22HNT5LT4/outs/all_contig.fasta"
csv <- "/20250326_scREPERTOIREseq/Sort-Target_22HNT5LT4/outs/all_contig_annotations.csv"
if (!file.exists(file.path(data, fasta)) & file.exists(file.path(path_cellranger, fasta))) {
  dir_create(path_dir(file.path(data, fasta)))
  file_copy(file.path(path_cellranger, fasta), file.path(data, fasta), overwrite = T)
}
if (!file.exists(file.path(data, csv)) & file.exists(file.path(path_cellranger, csv))) {
  dir_create(path_dir(file.path(data, csv)))
  file_copy(file.path(path_cellranger, csv), file.path(data, csv), overwrite = T)
}

# Sort-Target_22L2CTLT4 (Repertoire 3rd) ----
fasta <- "/20250502_scREPERTOIREseq/Sort-Target_22L2CTLT4/outs/all_contig.fasta"
csv <- "/20250502_scREPERTOIREseq/Sort-Target_22L2CTLT4/outs/all_contig_annotations.csv"
if (!file.exists(file.path(data, fasta)) & file.exists(file.path(path_cellranger, fasta))) {
  dir_create(path_dir(file.path(data, fasta)))
  file_copy(file.path(path_cellranger, fasta), file.path(data, fasta), overwrite = T)
}
if (!file.exists(file.path(data, csv)) & file.exists(file.path(path_cellranger, csv))) {
  dir_create(path_dir(file.path(data, csv)))
  file_copy(file.path(path_cellranger, csv), file.path(data, csv), overwrite = T)
}

# Sort-Target_22LWKVLT4 (Repertoire 3rd) ----
fasta <- "/20250602_scREPERTOIREseq/Sort-Target_22LWKVLT4/outs/all_contig.fasta"
csv <- "/20250602_scREPERTOIREseq/Sort-Target_22LWKVLT4/outs/all_contig_annotations.csv"
if (!file.exists(file.path(data, fasta)) & file.exists(file.path(path_cellranger, fasta))) {
  dir_create(path_dir(file.path(data, fasta)))
  file_copy(file.path(path_cellranger, fasta), file.path(data, fasta), overwrite = T)
}
if (!file.exists(file.path(data, csv)) & file.exists(file.path(path_cellranger, csv))) {
  dir_create(path_dir(file.path(data, csv)))
  file_copy(file.path(path_cellranger, csv), file.path(data, csv), overwrite = T)
}

# Data map ----
map <- "../01_seqData/Datamap.csv"
if (!file.exists(file.path(data, "Datamap.csv")) & file.exists(file.path(map))) {
  file_copy(file.path(map), file.path(data, "Datamap.csv"), overwrite = T)
}

# Remove values
rm(h5, fasta, map, csv, path_cellranger)