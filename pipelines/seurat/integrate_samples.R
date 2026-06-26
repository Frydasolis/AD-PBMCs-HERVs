#!/usr/bin/env Rscript

# 
library(Seurat)
library(Matrix)
library(tidyverse)

# 1. 
STELLARSCOPE_DIR <- "/storage/lemus_g/roldan/AD/stellarscope/results"
STARSOLO_DIR     <- "/storage/lemus_g/roldan/AD/results"
OUTPUT_DIR       <- "/storage/lemus_g/roldan/AD/stellarscope/results/final_analysis"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 2. 
samples_file <- file.path(STARSOLO_DIR, "samples.txt")
if (!file.exists(samples_file)) {
  stop(paste("Error: No se encontró el archivo de muestras en:", samples_file))
}
samples_list <- readLines(samples_file)
cat("Founded", length(samples_list), "samples to proccess.\n")

# 
seurat_list <- list()

# 3.
for (sample_id in samples_list) {
  cat("\n=========================================\n")
  cat("Procesando Muestra:", sample_id, "\n")
  cat("=========================================\n")
  
  #---
  sample_stellar_dir <- file.path(STELLARSCOPE_DIR, sample_id)
  te_mtx <- file.path(sample_stellar_dir, paste0(sample_id, "-TE_counts.mtx"))
  te_feats <- file.path(sample_stellar_dir, paste0(sample_id, "-features.tsv"))
  te_barcodes <- file.path(sample_stellar_dir, paste0(sample_id, "-barcodes.tsv"))
  
  if (!file.exists(te_mtx)) {
    cat("Warning: There is not counting TEs for ", sample_id, ". Saltando...\n")
    next
  }
  
  # 
  te_counts <- ReadMtx(mtx = te_mtx, features = te_feats, cells = te_barcodes, feature.column = 1)
  
  #--- 
  correction_path <- file.path(STELLARSCOPE_DIR, "correction_matrices", paste0(sample_id, "_duplicates.txt"))
  
  if (file.exists(correction_path)) {
    cat("Applying the doble count correction matrix...\n")
    correction_df <- read.delim(correction_path, stringsAsFactors = FALSE)
    
    if (nrow(correction_df) > 0) {
      match_i <- match(correction_df$Gene_TE, rownames(te_counts))
      match_j <- match(correction_df$Cell, colnames(te_counts))
      
      valid_idx <- !is.na(match_i) & !is.na(match_j)
      
      correction_matrix <- sparseMatrix(
        i = match_i[valid_idx],
        j = match_j[valid_idx],
        x = correction_df$Count_To_Subtract[valid_idx],
        dims = dim(te_counts),
        dimnames = dimnames(te_counts)
      )
      
      te_counts_corrected <- te_counts - correction_matrix
      te_counts_corrected[te_counts_corrected < 0] <- 0
    } else {
      cat("Empty correction file.\n")
      te_counts_corrected <- te_counts
    }
  } else {
    cat("Warning: There is not correction file", sample_id, ". Se usan TEs crudos.\n")
    te_counts_corrected <- te_counts
  }
  
  #--- 
  starsolo_path <- file.path(STARSOLO_DIR, paste0(sample_id, "_Solo.out"), "Gene", "filtered")
  
  if (!dir.exists(starsolo_path)) {
    cat("Warning: No se encontró directorio STARsolo para", sample_id, ". Saltando...\n")
    next
  }
  
  # 
  starsolo_counts <- tryCatch({
    Read10X(data.dir = starsolo_path)
  }, error = function(e) {
    cat("Detectados archivos planos o con prefijos de STARsolo, cargando manualmente...\n")
    mtx_f  <- list.files(starsolo_path, pattern = "matrix.mtx", full.names = TRUE)
    feat_f <- list.files(starsolo_path, pattern = "(features|genes).tsv", full.names = TRUE)
    barc_f <- list.files(starsolo_path, pattern = "barcodes.tsv", full.names = TRUE)
    
    if(length(mtx_f) == 1 && length(feat_f) == 1 && length(barc_f) == 1) {
      return(ReadMtx(mtx = mtx_f, features = feat_f, cells = barc_f))
    } else {
      stop("Error fatal There is not  STARsolo matrix.")
    }
  })
  
  #--- 
  common_cells <- intersect(colnames(starsolo_counts), colnames(te_counts_corrected))
  cat("Células en común encontradas:", length(common_cells), "\n")
  
  if (length(common_cells) == 0) {
    cat("Warning: 0 células en común. Saltando muestra.\n")
    next
  }
  
  starsolo_counts <- starsolo_counts[, common_cells]
  te_counts_corrected <- te_counts_corrected[, common_cells]
  
  #--- 
  merged_counts <- rbind(starsolo_counts, te_counts_corrected)
  
  #---
  seurat_obj <- CreateSeuratObject(counts = merged_counts, project = sample_id, min.cells = 3, min.features = 200)
  
  # 
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # 
  seurat_list[[sample_id]] <- seurat_obj
}

#
cat("\nallsalmples...\n")
if (length(seurat_list) < 2) {
  stop("Error: There is not samples to join.")
}

combined_seurat <- merge(
  x = seurat_list[[1]], 
  y = seurat_list[-1], 
  add.cell.ids = names(seurat_list), 
  project = "PBMC_AD_Stellarscope"
)

# 5. 
cat("Quality control Pipeline Downstream...\n")
combined_seurat <- subset(combined_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

combined_seurat <- NormalizeData(combined_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
combined_seurat <- FindVariableFeatures(combined_seurat, selection.method = "vst", nfeatures = 2000)
combined_seurat <- ScaleData(combined_seurat)
combined_seurat <- RunPCA(combined_seurat, features = VariableFeatures(object = combined_seurat), verbose = FALSE)
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:20, verbose = FALSE)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5, verbose = FALSE)
combined_seurat <- RunUMAP(combined_seurat, dims = 1:20, verbose = FALSE)

# 6. 
output_file <- file.path(OUTPUT_DIR, "seurat_final_genes_and_TEs.rds")
cat("Guardando objeto final de Seurat en:", output_file, "\n")
saveRDS(combined_seurat, file = output_file)

cat("\n¡Complete!\n")
