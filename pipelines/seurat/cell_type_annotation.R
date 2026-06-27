library(Seurat)
library(ggplot2)
library(dplyr)

# 1. 

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# 2. 
table(seurat_obj$seurat_clusters)

library(Seurat)
library(ggplot2)
library(dplyr)

# 
Idents(seurat_obj) <- "seurat_clusters"

# 
cluster_ids <- c(
  "CD4 T Cells",          # Cluster 0  (CD3D+)
  "CD4 T Cells",          # Cluster 1  (CD3D+)
  "Other T Cells",        # Cluster 2  (CD3D+ mixto/gamma-delta)
  "CD4 T Cells",          # Cluster 3  (CD3D+)
  "CD4 T Cells",          # Cluster 4  (CD3D+)
  "CD14 Monocytes",       # Cluster 5  (CD14+)
  "B Cells",              # Cluster 6  (CD19+)
  "B Cells",              # Cluster 7  (CD19+)
  "CD8 T Cells",          # Cluster 8  (CD3D+ / CD8A+)
  "CD4 T Cells",          # Cluster 9  (CD3D+)
  "CD8 T Cells",          # Cluster 10 (CD3D+ / CD8A+)
  "CD4 T Cells",          # Cluster 11 (CD3D+)
  "Natural Killer",       # Cluster 12 (NCAM1+)
  "CD8 T Cells",          # Cluster 13 (CD3D+ / CD8A+)
  "CD16 Monocytes",       # Cluster 14 (FCGR3A+)
  "CD4 T Cells",          # Cluster 15 (CD3D+)
  "CD8 T Cells",          # Cluster 16 (CD3D+ / CD8A+)
  "Other T Cells",        # Cluster 17 (CD3D+ )
  "B Cells",              # Cluster 18 (CD19+)
  "B Cells",              # Cluster 19 (CD19+)
  "Dendritic Cells",      # Cluster 20 (LILRA4+ / CD1C+ )
  "Other T Cells",        # Cluster 21 (CD3D+ transicional)
  "Other T Cells",        # Cluster 22 (CD3D+ residual)
  "Other T Cells"         # Cluster 23 (CD3D+ residual)
)

# 
names(cluster_ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, cluster_ids)

# 
seurat_obj$paper_celltype <- Idents(seurat_obj)
# 
plot_final_umap <- DimPlot(
  seurat_obj, 
  reduction = "umap", 
  label = TRUE, 
  label.size = 4, 
  repel = TRUE,
  cols = c(
    "CD4 T Cells" = "#7FC97F",
    "CD8 T Cells" = "#BEAED4",
    "Other T Cells" = "#FDC086",
    "B Cells" = "#FFFF99",
    "CD14 Monocytes" = "#386CB0",
    "CD16 Monocytes" = "#F0027F",
    "Natural Killer" = "#BF5B17",
    "Dendritic Cells" = "#666666"
  )
) +
  theme_classic(base_size = 14) +
  labs(title = "UMAP of De Novo Annotated PBMCs", x = "UMAP 1", y = "UMAP 2", color = "Cell Type") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text = element_text(color = "black"))

# 
ggsave(
