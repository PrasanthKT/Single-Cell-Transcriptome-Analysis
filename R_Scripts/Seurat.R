# Load necessary libraries
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(patchwork)
library(readr)
library(clusterProfiler)
library(org.Mm.eg.db) 
library(enrichplot)
library(SingleR)
library(celldex)
library(enrichplot)
library(RColorBrewer)

h5_files <- list.files(pattern = "*.h5")

print(h5_files)

seurat_list <- lapply(h5_files, function(file) {
  sample_name <- strsplit(file, "_")[[1]][2]  
  data <- Read10X_h5(file)
  seurat_obj <- CreateSeuratObject(counts = data, project = sample_name)
  seurat_obj$sample <- sample_name
  return(seurat_obj)
})
I 
print(seurat_list[[1]])

seurat_combined <- merge(x = seurat_list[[1]], 
                         y = seurat_list[-1], 
                         add.cell.ids = names(seurat_list), 
                         project = "SingleCell")


print(seurat_combined)

print(Layers(seurat_combined))


seurat_combined[["percent.mt"]] <- PercentageFeatureSet(seurat_combined, pattern = "^MT-")

summary(seurat_combined[["percent.mt"]])

mito_genes <- grep("^MT-", rownames(seurat_combined), value = TRUE)

print(mito_genes)

seurat_combined <- subset(seurat_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)

print(seurat_combined)

seurat_combined <- NormalizeData(seurat_combined, normalization.method = "LogNormalize", scale.factor = 10000)

print(seurat_combined)

seurat_combined <- FindVariableFeatures(seurat_combined, selection.method = "vst", nfeatures = 2000)


print(length(VariableFeatures(seurat_combined)))

print(head(VariableFeatures(seurat_combined), 50))

seurat_combined <- ScaleData(seurat_combined, features = VariableFeatures(seurat_combined))

print(seurat_combined)

seurat_combined <- RunPCA(seurat_combined, features = VariableFeatures(seurat_combined))

print(seurat_combined[["pca"]], dims = 1:5, nfeatures = 10)


ElbowPlot(seurat_combined, ndims = 50)


seurat_combined <- FindNeighbors(seurat_combined, dims = 1:15)

seurat_combined <- FindClusters(seurat_combined, resolution = 0.5)  

seurat_combined <- RunUMAP(seurat_combined, dims = 1:15)

DimPlot(seurat_combined, reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle("UMAP Clustering of Cells")


seurat_combined <- JoinLayers(seurat_combined)

cluster_markers <- FindAllMarkers(seurat_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

print(head(cluster_markers))
#write.csv(cluster_markers, "cluster_markers.csv")


custom_colors <- brewer.pal(n = 12, name = "Paired") 
custom_colors <- colorRampPalette(custom_colors)(length(unique(seurat_combined$seurat_clusters)))


DimPlot(seurat_combined, reduction = "umap", label = TRUE, repel = TRUE, 
        cols = custom_colors) + 
  ggtitle("UMAP Clustering of Cells") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# Load necessary libraries
library(ggplot2)
library(Seurat)

# Improved Violin Plot with better aesthetics
VlnPlot(seurat_combined, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, 
        pt.size = 0.1) + 
  theme_minimal() +
  ggtitle("Quality Control Metrics") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  scale_fill_brewer(palette = "Pastel1") +  # Using a visually appealing color scheme
  scale_y_continuous(labels = scales::comma)  # Formatting y-axis for better readability


top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

# Generate improved DotPlot
DotPlot(seurat_combined, features = unique(top_markers$gene)) + 
  RotatedAxis() +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  theme_minimal() +
  ggtitle("Top Marker Gene Expression Across Clusters") +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

head(seurat_combined$seurat_clusters)

table(seurat_combined$seurat_clusters)


cluster_18_cells <- subset(seurat_combined, idents = "18")

print(cluster_18_cells)

DimPlot(seurat_combined, reduction = "umap", label = TRUE, repel = TRUE, 
        cells.highlight = WhichCells(seurat_combined, idents = "18")) + 
  ggtitle("Highlighting Cluster 18 in UMAP")

cluster_18_markers <- FindMarkers(seurat_combined, ident.1 = "18", min.pct = 0.25, logfc.threshold = 0.25)

head(cluster_18_markers)

head(cluster_markers)

#write.csv(cluster_markers, "cluster_markers.csv")

unique(cluster_markers$cluster)

DefaultAssay(seurat_combined) <- "RNA"

ref.data <- MouseRNAseqData()

DefaultAssay(seurat_combined) <- "RNA"

seurat_combined_avg <- GetAssayData(seurat_combined, layer = "counts")

cell_annotations <- SingleR(test = seurat_combined_avg, ref = ref.data, labels = ref.data$label.main)

seurat_combined$SingleR_label <- cell_annotations$labels

table(seurat_combined$SingleR_label)

DimPlot(seurat_combined, group.by = "SingleR_label", label = TRUE, repel = TRUE) +
  ggtitle("Cell Type Annotation with SingleR")

saveRDS(seurat_combined, file = "seurat_annotated.rds")

print(dim(seurat_combined_avg))

DefaultAssay(seurat_combined) <- "RNA"

top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)  


gene_list <- unique(top_markers$gene)
entrez_ids <- mapIds(org.Mm.eg.db, 
                     keys = gene_list, 
                     column = "ENTREZID", 
                     keytype = "SYMBOL", 
                     multiVals = "first")

entrez_ids <- na.omit(entrez_ids)

go_results_bp <- enrichGO(gene = entrez_ids,
                          OrgDb = org.Mm.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          readable = TRUE)


go_results_mf <- enrichGO(gene = entrez_ids,
                          OrgDb = org.Mm.eg.db,
                          keyType = "ENTREZID",
                          ont = "MF",
                          pAdjustMethod = "BH",
                          readable = TRUE)

go_results_cc <- enrichGO(gene = entrez_ids,
                          OrgDb = org.Mm.eg.db,
                          keyType = "ENTREZID",
                          ont = "CC",
                          pAdjustMethod = "BH",
                          readable = TRUE)

head(go_results_bp)
head(go_results_mf)
head(go_results_cc)

# Save results to CSV
#write.csv(as.data.frame(go_results_bp), "GO_BP_results.csv")
#write.csv(as.data.frame(go_results_mf), "GO_MF_results.csv")
#write.csv(as.data.frame(go_results_cc), "GO_CC_results.csv")


# Visualization

barplot(go_results_bp, showCategory = 10) + 
  ggtitle("Top 10 GO Biological Processes")

barplot(go_results_mf, showCategory = 10) + 
  ggtitle("Top 10 GO Molecular Functions")

barplot(go_results_cc, showCategory = 10) + 
  ggtitle("Top 10 GO Cellular Components")

#  KEGG Pathway Enrichment

DefaultAssay(seurat_combined) <- "RNA"

top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

gene_list <- unique(top_markers$gene)
entrez_ids <- mapIds(org.Mm.eg.db, 
                     keys = gene_list, 
                     column = "ENTREZID", 
                     keytype = "SYMBOL", 
                     multiVals = "first")

entrez_ids <- na.omit(entrez_ids)

kegg_results <- enrichKEGG(gene = entrez_ids, 
                           organism = "mmu",
                           keyType = "kegg",
                           pAdjustMethod = "BH")

head(kegg_results)

#write.csv(as.data.frame(kegg_results), "KEGG_pathway_results.csv")

# KEGG Visualization
dotplot(kegg_results, showCategory = 10) + 
  ggtitle("KEGG Pathway Enrichment Analysis")

cnetplot(kegg_results, showCategory = 5) + 
  ggtitle("KEGG Pathway Network")


