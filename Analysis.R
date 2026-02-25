#Load packages
library(dplyr)
library(Seurat)
library(ggplot2)
#Load dataset
SLE1.data <- Read10X(data.dir = "/Users/rachanakamath/Desktop/SLE_project/SLE1")
SLE2.data <- Read10X(data.dir = "/Users/rachanakamath/Desktop/SLE_project/SLE2")
CTRL.data <- Read10X(data.dir = "/Users/rachanakamath/Desktop/SLE_project/CRTL")

#Create Seurat objects
SLE1 <- CreateSeuratObject(SLE1.data, project = "SLE1", min.cells = 3, min.features = 200)
SLE2 <- CreateSeuratObject(SLE2.data, project = "SLE2", min.cells = 3, min.features = 200)
CTRL <- CreateSeuratObject(CTRL.data, project = "CTRL", min.cells = 3, min.features = 200)

#Add metadeta
SLE1$group <- "SLE"
SLE2$group <- "SLE"
CTRL$group <- "Control"

#merge metadata and samples
Data_final <- merge(CTRL, 
                    y = c(SLE1, SLE2),
                    add.cells.ids = c("CTRL", "SLE1", "SLE2"),
                    project = "SLE_PBMC") 
dim(Data_final)

#Calculate QC Metrics
#ad. mitochondrial percentage
Data_final[["percent.mt"]] <- PercentageFeatureSet(Data_final, pattern = "^MT-")

#Violin plot (before filtering)
VlnPlot(Data_final,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

#define threshold and apply filter
Data_filtered <- subset(Data_final,
                        subset = nFeature_RNA > 200 &
                          nFeature_RNA < 6000 &
                          percent.mt < 15)
dim(Data_filtered)

#Violin plot post filtering 
VlnPlot(Data_filtered,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

#Normalisation
Data_filtered <- NormalizeData(Data_filtered,
                               normalization.method = "LogNormalize",
                               scale.factor = 10000)

#finding highly variable features:
Data_filtered <- FindVariableFeatures(Data_filtered,
                                      selection.method = "vst",
                                      nfeatures = 2000)



VariableFeaturePlot(Data_filtered)

#pca 
all.genes <- rownames(Data_filtered)

Data_filtered <- ScaleData(Data_filtered,
                           features = all.genes)
Data_filtered <- RunPCA(Data_filtered,
                        features = VariableFeatures(object = Data_filtered))

ElbowPlot(Data_filtered) +
  ggtitle("PCA Diagnostics – Elbow Plot") +
  theme(plot.title = element_text(hjust = 0.5))


#clustering (using PC1-PC12)
Data_filtered <- FindNeighbors(Data_filtered, dims = 1:12)

Data_filtered <- FindClusters(Data_filtered, resolution = 0.5)

Data_filtered <- RunUMAP(Data_filtered, dims = 1:12)
#Plot UMAP
DimPlot(Data_filtered, label = TRUE) +
  ggtitle("UMAP Colored by Clusters")
table(Idents(Data_filtered))

#Marker Heatmap (cell annotation)
Data_filtered <- JoinLayers(Data_filtered)
markers <- FindAllMarkers(Data_filtered,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)
#select top markers

top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

#Selecting fewer cluster (top 8 largest)
cluster_sizes <- table(Idents(Data_filtered))
major_clusters <- names(sort(cluster_sizes, decreasing = TRUE))[1:8]

Data_major <- subset(Data_filtered, idents = major_clusters)

#Top 3 Marker/cluster
top3 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3)

Data_major <- subset(Data_major, downsample = 100)

#Marker heatmap
DoHeatmap(Data_major,
          features = unique(top3$gene),
          size = 3) +
  ggtitle("Top Marker Genes (Major Clusters)") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 8))

#Dotplot
DotPlot(Data_filtered,
        features = unique(top3$gene)) +
  RotatedAxis() +
  ggtitle("Marker Gene Dot Plot") +
  theme(plot.title = element_text(hjust = 0.5))

#Differential expression
Data_filtered$disease_group <- ifelse(
  Data_filtered$orig.ident == "CTRL",
  "CTRL",
  "SLE"
)
table(Data_filtered$disease_group)
#perfrom DE
Idents(Data_filtered) <- Data_filtered$disease_group
DE_results <- FindMarkers(Data_filtered,
                          ident.1 = "SLE",
                          ident.2 = "CTRL",
                          logfc.threshold = 0.25,
                          min.pct = 0.1)
#inspect top genes
head(DE_results[order(DE_results$p_val_adj), ])

#volcano plot
DE_results$gene <- rownames(DE_results)

ggplot(DE_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  ggtitle("Differential Expression (SLE vs CTRL)") +
  theme_minimal()

#save files
write.csv(DE_results,
          file = "DE_results_SLE_vs_CTRL_full.csv",
          row.names = TRUE)

top20_DE <- DE_results[order(DE_results$p_val_adj), ][1:20, ]

write.csv(top20_DE,
          file = "DE_results_SLE_vs_CTRL_top20.csv",
          row.names = TRUE)

DE_significant <- DE_results %>%
  dplyr::filter(p_val_adj < 0.05 &
                  abs(avg_log2FC) > 0.5)

write.csv(DE_significant,
          file = "DE_results_SLE_vs_CTRL_significant.csv",
          row.names = TRUE)

DE_up_SLE <- DE_results %>%
  dplyr::filter(p_val_adj < 0.05 &
                  avg_log2FC > 0.5)

write.csv(DE_up_SLE,
          file = "DE_results_upregulated_in_SLE.csv",
          row.names = TRUE)

top10_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(top10_markers,
          file = "Top10_Markers_Per_Cluster.csv",
          row.names = FALSE)
