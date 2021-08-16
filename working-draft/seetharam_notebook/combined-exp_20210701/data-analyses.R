setwd("~/TutejaLab/combined-exp_20210701")
# load the modules
library(Seurat)
library(SeuratWrappers)
library(knitr)
library(kableExtra)
library(ggplot2)
library(cowplot)
library(patchwork)
library(metap)
library(multtest)
library(gridExtra)
library(dplyr)
library(stringr)
library(TissueEnrich)
library(gprofiler2)
library(tidyverse)
library(enhancedDimPlot)
library(calibrate)
library(ggrepel)
library(dittoSeq)
library(ComplexHeatmap)
library(scales)
library(VennDiagram)

# create seurat object
experiment_name = "BAP"
dataset_loc <- "~/TutejaLab/expression-data"
ids <- c("5pcO2_r1", "5pcO2_r2", "20pcO2_r1", "20pcO2_r2", "nCT_D5", "nCT_D10", "nTE_D2", "nTE_D3")
# function d10x.data
d10x.data <- sapply(ids, function(i){
  d10x <- Read10X(file.path(dataset_loc,i,"filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"), '[[' , 1L ), i, sep="-")
  d10x
})

experiment.data <- do.call("cbind", d10x.data)
bapd8.combined <- CreateSeuratObject(
  experiment.data,
  project = "BAPd8",
  min.cells = 10,
  min.genes = 200,
  names.field = 2,
  names.delim = "\\-")

# add the replicate information in the metadata table
# backup
bapd8.temp <- bapd8.combined

# metadata stats
bapd8.combined$log10GenesPerUMI <- log10(bapd8.combined$nFeature_RNA) / log10(bapd8.combined$nCount_RNA)
bapd8.combined$mitoRatio <- PercentageFeatureSet(object = bapd8.combined, pattern = "^MT-")
bapd8.combined$mitoRatio <- bapd8.combined@meta.data$mitoRatio / 100
metadata <- bapd8.combined@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA, 
                seq_folder = orig.ident)
mt <- ggplot(dat = metadata, aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point(alpha = 0.5) + 
  scale_colour_gradient(low = "gray90", high = "black") + labs(colour="MT ratio") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black")) + 
  xlab("RNA counts") + ylab("Gene counts") +
  stat_smooth(method=lm) +
  facet_wrap(~seq_folder, labeller = labeller(seq_folder = 
                                                c("20pcO2_r1" = "20% Oxygen (rep1)",
                                                  "20pcO2_r2" = "20% Oxygen (rep2)",
                                                  "5pcO2_r1" = "5% Oxygen (rep1)",
                                                  "5pcO2_r2" = "5% Oxygen (rep2)",
                                                  "nCT_D5" = "nCT day 5",
                                                  "nCT_D10" = "nCT day 10",
                                                  "nTE_D2" = "nTE day 2",
                                                  "nTE_D3" = "nTE day 3"))) +
  scale_y_continuous(label=comma) +
  scale_x_continuous(label=comma) 

ggsave("Figure_S2.svg", plot = mt, dpi=900, width = 7, height = 8)
ggsave("Figure_S2.png", plot = mt, dpi=900, width = 7, height = 8)
# other qc plots
# number of cells per sample
ncells <- ggplot(metadata, aes(x=seq_folder, fill=seq_folder)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells")
ggsave("Figure_Sx_number-of-cells-per-sample.png", plot = ncells, dpi=900, width = 5, height = 6)
# density of cells per sample
dcells <- ggplot(metadata, aes(color=seq_folder, x=nUMI, fill= seq_folder)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
ggsave("Figure_Sx_density-of-cells-per-sample.png", plot = dcells, dpi=900, width = 5, height = 6)
# nubmer of cells vs. genes
ngenes <- ggplot(metadata, aes(x=seq_folder, y=log10(nGene), fill=seq_folder)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
ggsave("Figure_Sx_genes-vs-cells-per-sample.png", plot = dcells, dpi=900, width = 5, height = 6)
# transcritps
dtranscripts <- ggplot(metadata, aes(x=log10GenesPerUMI, color = seq_folder, fill=seq_folder)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
ggsave("Figure_Sx_transcripts-density-per-sample.png", plot = dtranscripts, dpi=900, width = 5, height = 6)
# mito ratio
mtratio <- ggplot(metadata, aes(color=seq_folder, x=mitoRatio, fill=seq_folder)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
ggsave("Figure_Sx_mt-ratio-per-sample.png", plot = mtratio, dpi=900, width = 5, height = 6)
# restore backup
bapd8.combined <- bapd8.temp


# restore backup, resume analyses
head(bapd8.combined@meta.data)
df <- bapd8.combined@meta.data
df$replicate <- NA
df$replicate[which(str_detect(df$orig.ident, "5pcO2"))] <- "5pcO2"
df$replicate[which(str_detect(df$orig.ident, "20pcO2"))] <- "20pcO2"
df$replicate[which(str_detect(df$orig.ident, "nCT_"))] <- "nCT"
df$replicate[which(str_detect(df$orig.ident, "nTE_"))] <- "nTE"
bapd8.combined@meta.data <- df
bapd8.combined[["percent.mt"]] <- PercentageFeatureSet(bapd8.combined, pattern = "^MT-")

# Figure S2: A, B and C
p <- VlnPlot(bapd8.combined, features = "nFeature_RNA", pt.size = 1) + 
  geom_hline(yintercept=200, color = "red", size=1) +
  geom_hline(yintercept=10000, color = "red", size=1) +
  theme(legend.position = "none")

q <- VlnPlot(bapd8.combined, features = "nCount_RNA", pt.size = 1) +
  theme(legend.position = "none")

r <- VlnPlot(bapd8.combined, features = "percent.mt", pt.size = 1) +
  geom_hline(yintercept=15, color = "red", size=1) +
  theme(legend.position = "none")
h =(p | q | r)
#ggsave("Figure_S2_C.svg", plot=  r, dpi=900, width = 3, height = 4)
#ggsave("Figure_S2_A.svg", plot = p, dpi=900, width = 3, height = 4)
#ggsave("Figure_S2_B.svg", plot = q, dpi=900, width = 3, height = 4)
#ggsave("Figure_S2_ABC.svg", plot = h, dpi=900, width = 12, height = 6)
ggsave("Figure_S2_C.png", plot=  r, dpi=900, width = 3, height = 4)
ggsave("Figure_S2_A.png", plot = p, dpi=900, width = 3, height = 4)
ggsave("Figure_S2_B.png", plot = q, dpi=900, width = 3, height = 4)
ggsave("Figure_S2_ABC.png", plot = h, dpi=900, width = 12, height = 6)


B1 <- FeatureScatter(bapd8.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
B2 <- FeatureScatter(bapd8.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
bapd8.combined <- subset(bapd8.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 25)
I1 <- FeatureScatter(bapd8.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
I2 <- FeatureScatter(bapd8.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
bapd8.combined <- subset(bapd8.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 15)
A1 <- FeatureScatter(bapd8.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
A2 <- FeatureScatter(bapd8.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

B <- B1 | B2
I <- I1 | I2
A <- A1 | A2
ggsave("Figure_before.png", plot = B, dpi=900, width = 12, height = 6)
ggsave("Figure_inter.png", plot = I, dpi=900, width = 12, height = 6)
ggsave("Figure_after.png", plot = A, dpi=900, width = 12, height = 6)



counts <- GetAssayData(object = bapd8.combined, slot = "counts")
counts <- counts[grep(pattern = "^MT", x = rownames(counts), invert = TRUE),]
counts <- counts[grep(pattern = "^MT", x = rownames(counts), invert = TRUE),]
counts <- counts[grep(pattern = "^RPL", x = rownames(counts), invert = TRUE),]
counts <- counts[grep(pattern = "^RPS", x = rownames(counts), invert = TRUE),]
counts <- counts[grep(pattern = "^MRPS", x = rownames(counts), invert = TRUE),]
counts <- counts[grep(pattern = "^MRPL", x = rownames(counts), invert = TRUE),]

keep_genes <- Matrix::rowSums(counts) >= 10
filtered_counts <- counts[keep_genes, ]
bapd8.fcombined <- CreateSeuratObject(filtered_counts, meta.data = bapd8.combined@meta.data)
bapd8.fcombined@meta.data <- bapd8.fcombined@meta.data[1:4]

bapd8.combined <- bapd8.fcombined


bapd8.list <- SplitObject(bapd8.combined, split.by = "orig.ident")
bapd8.list <- lapply(X = bapd8.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
bapd8.anchors <- FindIntegrationAnchors(object.list = bapd8.list, dims = 1:20)
bapd8.integrated <- IntegrateData(anchorset = bapd8.anchors, dims = 1:20)
DefaultAssay(bapd8.integrated) <- "integrated"
bapd8.integrated <- ScaleData(bapd8.integrated, verbose = FALSE)
bapd8.integrated <- RunPCA(bapd8.integrated, npcs = 30, verbose = FALSE)
bapd8.integrated <- RunUMAP(bapd8.integrated, reduction = "pca", dims = 1:20)
bapd8.integrated <- FindNeighbors(bapd8.integrated, reduction = "pca", dims = 1:20)
bapd8.integrated <- FindClusters(bapd8.integrated, resolution = 0.5)
num.clusters <- nlevels(bapd8.integrated$seurat_clusters)
num.clusters

df <- bapd8.integrated@meta.data
df$new_clusters <- as.factor(as.numeric(df$seurat_clusters))
bapd8.integrated@meta.data <- df
Idents(bapd8.integrated) <- "new_clusters"

# Figure 2: A, B and C
A = enhancedDimPlot(object = bapd8.integrated,
                    grouping_var = 'ident',
                    reduction = "umap",
                    label = TRUE,
                    pt.size = 1,
                    alpha = 0.5) +
  ggtitle("A") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_classic() + 
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold"))

B <- enhancedDimPlot(object = bapd8.integrated,
                     grouping_var = 'replicate',
                     reduction = "umap",
                     label = FALSE,
                     pt.size = 1,
                     alpha = 0.4) +
  ggtitle("B") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_classic() + 
  theme(legend.justification = c(1, 1),
        legend.position = c(1, 1),
        plot.title = element_text(face = "bold")) +
  scale_colour_manual(name = "Conditions", 
                      labels = c(expression(paste('20% ', 'O'[2])),
                                 expression(paste('5% ', 'O'[2])),
                                 'nCT',
                                 'nTE'), 
                      values = c("20pcO2" = "#0571b0",
                                 "5pcO2" = "#ca0020",
                                 "nCT" = "#e333ff",
                                 "nTE" = "#3fff33" )) +
  scale_fill_manual(name = "Conditions", 
                    labels = c(expression(paste('20% ', 'O'[2])),
                               expression(paste('5% ', 'O'[2])),
                               'nCT',
                               'nTE'), 
                    values = c("20pcO2" = "#0571b0",
                               "5pcO2" = "#ca0020",
                               "nCT" = "#e333ff",
                               "nTE" = "#3fff33" )) +
  scale_linetype_manual(values = "blank")

C <- enhancedDimPlot(object = bapd8.integrated,
                     grouping_var = 'orig.ident',
                     reduction = "umap",
                     label = FALSE,
                     pt.size = 1,
                     alpha = 0.4) +
  ggtitle("C") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_classic() + 
  theme(legend.justification = c(1, 1),
        legend.position = c(1, 1),
        plot.title = element_text(face = "bold")) +
  scale_colour_manual(name = "Replicates", 
                      labels = c(expression(paste('20% ', 'O'[2], ' rep1')),
                                 expression(paste('20% ', 'O'[2], ' rep2')),
                                 expression(paste('5% ', 'O'[2], ' rep1')),
                                 expression(paste('5% ', 'O'[2], ' rep1')),
                                 "nCT day 5",
                                 "nCT day 10",
                                 "nTE day 3",
                                 "nTE day 5"), 
                      values = c("20pcO2_r1" = "#0571b0",
                                 "20pcO2_r2" = "#92c5de",
                                 "5pcO2_r1" = "#ca0020",
                                 "5pcO2_r2" = "#f4a582",
                                 "nCT_D5" = "#d133ff",
                                 "nCT_D10" = "#ff33f6",
                                 "nTE_D2" = "#33ffa2",
                                 "nTE_D3" = "#5bff33")) +
  scale_fill_manual(name = "Replicates", 
                    labels = c(expression(paste('20% ', 'O'[2], ' rep1')),
                               expression(paste('20% ', 'O'[2], ' rep2')),
                               expression(paste('5% ', 'O'[2], ' rep1')),
                               expression(paste('5% ', 'O'[2], ' rep1')),
                               "nCT day 5",
                               "nCT day 10",
                               "nTE day 3",
                               "nTE day 5"), 
                    values = c("20pcO2_r1" = "#0571b0",
                               "20pcO2_r2" = "#92c5de",
                               "5pcO2_r1" = "#ca0020",
                               "5pcO2_r2" = "#f4a582",
                               "nCT_D5" = "#d133ff",
                               "nCT_D10" = "#ff33f6",
                               "nTE_D2" = "#33ffa2",
                               "nTE_D3" = "#5bff33")) +
  scale_linetype_manual(values = "blank")

ggsave("Figure_2_A.png", plot = A, dpi=900, width = 5, height = 5)
ggsave("Figure_2_B.png", plot = B, dpi=900, width = 5, height = 5)
ggsave("Figure_2_C.png", plot = C, dpi=900, width = 5, height = 5)
ABC = (A | B | C)
ggsave("Figure_2_ABC.png", plot = ABC, dpi=900, width = 16, height = 6)


#multi_dittoPlot(cluster2356, vars = figs7, group.by = "new_clusters", split.by = "replicate",
#                vlnplot.lineweight = 0.2, jitter.size = 0.3, ncol = 1)

# finding conserved markers
DefaultAssay(bapd8.integrated) <- "RNA"
#for (i in 0:num.clusters){
#  try({
#    cluster.markers.all <- FindMarkers(bapd8.integrated, ident.1 = i)
#    cluster.markers.all$avg_FC <- exp(cluster.markers.all$avg_logFC)
#    cluster.markers <- cluster.markers.all %>% filter(exp.avg_FC >= 1.5) %>% filter(p_val_adj <= 0.05)
#    cluster=paste0("cluster", i, "names")
#    cluster <- rownames(cluster.markers)
#    write.csv(cluster.markers.pooled, file=paste0("markers.pooled.cluster_",i,".csv"))
#    write.csv(cluster.markers.pooled.roc, file=paste0("markers.pooled.roc.cluster_",i,".csv"))
#  })
#}
Idents(bapd8.integrated) <- "new_clusters"
markers.all.1 <- FindMarkers(bapd8.integrated, ident.1 = 1)
markers.all.2 <- FindMarkers(bapd8.integrated, ident.1 = 2)
markers.all.3 <- FindMarkers(bapd8.integrated, ident.1 = 3)
markers.all.4 <- FindMarkers(bapd8.integrated, ident.1 = 4)
markers.all.5 <- FindMarkers(bapd8.integrated, ident.1 = 5)
markers.all.6 <- FindMarkers(bapd8.integrated, ident.1 = 6)
markers.all.7 <- FindMarkers(bapd8.integrated, ident.1 = 7)
markers.all.8 <- FindMarkers(bapd8.integrated, ident.1 = 8)
markers.all.9 <- FindMarkers(bapd8.integrated, ident.1 = 9)
markers.all.10 <- FindMarkers(bapd8.integrated, ident.1 = 10)
markers.all.11 <- FindMarkers(bapd8.integrated, ident.1 = 11)
markers.all.12 <- FindMarkers(bapd8.integrated, ident.1 = 12)
markers.all.13 <- FindMarkers(bapd8.integrated, ident.1 = 13)

markers.filtered.1 <- markers.all.1 %>% filter(avg_log2FC >= 0.584962501) %>% filter(p_val_adj <= 0.05) %>% arrange(desc(avg_log2FC))
markers.filtered.2 <- markers.all.2 %>% filter(avg_log2FC >= 0.584962501) %>% filter(p_val_adj <= 0.05) %>% arrange(desc(avg_log2FC))
markers.filtered.3 <- markers.all.3 %>% filter(avg_log2FC >= 0.584962501) %>% filter(p_val_adj <= 0.05) %>% arrange(desc(avg_log2FC))
markers.filtered.4 <- markers.all.4 %>% filter(avg_log2FC >= 0.584962501) %>% filter(p_val_adj <= 0.05) %>% arrange(desc(avg_log2FC))
markers.filtered.5 <- markers.all.5 %>% filter(avg_log2FC >= 0.584962501) %>% filter(p_val_adj <= 0.05) %>% arrange(desc(avg_log2FC))
markers.filtered.6 <- markers.all.6 %>% filter(avg_log2FC >= 0.584962501) %>% filter(p_val_adj <= 0.05) %>% arrange(desc(avg_log2FC))
markers.filtered.7 <- markers.all.7 %>% filter(avg_log2FC >= 0.584962501) %>% filter(p_val_adj <= 0.05) %>% arrange(desc(avg_log2FC))
markers.filtered.8 <- markers.all.8 %>% filter(avg_log2FC >= 0.584962501) %>% filter(p_val_adj <= 0.05) %>% arrange(desc(avg_log2FC))
markers.filtered.9 <- markers.all.9 %>% filter(avg_log2FC >= 0.584962501) %>% filter(p_val_adj <= 0.05) %>% arrange(desc(avg_log2FC))
markers.filtered.10 <- markers.all.10 %>% filter(avg_log2FC >= 0.584962501) %>% filter(p_val_adj <= 0.05) %>% arrange(desc(avg_log2FC))
markers.filtered.11 <- markers.all.11 %>% filter(avg_log2FC >= 0.584962501) %>% filter(p_val_adj <= 0.05) %>% arrange(desc(avg_log2FC))
markers.filtered.12 <- markers.all.12 %>% filter(avg_log2FC >= 0.584962501) %>% filter(p_val_adj <= 0.05) %>% arrange(desc(avg_log2FC))
markers.filtered.13 <- markers.all.13 %>% filter(avg_log2FC >= 0.584962501) %>% filter(p_val_adj <= 0.05) %>% arrange(desc(avg_log2FC))



markers.filtered.names.1 <- rownames(markers.filtered.1)
markers.filtered.names.2 <- rownames(markers.filtered.2)
markers.filtered.names.3 <- rownames(markers.filtered.3)
markers.filtered.names.4 <- rownames(markers.filtered.4)
markers.filtered.names.5 <- rownames(markers.filtered.5)
markers.filtered.names.6 <- rownames(markers.filtered.6)
markers.filtered.names.7 <- rownames(markers.filtered.7)
markers.filtered.names.8 <- rownames(markers.filtered.8)
markers.filtered.names.9 <- rownames(markers.filtered.9)
markers.filtered.names.10 <- rownames(markers.filtered.10)
markers.filtered.names.11 <- rownames(markers.filtered.11)
markers.filtered.names.12 <- rownames(markers.filtered.12)
markers.filtered.names.13 <- rownames(markers.filtered.13)



length(markers.filtered.names.1)
length(markers.filtered.names.2)
length(markers.filtered.names.3)
length(markers.filtered.names.4)
length(markers.filtered.names.5)
length(markers.filtered.names.6)
length(markers.filtered.names.7)
length(markers.filtered.names.8)
length(markers.filtered.names.9)
length(markers.filtered.names.10)
length(markers.filtered.names.11)
length(markers.filtered.names.12)
length(markers.filtered.names.13)

# pie chart for # of cells
piedata <- tibble(
  cluster = bapd8.integrated@meta.data$new_clusters,
  cell_type = bapd8.integrated@meta.data$orig.ident
) %>%
  dplyr::group_by(cluster,cell_type) %>%
  dplyr::count() %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(
    percent=(100*n)/sum(n)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    cluster=paste("Cluster",cluster)) 

pie.clust1 <- subset(piedata, cluster %in% "Cluster 1")
pie.clust2 <- subset(piedata, cluster %in% "Cluster 2")
pie.clust3 <- subset(piedata, cluster %in% "Cluster 3")
pie.clust4 <- subset(piedata, cluster %in% "Cluster 4")
pie.clust5 <- subset(piedata, cluster %in% "Cluster 5")
pie.clust6 <- subset(piedata, cluster %in% "Cluster 6")
pie.clust7 <- subset(piedata, cluster %in% "Cluster 7")
pie.clust8 <- subset(piedata, cluster %in% "Cluster 8")
pie.clust9 <- subset(piedata, cluster %in% "Cluster 9")
pie.clust10 <- subset(piedata, cluster %in% "Cluster 10")
pie.clust11 <- subset(piedata, cluster %in% "Cluster 11")
pie.clust12 <- subset(piedata, cluster %in% "Cluster 12")
pie.clust13 <- subset(piedata, cluster %in% "Cluster 13")

pie.clust1  <- pie.clust1  %>% arrange(desc(cluster)) %>% mutate(lab.ypos = cumsum(percent) - 0.5*percent)
pie.clust2  <- pie.clust2  %>% arrange(desc(cluster)) %>% mutate(lab.ypos = cumsum(percent) - 0.5*percent)
pie.clust3  <- pie.clust3  %>% arrange(desc(cluster)) %>% mutate(lab.ypos = cumsum(percent) - 0.5*percent)
pie.clust4  <- pie.clust4  %>% arrange(desc(cluster)) %>% mutate(lab.ypos = cumsum(percent) - 0.5*percent)
pie.clust5  <- pie.clust5  %>% arrange(desc(cluster)) %>% mutate(lab.ypos = cumsum(percent) - 0.5*percent)
pie.clust6  <- pie.clust6  %>% arrange(desc(cluster)) %>% mutate(lab.ypos = cumsum(percent) - 0.5*percent)
pie.clust7  <- pie.clust7  %>% arrange(desc(cluster)) %>% mutate(lab.ypos = cumsum(percent) - 0.5*percent)
pie.clust8  <- pie.clust8  %>% arrange(desc(cluster)) %>% mutate(lab.ypos = cumsum(percent) - 0.5*percent)
pie.clust9  <- pie.clust9  %>% arrange(desc(cluster)) %>% mutate(lab.ypos = cumsum(percent) - 0.5*percent)
pie.clust10 <- pie.clust10 %>% arrange(desc(cluster)) %>% mutate(lab.ypos = cumsum(percent) - 0.5*percent)
pie.clust11 <- pie.clust11 %>% arrange(desc(cluster)) %>% mutate(lab.ypos = cumsum(percent) - 0.5*percent)
pie.clust12 <- pie.clust12 %>% arrange(desc(cluster)) %>% mutate(lab.ypos = cumsum(percent) - 0.5*percent)
pie.clust13 <- pie.clust13 %>% arrange(desc(cluster)) %>% mutate(lab.ypos = cumsum(percent) - 0.5*percent)

mybarplot <- function(pdata, i) {
  ggplot(data=pdata, aes(x=cell_type, y=n, fill=as.factor(cell_type))) +
    geom_bar(stat="identity") +
    theme_classic() +
    theme(legend.position="none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle(paste("Cluster", i)) + ylab("number of nuclei") + xlab("experiment")
  ggsave(paste("nuclei-clus", i, ".png", sep = "" ), dpi=900, width = 6, height = 4)
}
mybarplot(pie.clust1, 1)
mybarplot(pie.clust2, 2)
mybarplot(pie.clust3, 3)
mybarplot(pie.clust4, 4)
mybarplot(pie.clust5, 5)
mybarplot(pie.clust6, 6)
mybarplot(pie.clust7, 7)
mybarplot(pie.clust8, 8)
mybarplot(pie.clust9, 9)
mybarplot(pie.clust10, 10)
mybarplot(pie.clust11, 11)
mybarplot(pie.clust12, 12)
mybarplot(pie.clust13, 13)



# experiments
topmarkers <- function(number) {
  top.num.markers <- c(head(markers.filtered.names.1, number), head(markers.filtered.names.2, number), 
                       head(markers.filtered.names.3, number), head(markers.filtered.names.4, number), 
                       head(markers.filtered.names.5, number), head(markers.filtered.names.6, number), 
                       head(markers.filtered.names.7, number), head(markers.filtered.names.8, number), 
                       head(markers.filtered.names.9, number), head(markers.filtered.names.10, number),
                       head(markers.filtered.names.11, number), head(markers.filtered.names.12, number),
                       head(markers.filtered.names.13, number))
  return(top.num.markers)
}
genes <- topmarkers(5)
h = dittoHeatmap(bapd8.integrated, genes, annot.by = c("new_clusters", "replicate"))
ggsave("heatmap.svg", plot = h, dpi=900, width = 10, height = 10)

# colors for each clusters defined
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
color_list <- ggplotColours(n=13)

grouped_violinPlots <- function(markersfile, clusternumber, seuratobject = bapd8.integrated) {
  dittoPlotVarsAcrossGroups(seuratobject, markersfile, 
                            group.by = "new_clusters", main = paste("Cluster ", clusternumber, " markers"), 
                            xlab = "Clusters", 
                            ylab = "Mean z-score expression", 
                            x.labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4",
                                         "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8",
                                         "Cluster 9", "Cluster 10", "Cluster 11", "Cluster 12",
                                         "Cluster 13"), 
                            vlnplot.lineweight = 0.5, 
                            legend.show = FALSE, 
                            jitter.size = 0.5, 
                            color.panel = color_list)
  #ggsave(paste("cluster-markers_", clusternumber, ".svg"), dpi=900, width = 6, height = 5)
  ggsave(paste("cluster-markers_", clusternumber, ".png"), dpi=900, width = 6, height = 5)
}

grouped_violinPlots(markers.filtered.names.1, 1)
grouped_violinPlots(markers.filtered.names.2, 2)
grouped_violinPlots(markers.filtered.names.3, 3)
grouped_violinPlots(markers.filtered.names.4, 4)
grouped_violinPlots(markers.filtered.names.5, 5)
grouped_violinPlots(markers.filtered.names.6, 6)
grouped_violinPlots(markers.filtered.names.7, 7)
grouped_violinPlots(markers.filtered.names.8, 8)
grouped_violinPlots(markers.filtered.names.9, 9)
grouped_violinPlots(markers.filtered.names.10, 10)
grouped_violinPlots(markers.filtered.names.11, 11)
grouped_violinPlots(markers.filtered.names.12, 12)
grouped_violinPlots(markers.filtered.names.13, 13)


# load placenta cell enrcih dataset 
l <- load(file = "~/TutejaLab/PlacentaEnrich/combine-test-expression1.Rdata")
humanGeneMapping <- dataset$GRCH38$humanGeneMapping
d <- dataset$PlacentaDeciduaBloodData
data <- d$expressionData
cellDetails <- d$cellDetails

runpce <- function(inputgenelist, clusternumber) {
  inputGenes<-toupper(inputgenelist)
  humanGene<-humanGeneMapping[humanGeneMapping$Gene.name %in% inputGenes,]
  inputGenes<-humanGene$Gene
  expressionData<-data[intersect(row.names(data),humanGeneMapping$Gene),]
  se<-SummarizedExperiment(assays = SimpleList(as.matrix(expressionData)),rowData = row.names(expressionData),colData = colnames(expressionData))
  cellSpecificGenesExp<-teGeneRetrieval(se,expressedGeneThreshold = 1)
  print(length(inputGenes))
  gs<-GeneSet(geneIds=toupper(inputGenes))
  output2<-teEnrichmentCustom(gs,cellSpecificGenesExp)
  enrichmentOutput<-setNames(data.frame(assay(output2[[1]]),row.names = rowData(output2[[1]])[,1]),colData(output2[[1]])[,1])
  row.names(cellDetails)<-cellDetails$RName
  enrichmentOutput$Tissue<- cellDetails[row.names(enrichmentOutput),"CellName"]
  p <- ggplot(data = enrichmentOutput, mapping = aes(x = reorder (Tissue, -Log10PValue), Log10PValue, fill= Tissue)) +
    geom_bar(stat = "identity") + theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
          axis.text.y = element_text(size = 12),
          plot.margin = margin(10, 10, 10, 100), legend.position = "none",
          plot.title = element_text(color = "black", size=18, face="bold.italic"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color="black", size=14, face="bold")) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    ggtitle(paste("Cluster ", clusternumber )) + ylab("-log10 p-value")
  ggsave(paste("Cluster_", clusternumber, ".png"), dpi=700, width = 10, height = 8)
}

runpce(markers.filtered.names.1, 1)
runpce(markers.filtered.names.2, 2)
runpce(markers.filtered.names.3, 3)
runpce(markers.filtered.names.4, 4)
runpce(markers.filtered.names.5, 5)
runpce(markers.filtered.names.6, 6)
runpce(markers.filtered.names.7, 7)
runpce(markers.filtered.names.8, 8)
runpce(markers.filtered.names.9, 9)
runpce(markers.filtered.names.10, 10)
runpce(markers.filtered.names.11, 11)
runpce(markers.filtered.names.12, 12)
runpce(markers.filtered.names.13, 13)

# https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

features<- c("KRT8", "S100P", "XAGE2", "ERVV-1", "TBX3", "GRHL1") 
StackedVlnPlot(obj = bapd8.integrated, features = features)


ggsave("Figure_5_B.svg", dpi=700, width = 8, height = 12)


figs8a <- FeaturePlot(bapd8.integrated, features = "SOX2", min.cutoff = "q9")
figs8b <- VlnPlot(bapd8.integrated, "SOX2", group.by = "new_clusters")
figs8ab <- (figs8a | figs8b)
ggsave("Figure_S8_A.svg", plot = figs8a, dpi=900, width = 5, height = 5)
ggsave("Figure_S8_B.svg", plot = figs8b, dpi=900, width = 5, height = 5)
ggsave("Figure_S8_AB.svg", plot = figs8ab, dpi=900, width = 10, height = 5)


#mesoderm <- c("DLL3", "FOXC1", "RIPPLY", "TBX3", "TWIST2", "FOXA2", "MIXL1")
#endoderm <- c("AFP", "GATA4", "GATA6", "GDF1", "GDF3", "MIXL2")
#ectoderm <- c("NES", "FGF5", "OTX2", "SOX1", "PAX6")
#StackedVlnPlot(obj = bapd8.integrated, features = mesoderm)
#StackedVlnPlot(obj = bapd8.integrated, features = endoderm)
#StackedVlnPlot(obj = bapd8.integrated, features = ectoderm)

bapd8.markers <- FindAllMarkers(bapd8.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
bapd8.markers.ranked.2.percluster <- bapd8.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
bapd8.markers.ranked.4.percluster <- bapd8.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_logFC)
write.csv(bapd8.markers.ranked.2.percluster,"bapd8.markers.ranked.2.percluster.csv", row.names = TRUE)
bapd8.markers.ranked.20.percluster <- bapd8.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.csv(bapd8.markers.ranked.20.percluster,"bapd8.markers.ranked.20.percluster.csv", row.names = TRUE)


markers.conserved.1 <- FindConservedMarkers(bapd8.integrated, ident.1 = 1, grouping.var = "replicate", verbose = FALSE)
markers.conserved.2 <- FindConservedMarkers(bapd8.integrated, ident.1 = 2, grouping.var = "replicate", verbose = FALSE)
markers.conserved.3 <- FindConservedMarkers(bapd8.integrated, ident.1 = 3, grouping.var = "replicate", verbose = FALSE)
markers.conserved.4 <- FindConservedMarkers(bapd8.integrated, ident.1 = 4, grouping.var = "replicate", verbose = FALSE)
markers.conserved.5 <- FindConservedMarkers(bapd8.integrated, ident.1 = 5, grouping.var = "replicate", verbose = FALSE)
markers.conserved.6 <- FindConservedMarkers(bapd8.integrated, ident.1 = 6, grouping.var = "replicate", verbose = FALSE)
markers.conserved.7 <- FindConservedMarkers(bapd8.integrated, ident.1 = 7, grouping.var = "replicate", verbose = FALSE)
markers.conserved.8 <- FindConservedMarkers(bapd8.integrated, ident.1 = 8, grouping.var = "replicate", verbose = FALSE)
markers.conserved.9 <- FindConservedMarkers(bapd8.integrated, ident.1 = 9, grouping.var = "replicate", verbose = FALSE)

#write.csv(cluster.markers.conserved, file=paste0("markers.conserved.cluster_",i,".csv"))



markers.to.plot <- bapd8.markers.ranked.4.percluster[7]
pdf("conserved-markers.pdf", width=12, height=10)
DotPlot(bapd8.integrated, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "replicate") + RotatedAxis()
dev.off()


cluster.averages.data <- AverageExpression(bapd8.integrated, slot = "data", assays = "RNA")
condition.averages.data <- AverageExpression(bapd8.integrated, slot = "data", add.ident = "orig.ident", assays = "RNA")
replicate.averages.data <- AverageExpression(bapd8.integrated, slot = "data", add.ident = "replicate", assays = "RNA")

avg.data <- cluster.averages.data[["RNA"]]
avg.condition <- condition.averages.data[["RNA"]]
avg.replicate <- replicate.averages.data[["RNA"]]

write.table(avg.data, file="snn-average-data.tsv", sep= "\t")
write.table(avg.condition, file="snn-average-condition.tsv", sep= "\t")
write.table(avg.replicate, file="snn-average-replicate.tsv", sep= "\t")

cells <- bapd8.integrated@meta.data	%>%
  dplyr::group_by(orig.ident,	seurat_clusters)	%>%	
  dplyr::summarise(length(seurat_clusters))
write.table(cells , file = "cells-number-stats.tsv", sep="\t")

head(bapd8.integrated@meta.data)
df <- bapd8.integrated@meta.data
df$stim <- (paste(df$replicate,df$new_clusters, sep = "."))
df$stim <- gsub('5pcO2', 'FIVE', df$stim)
df$stim <- gsub('20pcO2', 'TWENTY', df$stim)
bapd8.integrated@meta.data <- df

Idents(bapd8.integrated) <- "stim"
clus1.five.twenty <- FindMarkers(bapd8.integrated, ident.1 = "FIVE.1", ident.2 = "TWENTY.1", verbose = FALSE, logfc.threshold = 0)
clus1.five.twenty$log2fc <- log2(exp(clus1.five.twenty$avg_logFC))
clus1.five.twenty$Gene <- row.names(clus1.five.twenty)
clus1.five.twenty$diffexpressed <- "N0"
clus1.five.twenty$diffexpressed[clus1.five.twenty$log2fc >= 1 & clus1.five.twenty$p_val < 0.05] <- "UP"
clus1.five.twenty$diffexpressed[clus1.five.twenty$log2fc < -1 & clus1.five.twenty$p_val < 0.05] <- "DOWN"
clus1.five.twenty$delabel <- ""
clus1.five.twenty$delabel[clus1.five.twenty$log2fc >= 1 & clus1.five.twenty$p_val < 0.05] <- clus1.five.twenty$Gene[clus1.five.twenty$log2fc >= 1 & clus1.five.twenty$p_val < 0.05]
clus1.five.twenty$delabel[clus1.five.twenty$log2fc < -1 & clus1.five.twenty$p_val < 0.05] <- clus1.five.twenty$Gene[clus1.five.twenty$log2fc < -1 & clus1.five.twenty$p_val < 0.05]
v.clus1 <- ggplot(data=clus1.five.twenty, aes(x=log2fc, y=-log10(p_val), col=diffexpressed, label=delabel)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  scale_color_manual(name = "Expression", values=c("#0571b0", "#4d4d4d", "#ca0020"), labels = c(expression(paste("Higher in 20% ", "O"[2])), "All other genes", expression(paste("Higher in 5% ", "O"[2])) )) +
  geom_text_repel(show.legend = F) +
  ggtitle(expression(paste("Cluster 1: 20% ", "O"[2], " vs. 5% ", "O"[2]))) +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " pvalue"))) +
  theme(legend.text.align = 0)

clus2.five.twenty <- FindMarkers(bapd8.integrated, ident.1 = "FIVE.2", ident.2 = "TWENTY.2", verbose = FALSE, logfc.threshold = 0)
clus2.five.twenty$log2fc <- log2(exp(clus2.five.twenty$avg_logFC))
clus2.five.twenty$Gene <- row.names(clus2.five.twenty)
clus2.five.twenty$diffexpressed <- "N0"
clus2.five.twenty$diffexpressed[clus2.five.twenty$log2fc >= 1 & clus2.five.twenty$p_val < 0.05] <- "UP"
clus2.five.twenty$diffexpressed[clus2.five.twenty$log2fc < -1 & clus2.five.twenty$p_val < 0.05] <- "DOWN"
clus2.five.twenty$delabel <- ""
clus2.five.twenty$delabel[clus2.five.twenty$log2fc >= 1 & clus2.five.twenty$p_val < 0.05] <- clus2.five.twenty$Gene[clus2.five.twenty$log2fc >= 1 & clus2.five.twenty$p_val < 0.05]
clus2.five.twenty$delabel[clus2.five.twenty$log2fc < -1 & clus2.five.twenty$p_val < 0.05] <- clus2.five.twenty$Gene[clus2.five.twenty$log2fc < -1 & clus2.five.twenty$p_val < 0.05]
v.clus2 <- ggplot(data=clus2.five.twenty, aes(x=log2fc, y=-log10(p_val), col=diffexpressed, label=delabel)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  scale_color_manual(name = "Expression", values=c("#0571b0", "#4d4d4d", "#ca0020"), labels = c(expression(paste("Higher in 20% ", "O"[2])), "All other genes", expression(paste("Higher in 5% ", "O"[2])) )) +
  geom_text_repel(show.legend = F) +
  ggtitle(expression(paste("Cluster 2: 20% ", "O"[2], " vs. 5% ", "O"[2]))) +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " pvalue"))) +
  theme(legend.text.align = 0)

clus3.five.twenty <- FindMarkers(bapd8.integrated, ident.1 = "FIVE.3", ident.2 = "TWENTY.3", verbose = FALSE, logfc.threshold = 0)
clus3.five.twenty$log2fc <- log2(exp(clus3.five.twenty$avg_logFC))
clus3.five.twenty$Gene <- row.names(clus3.five.twenty)
clus3.five.twenty$diffexpressed <- "N0"
clus3.five.twenty$diffexpressed[clus3.five.twenty$log2fc >= 1 & clus3.five.twenty$p_val < 0.05] <- "UP"
clus3.five.twenty$diffexpressed[clus3.five.twenty$log2fc < -1 & clus3.five.twenty$p_val < 0.05] <- "DOWN"
clus3.five.twenty$delabel <- ""
clus3.five.twenty$delabel[clus3.five.twenty$log2fc >= 1 & clus3.five.twenty$p_val < 0.05] <- clus3.five.twenty$Gene[clus3.five.twenty$log2fc >= 1 & clus3.five.twenty$p_val < 0.05]
clus3.five.twenty$delabel[clus3.five.twenty$log2fc < -1 & clus3.five.twenty$p_val < 0.05] <- clus3.five.twenty$Gene[clus3.five.twenty$log2fc < -1 & clus3.five.twenty$p_val < 0.05]
v.clus3 <- ggplot(data=clus3.five.twenty, aes(x=log2fc, y=-log10(p_val), col=diffexpressed, label=delabel)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  scale_color_manual(name = "Expression", values=c("#0571b0", "#4d4d4d", "#ca0020"), labels = c(expression(paste("Higher in 20% ", "O"[2])), "All other genes", expression(paste("Higher in 5% ", "O"[2])) )) +
  geom_text_repel(show.legend = F) +
  ggtitle(expression(paste("Cluster 3: 20% ", "O"[2], " vs. 5% ", "O"[2]))) +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " pvalue"))) +
  theme(legend.text.align = 0)

clus4.five.twenty <- FindMarkers(bapd8.integrated, ident.1 = "FIVE.4", ident.2 = "TWENTY.4", verbose = FALSE, logfc.threshold = 0)
clus4.five.twenty$log2fc <- log2(exp(clus4.five.twenty$avg_logFC))
clus4.five.twenty$Gene <- row.names(clus4.five.twenty)
clus4.five.twenty$diffexpressed <- "N0"
clus4.five.twenty$diffexpressed[clus4.five.twenty$log2fc >= 1 & clus4.five.twenty$p_val < 0.05] <- "UP"
clus4.five.twenty$diffexpressed[clus4.five.twenty$log2fc < -1 & clus4.five.twenty$p_val < 0.05] <- "DOWN"
clus4.five.twenty$delabel <- ""
clus4.five.twenty$delabel[clus4.five.twenty$log2fc >= 1 & clus4.five.twenty$p_val < 0.05] <- clus4.five.twenty$Gene[clus4.five.twenty$log2fc >= 1 & clus4.five.twenty$p_val < 0.05]
clus4.five.twenty$delabel[clus4.five.twenty$log2fc < -1 & clus4.five.twenty$p_val < 0.05] <- clus4.five.twenty$Gene[clus4.five.twenty$log2fc < -1 & clus4.five.twenty$p_val < 0.05]
v.clus4 <- ggplot(data=clus4.five.twenty, aes(x=log2fc, y=-log10(p_val), col=diffexpressed, label=delabel)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  scale_color_manual(name = "Expression", values=c("#0571b0", "#4d4d4d", "#ca0020"), labels = c(expression(paste("Higher in 20% ", "O"[2])), "All other genes", expression(paste("Higher in 5% ", "O"[2])) )) +
  geom_text_repel(show.legend = F) +
  ggtitle(expression(paste("Cluster 4: 20% ", "O"[2], " vs. 5% ", "O"[2]))) +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " pvalue"))) +
  theme(legend.text.align = 0)

clus5.five.twenty <- FindMarkers(bapd8.integrated, ident.1 = "FIVE.5", ident.2 = "TWENTY.5", verbose = FALSE, logfc.threshold = 0)
clus5.five.twenty$log2fc <- log2(exp(clus5.five.twenty$avg_logFC))
clus5.five.twenty$Gene <- row.names(clus5.five.twenty)
clus5.five.twenty$diffexpressed <- "N0"
clus5.five.twenty$diffexpressed[clus5.five.twenty$log2fc >= 1 & clus5.five.twenty$p_val < 0.05] <- "UP"
clus5.five.twenty$diffexpressed[clus5.five.twenty$log2fc < -1 & clus5.five.twenty$p_val < 0.05] <- "DOWN"
clus5.five.twenty$delabel <- ""
clus5.five.twenty$delabel[clus5.five.twenty$log2fc >= 1 & clus5.five.twenty$p_val < 0.05] <- clus5.five.twenty$Gene[clus5.five.twenty$log2fc >= 1 & clus5.five.twenty$p_val < 0.05]
clus5.five.twenty$delabel[clus5.five.twenty$log2fc < -1 & clus5.five.twenty$p_val < 0.05] <- clus5.five.twenty$Gene[clus5.five.twenty$log2fc < -1 & clus5.five.twenty$p_val < 0.05]
v.clus5 <- ggplot(data=clus5.five.twenty, aes(x=log2fc, y=-log10(p_val), col=diffexpressed, label=delabel)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  scale_color_manual(name = "Expression", values=c("#0571b0", "#4d4d4d", "#ca0020"), labels = c(expression(paste("Higher in 20% ", "O"[2])), "All other genes", expression(paste("Higher in 5% ", "O"[2])) )) +
  geom_text_repel(show.legend = F) +
  ggtitle(expression(paste("Cluster 5: 20% ", "O"[2], " vs. 5% ", "O"[2]))) +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " pvalue"))) +
  theme(legend.text.align = 0)

clus6.five.twenty <- FindMarkers(bapd8.integrated, ident.1 = "FIVE.6", ident.2 = "TWENTY.6", verbose = FALSE, logfc.threshold = 0)
clus6.five.twenty$log2fc <- log2(exp(clus6.five.twenty$avg_logFC))
clus6.five.twenty$Gene <- row.names(clus6.five.twenty)
clus6.five.twenty$diffexpressed <- "N0"
clus6.five.twenty$diffexpressed[clus6.five.twenty$log2fc >= 1 & clus6.five.twenty$p_val < 0.05] <- "UP"
clus6.five.twenty$diffexpressed[clus6.five.twenty$log2fc < -1 & clus6.five.twenty$p_val < 0.05] <- "DOWN"
clus6.five.twenty$delabel <- ""
clus6.five.twenty$delabel[clus6.five.twenty$log2fc >= 1 & clus6.five.twenty$p_val < 0.05] <- clus6.five.twenty$Gene[clus6.five.twenty$log2fc >= 1 & clus6.five.twenty$p_val < 0.05]
clus6.five.twenty$delabel[clus6.five.twenty$log2fc < -1 & clus6.five.twenty$p_val < 0.05] <- clus6.five.twenty$Gene[clus6.five.twenty$log2fc < -1 & clus6.five.twenty$p_val < 0.05]
v.clus6 <- ggplot(data=clus6.five.twenty, aes(x=log2fc, y=-log10(p_val), col=diffexpressed, label=delabel)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  scale_color_manual(name = "Expression", values=c("#0571b0", "#4d4d4d", "#ca0020"), labels = c(expression(paste("Higher in 20% ", "O"[2])), "All other genes", expression(paste("Higher in 5% ", "O"[2])) )) +
  geom_text_repel(show.legend = F) +
  ggtitle(expression(paste("Cluster 6: 20% ", "O"[2], " vs. 5% ", "O"[2]))) +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " pvalue"))) +
  theme(legend.text.align = 0)

clus7.five.twenty <- FindMarkers(bapd8.integrated, ident.1 = "FIVE.7", ident.2 = "TWENTY.7", verbose = FALSE, logfc.threshold = 0)
clus7.five.twenty$log2fc <- log2(exp(clus7.five.twenty$avg_logFC))
clus7.five.twenty$Gene <- row.names(clus7.five.twenty)
clus7.five.twenty$diffexpressed <- "N0"
clus7.five.twenty$diffexpressed[clus7.five.twenty$log2fc >= 1 & clus7.five.twenty$p_val < 0.05] <- "UP"
clus7.five.twenty$diffexpressed[clus7.five.twenty$log2fc < -1 & clus7.five.twenty$p_val < 0.05] <- "DOWN"
clus7.five.twenty$delabel <- ""
clus7.five.twenty$delabel[clus7.five.twenty$log2fc >= 1 & clus7.five.twenty$p_val < 0.05] <- clus7.five.twenty$Gene[clus7.five.twenty$log2fc >= 1 & clus7.five.twenty$p_val < 0.05]
clus7.five.twenty$delabel[clus7.five.twenty$log2fc < -1 & clus7.five.twenty$p_val < 0.05] <- clus7.five.twenty$Gene[clus7.five.twenty$log2fc < -1 & clus7.five.twenty$p_val < 0.05]
v.clus7 <- ggplot(data=clus7.five.twenty, aes(x=log2fc, y=-log10(p_val), col=diffexpressed, label=delabel)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  scale_color_manual(name = "Expression", values=c("#0571b0", "#4d4d4d", "#ca0020"), labels = c(expression(paste("Higher in 20% ", "O"[2])), "All other genes", expression(paste("Higher in 5% ", "O"[2])) )) +
  geom_text_repel(show.legend = F) +
  ggtitle(expression(paste("Cluster 7: 20% ", "O"[2], " vs. 5% ", "O"[2]))) +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " pvalue"))) +
  theme(legend.text.align = 0)

clus8.five.twenty <- FindMarkers(bapd8.integrated, ident.1 = "FIVE.8", ident.2 = "TWENTY.8", verbose = FALSE, logfc.threshold = 0)
clus8.five.twenty$log2fc <- log2(exp(clus8.five.twenty$avg_logFC))
clus8.five.twenty$Gene <- row.names(clus8.five.twenty)
clus8.five.twenty$diffexpressed <- "N0"
clus8.five.twenty$diffexpressed[clus8.five.twenty$log2fc >= 1 & clus8.five.twenty$p_val < 0.05] <- "UP"
clus8.five.twenty$diffexpressed[clus8.five.twenty$log2fc < -1 & clus8.five.twenty$p_val < 0.05] <- "DOWN"
clus8.five.twenty$delabel <- ""
clus8.five.twenty$delabel[clus8.five.twenty$log2fc >= 1 & clus8.five.twenty$p_val < 0.05] <- clus8.five.twenty$Gene[clus8.five.twenty$log2fc >= 1 & clus8.five.twenty$p_val < 0.05]
clus8.five.twenty$delabel[clus8.five.twenty$log2fc < -1 & clus8.five.twenty$p_val < 0.05] <- clus8.five.twenty$Gene[clus8.five.twenty$log2fc < -1 & clus8.five.twenty$p_val < 0.05]
v.clus8 <- ggplot(data=clus8.five.twenty, aes(x=log2fc, y=-log10(p_val), col=diffexpressed, label=delabel)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  scale_color_manual(name = "Expression", values=c("#0571b0", "#4d4d4d", "#ca0020"), labels = c(expression(paste("Higher in 20% ", "O"[2])), "All other genes", expression(paste("Higher in 5% ", "O"[2])) )) +
  geom_text_repel(show.legend = F) +
  ggtitle(expression(paste("Cluster 8: 20% ", "O"[2], " vs. 5% ", "O"[2]))) +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " pvalue"))) +
  theme(legend.text.align = 0)

clus9.five.twenty <- FindMarkers(bapd8.integrated, ident.1 = "FIVE.9", ident.2 = "TWENTY.9", verbose = FALSE, logfc.threshold = 0)
clus9.five.twenty$log2fc <- log2(exp(clus9.five.twenty$avg_logFC))
clus9.five.twenty$Gene <- row.names(clus9.five.twenty)
clus9.five.twenty$diffexpressed <- "N0"
clus9.five.twenty$diffexpressed[clus9.five.twenty$log2fc >= 1 & clus9.five.twenty$p_val < 0.05] <- "UP"
clus9.five.twenty$diffexpressed[clus9.five.twenty$log2fc < -1 & clus9.five.twenty$p_val < 0.05] <- "DOWN"
clus9.five.twenty$delabel <- ""
clus9.five.twenty$delabel[clus9.five.twenty$log2fc >= 1 & clus9.five.twenty$p_val < 0.05] <- clus9.five.twenty$Gene[clus9.five.twenty$log2fc >= 1 & clus9.five.twenty$p_val < 0.05]
clus9.five.twenty$delabel[clus9.five.twenty$log2fc < -1 & clus9.five.twenty$p_val < 0.05] <- clus9.five.twenty$Gene[clus9.five.twenty$log2fc < -1 & clus9.five.twenty$p_val < 0.05]
v.clus9 <- ggplot(data=clus9.five.twenty, aes(x=log2fc, y=-log10(p_val), col=diffexpressed, label=delabel)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  scale_color_manual(name = "Expression", values=c("#0571b0", "#4d4d4d", "#ca0020"), labels = c(expression(paste("Higher in 20% ", "O"[2])), "All other genes", expression(paste("Higher in 5% ", "O"[2])) )) +
  geom_text_repel(show.legend = F) +
  ggtitle(expression(paste("Cluster 9: 20% ", "O"[2], " vs. 5% ", "O"[2]))) +
  xlab(expression(paste("log"[2], " fold change"))) +
  ylab(expression(paste("-log"[10], " pvalue"))) +
  theme(legend.text.align = 0)

ggsave("Figure_6_E.svg", plot = v.clus1, dpi=900, width = 8, height = 6)
ggsave("Figure_6_A.svg", plot = v.clus2, dpi=900, width = 8, height = 6)
ggsave("Figure_6_B.svg", plot = v.clus3, dpi=900, width = 8, height = 6)
ggsave("Figure_6_F.svg", plot = v.clus4, dpi=900, width = 8, height = 6)
ggsave("Figure_6_C.svg", plot = v.clus5, dpi=900, width = 8, height = 6)
ggsave("Figure_6_D.svg", plot = v.clus6, dpi=900, width = 8, height = 6)
ggsave("Figure_6_G.svg", plot = v.clus7, dpi=900, width = 8, height = 6)
ggsave("Figure_6_H.svg", plot = v.clus8, dpi=900, width = 8, height = 6)
ggsave("Figure_6_I.svg", plot = v.clus9, dpi=900, width = 8, height = 6)


Idents(bapd8.integrated) <- "new_clusters"

cluster5n6 <- subset(bapd8.integrated, idents = c("5", "6"))
cluster2356 <- subset(bapd8.integrated, idents = c("2", "3", "5", "6"))

mesoderm <- c("FOXC1", "TWIST2")
endoderm <- c("AFP", "GATA6")
ectoderm <- c("NES", "PAX6")

warburg <- c("GAPDH", "LDHA", "CLIC3", "SLC2A3")
nonwarburg <- c("FN1", "COL3A1", "LUM")
cpm <- c("CCNB1", "MKI67", "PCNA", "CENPF")
a = VlnPlot(object = cluster2356, features = warburg, split.by = "replicate", assay="RNA", pt.size=0.5) 
b = VlnPlot(object = cluster2356, features = nonwarburg, split.by = "replicate", assay="RNA", pt.size=0.5) 
c <- StackedVlnPlot(obj = cluster2356, features = c(nonwarburg, warburg), split.by = "replicate")
ggsave("Figure_nonwarburg.png", plot = b, dpi=900, width = 12, height = 5)
ggsave("Figure_warburg.png", plot = a, dpi=900, width = 12, height = 10)
ggsave("Figure_warburg-nonwarburg.png", plot = c, dpi=900, width = 5, height = 8)

d = VlnPlot(object = bapd8.integrated, features = cpm, split.by = "replicate", assay="RNA", pt.size=0.5)
ggsave("Figure_cpm.png", plot = d, dpi=900, width = 12, height = 10)
figs6 <- c("TLE4", "PCDH9", "MAML2", "ACTB", "ACTG1", "MYL6", "TPM1", "TMSB4", "TMSB10", "TAGLN", "S100A11", "S100A6", "S100A10", "IL32", "CA3")
figs7 <- c("SLC2A3", "CLIC3", "FN1", "APOE",  "COL3A1", "LUM")
figs6a <- c("TLE4", "PCDH9", "MAML2", "TMSB4Y", "CA3")
figs6b <- c("TAGLN", "TMSB10", "S100A11", "S100A6", "S100A10")
figs6c <- c("ACTB", "ACTG1", "IL32", "MYL6", "TPM1")
s6a <- StackedVlnPlot(obj = bapd8.integrated, features = figs6a)
s6b <- StackedVlnPlot(obj = bapd8.integrated, features = figs6b)
s6c <- StackedVlnPlot(obj = bapd8.integrated, features = figs6c)
(s6a | s6b | s6c)

fig3 <- c("GATA3", "TFAP2A", "KRT7", "KRT23", "CGA", "PGF", "SLC40A1", "XAGE2", "CYP11A1", "HSD3B1", "MALAT1", "NEAT1")

colors.samples <- c("#0571b0", "#ca0020")

DefaultAssay(cluster2356) <- "RNA"
s7 <- multi_dittoPlot(cluster2356, vars = figs7, group.by = "new_clusters", split.by = "replicate",
                      vlnplot.lineweight = 0.2, jitter.size = 0.3, ncol = 1)

ggsave("Figure_S7.svg", plot = s7, dpi=900, width = 8, height = 12)
?multi_dittoPlot

cluster1 <- subset(bapd8.integrated, idents = "1")
cluster2 <- subset(bapd8.integrated, idents = "2")
cluster3 <- subset(bapd8.integrated, idents = "3")
cluster4 <- subset(bapd8.integrated, idents = "4")
cluster5 <- subset(bapd8.integrated, idents = "5")
cluster6 <- subset(bapd8.integrated, idents = "6")
cluster7 <- subset(bapd8.integrated, idents = "7")
cluster8 <- subset(bapd8.integrated, idents = "8")
cluster9 <- subset(bapd8.integrated, idents = "9")

Idents(cluster1) <- "replicate"
Idents(cluster2) <- "replicate"
Idents(cluster3) <- "replicate"
Idents(cluster4) <- "replicate"
Idents(cluster5) <- "replicate"
Idents(cluster6) <- "replicate"
Idents(cluster7) <- "replicate"
Idents(cluster8) <- "replicate"
Idents(cluster9) <- "replicate"
avg.cluster1 <- log1p(AverageExpression(cluster1, verbose = FALSE)$RNA)
avg.cluster2 <- log1p(AverageExpression(cluster2, verbose = FALSE)$RNA)
avg.cluster3 <- log1p(AverageExpression(cluster3, verbose = FALSE)$RNA)
avg.cluster4 <- log1p(AverageExpression(cluster4, verbose = FALSE)$RNA)
avg.cluster5 <- log1p(AverageExpression(cluster5, verbose = FALSE)$RNA)
avg.cluster6 <- log1p(AverageExpression(cluster6, verbose = FALSE)$RNA)
avg.cluster7 <- log1p(AverageExpression(cluster7, verbose = FALSE)$RNA)
avg.cluster8 <- log1p(AverageExpression(cluster8, verbose = FALSE)$RNA)
avg.cluster9 <- log1p(AverageExpression(cluster9, verbose = FALSE)$RNA)

avg.cluster1$gene <- rownames(cluster1)
avg.cluster2$gene <- rownames(cluster2)
avg.cluster3$gene <- rownames(cluster3)
avg.cluster4$gene <- rownames(cluster4)
avg.cluster5$gene <- rownames(cluster5)
avg.cluster6$gene <- rownames(cluster6)
avg.cluster7$gene <- rownames(cluster7)
avg.cluster8$gene <- rownames(cluster8)
avg.cluster9$gene <- rownames(cluster9)



label.clus1 <- rownames(head(clus1.five.twenty, n=20))
label.clus2 <- rownames(head(clus2.five.twenty, n=20))
label.clus3 <- rownames(head(clus3.five.twenty, n=20))
label.clus4 <- rownames(head(clus4.five.twenty, n=20))
label.clus5 <- rownames(head(clus5.five.twenty, n=20))
label.clus6 <- rownames(head(clus6.five.twenty, n=20))
label.clus7 <- rownames(head(clus7.five.twenty, n=20))
label.clus8 <- rownames(head(clus8.five.twenty, n=20))
label.clus9 <- rownames(head(clus9.five.twenty, n=20))
colnames(avg.cluster9) <- c("Five.pc.O2", "Twenty.pc.O2", "gene" )
colnames(avg.cluster1) <- c("Five.pc.O2", "Twenty.pc.O2", "gene" )
colnames(avg.cluster2) <- c("Five.pc.O2", "Twenty.pc.O2", "gene" )
colnames(avg.cluster3) <- c("Five.pc.O2", "Twenty.pc.O2", "gene" )
colnames(avg.cluster4) <- c("Five.pc.O2", "Twenty.pc.O2", "gene" )
colnames(avg.cluster5) <- c("Five.pc.O2", "Twenty.pc.O2", "gene" )
colnames(avg.cluster6) <- c("Five.pc.O2", "Twenty.pc.O2", "gene" )
colnames(avg.cluster7) <- c("Five.pc.O2", "Twenty.pc.O2", "gene" )
colnames(avg.cluster8) <- c("Five.pc.O2", "Twenty.pc.O2", "gene" )
colnames(avg.cluster9) <- c("Five.pc.O2", "Twenty.pc.O2", "gene" )
p9 <- ggplot(avg.cluster9, aes(Twenty.pc.O2, Five.pc.O2)) + geom_point() + ggtitle("Cluster 9")
p1 <- ggplot(avg.cluster1, aes(Twenty.pc.O2, Five.pc.O2)) + geom_point() + ggtitle("Cluster 1")
p2 <- ggplot(avg.cluster2, aes(Twenty.pc.O2, Five.pc.O2)) + geom_point() + ggtitle("Cluster 2")
p3 <- ggplot(avg.cluster3, aes(Twenty.pc.O2, Five.pc.O2)) + geom_point() + ggtitle("Cluster 3")
p4 <- ggplot(avg.cluster4, aes(Twenty.pc.O2, Five.pc.O2)) + geom_point() + ggtitle("Cluster 4")
p5 <- ggplot(avg.cluster5, aes(Twenty.pc.O2, Five.pc.O2)) + geom_point() + ggtitle("Cluster 5")
p6 <- ggplot(avg.cluster6, aes(Twenty.pc.O2, Five.pc.O2)) + geom_point() + ggtitle("Cluster 6")
p7 <- ggplot(avg.cluster7, aes(Twenty.pc.O2, Five.pc.O2)) + geom_point() + ggtitle("Cluster 7")
p8 <- ggplot(avg.cluster8, aes(Twenty.pc.O2, Five.pc.O2)) + geom_point() + ggtitle("Cluster 8")


p0 <- LabelPoints(plot = p0, points = label.clus0, repel = TRUE)
p1 <- LabelPoints(plot = p1, points = label.clus1, repel = TRUE)
p2 <- LabelPoints(plot = p2, points = label.clus2, repel = TRUE)
p3 <- LabelPoints(plot = p3, points = label.clus3, repel = TRUE)
p4 <- LabelPoints(plot = p4, points = label.clus4, repel = TRUE)
p5 <- LabelPoints(plot = p5, points = label.clus5, repel = TRUE)
p6 <- LabelPoints(plot = p6, points = label.clus6, repel = TRUE)
p7 <- LabelPoints(plot = p7, points = label.clus7, repel = TRUE)
p8 <- LabelPoints(plot = p8, points = label.clus8, repel = TRUE)
p1
plot_grid(p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)
pdf("dge_cluster_condition.pdf")
p0
p1
p2
p3
p4
p5
p6
p7
p8
dev.off()
write.csv(clus0.five.twenty,"clus0.five.twenty_DGE.csv", row.names = TRUE)
write.csv(clus1.five.twenty,"clus1.five.twenty_DGE.csv", row.names = TRUE)
write.csv(clus2.five.twenty,"clus2.five.twenty_DGE.csv", row.names = TRUE)
write.csv(clus3.five.twenty,"clus3.five.twenty_DGE.csv", row.names = TRUE)
write.csv(clus4.five.twenty,"clus4.five.twenty_DGE.csv", row.names = TRUE)
write.csv(clus5.five.twenty,"clus5.five.twenty_DGE.csv", row.names = TRUE)
write.csv(clus6.five.twenty,"clus6.five.twenty_DGE.csv", row.names = TRUE)
write.csv(clus7.five.twenty,"clus7.five.twenty_DGE.csv", row.names = TRUE)
write.csv(clus8.five.twenty,"clus8.five.twenty_DGE.csv", row.names = TRUE)

bapd8.integrated <- BuildClusterTree(
  bapd8.integrated,
  dims = 1:30,
  reorder = TRUE,
  reorder.numeric = TRUE
)
bapd8.integrated[['cluster']] <- factor(
  as.character(bapd8.integrated@meta.data$tree.ident),
  levels = sort(unique(bapd8.integrated@meta.data$tree.ident))
)
bapd8.integrated@meta.data$seurat_clusters <- NULL
bapd8.integrated@meta.data$RNA_snn_res.0.5 <- NULL
bapd8.integrated@meta.data$tree.ident <- NULL
bapd8.integrated <- RunTSNE(
  bapd8.integrated,
  reduction.name = 'tSNE',
  reduction.key = 'tSNE_',
  dims = 1:30,
  dim.embed = 2,
  perplexity = 30,
  seed.use = 100
)
bapd8.integrated <- RunTSNE(
  bapd8.integrated,
  reduction.name = 'tSNE_3D',
  reduction.key = 'tSNE3D_',
  dims = 1:30,
  dim.embed = 3,
  perplexity = 30,
  seed.use = 100
)
bapd8.integrated <- RunUMAP(
  bapd8.integrated,
  reduction.name = 'UMAP_3D',
  reduction.key = 'UMAP3D_',
  dims = 1:30,
  n.components = 3,
  seed.use = 100
)

bapd8.integrated@meta.data$sample <- factor('BAPd8_O2Level', levels = 'BAPd8_O2Level')
bapd8.integrated@misc$experiment <- list(
  experiment_name = 'BAPd8_O2Level',
  organism = 'hg',
  date_of_analysis = Sys.Date()
)

bapd8.integrated@misc$technical_info <- list(
  'R' = capture.output(devtools::session_info())
)

bapd8.integrated@misc$parameters <- list(
  gene_nomenclature = 'gene_name',
  discard_genes_expressed_in_fewer_cells_than = 10,
  keep_mitochondrial_genes = TRUE,
  variables_to_regress_out = 'nUMI',
  number_PCs = 30,
  tSNE_perplexity = 30,
  cluster_resolution = 0.5
)

bapd8.integrated@misc$parameters$filtering <- list(
  UMI_min = 100,
  UMI_max = Inf,
  genes_min = 200,
  genes_max = Inf
)


bapd8.integrated <- cerebroApp::addPercentMtRibo(
  bapd8.integrated,
  organism = 'hg',
  gene_nomenclature = 'name'
)



bapd8.integrated <- cerebroApp::getMostExpressedGenes(
  bapd8.integrated,
  groups = c('replicate', 'orig.ident', 'cluster' )
)

bapd8.integrated <- cerebroApp::getMarkerGenes(
  bapd8.integrated,
  organism = 'hg',
  groups = c('replicate', 'orig.ident', 'cluster' )
)

bapd8.integrated <- cerebroApp::getEnrichedPathways(
  bapd8.integrated,
  databases = c("GO_Biological_Process_2018", "GO_Cellular_Component_2018",
                "GO_Molecular_Function_2018", "KEGG_2016", "WikiPathways_2016", "Reactome_2016",
                "Panther_2016", "Human_Gene_Atlas", "Mouse_Gene_Atlas"),
  adj_p_cutoff = 0.05,
  max_terms = 100,
  URL_API = "http://amp.pharm.mssm.edu/Enrichr/enrich"
)



cerebroApp::exportFromSeurat(
  bapd8.integrated,
  assay = "RNA",
  slot = "data",
  file ="Cerebro_BAPd8_O2Level.crb",
  experiment_name = 'BAPd8_O2Level',
  organism = 'hg',
  groups = c("orig.ident", "replicate", "cluster"),
  cell_cycle = NULL,
  nUMI = 'nCount_RNA',
  nGene = 'nFeature_RNA',
  add_all_meta_data = TRUE,
  use_delayed_array = FALSE,
  verbose = FALSE
)

saveRDS(bapd8.integrated, 'bapd8.integrated.rds')


expdata <- as.matrix(GetAssayData(object = bapd8.integrated, slot = "data"))
df <- bapd8.integrated@meta.data
write.table(df, file="metadata.tsv", sep= "\t")
write.table(expdata, file="normalized-counts.tsv", sep= "\t")
expdata <- as.matrix(GetAssayData(object = bapd8.combined, slot = "counts"))
write.table(expdata, file="raw-counts.tsv", sep= "\t")

input <- read.csv("violin_plot-data.tsv", sep="\t", quote='', stringsAsFactors=TRUE, header=TRUE)
h<- ggplot(input, aes(x=info, y=Expression, fill=info))  + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_classic(base_size = 14) + theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5, face="bold")) + ggtitle("CGB3 + CGB5 + CGB7 + CGB8") + labs(y = "Expression Level", x ="")

