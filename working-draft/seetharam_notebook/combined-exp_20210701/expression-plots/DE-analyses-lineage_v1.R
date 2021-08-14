setwd("~/TutejaLab/combined-exp_20210701/expression-plots")
# load packages
library(sva)
library(tidyverse)
library(DESeq2)
library(vsn)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(reshape2)
require(biomaRt)
library(EnhancedVolcano)
library(TissueEnrich)
# read counts and batch files
counts = 'lineage-counts-genes.txt'
groupFile = 'batch.txt'
coldata <- read.csv(groupFile, row.names=1, sep="\t", stringsAsFactors = TRUE)
cts <- as.matrix(read.csv(counts,sep="\t",row.names="geneids"))
# inspect
head(cts)
colnames(cts)
# check to make sure all rows of counts and batch have the names
all(rownames(coldata) %in% colnames(cts))
# and arrange in the same order
cts <- cts[, rownames(coldata)]
# read in the project info for libraries

# create DESeq2 object with design as condition
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ group)
# vst for pca (not needed here actually)
vsd <- vst(dds, blind=FALSE)
# keep genes with counts >= 10
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, maxit=500)
# inspect
dds
normalized_counts <- counts(dds, normalized=T) %>%
  data.frame() %>%
  rownames_to_column(var="gene")
normalized_counts$gene <- str_replace_all(normalized_counts$gene, pattern=" ", repl="")
# convert to tibble
normalized_counts <- normalized_counts %>%
  as_tibble()
file1 <- "cluster-6-and-8-SCT.txt"
file2 <- "cluster-6-only-SCT.txt"
file3 <- "cluster-8-only-SCT.txt"
c68.sct <- read.csv(file1, sep="\t", header = TRUE, stringsAsFactors = FALSE)
c06.sct <- read.csv(file2, sep="\t", header = TRUE, stringsAsFactors = FALSE)
c08.sct <- read.csv(file3, sep="\t", header = TRUE, stringsAsFactors = FALSE)
genes.c68.sct  = c68.sct$genes
genes.c06.sct  = c06.sct$genes
genes.c08.sct  = c08.sct$genes
# biomart table
attributes=c('ensembl_gene_id_version','hgnc_symbol')
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# genes.c68.sct
G_list<- getBM(attributes=attributes, filters="hgnc_symbol",values=genes.c68.sct, mart=mart, uniqueRows=T)
head(G_list)
norm_selected <- normalized_counts %>%
  dplyr::filter(gene %in% G_list$ensembl_gene_id_version)

norm.selected.long <- melt(norm_selected, id.vars="gene")

norm.selected.long <- merge(norm.selected.long, G_list,
                               by.x = "gene",
                               by.y = "ensembl_gene_id_version",
                               all.x = TRUE,
                               all.y = FALSE)
head(norm.selected.long)
norm.selected.long$days <- NA
norm.selected.long$days[which(str_detect(norm.selected.long$variable, "D10"))] <- "D10"
norm.selected.long$days[which(str_detect(norm.selected.long$variable, "D12"))] <- "D12"
norm.selected.long$days[which(str_detect(norm.selected.long$variable, "D8"))] <- "D8"
norm.selected.long$cell <- NA
norm.selected.long$cell[which(str_detect(norm.selected.long$variable, "_L"))] <- "STB"
norm.selected.long$cell[which(str_detect(norm.selected.long$variable, "_S"))] <- "CTB"
norm.selected.long$cell[which(str_detect(norm.selected.long$variable, "_EVT"))] <- "EVT"
norm.selected.long$cell[which(is.na(norm.selected.long$cell))] <- "STB"
norm.selected.long$type <- paste(norm.selected.long$days, ".", norm.selected.long$cell, sep = "")
goi <- unique(norm.selected.long$hgnc_symbol)
pdf("genes_in_both_cluster_6-and-8.pdf")
for (i in goi) {
  p <- ggplot(filter(norm.selected.long, hgnc_symbol== i), aes(type, value, fill=cell)) + 
    geom_violin() + geom_jitter(shape=16, position=position_jitter(0.1)) +
    ggtitle(i) +
    xlab("cell types") +
    ylab("normalized expression") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
          axis.text.y = element_text(size = 12),
          plot.margin = margin(10, 10, 10, 100), legend.position = "none",
          plot.title = element_text(color = "black", size=18, face="bold.italic"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color="black", size=14))
  print(p)
}
dev.off()
norm.selected.wide <- dcast(norm.selected.long, hgnc_symbol ~ variable, value.var="value")
write_delim(norm.selected.wide, file="genes_in_both_cluster_6-and-8-SCT-normalized.tsv", delim = "\t")
hmap <- norm.selected.wide[,-1]
row.names(hmap) <- norm.selected.wide[,1]
annotation <- read.csv("annotation.tsv", sep="\t", stringsAsFactors = TRUE, row.names="names")
heat_colors <- brewer.pal(9, "YlOrRd")
all(rownames(annotation) %in% colnames(hmap))
hmap <- hmap[, rownames(annotation)]
g <- pheatmap(hmap,
              color = heat_colors,
              cluster_rows = T,
              show_rownames = T,
              annotation_col = annotation,
              border_color = NA,
              fontsize = 14,
              scale = "row",
              fontsize_row = 10,
              height = 20)
ggsave("genes_in_both_cluster_6-and-8_heatmap.png", g, dpi=900, width = 28, height = 12)
# genes.c06.sct
G_list<- getBM(attributes=attributes, filters="hgnc_symbol",values=genes.c06.sct, mart=mart, uniqueRows=T)
head(G_list)
norm_selected <- normalized_counts %>%
  dplyr::filter(gene %in% G_list$ensembl_gene_id_version)

norm.selected.long <- melt(norm_selected, id.vars="gene")

norm.selected.long <- merge(norm.selected.long, G_list,
                            by.x = "gene",
                            by.y = "ensembl_gene_id_version",
                            all.x = TRUE,
                            all.y = FALSE)
head(norm.selected.long)
norm.selected.long$days <- NA
norm.selected.long$days[which(str_detect(norm.selected.long$variable, "D10"))] <- "D10"
norm.selected.long$days[which(str_detect(norm.selected.long$variable, "D12"))] <- "D12"
norm.selected.long$days[which(str_detect(norm.selected.long$variable, "D8"))] <- "D8"
norm.selected.long$cell <- NA
norm.selected.long$cell[which(str_detect(norm.selected.long$variable, "_L"))] <- "STB"
norm.selected.long$cell[which(str_detect(norm.selected.long$variable, "_S"))] <- "CTB"
norm.selected.long$cell[which(str_detect(norm.selected.long$variable, "_EVT"))] <- "EVT"
norm.selected.long$cell[which(is.na(norm.selected.long$cell))] <- "STB"
norm.selected.long$type <- paste(norm.selected.long$days, ".", norm.selected.long$cell, sep = "")
goi <- unique(norm.selected.long$hgnc_symbol)
pdf("genes_unique_to_cluster_6-SCT.pdf")
for (i in goi) {
  p <- ggplot(filter(norm.selected.long, hgnc_symbol== i), aes(type, value, fill=cell)) + 
    geom_violin() + geom_jitter(shape=16, position=position_jitter(0.1)) +
    ggtitle(i) +
    xlab("cell types") +
    ylab("normalized expression") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
          axis.text.y = element_text(size = 12),
          plot.margin = margin(10, 10, 10, 100), legend.position = "none",
          plot.title = element_text(color = "black", size=18, face="bold.italic"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color="black", size=14))
  print(p)
}
dev.off()
norm.selected.wide <- dcast(norm.selected.long, hgnc_symbol ~ variable, value.var="value")
write_delim(norm.selected.wide, file="genes_unique_to_cluster_6-SCT-normalized.tsv", delim = "\t")
hmap <- norm.selected.wide[,-1]
row.names(hmap) <- norm.selected.wide[,1]
annotation <- read.csv("annotation.tsv", sep="\t", stringsAsFactors = TRUE, row.names="names")
heat_colors <- brewer.pal(9, "YlOrRd")
all(rownames(annotation) %in% colnames(hmap))
hmap <- hmap[, rownames(annotation)]
g <- pheatmap(hmap,
              color = heat_colors,
              cluster_rows = T,
              show_rownames = T,
              annotation_col = annotation,
              border_color = NA,
              fontsize = 14,
              scale = "row",
              fontsize_row = 10,
              height = 20)
ggsave("genes_unique_to_cluster_6-SCT_heatmap.png", g, dpi=900, width = 28, height = 12)

# genes.c08.sct
G_list<- getBM(attributes=attributes, filters="hgnc_symbol",values=genes.c08.sct, mart=mart, uniqueRows=T)
head(G_list)
norm_selected <- normalized_counts %>%
  dplyr::filter(gene %in% G_list$ensembl_gene_id_version)
norm.selected.long <- melt(norm_selected, id.vars="gene")
norm.selected.long <- merge(norm.selected.long, G_list,
                            by.x = "gene",
                            by.y = "ensembl_gene_id_version",
                            all.x = TRUE,
                            all.y = FALSE)
head(norm.selected.long)
norm.selected.long$days <- NA
norm.selected.long$days[which(str_detect(norm.selected.long$variable, "D10"))] <- "D10"
norm.selected.long$days[which(str_detect(norm.selected.long$variable, "D12"))] <- "D12"
norm.selected.long$days[which(str_detect(norm.selected.long$variable, "D8"))] <- "D8"
norm.selected.long$cell <- NA
norm.selected.long$cell[which(str_detect(norm.selected.long$variable, "_L"))] <- "STB"
norm.selected.long$cell[which(str_detect(norm.selected.long$variable, "_S"))] <- "CTB"
norm.selected.long$cell[which(str_detect(norm.selected.long$variable, "_EVT"))] <- "EVT"
norm.selected.long$cell[which(is.na(norm.selected.long$cell))] <- "STB"
norm.selected.long$type <- paste(norm.selected.long$days, ".", norm.selected.long$cell, sep = "")
goi <- unique(norm.selected.long$hgnc_symbol)
pdf("genes_unique_to_cluster_8-SCT.pdf")
for (i in goi) {
  p <- ggplot(filter(norm.selected.long, hgnc_symbol== i), aes(type, value, fill=cell)) + 
    geom_violin() + geom_jitter(shape=16, position=position_jitter(0.1)) +
    ggtitle(i) +
    xlab("cell types") +
    ylab("normalized expression") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
          axis.text.y = element_text(size = 12),
          plot.margin = margin(10, 10, 10, 100), legend.position = "none",
          plot.title = element_text(color = "black", size=18, face="bold.italic"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color="black", size=14))
  print(p)
}
dev.off()
norm.selected.wide <- dcast(norm.selected.long, hgnc_symbol ~ variable, value.var="value")
write_delim(norm.selected.wide, file="genes_unique_to_cluster_8-SCT-normalized.tsv", delim = "\t")
hmap <- norm.selected.wide[,-1]
row.names(hmap) <- norm.selected.wide[,1]
annotation <- read.csv("annotation.tsv", sep="\t", stringsAsFactors = TRUE, row.names="names")
heat_colors <- brewer.pal(9, "YlOrRd")
all(rownames(annotation) %in% colnames(hmap))
hmap <- hmap[, rownames(annotation)]
g <- pheatmap(hmap,
              color = heat_colors,
              cluster_rows = T,
              show_rownames = T,
              annotation_col = annotation,
              border_color = NA,
              fontsize = 14,
              scale = "row",
              fontsize_row = 10,
              height = 20)
ggsave("genes_unique_to_cluster_8-SCT_heatmap.png", g, dpi=900, width = 28, height = 12)
