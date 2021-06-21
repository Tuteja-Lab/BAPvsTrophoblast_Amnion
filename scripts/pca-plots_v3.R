setwd("~/TutejaLab/lo-etal/counts-v4/Io-data_v3/")
library(sva)
library(tidyverse)
library(DESeq2)
library(vsn)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
# read counts and metadata
counts = 'assets/genes-only-selected-runs-counts-renamed.txt'
groupFile = 'assets/batch.txt'
coldata <- read.csv(groupFile, row.names=1, sep="\t", stringsAsFactors = TRUE)
cts <- as.matrix(read.csv(counts,sep="\t",row.names="gene.ids"))
# check
head(cts)
colnames(cts)
# make sure the columns match
all(rownames(coldata) %in% colnames(cts))
cts <- cts[, rownames(coldata)]
# create DEseq object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
# vst
vsd <- vst(dds, blind=FALSE)
# PCA plot before batch correction
pcaData <- plotPCA(vsd, intgroup=c("condition", "BioProject"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
B1 <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=BioProject)) +
  geom_label_repel(aes(PC1, PC2, label = condition),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.5, "lines"),
                   segment.color = 'grey50', max.overlaps = Inf) +
  scale_shape_manual(values=1:8) + theme(legend.title = element_blank()) +
  geom_point(size=2, stroke = 2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
ggsave("pca-before.png", dpi=900, width = 24, height = 16)


# batch correction
# read the bioproject column
cov1 <- as.factor(coldata$BioProject)
# run batch correction
adjusted_counts <- ComBat_seq(cts, batch=cov1, group=NULL)
all(rownames(coldata) %in% colnames(cts))
cts <- cts[, rownames(coldata)]
# create deseq object
dds <- DESeqDataSetFromMatrix(countData = adjusted_counts,
                              colData = coldata,
                              design = ~ condition)
# vst
vsd <- vst(dds, blind=FALSE)
# PCA plot after
pcaData <- plotPCA(vsd, intgroup=c("condition", "BioProject"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
A1 <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=BioProject)) +
  geom_label_repel(aes(PC1, PC2, label = condition),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.5, "lines"),
                   segment.color = 'grey50', max.overlaps = Inf) +
  scale_shape_manual(values=1:8) + theme(legend.title = element_blank()) +
  geom_point(size=2, stroke = 2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
ggsave("pca-after.png", dpi=900, width = 24, height = 16)






