setwd("~/TutejaLab/lo-etal/counts-v4/Io-data_v3/")
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
library(reshape2)

# read counts and batch files
counts = 'assets/genes-only-selected-runs-counts-renamed.txt'
groupFile = 'assets/batch.txt'
coldata <- read.csv(groupFile, row.names=1, sep="\t", stringsAsFactors = TRUE)
cts <- as.matrix(read.csv(counts,sep="\t",row.names="gene.ids"))
# check files
head(cts)
colnames(cts)
# ensure columns match up between 2 files
all(rownames(coldata) %in% colnames(cts))
cts <- cts[, rownames(coldata)]
# read the projects for each sample
cov1 <- as.factor(coldata$BioProject)
# batch correct based on projects
adjusted_counts <- ComBat_seq(cts, batch=cov1, group=NULL)
# double check
all(rownames(coldata) %in% colnames(cts))
cts <- cts[, rownames(coldata)]
# create DESeq object
dds <- DESeqDataSetFromMatrix(countData = adjusted_counts,
                              colData = coldata,
                              design = ~ condition)
# vst
vsd <- vst(dds, blind=FALSE)
# remove low counts 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# run DESeq
dds <- DESeq(dds)
# check names in the dds
resultsNames(dds)
# create table for normalized counts
normalized_counts <- counts(dds, normalized=T) %>%
  data.frame() %>%
  rownames_to_column(var="gene")
# convert to tibble
normalized_counts <- normalized_counts %>%
  as_tibble()
file1 <- "assets/supp-table-s6-filtered.txt"
file2 <- "assets/roost-9w-human-amnion-genes.txt"
file3 <- "assets/monkey-amnion-genes.txt"
table.s6 <- read.csv(file1, sep="\t", header = TRUE, stringsAsFactors = FALSE)
roost.9w <- read.csv(file2, sep="\t", header = TRUE, stringsAsFactors = FALSE)
monkey.a <- read.csv(file3, sep="\t", header = TRUE, stringsAsFactors = FALSE)

genes.tableS6 = table.s6$genes
genes.monkey.a = monkey.a$genes
genes.roost.9w = roost.9w$genes

# biomart table
attributes=c('ensembl_gene_id_version','hgnc_symbol')
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# table S6
G_list<- getBM(attributes=attributes, filters="hgnc_symbol",values=genes.tableS6, mart=mart, uniqueRows=T)
head(G_list)
norm_selected <- normalized_counts %>%
  filter(gene %in% G_list$ensembl_gene_id_version)

norm.selected.long <- melt(norm_selected, id.vars="gene")
defile= "assets/selected-datasets-for-heatmap.txt"
decoldata <- read.csv(defile, sep="\t", stringsAsFactors = TRUE)
de_names <- decoldata$name
de.norm.selected.long <- norm.selected.long %>% filter (variable %in% de_names)
colnames(de.norm.selected.long) <- c("gene", "condition", "norm.expression")
de.norm.selected.long <- merge(de.norm.selected.long, G_list,
                               by.x = "gene",
                               by.y = "ensembl_gene_id_version",
                               all.x = TRUE,
                               all.y = FALSE)
head(de.norm.selected.long)
de.norm.selected.wide <- dcast(de.norm.selected.long, hgnc_symbol ~ condition, value.var="norm.expression")
head(de.norm.selected.wide)
write_delim(de.norm.selected.wide, file="table-S6_normalized-table.tsv", delim = "\t")
temp <- de.norm.selected.wide[,-1]
row.names(temp) <- de.norm.selected.wide[,1]
pheatmap(temp)
colnames(temp)
# clean labels
annotation <- data.frame(Condition = factor(1:42, labels = c("CytoTB11", "CytoTB7", "CytoTB7","CytoTB9","STB.gt70um","STB.gt70um","STB.gt70um","STB.lt40um","STB.lt40um","STB.lt40um","nCT.P10","nCT.P10","nTE.D3","nTE.D3","pBAP.D3","pBAP.D3","CT.EVT","CT.STB","CT.TSC","CT.EVT","CT.STB","CT.TSC","CT.EVT","CT.STB","CT.TSC","Amnion.female","Amnion.male","Amnion.female","Amnion.female","Amnion.female","Amnion.male","Amnion.female","Amnion.male","Amnion.male","Amnion.female","Amnion.male","Amnion.female","amnion.g18","amnion.g16","amnion.g16","amnion.g22","amnion.g22")))
annotation
heat_colors <- brewer.pal(9, "YlOrRd")
rownames(annotation) <- colnames(temp)
annotation
g <- pheatmap(temp,
              color = heat_colors,
              cluster_rows = T,
              show_rownames = T,
              annotation_col = annotation,
              border_color = NA,
              fontsize = 14,
              scale = "row",
              fontsize_row = 10,
              height = 20)
ggsave("heatmap-tables6.png", g, dpi=900, width = 16, height = 30)

# monkey genes heatmap
attributes=c('ensembl_gene_id_version','hgnc_symbol')
G_list<- getBM(attributes=attributes, filters="hgnc_symbol",values=genes.monkey.a, mart=mart, uniqueRows=T)
norm_selected <- normalized_counts %>%
  filter(gene %in% G_list$ensembl_gene_id_version)
norm.selected.long <- melt(norm_selected, id.vars="gene")
defile= "selected-datasets-for-heatmap.txt"
decoldata <- read.csv(defile, sep="\t", stringsAsFactors = TRUE)
de_names <- decoldata$name
de.norm.selected.long <- norm.selected.long %>% filter (variable %in% de_names)
colnames(de.norm.selected.long) <- c("gene", "condition", "norm.expression")
de.norm.selected.long <- merge(de.norm.selected.long, G_list,
                               by.x = "gene",
                               by.y = "ensembl_gene_id_version",
                               all.x = TRUE,
                               all.y = FALSE)
head(de.norm.selected.long)
de.norm.selected.wide <- dcast(de.norm.selected.long, hgnc_symbol ~ condition, value.var="norm.expression")
head(de.norm.selected.wide)
write_delim(de.norm.selected.wide, file="table-S6_normalized-table.tsv", delim = "\t")

temp <- de.norm.selected.wide[,-1]
row.names(temp) <- de.norm.selected.wide[,1]
pheatmap(temp)
colnames(temp)
annotation <- data.frame(Condition = factor(1:42, labels = c("CytoTB11", "CytoTB7", "CytoTB7","CytoTB9","STB.gt70um","STB.gt70um","STB.gt70um","STB.lt40um","STB.lt40um","STB.lt40um","nCT.P10","nCT.P10","nTE.D3","nTE.D3","pBAP.D3","pBAP.D3","CT.EVT","CT.STB","CT.TSC","CT.EVT","CT.STB","CT.TSC","CT.EVT","CT.STB","CT.TSC","Amnion.female","Amnion.male","Amnion.female","Amnion.female","Amnion.female","Amnion.male","Amnion.female","Amnion.male","Amnion.male","Amnion.female","Amnion.male","Amnion.female","amnion.g18","amnion.g16","amnion.g16","amnion.g22","amnion.g22")))
annotation
heat_colors <- brewer.pal(9, "YlOrRd")
rownames(annotation) <- colnames(temp)
annotation
g <- pheatmap(temp,
              color = heat_colors,
              cluster_rows = T,
              show_rownames = T,
              annotation_col = annotation,
              border_color = NA,
              fontsize = 14,
              scale = "row",
              fontsize_row = 10,
              height = 20)
ggsave("heatmap-monkey.png", g, dpi=900, width = 12, height = 8)

# Roost 9w human amnion genes heatmap
attributes=c('ensembl_gene_id_version','hgnc_symbol')
G_list<- getBM(attributes=attributes, filters="hgnc_symbol",values=genes.roost.9w, mart=mart, uniqueRows=T)
norm_selected <- normalized_counts %>%
  filter(gene %in% G_list$ensembl_gene_id_version)
norm.selected.long <- melt(norm_selected, id.vars="gene")
defile= "selected-datasets-for-heatmap.txt"
decoldata <- read.csv(defile, sep="\t", stringsAsFactors = TRUE)
de_names <- decoldata$name
de.norm.selected.long <- norm.selected.long %>% filter (variable %in% de_names)
colnames(de.norm.selected.long) <- c("gene", "condition", "norm.expression")
de.norm.selected.long <- merge(de.norm.selected.long, G_list,
                               by.x = "gene",
                               by.y = "ensembl_gene_id_version",
                               all.x = TRUE,
                               all.y = FALSE)
head(de.norm.selected.long)
de.norm.selected.wide <- dcast(de.norm.selected.long, hgnc_symbol ~ condition, value.var="norm.expression")
head(de.norm.selected.wide)
write_delim(de.norm.selected.wide, file="table-S6_normalized-table.tsv", delim = "\t")

temp <- de.norm.selected.wide[,-1]
row.names(temp) <- de.norm.selected.wide[,1]
pheatmap(temp)
colnames(temp)
annotation <- data.frame(Condition = factor(1:42, labels = c("CytoTB11", "CytoTB7", "CytoTB7","CytoTB9","STB.gt70um","STB.gt70um","STB.gt70um","STB.lt40um","STB.lt40um","STB.lt40um","nCT.P10","nCT.P10","nTE.D3","nTE.D3","pBAP.D3","pBAP.D3","CT.EVT","CT.STB","CT.TSC","CT.EVT","CT.STB","CT.TSC","CT.EVT","CT.STB","CT.TSC","Amnion.female","Amnion.male","Amnion.female","Amnion.female","Amnion.female","Amnion.male","Amnion.female","Amnion.male","Amnion.male","Amnion.female","Amnion.male","Amnion.female","amnion.g18","amnion.g16","amnion.g16","amnion.g22","amnion.g22")))
annotation
heat_colors <- brewer.pal(9, "YlOrRd")
rownames(annotation) <- colnames(temp)
annotation
g <- pheatmap(temp,
              color = heat_colors,
              cluster_rows = T,
              show_rownames = T,
              annotation_col = annotation,
              border_color = NA,
              fontsize = 14,
              scale = "row",
              fontsize_row = 10,
              height = 20)
ggsave("heatmap-roost9w.png", g, dpi=900, width = 16, height = 30)
