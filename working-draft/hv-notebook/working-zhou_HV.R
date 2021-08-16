setwd("Z:/hhvu/Project5_amnion/amnion.vs.other_RNASeq/seetharam_notebook/zhou-analyses-v1/")
# load libraries
library(tidyverse)
library(RColorBrewer)
# read-in the counts matrix
# counts were downloaded from NCBI
# link: ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109555/suppl/GSE109555%5FAll%5FEmbryo%5FTPM%2Etxt%2Egz
cts <- as.matrix(read.csv("GSE109555_All_Embryo_TPM.txt", sep="\t"))
#,row.names="genes" threw me an error because the count file I got from the link above doesn't have a "genes" column
#Error in data[[rowvar]] : 
# attempt to select less than one element in get1index
#It still got read in fine (genes are still row names) because of the default behavior:
#If there is a header and the first row contains one fewer field than the number of columns, the first column in the input is used for the row names. 

colnames(cts)
# create annotation information for the cells
df <- as.data.frame(colnames(cts))
colnames(df) <- "cells"
df$days <- NA
df$days[which(str_detect(df$cells, "D14"))] <- "D14"
df$days[which(str_detect(df$cells, "D12"))] <- "D12"
df$days[which(str_detect(df$cells, "D10"))] <- "D10"
df$days[which(str_detect(df$cells, "D8"))] <- "D8"
df$days[which(str_detect(df$cells, "D6"))] <- "D6"
df$days <- as.factor(df$days)
df <- df %>% remove_rownames %>% column_to_rownames(var="cells")

# SCT genes obtained from snRNAseq analyses
file1 <- "cluster-6-and-8-SCT.txt"
file2 <- "cluster-8-only-SCT.txt"
file3 <- "cluster-6-only-SCT.txt"
file4 <- "cluster8-SCT.txt"
file5 <- "cluster6-SCT.txt"
# read them
c6n8.sct <- read.csv(file1, sep="\t", header = TRUE, stringsAsFactors = FALSE)
c8ex.sct <- read.csv(file2, sep="\t", header = TRUE, stringsAsFactors = FALSE)
c6ex.sct <- read.csv(file3, sep="\t", header = TRUE, stringsAsFactors = FALSE)
c8in.sct <- read.csv(file4, sep="\t", header = TRUE, stringsAsFactors = FALSE)
c6in.sct <- read.csv(file5, sep="\t", header = TRUE, stringsAsFactors = FALSE)

# save them as meaningful name
genes.c6n8.sct  = c6n8.sct$genes
genes.c8ex.sct  = c8ex.sct$genes
genes.c6ex.sct  = c6ex.sct$genes
genes.c8in.sct  = c8in.sct$genes
genes.c6in.sct  = c6in.sct$genes
# subset the counts table (to keep only sct genes for plotting heatmap)
genes.c6n8.sct.cts <- cts[rownames(cts) %in% genes.c6n8.sct, ]
genes.c8ex.sct.cts <- cts[rownames(cts) %in% genes.c8ex.sct, ]
genes.c6ex.sct.cts <- cts[rownames(cts) %in% genes.c6ex.sct, ]
genes.c8in.sct.cts <- cts[rownames(cts) %in% genes.c8in.sct, ]
genes.c6in.sct.cts <- cts[rownames(cts) %in% genes.c6in.sct, ]

# function for heatmap
# subset dataset is grouped based on days and the average is computed to plot as single raster in heatmap
mydatahm <- function(mycts, name) {
  mycts.t <- as.data.frame(t(mycts))
  mycts.t$days <- NA
  mycts.t$days[which(str_detect(rownames(mycts.t), "D14"))] <- "D14"
  mycts.t$days[which(str_detect(rownames(mycts.t), "D12"))] <- "D12"
  mycts.t$days[which(str_detect(rownames(mycts.t), "D10"))] <- "D10"
  mycts.t$days[which(str_detect(rownames(mycts.t), "D8"))] <- "D8"
  mycts.t$days[which(str_detect(rownames(mycts.t), "D6"))] <- "D6"
  mycts.mean <- aggregate(. ~ days, mycts.t[-2], mean) #why do we need to remove column 2 here?
  ordered <- as.data.frame(t(mycts.mean %>% remove_rownames %>% column_to_rownames(var="days")))
  ordered <- ordered[c("D6", "D8", "D10", "D12", "D14")]
  hmap <- as.matrix(ordered)
  g <- pheatmap(hmap,
                color = heat_colors,
                cluster_rows = T,
                cluster_cols = F,
                show_rownames = T,
                border_color = NA,
                fontsize = 12,
                scale = "row",
                fontsize_row = 10)
  ggsave(paste("heatmap_", name, ".png", sep = ""), g, dpi=900, width = 8, height = 10)
}

# plot heatmaps
mydatahm(genes.c6n8.sct.cts, "cluster_both-6-and-8-SCT")
mydatahm(genes.c8ex.sct.cts, "cluster_8-exclusive-SCT")
mydatahm(genes.c6ex.sct.cts, "cluster_6-exclusive-SCT")
mydatahm(genes.c8in.sct.cts, "cluster_8-inclusive-SCT")
mydatahm(genes.c6in.sct.cts, "cluster_6-inclusive-SCT")

#just testing
#mycts <- genes.c6n8.sct.cts
#mycts.t <- as.data.frame(t(mycts))
#mycts.t$days <- NA
#mycts.t$days[which(str_detect(rownames(mycts.t), "D14"))] <- "D14"
#mycts.t$days[which(str_detect(rownames(mycts.t), "D12"))] <- "D12"
#mycts.t$days[which(str_detect(rownames(mycts.t), "D10"))] <- "D10"
#mycts.t$days[which(str_detect(rownames(mycts.t), "D8"))] <- "D8"
#mycts.t$days[which(str_detect(rownames(mycts.t), "D6"))] <- "D6"
#mycts.mean <- aggregate(. ~ days, mycts.t[-2], mean)
#ordered <- as.data.frame(t(mycts.mean %>% remove_rownames %>% column_to_rownames(var="days")))
#ordered <- ordered[c("D6", "D8", "D10", "D12", "D14")]
#hmap <- as.matrix(ordered)
