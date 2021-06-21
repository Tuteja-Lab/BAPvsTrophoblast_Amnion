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
library(EnhancedVolcano)
library(TissueEnrich)
# read counts and batch files
counts = 'genes-only-selected-runs-counts-renamed.txt'
groupFile = 'batch.txt'
coldata <- read.csv(groupFile, row.names=1, sep="\t", stringsAsFactors = TRUE)
cts <- as.matrix(read.csv(counts,sep="\t",row.names="gene.ids"))
# inspect
head(cts)
colnames(cts)
# check to make sure all rows of counts and batch have the names 
all(rownames(coldata) %in% colnames(cts))
# and arrange in the same order
cts <- cts[, rownames(coldata)]
# read in the project info for libraries
cov1 <- as.factor(coldata$BioProject)
# batch correction
adjusted_counts <- ComBat_seq(cts, batch=cov1, group=NULL)
# check again
all(rownames(coldata) %in% colnames(cts))
cts <- cts[, rownames(coldata)]
# create DESeq2 object with design as condition
dds <- DESeqDataSetFromMatrix(countData = adjusted_counts,
                              colData = coldata,
                              design = ~ condition)
# vst for pca (not needed here actually)
vsd <- vst(dds, blind=FALSE)
# keep genes with counts >= 10
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
# inspect
dds
# see the names for conditions
resultsNames(dds)
# run DE for comparisons we are interested
res.nCTvsBAP <- results(dds, contrast=c("condition","nCT_P10.Naive_H9_hESCs","pBAP_D3.Primed_H9_hESCs"))
res.nTEvsSTB <- results(dds, contrast=c("condition","nTE_D3.Naive_H9_hESCs","hESC_H1_STB_gt70um_D8_BAP"))
res.nCTvsSTB <- results(dds, contrast=c("condition","nCT_P10.Naive_H9_hESCs","hESC_H1_STB_gt70um_D8_BAP"))
res.nTEvsBAP <- results(dds, contrast=c("condition","nTE_D3.Naive_H9_hESCs","pBAP_D3.Primed_H9_hESCs"))
res.nTEvsL40 <- results(dds, contrast=c("condition","nTE_D3.Naive_H9_hESCs","hESC_H1_STB_lt40um_D8_BAP"))
res.nCTvsL40 <- results(dds, contrast=c("condition","nCT_P10.Naive_H9_hESCs","hESC_H1_STB_lt40um_D8_BAP"))
res.BAPvsL40 <- results(dds, contrast=c("condition","pBAP_D3.Primed_H9_hESCs","hESC_H1_STB_lt40um_D8_BAP"))
res.BAPvsSTB <- results(dds, contrast=c("condition","pBAP_D3.Primed_H9_hESCs","hESC_H1_STB_gt70um_D8_BAP"))
res.STBvsL40 <- results(dds, contrast=c("condition","hESC_H1_STB_gt70um_D8_BAP","hESC_H1_STB_lt40um_D8_BAP"))
# check the number of significant genes (regardless of FC)
table(res.nCTvsBAP$padj<0.05)
table(res.nCTvsSTB$padj<0.05)
table(res.nTEvsSTB$padj<0.05)
table(res.nTEvsBAP$padj<0.05)
table(res.nTEvsL40$padj<0.05)
table(res.nCTvsL40$padj<0.05)
table(res.BAPvsL40$padj<0.05)
table(res.STBvsL40$padj<0.05)
table(res.BAPvsSTB$padj<0.05)
# order them based on p_adj
res.nCTvsBAP <- res.nCTvsBAP[order(res.nCTvsBAP$padj), ]
res.nCTvsSTB <- res.nCTvsSTB[order(res.nCTvsSTB$padj), ]
res.nTEvsSTB <- res.nTEvsSTB[order(res.nTEvsSTB$padj), ]
res.nTEvsBAP <- res.nTEvsBAP[order(res.nCTvsBAP$padj), ]
res.nTEvsL40 <- res.nTEvsL40[order(res.nTEvsL40$padj), ]
res.nCTvsL40 <- res.nCTvsL40[order(res.nCTvsL40$padj), ]
res.BAPvsL40 <- res.BAPvsL40[order(res.BAPvsL40$padj), ]
res.STBvsL40 <- res.STBvsL40[order(res.STBvsL40$padj), ]
res.BAPvsSTB <- res.BAPvsSTB[order(res.BAPvsSTB$padj), ]
# add in the normalized coutns for the results file
res.nCTvsBAPdata <- merge(as.data.frame(res.nCTvsBAP), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
res.nCTvsSTBdata <- merge(as.data.frame(res.nCTvsSTB), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
res.nTEvsSTBdata <- merge(as.data.frame(res.nTEvsSTB), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
res.nTEvsBAPdata <- merge(as.data.frame(res.nTEvsBAP), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
res.nTEvsL40data <- merge(as.data.frame(res.nTEvsL40), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
res.nCTvsL40data <- merge(as.data.frame(res.nCTvsL40), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
res.BAPvsL40data <- merge(as.data.frame(res.BAPvsL40), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
res.STBvsL40data <- merge(as.data.frame(res.STBvsL40), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
res.BAPvsSTBdata <- merge(as.data.frame(res.BAPvsSTB), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
# fix the header
names(res.nCTvsBAPdata)[1] <- "Gene"
names(res.nCTvsSTBdata)[1] <- "Gene"
names(res.nTEvsSTBdata)[1] <- "Gene"
names(res.nTEvsBAPdata)[1] <- "Gene"
names(res.nTEvsL40data)[1] <- "Gene"
names(res.nCTvsL40data)[1] <- "Gene"
names(res.BAPvsL40data)[1] <- "Gene"
names(res.STBvsL40data)[1] <- "Gene"
names(res.BAPvsSTBdata)[1] <- "Gene"
# save the results
write_delim(res.nCTvsBAPdata, file="DESeq2results-nCTvsBAP_fc.tsv", delim = "\t")
write_delim(res.nCTvsSTBdata, file="DESeq2results-nCTvsSTB_fc.tsv", delim = "\t")
write_delim(res.nTEvsSTBdata, file="DESeq2results-nTEvsSTB_fc.tsv", delim = "\t")
write_delim(res.nTEvsBAPdata, file="DESeq2results-nTEvsBAP_fc.tsv", delim = "\t")
write_delim(res.nTEvsL40data, file="DESeq2results-nTEvsL40_fc.tsv", delim = "\t")
write_delim(res.nCTvsL40data, file="DESeq2results-nCTvsL40_fc.tsv", delim = "\t")
write_delim(res.BAPvsL40data, file="DESeq2results-BAPvsL40_fc.tsv", delim = "\t")
write_delim(res.STBvsL40data, file="DESeq2results-STBvsL40_fc.tsv", delim = "\t")
write_delim(res.BAPvsSTBdata, file="DESeq2results-BAPvsSTB_fc.tsv", delim = "\t")
# for PCE, create list of genes (up and down)
# for nCTvsBAP
nCTvsBAP.up <- res.nCTvsBAPdata %>%
  filter(log2FoldChange >= 1) %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene)
nCTvsBAP.down <- res.nCTvsBAPdata %>%
  filter(log2FoldChange <= -1) %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene)
# for nCTvsSTB
nCTvsSTB.up <- res.nCTvsSTBdata %>%
  filter(log2FoldChange >= 1) %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene)
nCTvsSTB.down <- res.nCTvsSTBdata %>%
  filter(log2FoldChange <= -1) %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene)
# for nTEvsSTB
nTEvsSTB.up <- res.nTEvsSTBdata %>%
  filter(log2FoldChange >= 1) %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene)
nTEvsSTB.down <- res.nTEvsSTBdata %>%
  filter(log2FoldChange <= -1) %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene)
# for nTEvsBAP
nTEvsBAP.up <- res.nTEvsBAPdata %>%
  filter(log2FoldChange >= 1) %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene)
nTEvsBAP.down <- res.nTEvsBAPdata %>%
  filter(log2FoldChange <= -1) %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene)
# for nTEvsL40
nTEvsL40.up <- res.nTEvsL40data %>%
  filter(log2FoldChange >= 1) %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene)
nTEvsL40.down <- res.nTEvsL40data %>%
  filter(log2FoldChange <= -1) %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene)
# for nCTvsL40
nCTvsL40.up <- res.nCTvsL40data %>%
  filter(log2FoldChange >= 1) %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene)
nCTvsL40.down <- res.nCTvsL40data %>%
  filter(log2FoldChange <= -1) %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene)
# for BAPvsL40
BAPvsL40.up <- res.BAPvsL40data %>%
  filter(log2FoldChange >= 1) %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene)
BAPvsL40.down <- res.BAPvsL40data %>%
  filter(log2FoldChange <= -1) %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene)
# for STBvsL40
STBvsL40.up <- res.STBvsL40data %>%
  filter(log2FoldChange >= 1) %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene)
STBvsL40.down <- res.STBvsL40data %>%
  filter(log2FoldChange <= -1) %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene)
# for BAPvsSTB
BAPvsSTB.up <- res.BAPvsSTBdata %>%
  filter(log2FoldChange >= 1) %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene)
BAPvsSTB.down <- res.BAPvsSTBdata %>%
  filter(log2FoldChange <= -1) %>%
  filter(padj <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene)
# for converting ensembl gene id and version -> ensembl gene ids, get the 
# details from BioMart
yay
# remove duplicates, if any
isDup <- duplicated(annot$ensembl_gene_id)
dup <- annot$ensembl_gene_id[isDup]
annot <- annot[annot$ensembl_gene_id%in%dup,]
# add in the ids for the PCE up/down list files
nCTvsBAP.up.new <- annot[annot$ensembl_gene_id_version %in% nCTvsBAP.up$Gene,]
nCTvsBAP.down.new <- annot[annot$ensembl_gene_id_version %in% nCTvsBAP.down$Gene,]
nCTvsSTB.up.new <- annot[annot$ensembl_gene_id_version %in% nCTvsSTB.up$Gene,]
nCTvsSTB.down.new <- annot[annot$ensembl_gene_id_version %in% nCTvsSTB.down$Gene,]
nTEvsSTB.up.new <- annot[annot$ensembl_gene_id_version %in% nTEvsSTB.up$Gene,]
nTEvsSTB.down.new <- annot[annot$ensembl_gene_id_version %in% nTEvsSTB.down$Gene,]
nTEvsBAP.up.new <- annot[annot$ensembl_gene_id_version %in% nTEvsBAP.up$Gene,]
nTEvsBAP.down.new <- annot[annot$ensembl_gene_id_version %in% nTEvsBAP.down$Gene,]
nTEvsL40.up.new <- annot[annot$ensembl_gene_id_version %in% nTEvsL40.up$Gene,]
nTEvsL40.down.new <- annot[annot$ensembl_gene_id_version %in% nTEvsL40.down$Gene,]
nCTvsL40.up.new <- annot[annot$ensembl_gene_id_version %in% nCTvsL40.up$Gene,]
nCTvsL40.down.new <- annot[annot$ensembl_gene_id_version %in% nCTvsL40.down$Gene,]
BAPvsL40.up.new <- annot[annot$ensembl_gene_id_version %in% BAPvsL40.up$Gene,]
BAPvsL40.down.new <- annot[annot$ensembl_gene_id_version %in% BAPvsL40.down$Gene,]
STBvsL40.up.new <- annot[annot$ensembl_gene_id_version %in% STBvsL40.up$Gene,]
STBvsL40.down.new <- annot[annot$ensembl_gene_id_version %in% STBvsL40.down$Gene,]
BAPvsSTB.up.new <- annot[annot$ensembl_gene_id_version %in% BAPvsSTB.up$Gene,]
BAPvsSTB.down.new <- annot[annot$ensembl_gene_id_version %in% BAPvsSTB.down$Gene,]
# just retain the ids and nothing else
# since PCE natively uses ensembl ids, we just used that here
nCTvsBAP.up.pce <- nCTvsBAP.up.new$ensembl_gene_id
nCTvsBAP.down.pce <- nCTvsBAP.down.new$ensembl_gene_id
nCTvsSTB.up.pce <- nCTvsSTB.up.new$ensembl_gene_id
nCTvsSTB.down.pce <- nCTvsSTB.down.new$ensembl_gene_id
nTEvsSTB.up.pce <- nTEvsSTB.up.new$ensembl_gene_id
nTEvsSTB.down.pce <- nTEvsSTB.down.new$ensembl_gene_id
nTEvsBAP.up.pce <- nTEvsBAP.up.new$ensembl_gene_id
nTEvsBAP.down.pce <- nTEvsBAP.down.new$ensembl_gene_id
nTEvsL40.up.pce <- nTEvsL40.up.new$ensembl_gene_id
nTEvsL40.down.pce <- nTEvsL40.down.new$ensembl_gene_id
nCTvsL40.up.pce <- nCTvsL40.up.new$ensembl_gene_id
nCTvsL40.down.pce <- nCTvsL40.down.new$ensembl_gene_id
BAPvsL40.up.pce <- BAPvsL40.up.new$ensembl_gene_id
BAPvsL40.down.pce <- BAPvsL40.down.new$ensembl_gene_id
STBvsL40.up.pce <- STBvsL40.up.new$ensembl_gene_id
STBvsL40.down.pce <- STBvsL40.down.new$ensembl_gene_id
BAPvsSTB.up.pce <- BAPvsSTB.up.new$ensembl_gene_id
BAPvsSTB.down.pce <- BAPvsSTB.down.new$ensembl_gene_id
# load the PCE data
l <- load(file = "~/TutejaLab/PlacentaEnrich/combine-test-expression1.Rdata")
humanGeneMapping <- dataset$GRCH38$humanGeneMapping
d <- dataset$PlacentaDeciduaBloodData
data <- d$expressionData
cellDetails <- d$cellDetails
# create a run PCE function
runpce <- function(inputgenelist, title) {
  inputGenes<-toupper(inputgenelist)
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
    ggtitle(title) + ylab("-log10 p-value")
}
# run PCE and save results as png plots
runpce(nCTvsBAP.up.pce, "up regulated in nCT (compared to BAP.d3)"                     )
ggsave("nCTvsBAP.up.pce.png", dpi=900, width = 12, height = 10)
runpce(nCTvsBAP.up.pce , "up regulated in nCT (compared to BAP.d3)")
ggsave("nCTvsBAP.up.pce.png", dpi=900, width = 12, height = 10)
runpce(nCTvsBAP.down.pce , "up regulated in BAP.d3 (compared to nCT)")
ggsave("nCTvsBAP.down.pce.png", dpi=900, width = 12, height = 10)
runpce(nCTvsSTB.up.pce , "up regulated in nCT (compared to Yabe.BAP>70um)")
ggsave("nCTvsSTB.up.pce.png", dpi=900, width = 12, height = 10)
runpce(nCTvsSTB.down.pce , "up regulated in Yabe.BAP>70um (compared to nCT)")
ggsave("nCTvsSTB.down.pce.png", dpi=900, width = 12, height = 10)
runpce(nTEvsSTB.up.pce , "up regulated in nTE (compared to Yabe.BAP>70um)")
ggsave("nTEvsSTB.up.pce.png", dpi=900, width = 12, height = 10)
runpce(nTEvsSTB.down.pce , "up regulated in Yabe.BAP>70um (compared to nTE)")
ggsave("nTEvsSTB.down.pce.png", dpi=900, width = 12, height = 10)
runpce(nTEvsBAP.up.pce , "up regulated in nTE (compared to BAP.d3)")
ggsave("nTEvsBAP.up.pce.png", dpi=900, width = 12, height = 10)
runpce(nTEvsBAP.down.pce , "up regulated in BAP.d3 (compared to nTE)")
ggsave("nTEvsBAP.down.pce.png", dpi=900, width = 12, height = 10)
runpce(nTEvsL40.up.pce , "up regulated in nTE (compared to Yabe.BAP<40um)")
ggsave("nTEvsL40.up.pce.png", dpi=900, width = 12, height = 10)
runpce(nTEvsL40.down.pce , "up regulated in Yabe.BAP<40um (compared to nTE)")
ggsave("nTEvsL40.down.pce.png", dpi=900, width = 12, height = 10)
runpce(nCTvsL40.up.pce , "up regulated in nCT (compared to Yabe.BAP<40um)")
ggsave("nCTvsL40.up.pce.png", dpi=900, width = 12, height = 10)
runpce(nCTvsL40.down.pce , "up regulated in Yabe.BAP<40um (compared to nCT)")
ggsave("nCTvsL40.down.pce.png", dpi=900, width = 12, height = 10)
runpce(BAPvsL40.up.pce , "up regulated in BAP.d3 (compared to Yabe.BAP<40um)")
ggsave("BAPvsL40.up.pce.png", dpi=900, width = 12, height = 10)
runpce(BAPvsL40.down.pce , "up regulated in Yabe.BAP<40um (compared to BAP.d3)")
ggsave("BAPvsL40.down.pce.png", dpi=900, width = 12, height = 10)
runpce(STBvsL40.up.pce , "up regulated in Yabe.BAP>70um (compared to Yabe.BAP<40um)")
ggsave("STBvsL40.up.pce.png", dpi=900, width = 12, height = 10)
runpce(STBvsL40.down.pce , "up regulated in Yabe.BAP<40um (compared to Yabe.BAP>70um)")
ggsave("STBvsL40.down.pce.png", dpi=900, width = 12, height = 10)
runpce(BAPvsSTB.up.pce , "up regulated in BAP.d3 (compared to Yabe.BAP>70um)")
ggsave("BAPvsSTB.up.pce.png", dpi=900, width = 12, height = 10)
runpce(BAPvsSTB.down.pce , "up regulated in Yabe.BAP>70um (compared to BAP.d3)")
ggsave("BAPvsSTB.down.pce.png", dpi=900, width = 12, height = 10)

# I had the handy file for conversion, so loaded it up here
# but you can actually use the annot table create above as well
mart <- read.csv("mart-genes.tsv", sep="\t", stringsAsFactors = TRUE, header = TRUE)
# for volcano plots, remove the normalized counts from the results
res.nCTvsBAPdata.clean <- res.nCTvsBAPdata[, c(1:7)]
res.nCTvsSTBdata.clean <- res.nCTvsSTBdata[, c(1:7)]
res.nTEvsSTBdata.clean <- res.nTEvsSTBdata[, c(1:7)]
res.nTEvsBAPdata.clean <- res.nTEvsBAPdata[, c(1:7)]
res.nTEvsL40data.clean <- res.nTEvsL40data[, c(1:7)]
res.nCTvsL40data.clean <- res.nCTvsL40data[, c(1:7)]
res.BAPvsL40data.clean <- res.BAPvsL40data[, c(1:7)]
res.STBvsL40data.clean <- res.STBvsL40data[, c(1:7)]
res.BAPvsSTBdata.clean <- res.BAPvsSTBdata[, c(1:7)]
# merge with the annotation to add-in gene symbols (instead of ensembl ids)
res.nCTvsBAPdata.merged <- merge(res.nCTvsBAPdata.clean, mart, by.x=c("Gene"), by.y = c("ensembl_gene_id_version"))
res.nCTvsSTBdata.merged <- merge(res.nCTvsSTBdata.clean, mart, by.x=c("Gene"), by.y = c("ensembl_gene_id_version"))
res.nTEvsSTBdata.merged <- merge(res.nTEvsSTBdata.clean, mart, by.x=c("Gene"), by.y = c("ensembl_gene_id_version"))
res.nTEvsBAPdata.merged <- merge(res.nTEvsBAPdata.clean, mart, by.x=c("Gene"), by.y = c("ensembl_gene_id_version"))
res.nTEvsL40data.merged <- merge(res.nTEvsL40data.clean, mart, by.x=c("Gene"), by.y = c("ensembl_gene_id_version"))
res.nCTvsL40data.merged <- merge(res.nCTvsL40data.clean, mart, by.x=c("Gene"), by.y = c("ensembl_gene_id_version"))
res.BAPvsL40data.merged <- merge(res.BAPvsL40data.clean, mart, by.x=c("Gene"), by.y = c("ensembl_gene_id_version"))
res.STBvsL40data.merged <- merge(res.STBvsL40data.clean, mart, by.x=c("Gene"), by.y = c("ensembl_gene_id_version"))
res.BAPvsSTBdata.merged <- merge(res.BAPvsSTBdata.clean, mart, by.x=c("Gene"), by.y = c("ensembl_gene_id_version"))
# create volcano plots
nCTvsBAPvolPlot <- EnhancedVolcano(res.nCTvsBAPdata.merged, lab = res.nCTvsBAPdata.merged$gene_symbol, x = 'log2FoldChange', y = 'padj', title = "nCT compared to BAP.d3", subtitle= "", pCutoff = 0.05) + ylab(expression(paste("-log"[10], " pvalue (adj)")))
nCTvsSTBvolPlot <- EnhancedVolcano(res.nCTvsSTBdata.merged, lab = res.nCTvsSTBdata.merged$gene_symbol, x = 'log2FoldChange', y = 'padj', title = "nCT compared to Yabe.BAP>70um", subtitle= "", pCutoff = 0.05) + ylab(expression(paste("-log"[10], " pvalue (adj)")))
nTEvsSTBvolPlot <- EnhancedVolcano(res.nTEvsSTBdata.merged, lab = res.nTEvsSTBdata.merged$gene_symbol, x = 'log2FoldChange', y = 'padj', title = "nTE compared to Yabe.BAP>70um", subtitle= "", pCutoff = 0.05) + ylab(expression(paste("-log"[10], " pvalue (adj)")))
nTEvsBAPvolPlot <- EnhancedVolcano(res.nTEvsBAPdata.merged, lab = res.nTEvsBAPdata.merged$gene_symbol, x = 'log2FoldChange', y = 'padj', title = "nTE compared to BAP.d3", subtitle= "", pCutoff = 0.05) + ylab(expression(paste("-log"[10], " pvalue (adj)")))
nTEvsL40volPlot <- EnhancedVolcano(res.nTEvsL40data.merged, lab = res.nTEvsL40data.merged$gene_symbol, x = 'log2FoldChange', y = 'padj', title = "nTE compared to Yabe.BAP<40um", subtitle= "", pCutoff = 0.05) + ylab(expression(paste("-log"[10], " pvalue (adj)")))
nCTvsL40volPlot <- EnhancedVolcano(res.nCTvsL40data.merged, lab = res.nCTvsL40data.merged$gene_symbol, x = 'log2FoldChange', y = 'padj', title = "nCT compared to Yabe.BAP<40um", subtitle= "", pCutoff = 0.05) + ylab(expression(paste("-log"[10], " pvalue (adj)")))
BAPvsL40volPlot <- EnhancedVolcano(res.BAPvsL40data.merged, lab = res.BAPvsL40data.merged$gene_symbol, x = 'log2FoldChange', y = 'padj', title = "BAP.d3 compared to Yabe.BAP<40um", subtitle= "", pCutoff = 0.05) + ylab(expression(paste("-log"[10], " pvalue (adj)")))
STBvsL40volPlot <- EnhancedVolcano(res.STBvsL40data.merged, lab = res.STBvsL40data.merged$gene_symbol, x = 'log2FoldChange', y = 'padj', title = "Yabe.BAP>70um compared to Yabe.BAP<40um", subtitle= "", pCutoff = 0.05) + ylab(expression(paste("-log"[10], " pvalue (adj)")))
BAPvsSTBvolPlot <- EnhancedVolcano(res.BAPvsSTBdata.merged, lab = res.BAPvsSTBdata.merged$gene_symbol, x = 'log2FoldChange', y = 'padj', title = "BAP.d3 compared to Yabe.BAP>70um", subtitle= "", pCutoff = 0.05) + ylab(expression(paste("-log"[10], " pvalue (adj)")))
# save plots
ggsave("nCTvsBAPvolPlot.png", nCTvsBAPvolPlot,  dpi=900, width = 12, height = 10)
ggsave("nCTvsSTBvolPlot.png", nCTvsSTBvolPlot,  dpi=900, width = 12, height = 10)
ggsave("nTEvsSTBvolPlot.png", nTEvsSTBvolPlot,  dpi=900, width = 12, height = 10)
ggsave("nTEvsBAPvolPlot.png", nTEvsBAPvolPlot,  dpi=900, width = 12, height = 10)
ggsave("nTEvsL40volPlot.png", nTEvsL40volPlot,  dpi=900, width = 12, height = 10)
ggsave("nCTvsL40volPlot.png", nCTvsL40volPlot,  dpi=900, width = 12, height = 10)
ggsave("BAPvsL40volPlot.png", BAPvsL40volPlot,  dpi=900, width = 12, height = 10)
ggsave("STBvsL40volPlot.png", STBvsL40volPlot,  dpi=900, width = 12, height = 10)
ggsave("BAPvsSTBvolPlot.png", BAPvsSTBvolPlot,  dpi=900, width = 12, height = 10)
