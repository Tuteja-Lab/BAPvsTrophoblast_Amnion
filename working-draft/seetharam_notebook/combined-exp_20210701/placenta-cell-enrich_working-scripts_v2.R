library(TissueEnrich)
library(gprofiler2)
library(tidyverse)

setwd("~/TutejaLab/PlacentaEnrich")
l <- load(file = "combine-test-expression1.Rdata")
humanGeneMapping <- dataset$GRCH38$humanGeneMapping
d <- dataset$PlacentaDeciduaBloodData
data <- d$expressionData
cellDetails <- d$cellDetails

setwd("~/TutejaLab/combined-exp_20210701/pce-results")

inputGenes<-toupper(markers.filtered.names.1)
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
p1<- ggplot(data = enrichmentOutput, mapping = aes(x = reorder (Tissue, -Log10PValue), Log10PValue, fill= Tissue)) +
  geom_bar(stat = "identity") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 100), legend.position = "none",
        plot.title = element_text(color = "black", size=18, face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggtitle("Cluster 1") + ylab("-log10 p-value")
ggsave("Cluster_1.png", dpi=700, width = 10, height = 8)
write.table(enrichmentOutput[order(-enrichmentOutput$Log10PValue),], file = "cluster1-cellEnrich.csv", quote = T, sep = ",")
enriched <- enrichmentOutput %>% filter(Log10PValue != 0) %>% rownames()
for (p in enriched) {
  new <- as.data.frame(assay(output2[[2]][[p]]))[1]
  conv <- humanGeneMapping[humanGeneMapping$Gene %in% as.list(new$Gene),]
  tegenes <- paste( unlist(conv$Gene.name), collapse=';')
  df <- data.frame(row.names = p, tegenes)
  write.table(df, file = paste0("cluster1-",p,"-cellEnrich.txt"), quote = T, sep = ",")
}
inputGenes<-toupper(markers.filtered.names.2)
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
p2 <- ggplot(data = enrichmentOutput, mapping = aes(x = reorder (Tissue, -Log10PValue), Log10PValue, fill= Tissue)) +
  geom_bar(stat = "identity") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 100), legend.position = "none",
        plot.title = element_text(color = "black", size=18, face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggtitle("Cluster 2") + ylab("-log10 p-value")
ggsave("Cluster_2.png", dpi=700, width = 10, height = 8)
write.table(enrichmentOutput[order(-enrichmentOutput$Log10PValue),], file = "cluster2-cellEnrich.csv", quote = T, sep = ",")
enriched <- enrichmentOutput %>% filter(Log10PValue != 0) %>% rownames()
for (p in enriched) {
  new <- as.data.frame(assay(output2[[2]][[p]]))[1]
  conv <- humanGeneMapping[humanGeneMapping$Gene %in% as.list(new$Gene),]
  tegenes <- paste( unlist(conv$Gene.name), collapse=';')
  df <- data.frame(row.names = p, tegenes)
  write.table(df, file = paste0("cluster2-",p,"-cellEnrich.txt"), quote = T, sep = ",")
}
inputGenes<-toupper(markers.filtered.names.3)
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
p3 <- ggplot(data = enrichmentOutput, mapping = aes(x = reorder (Tissue, -Log10PValue), Log10PValue, fill= Tissue)) +
  geom_bar(stat = "identity") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 100), legend.position = "none",
        plot.title = element_text(color = "black", size=18, face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggtitle("Cluster 3") + ylab("-log10 p-value")
ggsave("Cluster_3.png", dpi=700, width = 10, height = 8)
write.table(enrichmentOutput[order(-enrichmentOutput$Log10PValue),], file = "cluster3-cellEnrich.csv", quote = T, sep = ",")
enriched <- enrichmentOutput %>% filter(Log10PValue != 0) %>% rownames()
for (p in enriched) {
  new <- as.data.frame(assay(output2[[2]][[p]]))[1]
  conv <- humanGeneMapping[humanGeneMapping$Gene %in% as.list(new$Gene),]
  tegenes <- paste( unlist(conv$Gene.name), collapse=';')
  df <- data.frame(row.names = p, tegenes)
  write.table(df, file = paste0("cluster3-",p,"-cellEnrich.txt"), quote = T, sep = ",")
}
inputGenes<-toupper(markers.filtered.names.4)
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
p4 <- ggplot(data = enrichmentOutput, mapping = aes(x = reorder (Tissue, -Log10PValue), Log10PValue, fill= Tissue)) +
  geom_bar(stat = "identity") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 100), legend.position = "none",
        plot.title = element_text(color = "black", size=18, face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggtitle("Cluster 4") + ylab("-log10 p-value")
ggsave("Cluster_4.png", dpi=700, width = 10, height = 8)
write.table(enrichmentOutput[order(-enrichmentOutput$Log10PValue),], file = "cluster4-cellEnrich.csv", quote = T, sep = ",")
enriched <- enrichmentOutput %>% filter(Log10PValue != 0) %>% rownames()
for (p in enriched) {
  new <- as.data.frame(assay(output2[[2]][[p]]))[1]
  conv <- humanGeneMapping[humanGeneMapping$Gene %in% as.list(new$Gene),]
  tegenes <- paste( unlist(conv$Gene.name), collapse=';')
  df <- data.frame(row.names = p, tegenes)
  write.table(df, file = paste0("cluster4-",p,"-cellEnrich.txt"), quote = T, sep = ",")
}
inputGenes<-toupper(markers.filtered.names.5)
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
p5 <- ggplot(data = enrichmentOutput, mapping = aes(x = reorder (Tissue, -Log10PValue), Log10PValue, fill= Tissue)) +
  geom_bar(stat = "identity") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 100), legend.position = "none",
        plot.title = element_text(color = "black", size=18, face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggtitle("Cluster 5") + ylab("-log10 p-value")
ggsave("Cluster_5.png", dpi=700, width = 10, height = 8)
write.table(enrichmentOutput[order(-enrichmentOutput$Log10PValue),], file = "cluster5-cellEnrich.csv", quote = T, sep = ",")
enriched <- enrichmentOutput %>% filter(Log10PValue != 0) %>% rownames()
for (p in enriched) {
  new <- as.data.frame(assay(output2[[2]][[p]]))[1]
  conv <- humanGeneMapping[humanGeneMapping$Gene %in% as.list(new$Gene),]
  tegenes <- paste( unlist(conv$Gene.name), collapse=';')
  df <- data.frame(row.names = p, tegenes)
  write.table(df, file = paste0("cluster5-",p,"-cellEnrich.txt"), quote = T, sep = ",")
}
inputGenes<-toupper(markers.filtered.names.6)
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
p6 <- ggplot(data = enrichmentOutput, mapping = aes(x = reorder (Tissue, -Log10PValue), Log10PValue, fill= Tissue)) +
  geom_bar(stat = "identity") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 100), legend.position = "none",
        plot.title = element_text(color = "black", size=18, face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggtitle("Cluster 6") + ylab("-log10 p-value")
ggsave("Cluster_6.png", dpi=700, width = 10, height = 8)
write.table(enrichmentOutput[order(-enrichmentOutput$Log10PValue),], file = "cluster6-cellEnrich.csv", quote = T, sep = ",")
enriched <- enrichmentOutput %>% filter(Log10PValue != 0) %>% rownames()
for (p in enriched) {
  new <- as.data.frame(assay(output2[[2]][[p]]))[1]
  conv <- humanGeneMapping[humanGeneMapping$Gene %in% as.list(new$Gene),]
  tegenes <- paste( unlist(conv$Gene.name), collapse=';')
  df <- data.frame(row.names = p, tegenes)
  write.table(df, file = paste0("cluster6-",p,"-cellEnrich.txt"), quote = T, sep = ",")
}

inputGenes<-toupper(markers.filtered.names.7)
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
p7 <- ggplot(data = enrichmentOutput, mapping = aes(x = reorder (Tissue, -Log10PValue), Log10PValue, fill= Tissue)) +
  geom_bar(stat = "identity") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 100), legend.position = "none",
        plot.title = element_text(color = "black", size=18, face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggtitle("Cluster 7") + ylab("-log10 p-value")
ggsave("Cluster_7.png", dpi=700, width = 10, height = 8)
write.table(enrichmentOutput[order(-enrichmentOutput$Log10PValue),], file = "cluster7-cellEnrich.csv", quote = T, sep = ",")
enriched <- enrichmentOutput %>% filter(Log10PValue != 0) %>% rownames()
for (p in enriched) {
  new <- as.data.frame(assay(output2[[2]][[p]]))[1]
  conv <- humanGeneMapping[humanGeneMapping$Gene %in% as.list(new$Gene),]
  tegenes <- paste( unlist(conv$Gene.name), collapse=';')
  df <- data.frame(row.names = p, tegenes)
  write.table(df, file = paste0("cluster7-",p,"-cellEnrich.txt"), quote = T, sep = ",")
}
inputGenes<-toupper(markers.filtered.names.8)
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
p8 <- ggplot(data = enrichmentOutput, mapping = aes(x = reorder (Tissue, -Log10PValue), Log10PValue, fill= Tissue)) +
  geom_bar(stat = "identity") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 100), legend.position = "none",
        plot.title = element_text(color = "black", size=18, face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggtitle("Cluster_8") + ylab("-log10 p-value")
ggsave("Cluster_8.png", dpi=700, width = 10, height = 8)
write.table(enrichmentOutput[order(-enrichmentOutput$Log10PValue),], file = "cluster8-cellEnrich.csv", quote = T, sep = ",")
enriched <- enrichmentOutput %>% filter(Log10PValue != 0) %>% rownames()
for (p in enriched) {
  new <- as.data.frame(assay(output2[[2]][[p]]))[1]
  conv <- humanGeneMapping[humanGeneMapping$Gene %in% as.list(new$Gene),]
  tegenes <- paste( unlist(conv$Gene.name), collapse=';')
  df <- data.frame(row.names = p, tegenes)
  write.table(df, file = paste0("cluster8-",p,"-cellEnrich.txt"), quote = T, sep = ",")
}

inputGenes<-toupper(markers.filtered.names.9)
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
p9<- ggplot(data = enrichmentOutput, mapping = aes(x = reorder (Tissue, -Log10PValue), Log10PValue, fill= Tissue)) +
  geom_bar(stat = "identity") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 100), legend.position = "none",
        plot.title = element_text(color = "black", size=18, face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggtitle("Cluster 9") + ylab("-log10 p-value")
ggsave("Cluster_9.png", dpi=700, width = 10, height = 8)
write.table(enrichmentOutput[order(-enrichmentOutput$Log10PValue),], file = "cluster9-cellEnrich.csv", quote = T, sep = ",")
enriched <- enrichmentOutput %>% filter(Log10PValue != 0) %>% rownames()
for (p in enriched) {
  new <- as.data.frame(assay(output2[[2]][[p]]))[1]
  conv <- humanGeneMapping[humanGeneMapping$Gene %in% as.list(new$Gene),]
  tegenes <- paste( unlist(conv$Gene.name), collapse=';')
  df <- data.frame(row.names = p, tegenes)
  write.table(df, file = paste0("cluster9-",p,"-cellEnrich.txt"), quote = T, sep = ",")
}




inputGenes<-toupper(markers.filtered.names.10)
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
p9<- ggplot(data = enrichmentOutput, mapping = aes(x = reorder (Tissue, -Log10PValue), Log10PValue, fill= Tissue)) +
  geom_bar(stat = "identity") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 100), legend.position = "none",
        plot.title = element_text(color = "black", size=18, face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggtitle("Cluster 10") + ylab("-log10 p-value")
ggsave("Cluster_10.png", dpi=700, width = 10, height = 8)
write.table(enrichmentOutput[order(-enrichmentOutput$Log10PValue),], file = "cluster10-cellEnrich.csv", quote = T, sep = ",")
enriched <- enrichmentOutput %>% filter(Log10PValue != 0) %>% rownames()
for (p in enriched) {
  new <- as.data.frame(assay(output2[[2]][[p]]))[1]
  conv <- humanGeneMapping[humanGeneMapping$Gene %in% as.list(new$Gene),]
  tegenes <- paste( unlist(conv$Gene.name), collapse=';')
  df <- data.frame(row.names = p, tegenes)
  write.table(df, file = paste0("cluster10-",p,"-cellEnrich.txt"), quote = T, sep = ",")
}


inputGenes<-toupper(markers.filtered.names.11)
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
p9<- ggplot(data = enrichmentOutput, mapping = aes(x = reorder (Tissue, -Log10PValue), Log10PValue, fill= Tissue)) +
  geom_bar(stat = "identity") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 100), legend.position = "none",
        plot.title = element_text(color = "black", size=18, face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggtitle("Cluster 11") + ylab("-log10 p-value")
ggsave("Cluster_11.png", dpi=700, width = 10, height = 8)
write.table(enrichmentOutput[order(-enrichmentOutput$Log10PValue),], file = "cluster11-cellEnrich.csv", quote = T, sep = ",")
enriched <- enrichmentOutput %>% filter(Log10PValue != 0) %>% rownames()
for (p in enriched) {
  new <- as.data.frame(assay(output2[[2]][[p]]))[1]
  conv <- humanGeneMapping[humanGeneMapping$Gene %in% as.list(new$Gene),]
  tegenes <- paste( unlist(conv$Gene.name), collapse=';')
  df <- data.frame(row.names = p, tegenes)
  write.table(df, file = paste0("cluster11-",p,"-cellEnrich.txt"), quote = T, sep = ",")
}



inputGenes<-toupper(markers.filtered.names.12)
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
p9<- ggplot(data = enrichmentOutput, mapping = aes(x = reorder (Tissue, -Log10PValue), Log10PValue, fill= Tissue)) +
  geom_bar(stat = "identity") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 100), legend.position = "none",
        plot.title = element_text(color = "black", size=18, face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggtitle("Cluster 12") + ylab("-log10 p-value")
ggsave("Cluster_12.png", dpi=700, width = 10, height = 8)
write.table(enrichmentOutput[order(-enrichmentOutput$Log10PValue),], file = "cluster12-cellEnrich.csv", quote = T, sep = ",")
enriched <- enrichmentOutput %>% filter(Log10PValue != 0) %>% rownames()
for (p in enriched) {
  new <- as.data.frame(assay(output2[[2]][[p]]))[1]
  conv <- humanGeneMapping[humanGeneMapping$Gene %in% as.list(new$Gene),]
  tegenes <- paste( unlist(conv$Gene.name), collapse=';')
  df <- data.frame(row.names = p, tegenes)
  write.table(df, file = paste0("cluster12-",p,"-cellEnrich.txt"), quote = T, sep = ",")
}

inputGenes<-toupper(markers.filtered.names.13)
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
p9<- ggplot(data = enrichmentOutput, mapping = aes(x = reorder (Tissue, -Log10PValue), Log10PValue, fill= Tissue)) +
  geom_bar(stat = "identity") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 100), legend.position = "none",
        plot.title = element_text(color = "black", size=18, face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggtitle("Cluster 13") + ylab("-log10 p-value")
ggsave("Cluster_13.png", dpi=700, width = 10, height = 8)
write.table(enrichmentOutput[order(-enrichmentOutput$Log10PValue),], file = "cluster13-cellEnrich.csv", quote = T, sep = ",")
enriched <- enrichmentOutput %>% filter(Log10PValue != 0) %>% rownames()
for (p in enriched) {
  new <- as.data.frame(assay(output2[[2]][[p]]))[1]
  conv <- humanGeneMapping[humanGeneMapping$Gene %in% as.list(new$Gene),]
  tegenes <- paste( unlist(conv$Gene.name), collapse=';')
  df <- data.frame(row.names = p, tegenes)
  write.table(df, file = paste0("cluster13-",p,"-cellEnrich.txt"), quote = T, sep = ",")
}
