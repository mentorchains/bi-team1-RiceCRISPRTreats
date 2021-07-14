#By Arian Veyssi and Shreya Vora

#loading in the metadata
data <- read.delim("New_MetaData.txt", header = TRUE)
#data <- data[1:120,]
#importing the raw data using the simplyaffy function
library(simpleaffy)
raw.data = ReadAffy(compress = T, celfile.path = "C:/ShreyaVora/stemaway/R_Projects/Module5/GSE19804_RAW")
#aw.data

#normalizing the raw data
norData <- rma(raw.data, normalize = T, background = T)
norData <- exprs(norData)

library(hgu133plus2.db)
norData_rown <- rownames(norData)
anno <- select(hgu133plus2.db, keys = norData_rown, columns = 'SYMBOL')

#removes duplicates
anno <- anno[!duplicated(anno$PROBEID),]

#changing the rownames in the normaized data so that the rownames will be the same so that the norData and the anno can merge
rownames(anno) <- anno$PROBEID
#merge the anno (gives the probids), and the normalized data frame that has the symbols
merged_data <- merge(anno, norData, by = "row.names")

#removing duplicated symbols
merged_data <- merged_data[!duplicated(merged_data$SYMBOL),]

#removing NA values
merged_data <- na.omit(merged_data)
rownames(merged_data) <- merged_data$SYMBOL
#removing unecessary colunms from the merged data and storing it in the normData
norData <- merged_data[-c(1,2,3)]
rownames(norData) <- merged_data$SYMBOL

#getting the mean of each row in the dataset
r_mean <- rowMeans(norData)
quan <- quantile(r_mean , 0.02)
filtered <- norData[which(r_mean > quan),]

#limma analysis
rownames(data) <- data$name
data_set <- data.frame(tissue = data$tissue)
rownames(data_set) <- rownames(data)



#limma analysis
library(limma)
mm <- model.matrix(~tissue , data)
linear_model <- lmFit(filtered, mm)
ebays_mod <- eBayes(fit = linear_model)


#topTable
DEG <- topTable(ebays_mod, p.value = 0.05, adjust.method = "fdr", sort.by = "P", gene = rownames(norData), number = 10000)
DEG[c(1)] <- rownames(DEG)


#module5
library(ggplot2)
library(GEOquery)
library(pheatmap)
library(simpleaffy)
library(tidyverse)
library(hgu133plus2.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(limma)
library(pheatmap)
library(edgeR)
library(Biobase)
library(dplyr)
#module 5
library(clusterProfiler)
library(enrichplot)
library(magrittr)
library(tidyr)
library(ggnewscale)
library(msigdb)


#generating DEG Vector
#DEG 
logFC <- DEG %>% dplyr::select('logFC')
logFC_vec <- as.vector(logFC$logFC)
element_names <- rownames(DEG)
names(logFC_vec) <- element_names

#threshold + filtering (sorted, named, numeric vector)
filtered_FC <- logFC_vec[logFC_vec > 2]
arrange_FC <- sort(filtered_FC, decreasing = TRUE)

#selecting symbol and entrezid
Entrezid_symbol <- AnnotationDbi::select(hgu133plus2.db, keys = rownames(raw.data), 
                                         columns = c("SYMBOL", "ENTREZID"))
df_entrezid <- Entrezid_symbol %>% dplyr::select("SYMBOL", "ENTREZID")
names <- names(arrange_FC)
selected <- df_entrezid[df_entrezid$SYMBOL %in% names, ]
selected <- unique(selected)
rownames(selected) <- 1:nrow(selected)

#gene ontology 

BP <- enrichGO(selected$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", readable = T )
barPlot_BP <- barplot(BP)

MF <- enrichGO(selected$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF", readable = T )
barPlot_MF <- barplot(MF)

CC <- enrichGO(selected$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC", readable = T )
barPlot_CC <- barplot(CC)