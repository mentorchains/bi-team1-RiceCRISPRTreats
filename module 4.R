library(WGCNA)
library(EnhancedVolcano)
library(ggplot2)
library(pheatmap)
library(simpleaffy)
library(GEOquery)
library(tidyverse)
library(hgu133plus2.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(limma)
library(pheatmap)
library(edgeR)


gse <- ReadAffy(celfile.path = "GSE19804_RAW/")
SeriesMatrix <- getGEO(filename = "GSE19804_series_matrix.txt")
metadata19804 <- SeriesMatrix@phenoData@data
df_metaframe <- model.frame(metadata19804)

rma <- rma(gse) #normalisastion

#annotation

sampleNames(rma) = paste(rep(c("C", "N"), each = 60), rep(c(1:60, 1:60), each = 1), sep="")
sampleNames(rma)

#limma

interest <- factor(paste(rep(c("C", "N"), each = 60), sep = ""))
matrix <- model.matrix(~ 0 + interest)


fit <- limma::lmFit(rma, matrix)

make.names(c("interestlung cancer", "interestpaired normal adjacent"), unique = TRUE)
contrast.matrix <- makeContrasts(
  interestC-interestN, 
  levels = matrix)

fit.contrast = contrasts.fit(fit, contrast.matrix)

efit <- eBayes(fit.contrast)

genes=geneNames(gse)
limma_output <- topTable(efit, n = 54670, genelist = genes)

EnhancedVolcano( toptable = limma_output, 
                 lab = rownames(limma_output), 
                 x = "logFC", 
                 y = "P.Value")

#DEG 

genesymbols <- AnnotationDbi::select(hgu133plus2.db, probeID, columns = "SYMBOL")
top10 <- c(rownames(topTable(efit, n = 10)))
top10gene <- subset(genesymbols, genesymbols$PROBEID %in% top10)
top10gene[,2]

#still fixing