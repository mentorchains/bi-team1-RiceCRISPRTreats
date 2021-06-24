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
library(Biobase)

gse <- ReadAffy(celfile.path = "GSE19804_RAW/")

rma <- rma(gse) #normalisastion

#gene filtering

ID <- rownames(gse)

probe <- AnnotationDbi::select(hgu133plus2.db, keys = ID,  columns = "SYMBOL")
duplicate_probeID <- probe[!duplicated(probe$PROBEID),]
duplicate_probeID_symbol <- duplicate_probeID[!duplicated(duplicate_probeID$SYMBOL),]
filtered <- duplicate_probeID_symbol[!is.na(duplicate_probeID_symbol$SYMBOL),]

rma_exprs <- Biobase::exprs(rma)

table_merge <- merge(x = filtered, y = rma_exprs)

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

top10 <- c(rownames(topTable(efit, n = 10)))
top10gene <- subset(filtered, filtered$PROBEID %in% top10)
top10gene[,2]

#still fixing