library(WGCNA)
library(EnhancedVolcano)
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

#metadata

gse <- ReadAffy(celfile.path = "GSE19804_RAW/")
gse19804 <- getGEO(filename = "GSE19804_series_matrix.txt")
metadata <- gse19804@phenoData@data
CN <- metadata[c("tissue:ch1")]

#Normalisation

rma <- rma(gse) 

#Annotation

ID <- rownames(gse)

probe <- AnnotationDbi::select(hgu133plus2.db, keys = ID,  columns = "SYMBOL")
duplicate_probeID <- probe[!duplicated(probe$PROBEID),]

rma_exprs <- Biobase::exprs(rma)

table_merge <- merge(x = duplicate_probeID, y = rma_exprs, by.x = "PROBEID", by.y = "row.names")

table_merge <- table_merge[!duplicated(table_merge$SYMBOL),]
table_merge <- na.omit(table_merge)
rownames(table_merge) <- table_merge$SYMBOL
annotated <- table_merge[-c(1,2)]

#Gene filtering

mean <- rowMeans(annotated)
remove_lower_0.02 <- annotated[which(mean > 0.02),]


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

top10 <- topTable(efit, number = 10)
rownames(top10)

#still fixing