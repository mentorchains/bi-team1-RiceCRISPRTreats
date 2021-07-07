# Kelly Zhang and Roman Ramirez
# rr8rk@virginia.edu
# Level 1, Module 5, STEM-Away

### 1: LOAD LIBRARIES ###

# library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(magrittr)
library(tidyr)
library(ggnewscale)

### 2.1: LOADING VECTOR FROM MODULE 4 ###

# library(WGCNA)
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

gse <- ReadAffy(celfile.path = "data/GSE19804_RAW/")
gse19804 <- getGEO(filename = "data/GSE19804_series_matrix.txt")
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
#creating model
CN <- data.frame(Tissue = metadata$`tissue:ch1`)
rownames(CN) <- rownames(metadata)
CN$Tissue <- ifelse(CN$Tissue == 'lung cancer', 1, 0)
matrix <- model.matrix(~ Tissue, CN)
fit <- limma::lmFit(remove_lower_0.02, matrix)
efit <- eBayes(fit)
genes=geneNames(gse)
limma_output <- topTable(efit, p.value=0.05, adjust.method="fdr", sort.by=, n = 50000)
# EnhancedVolcano( toptable = limma_output, 
#                  lab = rownames(limma_output), 
#                  x = "logFC", 
#                  y = "P.Value")

#DEG 
top10 <- rownames(topTable(efit, n = 10))
top10_df <- as.data.frame(top10)
top10_merge <- merge(x = top10_df, y = probe, by.x = "top10", by.y = "PROBEID")
DEG <- top10_merge$SYMBOL

### 2: GENERATE MODIFIED DEG VECTOR ###

logFC <- limma_output$logFC
genes <- rownames(limma_output)
# AnnotationDbi::select(hgu133plus2.db, keys = topTableVector,  columns = rownames(topTableGED))

DEGvector <- data.frame(genes, logFC)

THRESHOLD = 1.5
# thresholdDEGVector <- DEGvector[which(DEGvector$logFC > THRESHOLD | DEGvector$logFC < -THRESHOLD),]
thresholdDEGVector <- DEGvector[which(DEGvector$logFC > THRESHOLD),]
sorted <- thresholdDEGVector[order(thresholdDEGVector$logFC, decreasing=TRUE),]

dfEntrezid <- AnnotationDbi::select(org.Hs.eg.db, keys=sorted$genes, columns="ENTREZID", keytype="SYMBOL")


### 3: GENE ONTOLOGY ###

goOutput <- enrichGO(dfEntrezid$ENTREZID, OrgDb = org.Hs.eg.db, readable=T)
barplot(goOutput)

### 4: KEGG (KYOTO ENCYCLOPEDIA OF GENES AND GENOMES) ANALYSIS ###

keggOutput <- enrichKEGG(dfEntrezid$ENTREZID)
dotplot(keggOutput)

### 5: GENE-CONCEPT NETWORK ###

enrichResultObject <- setReadable(keggOutput, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
cnetplot(enrichResultObject, foldChange = sorted$logFC, categorySize="pvalue")

### 6: GLOBAL/UNIVERSAL GENE SET ENRICHMENT ANALYSIS (GSEA) ###

# using limma_output to get an unfiltered data set
unfilteredDEGs <- data.frame(genes=rownames(limma_output), logFC=limma_output$logFC)
GSEAlogFC <- limma_output$logFC
names(GSEAlogFC) <- rownames(limma_output)
# read gmt downloaded from http://www.gsea-msigdb.org/gsea/downloads.jsp#msigdb
hallmarkGeneSets <- read.gmt("data/h.all.v7.4.symbols.gmt")
# obtain gene set
geneSet <- msigdbr(species = "Homo sapiens", category = "H")
# select desired information
h <- geneSet %>% dplyr::select(gs_name, entrez_gene)


# Generating Global DEG vector

# Convert Gene Symbols to EntrezID
converter <- AnnotationDbi::select(org.Hs.eg.db, keys=unfilteredDEGs$genes, columns="ENTREZID", keytype="SYMBOL")
converterNoDup <- converter[!duplicated(converter$SYMBOL),]
row.names(converterNoDup) <- converterNoDup$SYMBOL
mergedDEGs <- merge(converterNoDup, GSEAlogFC, by = "row.names")

### WIP?
globalDEGVector <- data.frame(symbol=unfilteredDEGs$genes, entrezid=converterNoDup$"ENTREZID", logFC=unfilteredDEGs$logFC)
globalDEGVectorSorted <- globalDEGVector[order(globalDEGVector$logFC, decreasing=TRUE),]
logFCVector <- data.frame(globalDEGVectorSorted$logFC)
rownames(logFCVector) <- globalDEGVectorSorted$entrezid

# Making a list of values: ENTREZID and logFC
geneVector <- mergedDEGs$y
names(geneVector) <- mergedDEGs$ENTREZID
geneVectorSorted <- sort(geneVector, decreasing=TRUE)

# Perform GSEA Analysis
gseaResultObject <- GSEA(geneVectorSorted,TERM2GENE = h)
gseaplot2(gseaResultObject,
          geneSetID = c("HALLMARK_E2F_TARGETS",
                        "HALLMARK_G2M_CHECKPOINT",
                        "HALLMARK_MTORC1_SIGNALING",
                        "HALLMARK_MYC_TARGETS_V1",
                        "HALLMARK_MYC_TARGETS_V2"),
         
          #,color = c("red", "yellow", "green", "blue", "purple")
    )
