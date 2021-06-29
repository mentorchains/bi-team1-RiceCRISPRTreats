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

#from module 4
#metadata
gse <- ReadAffy(celfile.path = "GSE19804_RAW/")
gse19804 <- getGEO(filename = "GSE19804_series_matrix.txt")
metadata <- gse19804@phenoData@data
CN <- metadata[c("tissue:ch1")]
#Normalisation
rma <- rma(gse) 
#Annotation and remove duplicates + NA
ID <- rownames(gse)
probe <- AnnotationDbi::select(hgu133plus2.db, keys = ID,  columns = "SYMBOL")
duplicate_probeID <- probe[!duplicated(probe$PROBEID),]
rma_exprs <- Biobase::exprs(rma)
table_merge <- merge(x = duplicate_probeID, y = rma_exprs, by.x = "PROBEID", by.y = "row.names")
table_merge <- table_merge[!duplicated(table_merge$SYMBOL),]
table_merge <- na.omit(table_merge)
rownames(table_merge) <- table_merge$SYMBOL
annotated <- table_merge[-c(1,2)]
#limma
CN <- data.frame(Tissue = metadata$`tissue:ch1`)
rownames(CN) <- rownames(metadata)
CN$Tissue <- ifelse(CN$Tissue == 'lung cancer', 1, 0)
matrix <- model.matrix(~ Tissue, CN)
fit <- limma::lmFit(annotated, matrix)
efit <- eBayes(fit)
genes=geneNames(gse)
limma_output <- topTable(efit, p.value=0.05, adjust.method="fdr", sort.by=, n = 50000)

#MODULE 5
#DEG 
logFC <- limma_output %>% dplyr::select('logFC')
SYMBOLS_FC <-rownames(logFC)
logFC_3 <- mutate(logFC, id = SYMBOLS_FC)
logFC_id_probe <- merge(x = logFC_3, y = probe, by.x = "id", by.y = "SYMBOL")

#threshold + filtering (sorted, named, numeric vector)
filtered_FC <- filter(logFC_id_probe, logFC < -2)
arrange_FC <- filtered_FC %>% arrange(logFC)
FC_vec <- as.vector(arrange_FC$logFC)
names(FC_vec) <- arrange_FC$id

#selecting symbol and entrezid
Entrezid_symbol <- AnnotationDbi::select(hgu133plus2.db, keys = ID, 
                                         columns = c("ENTREZID", "SYMBOL"))
df_entrezid <- Entrezid_symbol %>% dplyr::select("SYMBOL", "ENTREZID")

#gene ontology
CC <- enrichGO(df_entrezid$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC" , readable = T )
barplot_GO_CC <- barplot(CC)

MF <- enrichGO(df_entrezid$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF" , readable = T )
barplot_GO_MF <- barplot(MF)

BP <- enrichGO(df_entrezid$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP" , readable = T )
barplot_GO_BP <- barplot(BP)
