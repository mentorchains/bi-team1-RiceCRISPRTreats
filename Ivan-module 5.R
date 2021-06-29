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


#DEG 
top <- topTable(efit, n = 50000)

logFC <- top %>% dplyr::select('logFC')
rownames_FC <-rownames(logFC)
logFC_3 <- mutate(logFC, id = rownames_FC)
logFC_id_probe <- merge(x = logFC_3, y = probe, by.x = "id", by.y = "PROBEID")
logFC_symbol <- as.vector(logFC_id_probe$SYMBOL)
IDnames <- rownames(top)

#threshold + filtering
filtered_FC <- filter(logFC, logFC < -2)
arrange_FC <- filtered_FC %>% arrange(logFC)
FC_vec <- as.vector(arrange_FC$logFC)
names(FC_vec) <- rownames(arrange_FC)

#selecting symbol and entrezid
Entrezid_symbol <- AnnotationDbi::select(hgu133plus2.db, keys = ID, 
                                         columns = c("ENTREZID", "SYMBOL"))
df_entrezid <- Entrezid_symbol %>% dplyr::select("SYMBOL", "ENTREZID")
filter_FC <- filter(logFC_3, logFC < -2)
filtered_with_sym <-  merge(x = filter_FC, y = probe, by.x = "id", by.y = "PROBEID")
for_ego <- merge(x = filtered_with_sym, y = df_entrezid, by.x = "SYMBOL", by.y = "SYMBOL")
#gene ontology
enrichGO(df_entrezid, keyType = "SYMBOL" , OrgDb = org.Hs.eg.db, ont = "CC", "MF", "BP" , readable = T )
