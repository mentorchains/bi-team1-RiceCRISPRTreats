library(WGCNA)
library(EnhancedVolcano)
library(ggplot2)
library(affyQCReport)
library(GEOquery)
library(simpleaffy)
library(tidyverse)
library(hgu133plus2.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(limma)
library(pheatmap)
library(edgeR)
library(Biobase)
library(affyQCReport)
library(affy)
library(simpleaffy)
library(affyPLM)
library(GEOquery)
library(utilsIPEA)
library(clusterProfiler)
library(enrichplot)
library(magrittr)
library(tidyr)
library(ggnewscale)
library(msigdbr)
library(readxl)

#PCA
library(simpleaffy)
library(tidyverse)
library(ggfortify)

#metadata

gse <- ReadAffy(celfile.path = "gse4107/GSE4107_RAW")
gse19804 <- getGEO(filename = "gse4107/GSE4107_series_matrix.txt")
metadata <- gse19804@phenoData@data
CN <- metadata[c("title")]
CN_1 <- metadata[c("title")]

rma <- rma(gse) 

#Raw Data Plot
colour <- c(rep('cancer',12), rep('normal', 10))

rawGSE <- exprs(gse)
pca_OG <- prcomp(rawGSE, scale=F, center=F)
pca_OG <- as.data.frame(pca_OG$rotation)

pca_2 <- cbind(pca_OG, new_col = colour)

PCA_raw <- ggplot(pca_2, aes(x=PC1, y=PC2, color=colour))+ geom_point() + stat_ellipse() +
  labs(title = "Raw Data PCA Plot")

#Normalized Plot
normalized <- rma(gse)
normalizedexprs <- exprs(normalized)
pca_normalized <- prcomp(normalizedexprs, scale=F, center=F)
pca_normalized <- as.data.frame(pca_normalized$rotation)

pca_normalized2 <- cbind(pca_normalized, new_col = colour)

PCA_normalised <- ggplot(pca_normalized2, aes(x=PC1, y=PC2, color=colour))+ geom_point() + stat_ellipse() +
  labs(title = "Normalized Data PCA Plot") 


#rma normalisation boxplot 
norm <- rma(gse)
df <- exprs(norm)
rma_boxplot <- boxplot(df, main="Relative Signal BoxPlot map", ylab="Relative log expression signal-RMA", las=2)




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
geneFilt <- quantile(mean, p=0.02)
geneFilt <- annotated[which(mean>geneFilt),]

#limma
#creating model
CN <- data.frame(Tissue = metadata$`title`)
rownames(CN) <- rownames(metadata)
CN$Tissue <- ifelse(str_detect(CN$Tissue, regex("patient")), 1, 0)
matrix <- model.matrix(~ Tissue, CN)
fit <- limma::lmFit(remove_lower_0.02, matrix)
efit <- eBayes(fit)
genes=geneNames(gse)
limma_output <- topTable(efit, adjust.method="fdr", sort.by=, n = 50000)
EnhancedVolcano( toptable = limma_output, 
                 lab = rownames(limma_output), 
                 x = "logFC", 
                 y = "P.Value")

#heatmap
group <- data.frame(Tissue = metadata$title)
top50 <- topTable(efit, number=50) 
filt_50 <- data.frame(geneFilt[rownames(geneFilt) %in% rownames(top50),])
heatmap_4107 <- pheatmap(filt_50, main="Top 50 DEG", 
         labels_row = rownames(filt_50),  
         labels_col = colnames(filt_50),
         annotation_row=group, 
         annotation_col=group, 
         annotation_colors=list(Tissue=c(Cancer="red", Normal="blue")))

#DEG 

top10 <- rownames(topTable(efit, n = 10))

#MODULE 5
#DEG 
logFC <- limma_output %>% dplyr::select('logFC')
logFC_vec <- as.vector(logFC$logFC)
element_names <- rownames(limma_output)
names(logFC_vec) <- element_names

#threshold + filtering (sorted, named, numeric vector)
filtered_FC <- logFC_vec[logFC_vec > 1.5]
arrange_FC <- sort(filtered_FC, decreasing = TRUE)
topDEG <- data.frame(arrange_FC)

#selecting symbol and entrezid
Entrezid_symbol <- AnnotationDbi::select(hgu133plus2.db, keys = ID, 
                                         columns = c("SYMBOL", "ENTREZID"))
df_entrezid <- Entrezid_symbol %>% dplyr::select("SYMBOL", "ENTREZID")
names <- names(arrange_FC)
selected <- df_entrezid[df_entrezid$SYMBOL %in% names, ]
selected <- unique(selected)
rownames(selected) <- 1:nrow(selected)

#gene ontology + boxplots
CC <- enrichGO(selected$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC" , readable = T )
barplot_GO_CC <- barplot(CC)

MF <- enrichGO(selected$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF" , readable = T )
barplot_GO_MF <- barplot(MF)

BP <- enrichGO(selected$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP" , readable = T )
barplot_GO_BP <- barplot(BP)

#KEGG
KEGG <- enrichKEGG(selected$ENTREZID, pvalueCutoff = 0.05)
dotplot_KEGG <- dotplot(KEGG)

#Gene-concept network
convert <- setReadable(KEGG, OrgDb=org.Hs.eg.db, keyType = "ENTREZID")
cnetplot_GCN <- cnetplot(convert, foldChange = arrange_FC, categorySize = "pvalue")

#Global/universal gene set enrichment analysis (GSEA)
gene_set <- msigdbr(species = "Homo sapiens", category = "H")
h <- gene_set %>% select (gs_name, entrez_gene)

removed_duplicate <- Entrezid_symbol[!duplicated(Entrezid_symbol$SYMBOL),]
removed_duplicate <- na.omit(removed_duplicate)
rowname_symbol <- removed_duplicate$SYMBOL
rownames(removed_duplicate) <- rowname_symbol
GSEA_merge <- merge(x = removed_duplicate, y = limma_output, by = "row.names")
GSEA_logFC <- as.vector(GSEA_merge$logFC)
GSEA_entrezid <- as.vector(GSEA_merge$ENTREZID)
names(GSEA_logFC) <- GSEA_entrezid
sorted <- sort(GSEA_logFC, decreasing = TRUE)

GSEA_analysis <- GSEA(sorted, TERM2GENE = h)
GSEA_plot <- gseaplot(GSEA_analysis, geneSetID = 2)

for_module_6 <- names(filtered_FC)
write.table(arrange_FC, file = 'DEG.txt')



