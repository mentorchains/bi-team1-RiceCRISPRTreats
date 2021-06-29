# Module 4 Deliverables
# Level 1 Module 4 Group 4 (L1M4G4)
# Aditi Verma 
# Shreya Vora
# Kelly Zhang

# DOWNLOAD FOLLOWING PACKAGES
# BiocManager::install("hgu133plus2.db")
# BiocManager::install("WGCNA")
# BiocManager::install("limma")
# BiocManager::install('EnhancedVolcano')

# INSTALL PACKAGES
library("hgu133plus2.db")
library("WGCNA")
library("limma")
library("EnhancedVolcano")

# DATA
# remove outliers from metadata GSE19804_METADATA
GSE19804_METADATA <- read_excel("GSE19804_METADATA.xlsx")
GSE19804_RAW <- ReadAffy(compress = T, celfile.path = "GSE19804_RAW/")
GSE19804_METADATA <- na.omit(GSE19804_METADATA)
rownames(GSE19804_METADATA) <- GSE19804_METADATA$name


# NORMALIZE RAW DATA
normalization <- rma(GSE19804_RAW, background = T, normalize = T)
rma <- exprs(normalization)


# ANNOTATE DATA WITH hgu133plus.db
probeIDs <- rownames(normalization)
annotate <- select(hgu133plus2.db, keys=probeIDs, columns="SYMBOL")


# REMOVE DUPLICATED PROBES EXCEPT FOR FIRST AND LAST PROBE
annotate <- annotate[!duplicated(annotate $PROBEID), ]
rownames(annotate) <- annotate$PROBEID # change row names to PROBEIDs
merge_df <- merge(annotate, rma, by = "row.names") # merge exprs normalized and annotated data together by rownames

# remove duplicate symbols
merge_df <- merge_df[!duplicated(merge_df $SYMBOL), ]

# remove NA values
merge_df <- na.omit(merge_df)
rownames(merge_df) <- merge_df$SYMBOL # set row names of merge_df to SYMBOLs
annotated_df <- merge_df[-c(1, 2, 3)] # remove unnecessary columns


# GENE FILTERING
row_mean <- rowMeans(annotated_df) # calculate mean of each row in matrix
# filter out genes below 2nd percentile of expression distribution of data set
geneFilt <- quantile(row_mean, p=0.02)
geneFilt <- annotated_df[which(row_mean>geneFilt),]


# LIMMA ANALYSIS
# Build model preserving whether sample is cancer or normal
group <- data.frame(Tissue=GSE19804_METADATA$col2)
rownames(group) <- GSE19804_METADATA$name

# analysis models
matrix <- model.matrix(~ col2, GSE19804_METADATA)
linear_model <- lmFit(geneFilt, matrix)
eBayes_model <- eBayes(fit=linear_model)

# TOP TABLE OF TOP DEG
deg <- topTable(eBayes_model, p.value=0.05, adjust.method="fdr", sort.by="P",genelist=rownames(geneFilt), number=10000)
write.table(eBayes_model$p.value, file="DEG Sorted by P-Value.txt") # deliverable: differential expression results to a file sorted by adjusted p-value
write.table(deg, file="Significant Gene with Adjusted P-value5%.txt") # deliverable: number of significant genes with an adjusted p-value < 0.05

# TOP TABLE OF TOP 10 DEG
top10 <- topTable(eBayes_model, n=10)
filt_10 <- data.frame(geneFilt[rownames(geneFilt) %in% rownames(top10),])

# VOLCANO PLOT  
EnhancedVolcano(toptable=deg, lab=deg$ID, x="logFC",y="P.Value", xlab="Log 2 Fold Change", ylab="-Log10 b", title="Volcano Plot")

# HEAT MAP: TOP 50 DEG
pheatmap(filt_10, main="Top 10 DEG", labels_row = rownames(filt_10),  labels_col = colnames(filt_10), annotation_row=top10, annotation_col=group, annotation_colors=list(Tissue=c(Cancer="red", Normal="blue")))

top50 <- topTable(eBayes_model, number=50)
filt_50 <- data.frame(geneFilt[rownames(geneFilt) %in% rownames(top50),])
pheatmap(filt_50, main="Top 50 DEG", labels_row = rownames(filt_50),  labels_col = colnames(filt_50), annotation_row=group, annotation_col=group, annotation_colors=list(Tissue=c(Cancer="red", Normal="blue")))
