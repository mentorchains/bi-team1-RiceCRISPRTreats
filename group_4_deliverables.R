# Module 4 Deliverables
# Kelly Zhang, zhangk@brandeis.edu
# Aditi Verma, 
# Shreya Vora, 

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

# NORMALIZE RAW DATA
normalization <- rma(GSE19804_RAW, background = T, normalize = T)
rma <- exprs(normalization)

# ANNOTATE DATA WITH hgu133plus.db
probeIDs <- rownames(normalization)
annotate <- select(hgu133plus2.db, keys=probeIDs, columns="SYMBOL")

# check for 1st and last occurrence

# REMOVE DUPLICATED PROBES EXCEPT FOR FIRST AND LAST PROBE
annotate <- annotate[!duplicated(annotate $PROBEID), ]
rownames(annotate) <- annotate$PROBEID # change row names to PROBEIDs
merge_df <- merge(annotate, rma, by = "row.names") # merge exprs normalized and annotated data together by rownames

# remove duplicate symbols
merge_df <- merge_df[!duplicated(merge_df $SYMBOL), ]

# remove NA values
merge_df <- na.omit(merge_df)
rownames(merge_df) <- merge_df$SYMBOL # set row names of merge_df to SYMBOLs
normalization <- merge_df[-c(1, 2, 3)] # remove unnecessary columns


# GENE FILTERING
row_mean <- rowMeans(normalization) # calculate mean of each row in matrix
# filter out genes below 2nd percentile of expression distribution of data set
geneFilt <- quantile(row_mean, p=0.02)
geneFilt <- normalization[which(row_mean<geneFilt),]

# LIMMA ANALYSIS
# Build model preserving whether sample is cancer or normal
group <- data.frame(Tissue=GSE19804_METADATA$col2)
rownames(group) <- GSE19804_METADATA$
  
# set col names of filtered set to rownames
colnames(geneFilt) <- rownames(group)

# analysis models
matrix <- model.matrix(~ col2, GSE19804_METADATA)
linear_model <- lmFit(data.frame(geneFilt), matrix)
eBayes_model <- eBayes(fit=linear_model)

# TOP TABLE OF TOP TEN DEG
DEG <- topTable(eBayes_model, p.value=0.05, adjust.method="fdr", sort.by=)

# VOLCANO PLOT
volcanoPlot <- volcanoplot(linear_model)

# HEAT MAP
top50 <- topTable(eBayes_model, number=50)
hm_df <- data.frame(geneFilt[rownames(geneFilt) %in% rownames(top50),])
#annotateCol <- 
heatmap <- pheatmap(hm_df, annotation_col=annotateCol, cluster_rows=T)
