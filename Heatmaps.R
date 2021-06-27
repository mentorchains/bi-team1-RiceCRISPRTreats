# Module 3 Deliverables - Heat Map
# Kelly Zhang
# Leila Li

# Download the following packages
library("affy")
library("affyPLM")
library("simpleaffy")
library("arrayQualityMetrics")
library("affyQCReport")
library("sva")
library("ggplot2")
library("pheatmap")

gse <- ReadAffy(celfile.path = "GSE19804_RAW/")
df <- as.data.frame(exprs(gse))

# Normalize data
rma <- rma(gse)
rma_exprs <- exprs(rma)
saveRDS(rma, "./QC/data/norm/rma.rds")
par(main=c(3,5,1,1,1))
boxplot(exprs(rma), main="Boxplot rma Data", ylab="Probe Intensities", las=2, col=c("red","orange","yellow","green","blue","purple"))

# Create a heatmap with normalized and raw data
corrltn <- cor(df)
pheatmap(corrltn, annotation_row=df, annotation_col=df)
