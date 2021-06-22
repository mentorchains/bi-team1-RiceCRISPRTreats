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
mas5 <- mas5(gse)
saveRDS(mas5, "./QC/data/norm/mas5.rds")
par(mai=c(3,5,1,1,1))
boxplot(log2(exprs(mas5)), main="Boxplot mas5 Data", ylab="Probe Intensities", las=2, col=c("red","orange","yellow","green","blue","purple"))

rma <- rma(gse)
saveRDS(rma, "./QC/data/norm/rma.rds")
par(mai=c(3,5,1,1,1))
boxplot(exprs(rma), main="Boxplot rma Data", ylab="Probe Intensities", las=2, col=c("red","orange","yellow","green","blue","purple"))

gcrma <- gcrma(gse)
saveRDS(gcrma, "./QC/data/norm/gcrma.rds")
par(mai=c(3,5,1,1,1))
boxplot(exprs(gcrma), main="Boxplot gcrma Data", ylab="Probe Intensities", las=2, col=c("red","orange","yellow","green","blue","purple"))

# Create a heatmap with normalized and raw data
pheatmap <- pheatmap(cor(df))
