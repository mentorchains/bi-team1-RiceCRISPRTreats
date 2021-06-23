library(affyQCReport)
library(affy)
library(simpleaffy)
library(affyPLM)
library(GEOquery)
library(utilsIPEA)

gse <- ReadAffy(celfile.path = "GSE19804_RAW/")

# QCReport(gse, file = "GSE19804.pdf")

normalise <- fitPLM(gse, normalize = T, background = T)
boxplot <- boxplot(normalise)


rle_gse <- RLE(normalise, main = "RLE for GSE dataset")
stat_rle_gse <- RLE(normalise, type = "stats")
#Assuming that most genes are not changing in expression across arrays means ideally most of these RLE values will be near 0.

normalisation <- rma(gse)
rma <-exprs(normalisation)
# write.table(rma, file = "rma.txt", quote = FALSE, sep = "\t")

nuse(normalise, main = "NUSE for GSE dataset", las = 3, cex.lab = 0.1)
axis(1, cex.axis = 0.1)
