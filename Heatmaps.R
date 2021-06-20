gse <- ReadAffy(celfile.path = "GSE19804_RAW/")
df <- as.data.frame(exprs(gse))

pheatmap(cor(df))
