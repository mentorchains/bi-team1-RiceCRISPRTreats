#Module 4 Group 2
#Ananya Kaushik
# Leila Li


#Loading the necessary packages
library(affy) 
library(affyPLM) 
library(hgu133plus2.db) 
library(limma) 
library(EnhancedVolcano) 
library(ggplot2)

#Obtaining the data 
raw <- ReadAffy(celfile.path = "GSE19804_RAW")

normalized = rma(raw)         #Normalizing data 
exprs_data= exprs(normalized)

#Obtaining the probe ids and the corresponding gene symbols
ID <- rownames(raw)
probe <- select(hgu133plus2.db, keys = ID,  columns = "SYMBOL")
duplicate_probeID <- probe[!duplicated(probe$PROBEID),]

#Creating the annotated data frame with row names being gene symbols
table_merge <- merge(x = duplicate_probeID, y = exprs_data, by.x = "PROBEID", by.y = "row.names")
table_merge <- table_merge[!duplicated(table_merge$SYMBOL),]
table_merge <- na.omit(table_merge)
rownames(table_merge) <- table_merge$SYMBOL
annotated <- table_merge[-c(1,2)]

#Creating the design matrix 
interest <- factor(paste(rep(c("C", "N"), each = 60), sep = ""))
design.mat <- model.matrix(~ 0 + interest)


#Next section is limma analysis
fit <- lmFit(annotated, design.mat)

make.names(c("interestlung cancer", "interestpaired normal adjacent"), unique = TRUE)
contrast.matrix <- makeContrasts(
  interestC-interestN, 
  levels = design.mat)

fit.contrast = contrasts.fit(fit, contrast.matrix)

efit <- eBayes(fit.contrast)

genes=geneNames(raw)
limma_output <- topTable(efit, n = 54675, p.value=0.05, adjust.method="fdr")


#Creating the volcano plot
EnhancedVolcano( toptable = limma_output, 
                 lab =rownames(limma_output), 
                 x = "logFC", 
                 y = "P.Value")
