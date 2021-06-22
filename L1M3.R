# Roman Ramirez, rr8rk@virginia.edu
# Arian Veyssi, veyssi@mentorchains.com
# Level 1 Module 3, BoxPlot

# imports
library(affy)
library(simpleaffy)
library(arrayQualityMetrics)

# creating variable paths
directory <- getwd()
gseDir <- paste(directory, "/data/GSE19804_RAW", sep="")

# loading using ReadAffy into variable: myData
raw <- ReadAffy(celfile.path = gseDir)

# norm <- mas5(raw)
# df <- log2(exprs(norm))
norm <- rma(raw)
df <- exprs(norm)
boxplot(df, main="Relative Signal BoxPlot map", ylab="Relative log expression signal-RMA", las=2)