#loading in the metadata
data <- read.delim("the_metadata.txt", header = TRUE)
data <- data[1:120,]
#importing the raw data using the simplyaffy function
library(simpleaffy)
raw.data = ReadAffy()
raw.data


#raw data fixing
rawGSE <- exprs(raw.data)
#pca model
pc <- prcomp(rawGSE, scale = F, center = F)
pc <- as.data.frame(pc$rotation)

new_data <- cbind(pc, new_col = data['tissue'])
new_data

library(tidyverse)
Raw <- ggplot(new_data, aes(x = PC1, y = PC2, color = tissue)) + geom_point() + stat_ellipse()+ labs(title = "Raw Data PCA Plot")


#normalized plot
#background correctiona nd normalization
norData <- rma(raw.data)
#norData <- norData[1:120, ]
norData <- exprs(norData)
Normpc <- prcomp(norData,scale = F, center = F)
Normpc <- as.data.frame(Normpc$rotation)

new_data2 <- cbind(Normpc, new_col = data['tissue'])
new_data2

normPlot<- ggplot(new_data2, aes(x = PC1, y = PC2, color = tissue)) + geom_point() + stat_ellipse() + labs(title = "Normalized Data PCA Plot")

