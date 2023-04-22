BiocManager::install("affy")
BiocManager::install("GEOquery")
library(affy)
library(GEOquery)
library(tidyverse)
library(dplyr)

getGEOSuppFiles("GSE148537")
install.packages("utils")
library(utils)
untar("GSE148537/GSE148537_Raw.tar",exdir = "data/")
?untar

raw_data <- ReadAffy(celfile.path = "data/")

normalized_data <- rma(raw_data)

normalized_expression <- as.data.frame(exprs(normalized_data))

#map probe ids to gene symbols
gse <- getGEO('GSE148537', GSEMatrix = TRUE)

feature_Data <- gse$GSE148537_series_matrix.txt.gz@featureData@data

feature_Data <- feature_Data[,c(1,11)]
install.packages("tibble")
library(tibble)
nornmalized_expression <-normalized_expression %>% 
  rownames_to_column(var = "ID")%>%
  inner_join (.,feature_Data ,by ="ID")



