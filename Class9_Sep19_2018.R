# get lines 73-76, 80-86, 88-89, 112, 114-122

#setwd("C:/Users/Natalia/Documents/GitHub/Class_4/estrogen")
#targetsFile <- "estrogen.txt"

library(affy)
setwd("C:/Users/Natalia/Documents/GitHub/Class_4/estrogen")
targetsFile <- "estrogen.txt"
pd <- read.AnnotatedDataFrame(targetsFile,header=TRUE,sep="",row.names=1)

ER <- pData(pd)$estrogen
Time <- factor(pData(pd)$time.h)
design <- model.matrix(~ER+Time)
design

design2 <- model.matrix(~ER*Time)
design2

raw <-ReadAffy(celfile.path = "C:/Users/Natalia/Documents/GitHub/Class_4/estrogen", filenames=rownames(pData(pd)),phenoData = pd)
raw

eset <- rma(raw)

library(limma)
fit1 <- lmFit(eset, design) # fit of expression values on linerr line
fit1 <- eBayes(fit1) # Baysian algoritm to calculate P values, fold
topTable(fit1, coef=2)


fit2 <- lmFit(eset, design2)
fit2 <- eBayes(fit2)
topTable(fit2, coef=2)

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("genefilter")

# Annotation of probe Ids
library(genefilter)


## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("GEOquery")

library(GEOquery)
library(limma)
url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE33nnn/GSE33126/matrix/GSE33126_series_matrix.txt.gz"
filenm <- "GSE33126_series_matrix.txt.gz"
if(!file.exists(filenm)) download.file(url, destfile=filenm)
gse <- getGEO(filename=filenm)

# 3 Functions: varFilter, fData, anno
gse.expfilter <- varFilter(gse)
anno <- fData(gse.expfilter)
anno <- anno[,c("symbol", "Entrez_Gene_ID", "Chromosome", "Cytoband")]
fit2$genes <- anno
topTable(fit2)

# Make a volcanoplot on fit2 object....
# Volcano Plots
volcanoplot(fit2)# fit 2 is our p values

