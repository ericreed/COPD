require(ConsensusClusterPlus)
require(tidyr)
require(Biobase)

#Get list of DGE genes
baseDir<-"/Users/ericreed/Google Drive/COPD_project/FormattedData"
#Read in DGE file
resDir<-"/Users/ericreed/Google Drive/COPD_project/Results"
DGEout<-read.table(file.path(resDir, "GeneFiltering", "GeneFil_Results.txt"), header=T)


# Read in SCAN normalized expression matrix
expr<-read.table(file.path(baseDir, "NormalizedExpMat.txt"), stringsAsFactors = FALSE, header = T)

# Read in subject info
samples<-read.table(file.path(baseDir, "clinical.txt"), sep="\t", header = TRUE, stringsAsFactors = FALSE)
samples<-samples[samples$included=="yes",]

# Get list of subjects with COPD
COPDid<-samples$sampGEO[samples$copd=="yes"]

# Subset expression matrix for COPD subjects
exprCOPD<-expr[,colnames(expr)%in%COPDid]

# Subset expression matrix for COPD DGE genes
exprCOPD<-exprCOPD[row.names(exprCOPD)%in%DGEout$gene,]
exprMat<-as.matrix(exprCOPD)


# write table of genes
write.table(row.names(exprMat), file.path(baseDir, "genes4clustering.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)


#Create heatmap straight up
corrdist = function(x) as.dist(1-cor(t(x)))
hclust.avl = function(x) hclust(x, method = "ward.D2")

require(Heatplus)
require(RColorBrewer)

colfunc <- colorRampPalette(c("blue","white", "red"))

reg = regHeatmap(exprMat, legend=2, col = colfunc, dendrogram = list(clustfun=hclust.avl, distfun=corrdist,labels = list()))
reg$labels<-NULL
plot(reg)

setwd("~/Desktop/Challenge Project/COPD")

resultsSample = ConsensusClusterPlus(exprMat, maxK=10,reps=1000,pItem=0.8,pFeature=1,clusterAlg="hc",distance="pearson",seed=1, title = "sample-level", plot = "png", , innerLinkage="ward.D2", finalLinkage="ward.D2")

iclSample = calcICL(resultsSample)

ICsample<-iclSample$itemConsensus

write.table(ICsample, file.path(resDir, "Clustering", "SampleConClust_Results.txt"), col.names = T, row.names = F)

