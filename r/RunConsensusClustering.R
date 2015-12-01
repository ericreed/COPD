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

# Filter for genes with high variability
# mads <- apply(exprCOPD,1,mad)
# exprCOPD <- exprCOPD[rev(order(mads))[1:2000],]

exprMat<-as.matrix(exprCOPD)

# write table of genes
write.table(row.names(exprMat), file.path(baseDir, "genes4clustering.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)


#Create heatmap straight up
corrdist = function(x) as.dist(1-cor(t(x)))
hclust.avl = function(x) hclust(x, method="ward.D")

require(Heatplus)
require(RColorBrewer)

colfunc <- colorRampPalette(c("blue","white", "red"))

reg = regHeatmap(exprMat, legend=2,col=colfunc, dendrogram = list(clustfun=hclust.avl, distfun=corrdist,labels = list()))
reg$labels<-NULL
plot(reg)

setwd("~/Desktop/Challenge Project/COPD")

#resultsGene = ConsensusClusterPlus(t(exprMat), maxK=10,reps=1000,pItem=0.8,pFeature=1,clusterAlg="hc",distance="pearson",seed=1, title = "gene-level", plot = "png", , innerLinkage="ward",finalLinkage="ward")

#iclGene = calcICL(resultsGene)

#ICgene<-iclGene$itemConsensus

#write.table(ICgene, file.path(baseDir, "GeneConClust_Results.txt"), col.names = T, row.names = F)


resultsSample = ConsensusClusterPlus(exprMat, maxK=10,reps=1000,pItem=0.8,pFeature=1,clusterAlg="hc",distance="pearson",seed=1, title = "sample-level", plot = "png", , innerLinkage="ward.D",finalLinkage="ward.D")

iclSample = calcICL(resultsSample)

ICsample<-iclSample$itemConsensus

write.table(ICsample, file.path(resDir, "Clustering", "SampleConClust_Results.txt"), col.names = T, row.names = F)

#save(resultsGene, iclGene, resultsSample, iclSample, exprCOPD,exprMat, file = "ConClustResults.RData")

#Samples -- k=3
sampClust<-ICsample[ICsample$k==3,]
sampClust<-spread(sampClust, cluster, itemConsensus)
row.names(sampClust)<-sampClust$item
sampClust<-sampClust[, !colnames(sampClust)%in% c("k", "item") ]

#Assign cluster to samps
cluster<-apply(sampClust, 1, which.max)
maxes<-apply(sampClust, 1, max)
sums<-rowSums(sampClust)

sampClust$prop.max<-maxes/sums
sampClust$cluster<-cluster

write.table(sampClust, file.path(resDir, "Clustering", "sampK3_clusters.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)

#Samples -- k=2
sampClust<-ICsample[ICsample$k==2,]
sampClust<-spread(sampClust, cluster, itemConsensus)
row.names(sampClust)<-sampClust$item
sampClust<-sampClust[, !colnames(sampClust)%in% c("k", "item") ]

#Assign cluster to samps
cluster<-apply(sampClust, 1, which.max)
maxes<-apply(sampClust, 1, max)
sums<-rowSums(sampClust)

sampClust$prop.max<-maxes/sums
sampClust$cluster<-cluster

write.table(sampClust, file.path(resDir, "Clustering", "sampK2_clusters.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)


#Qualitative assessment 

#Genes -- k=5
#geneClust<-ICgene[ICgene$k==5,]
#geneClust<-spread(geneClust, cluster, itemConsensus)
#row.names(geneClust)<-geneClust$item
#geneClust<-geneClust[, !colnames(geneClust)%in% c("k", "item") ]

#Assign cluster to genes
#cluster<-apply(geneClust, 1, which.max)
#maxes<-apply(geneClust, 1, max)
#sums<-rowSums(geneClust)

#geneClust$cluster<-cluster
#geneClust$prop.max<-maxes/sums

#write.table(geneClust, file.path(baseDir, "GeneK5_clusters.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)
















