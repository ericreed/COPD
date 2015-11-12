require(ConsensusClusterPlus)
require(tidyr)

#Get list of DGE genes
baseDir<-"/Users/ericreed/Google Drive/COPD_project/FormattedData"


#Read in DGE file
DGEout<-read.table(file.path(baseDir, "COPD_DGEresults.txt"), header=T)

#Get FDR adjusted p-values
DGEout$FDR<-p.adjust(DGEout$pvalue, method = "BH")

#Get list of FDR < 0.05
DGEsig<-DGEout[DGEout$FDR<=0.05,]

#Subset list for fold change > 1.5
#log(1.25, 2)
#[1] 0.584962

DGEkeep<-DGEsig[abs(DGEsig$beta_copdyes)>0.3219281,]
#This remove all results but a few, instead we will keep genes with FDR < 0.0.05

DGEsig<-DGEout[DGEout$FDR<=0.05,]

#This keeps around 2331 genes for clustering. 

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
exprCOPD<-exprCOPD[row.names(exprCOPD)%in%DGEsig$gene,]

# Filter for genes with high variability
#mads <- apply(exprCOPD,1,mad)
#exprCOPD<-exprCOPD[rev(order(mads))[1:2000],]

exprMat<-as.matrix(exprCOPD)

#Create heatmap straight up
corrdist = function(x) as.dist(1-abs(cor(t(x))))
hclust.avl = function(x) hclust(x, method="average")

require(Heatplus)
require(RColorBrewer)

colfunc <- colorRampPalette(c("blue","white", "red"))

reg = regHeatmap(exprMat, legend=2,col=colfunc, dendrogram = list(clustfun=hclust.avl, distfun=corrdist,labels = list()))
plot(reg)

title=tempdir()

resultsGene = ConsensusClusterPlus(t(exprMat), maxK=10,reps=1000,pItem=0.8,pFeature=1,clusterAlg="hc",distance="pearson",seed=1, title = "gene-level", plot = "png")

iclGene = calcICL(resultsGene)

#Clearly two samples is best
#We can trim the gene list to highest consensus

IC<-iclGene$itemConsensus
IC2<-IC[IC$k==2,]

IC2<-spread(IC2, cluster, itemConsensus)

row.names(IC2)<-IC2$item
IC2<-IC2[,!colnames(IC2)%in%c("k", "item")]
colnames(IC2)<-c("C1", "C2")
IC2$max<-apply(IC2, 1, max)

hist(IC2$max, breaks = 100)

#Remove genes with max consensus less than 0.98
ICkeep<-IC2[IC2$max>0.98,]
exprTrim<-exprMat[row.names(exprMat)%in%row.names(ICkeep),]


resultsSample = ConsensusClusterPlus(exprTrim, maxK=10,reps=1000,pItem=0.8,pFeature=1,clusterAlg="hc",distance="pearson",seed=1, title = "sample-level", plot = "png")
iclSample = calcICL(resultsSample)

#Three cluster is the way to go.  The third is only two samples large

IC3<-iclSample$itemConsensus

IC3<-IC3[IC3$k==3,]
IC3<-spread(IC3, cluster, itemConsensus)

row.names(IC3)<-IC3$item
IC3<-IC3[,!colnames(IC3)%in%c("k", "item")]
colnames(IC3)<-c("C1", "C2", "C3")

IC3$max<-apply(IC3, 1, which.max)

save(resultsGene, iclGene, resultsSample, iclSample,IC3, IC2,exprMat, exprTrim, file = "ConClustResults.RData")



















