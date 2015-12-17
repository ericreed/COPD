require(tidyr)

# Set Base directory
baseDir<-"/Users/ericreed/Google Drive/COPD_project/FormattedData"

#Directory with clusters
ClusDir<-"/Users/ericreed/Google Drive/COPD_project/Results/Clustering"

#Directory to output gene signatures
SigDir<-"/Users/ericreed/Google Drive/COPD_project/Results/GeneSignatures"

# Read in the mother fucking expression matrices

#Training Data
eTrain<-read.table(file.path(baseDir, "NormalizedExpMat.txt"), stringsAsFactors = FALSE, header = T)

#Validation Data
eVal<-expr<-read.table(file.path(baseDir, "NormalizedExpMat_val.txt"), stringsAsFactors = FALSE, header = T)

#Subset the data

#List of relavent genes
geneList<-read.table(file.path(baseDir, "genes4clustering.txt"), stringsAsFactors = FALSE, header = T)

eTrain<-eTrain[rownames(eTrain)%in%geneList[,1],]
eVal<-eVal[rownames(eVal)%in%geneList[,1],]

#List of COPD subjects in training data

ClustRes<-read.table(file.path(ClusDir, "SampleConClust_Results.txt"), header = T, stringsAsFactors = FALSE)
eTrain<-eTrain[,colnames(eTrain)%in%ClustRes$item]

ClustResVal<-read.table(file.path(ClusDir, "SampleConClust_Results_val.txt"), header = T, stringsAsFactors = FALSE)

kGetFunc<-function(k, ClustRes){
  ClustSub<-ClustRes[ClustRes$k==k, c("item", "cluster", "itemConsensus")]
  ClustSub<-spread(ClustSub, cluster, itemConsensus)
  rownames(ClustSub)<-ClustSub$item
  ClustSub<-ClustSub[,!colnames(ClustSub)%in%"item"]
  ClustSub$cluster<-apply(ClustSub, 1, which.max)
  
  ClustOut<-data.frame(sample = row.names(ClustSub), cluster = ClustSub$cluster)
  return(ClustOut)}

getGeneSig<-function(kOut, expr){
  expr<-expr[,colnames(expr)%in%kOut$sample]
  expr<-expr[,order(match(colnames(expr), kOut$sample))]
  
  geneT<-function(genei, k = i){
    geneSub<-unlist(expr[row.names(expr)==genei,])
    tOut<-t.test(geneSub~outFac)
    pOut<-tOut$p.value
    tOut<-tOut$statistic
    return(data.frame(cluster = k, gene = genei, t = tOut, p = pOut))}
  
  kList<-list()
  for(i in unique(kOut$cluster)){
    outFac<-factor(kOut$cluster==i, levels = c(TRUE, FALSE))
    outSig<-lapply(rownames(expr), geneT, k=i)
    outSig<-do.call(rbind, outSig)
    kList[[i]]<-outSig}
  kList<-do.call(rbind, kList)
  return(kList)}


getGeneList<-function(k, kD){
  kSub<-kD[kD$cluster==k,]
  kSub$FDR<-p.adjust(kSub$p, "BH")
  kSub<-kSub[kSub$t>0,]
  kSub<-kSub[kSub$FDR<0.001,]
  return(kSub)}

fish<-function(list1, list2, geneAll){
  xIn<-geneAll%in%list1*1
  yIn<-geneAll%in%list2*1
  tab<-matrix(table(xIn, yIn), 2, 2)
  FT<-fisher.test(tab, alternative = "greater")
  pval<-FT$p.value
  return(pval)}

outFunc <- function(tK, vK){

init <- kGetFunc(tK, ClustRes)
initT <- getGeneSig(init, eTrain)
initG <- lapply(unique(initT$cluster), getGeneList, initT)
initG <- do.call(rbind, initG)
initG$gene<-as.character(initG$gene)

val <- kGetFunc(vK, ClustResVal)
valT <- getGeneSig(val, eVal)
valG <- lapply(unique(valT$cluster), getGeneList, valT)
valG <- do.call(rbind, valG)
valG$gene<-as.character(valG$gene)

pMat<-matrix(NA, length(unique(initG$cluster)), length(unique(valG$cluster)))

for(i in 1:length(unique(initG$cluster))){
  for(j in 1:length(unique(valG$cluster))){
    initSub<-initG[initG$cluster==i,]
    list1<-initSub$gene
    
    valSub<-valG[valG$cluster==j,]
    list2<-valSub$gene
    pMat[i,j]<-fish(list1, list2, rownames(eTrain))
    
    pMat<-matrix(p.adjust(pMat, method = "BH"), length(unique(initG$cluster)), length(unique(valG$cluster)))
    rownames(pMat)<-paste("t", 1:length(unique(initG$cluster)), sep="")
    colnames(pMat)<-paste("v", 1:length(unique(valG$cluster)), sep="")}}

return(list(pMat, initG, valG))}

p23<-outFunc(2, 3)

names(p23)<-c("pMatrix", "trainGeneSigs", "valGeneSigs")

write.csv(p23$trainGeneSigs, file.path(SigDir, "trainGeneSigs_k2.csv"), row.names=FALSE)
write.csv(p23$trainGeneSigs, file.path(SigDir, "valGeneSigs_k3.csv"), row.names=FALSE)



library(hugene10sttranscriptcluster.db)
annodb <- "hugene10sttranscriptcluster.db"

ens <- unlist(as.list(hugene10sttranscriptclusterENSEMBL2PROBE))
ens<-data.frame(AddID = ens, Ens = names(ens))

symbols <- unlist(as.list(hugene10sttranscriptclusterSYMBOL))
symbols <- data.frame(AddID = names(symbols), gene = symbols)


map<-read.table("/Users/ericreed/Google Drive/COPD_project/RawData/hugene10st_Hs_ENTREZG_19.0.0/hugene10st_Hs_ENTREZG_mapping.txt", sep="\t", header = T, stringsAsFactors = FALSE)
map<-map[, c(1,7)]
colnames(map)<-c("ID", "AddID")
map<-unique(map)

####p23

p23T<-p23[[2]]
p23V<-p23[[3]]

p23T2<-p23T$gene[p23T$cluster==2]
p23V1<-p23V$gene[p23V$cluster==1]
p23GS2<-p23T2[p23T2%in%p23V1]


p23T1<-p23T$gene[p23T$cluster==1]
p23V2<-p23V$gene[p23V$cluster==2]
p23GS1<-p23T1[p23T1%in%p23V2]

p23GS2<-map[map$ID%in%p23GS2,]
p23GS2<-p23GS2[!duplicated(p23GS2$ID),]


p23GS1<-map[map$ID%in%p23GS1,]
p23GS1<-p23GS1[!duplicated(p23GS1$ID),]

C2sig<-merge(p23GS2, ens, all.x=TRUE)
C2sig<-C2sig[!duplicated(C2sig$ID),]

C1sig<-merge(p23GS1, ens, all.x=TRUE)
C1sig<-C1sig[!duplicated(C1sig$ID),]

C2sig<-merge(C2sig, symbols, all.x=TRUE)
C2sig<-merge(C2sig, p23T, by.x = "ID", by.y = "gene")

C1sig<-merge(C1sig, symbols, all.x=TRUE)
C1sig<-merge(C1sig, p23T, by.x = "ID", by.y = "gene")

write.csv(C1sig, "Cluster1_sig_p23.csv", row.names = FALSE, na="")
write.csv(C2sig, "Cluster2_sig_p23.csv", row.names = FALSE,  na="")

####p33

p33T<-p33[[2]]
p33V<-p33[[3]]

p33T2<-p33T$gene[p33T$cluster==2]
p33V1<-p33V$gene[p33V$cluster==1]
p33GS2<-p33T2[p33T2%in%p33V1]


p33T3<-p33T$gene[p33T$cluster==3]
p33V2<-p33V$gene[p33V$cluster==2]
p33GS3<-p33T3[p33T3%in%p33V2]

p33GS2<-map[map$ID%in%p33GS2,]
p33GS2<-p33GS2[!duplicated(p33GS2$ID),]


p33GS3<-map[map$ID%in%p33GS3,]
p33GS3<-p33GS3[!duplicated(p33GS3$ID),]

C2sig<-merge(p33GS2, ens, all.x=TRUE)
C2sig<-C2sig[!duplicated(C2sig$ID),]
C2sig<-merge(C2sig, p23T, by.x = "ID", by.y = "gene")

C3sig<-merge(p33GS3, ens, all.x=TRUE)
C3sig<-C3sig[!duplicated(C3sig$ID),]
C3sig<-merge(C3sig, p23T, by.x = "ID", by.y = "gene")

symbols <- unlist(as.list(hugene10sttranscriptclusterSYMBOL))
symbols <- data.frame(AddID = names(symbols), gene = symbols)

C2sig<-merge(C2sig, symbols, all.x=TRUE)
C2sig<-merge(C2sig, p33T, by.x = "ID", by.y = "gene")

C3sig<-merge(C3sig, symbols, all.x=TRUE)
C3sig<-merge(C3sig, p33T, by.x = "ID", by.y = "gene")

write.csv(C2sig, "Cluster2_sig_p33.csv", row.names = FALSE, na="")
write.csv(C3sig, "Cluster3_sig_p33.csv", row.names = FALSE,  na="")






