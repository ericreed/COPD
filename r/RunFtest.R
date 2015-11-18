require(parallel)
# Run Goldfeld-Quandt Test using LM

# Set Base directory
baseDir<-"/Users/ericreed/Google Drive/COPD_project/FormattedData"

# Read in SCAN normalized expression matrix
expr<-read.table(file.path(baseDir, "NormalizedExpMat.txt"), stringsAsFactors = FALSE, header = T)

# Read in Sample File
samples<-read.table(file.path(baseDir, "clinical.txt"), sep="\t", header = TRUE, stringsAsFactors = FALSE)
samples<-samples[samples$included=="yes",]

# Adjust for these covariates
pheno<-samples[, c("sampGEO", "copd", "age", "sex", "smoking", "fev1", "fev1.fvc")]

pheno<-pheno[!(pheno$fev1<80 & pheno$fev1.fvc>70),]

print(table(pheno$copd))
# no yes 
# 105  87 

# Check if any missing values
print(sum(is.na(pheno)))
#No--Chill

#Format the columns that are categorical as factors
pheno$copd<-factor(pheno$copd, levels = c("no", "yes"))
pheno$sex<-factor(pheno$sex)
pheno$smoking<-factor(pheno$smoking)

str(pheno)

expGQ<-function(gene){
  exprGene<-as.data.frame(t(expr[row.names(expr)==gene,]))
  colnames(exprGene)<-"expr"
  exprGene$sampGEO<-row.names(exprGene)
  
  phenoSub<-merge(exprGene, pheno)
  
  lmCOPD<-lm(expr ~ age + sex + smoking, data = phenoSub[phenoSub$copd=="yes",])
  lmControl<-lm(expr ~ age + sex + smoking, data = phenoSub[phenoSub$copd=="no",])
  
  copdSSR<-sum(lmCOPD$residuals^2)
  contSSR<-sum(lmControl$residuals^2)
  
  d1<-(nrow(phenoSub[phenoSub$copd=="yes",])-lmCOPD$rank)
  d2<-(nrow(phenoSub[phenoSub$copd=="no",])-lmControl$rank)
  
  GQ<-(copdSSR/d1)/(contSSR/d2)
  
  GQp<-pf(GQ, d1, d2, lower.tail = FALSE)
    
  FtestOut<-data.frame(gene = gene, d1 = d1, d2 = d2, GQ = GQ, GQp = GQp)
    
  return(FtestOut)}

#Run in parallel of course
geneVec<-row.names(expr)
Fout<-mclapply(geneVec, expGQ, mc.cores = 4)
Ffram<-as.data.frame(do.call(rbind, Fout))
colnames(Ffram)<-c("gene", "d1", "d2", "GQ", "GQp")

#Write table of differential variance results
write.table(Ffram, file.path(baseDir, "COPD_DiffVarResults.txt"), row.names = FALSE, col.names = TRUE)
