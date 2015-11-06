require(parallel)
# Run DGE using LM

# Set Base directory
baseDir<-"/Users/ericreed/Google Drive/COPD_project/FormattedData"

# Read in SCAN normalized expression matrix
expr<-read.table(file.path(baseDir, "NormalizedExpMat.txt"), stringsAsFactors = FALSE, header = T)

# Read in Sample File
samples<-read.table(file.path(baseDir, "clinical.txt"), sep="\t", header = TRUE, stringsAsFactors = FALSE)
samples<-samples[samples$included=="yes",]

# Adjust for these covariates
pheno<-samples[, c("sampGEO", "copd", "age", "sex", "smoking", "asthma", "inhaledMedications")]

# Check if any missing values
print(sum(is.na(pheno)))
#No--Chill

#Format the columns that are categorical as factors
pheno$copd<-factor(pheno$copd, levels = c("no", "yes"))
pheno$sex<-factor(pheno$sex)
pheno$smoking<-factor(pheno$smoking)
pheno$asthma<-factor(pheno$asthma)
pheno$inhaledMedications<-factor(pheno$inhaledMedications)

str(pheno)

#Create function to create linear model
###gene expression = copd + age + sex + smoking + asthma + inhaledMedications

expLM<-function(gene){
  exprGene<-as.data.frame(t(expr[row.names(expr)==gene,]))
  colnames(exprGene)<-"expr"
  exprGene$sampGEO<-row.names(exprGene)
  
  phenoSub<-merge(exprGene, pheno)
  
  lmOut<-summary(lm(expr ~ copd + age + sex + smoking + asthma + inhaledMedications, data = phenoSub))

  coefOut<-lmOut$coefficients["copdyes",]
  return(coefOut)}

#Run in parallel of course
geneVec<-row.names(expr)
DGEout<-mclapply(geneVec, expLM, mc.cores = 4)

DGEfram<-as.data.frame(do.call(rbind, DGEout))
DGEfram$gene<-geneVec

DGEfram<-DGEfram[,c("gene", "Estimate", "Std. Error", "t value", "Pr(>|t|)")]
colnames(DGEfram)<-c("gene", "beta_copdyes","se", "testStat", "pvalue")

write.table(DGEfram, file.path(baseDir, "COPD_DGEresults.txt"), row.names = FALSE, col.names = TRUE)

  