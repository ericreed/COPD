#Get list of DGE genes
baseDir<-"/Users/ericreed/Google Drive/COPD_project/"
inDirDGE<-file.path(baseDir, "Results/DGE")
inDirF<-file.path(baseDir, "Results/Ftest")

outDir<-file.path(baseDir, "Results/GeneFiltering")

#Read in DGE file
DGEout<-read.table(file.path(inDirDGE, "COPD_DGEresults.txt"), header=T)
#Get FDR adjusted p-values
DGEout$FDR<-p.adjust(DGEout$pvalue, method = "BH")
#Get list of FDR < 0.05
DGEsig<-DGEout[DGEout$FDR<=0.05,]

print(nrow(DGEsig))
# 3165 genes are differenitally expressed FDR<0.05

DGEfil<-DGEsig[abs(DGEsig$beta_copdyes)>log(1.25, 2),]
# Only 5 gave FC greater than 1.25

# Read in DGEV file
Fout<-read.table(file.path(inDirF, "COPD_DiffVarResults.txt"), header=T)
# Get FDR adjusted p-values
Fout$FDR_F<-p.adjust(Fout$GQp, method = "BH")
# Get list of FDR < 0.05
Fsig<-Fout[Fout$GQp<=0.05,]
# 19697 genes raw p-value < 0.05

Genelist<-merge(DGEsig, Fsig)
# Only 145 are differenitally expressed and have larger 

# We will use just raw DGE FDR < 0.05. 
write.table(DGEsig, file.path(outDir, "GeneFil_Results.txt"), row.names=F)


# We will use Fstat and DGE p-values.
write.table(Genelist, file.path(outDir, "GeneFil_Results_HigherVar.txt"), row.names=F)

