require(reshape2)
require(ggplot2)

# Create Expression Matrix
baseDir<-"/Users/ericreed/Google Drive/COPD_project/FormattedData"
inDir<-file.path(baseDir, "GLUCOLD_SCAN_normalized")

# Subset for 138 files used in study
samples<-read.table(file.path(baseDir, "clinical_val.txt"), sep="\t", header = TRUE, stringsAsFactors = FALSE)
samples<-samples[samples$timePoint==0,]

# Create matrix

expMat<-NULL
for(i in unique(samples$sampGEO)){
  print(i)
  expSet<-read.table(file.path(inDir, i))
  colnames(expSet)<-i
  if(is.null(expMat)) expMat<-expSet else expMat<-cbind(expMat, expSet)}

write.table(expMat, file.path(baseDir, "NormalizedExpMat_val.txt"), col.names=TRUE, row.names=TRUE)

# Check out distributions of each sample
ggplot(data = melt(expMat), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable))

# Chill
  
  
  
  



