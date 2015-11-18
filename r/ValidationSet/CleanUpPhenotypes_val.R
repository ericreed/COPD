#Format phenotypes
dir<-"/Users/ericreed/Google Drive/COPD_project"

#Read in the data from "GSE37147_series_matrix.txt"
tab<-read.table(file.path(dir, "RawData","GLUCOLD", "GSE36221_series_matrix.txt"), skip=27, nrows = 35)
tab<-as.data.frame(t(tab), stringsAsFactors = FALSE)
colnames(tab)<-sub("!", "", tab[1,])
tab<-tab[-1,]
row.names(tab)<-NULL

#Keep necessary rows
tab<-tab[,c(2,10:14)]
print(head(tab))

#Clean up the data to include just the necessary values
cleanData<-function(col){
  colSub<-unlist(strsplit(col, ": "))
  if(length(colSub)==(2*length(col))){
  colVal<-colSub[seq(2,length(colSub), by = 2)]
  return(colVal)} else {
    return(col)}}
tabOut<-as.data.frame(apply(tab, 2, cleanData))

#Clean up column names
colnames(tabOut)<-c("sampGEO", "timePoint", "PID", "age", "treatment","sex")
tabOut$sex<-as.character(tabOut$sex)

tabOut$sex[tabOut$sex=="0"]<-"F"
tabOut$sex[tabOut$sex=="1"]<-"M"

#Write out clinical data file
write.table(tabOut, file.path(dir, "FormattedData", "clinical_val.txt"), row.names = FALSE, quote = FALSE, sep="\t")

tabIn<-read.table(file.path(dir, "FormattedData", "clinical_val.txt"), sep="\t", header = T)


