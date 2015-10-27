#Format phenotypes
dir<-"/Users/ericreed/Google Drive/COPD_project"

#Read in the data from "GSE37147_series_matrix.txt"
tab<-read.table(file.path(dir, "RawData", "GSE37147_series_matrix.txt"), skip=27, nrows = 43)
tab<-as.data.frame(t(tab), stringsAsFactors = FALSE)
colnames(tab)<-sub("!", "", tab[1,])
tab<-tab[-1,]
row.names(tab)<-NULL

#Keep necessary rows
tab<-tab[,c(1:2,10:20)]
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
colnames(tabOut)<-c("sampID", "sampGEO", "patientID", "included", "fev1", "fev1/fvc", "copd", "age", "smoking", "sex", "packYears", "asthma", "inhaledMedications")

#Write out clinical data file
write.table(tabOut, file.path(dir, "FormattedData", "clinical.txt"), row.names = FALSE, quote = FALSE, sep="\t")

tabIn<-read.table(file.path(dir, "FormattedData", "clinical.txt"), sep="\t", header = T)


