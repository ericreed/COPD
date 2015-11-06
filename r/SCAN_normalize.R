#Normalize samples using SCAN

dirRaw<-"/Users/ericreed/Google Drive/COPD_project/RawData/GSE37147_RAW"
annotDir<-"/Users/ericreed/Google Drive/COPD_project/RawData"
dirFormat<-"/Users/ericreed/Google Drive/COPD_project/FormattedData"

require(SCAN.UPC)

#Get list of 238 subjects used in original study.  This had 10 subjects removed with lung cancer

###A total of 269 arrays from 267 samples including two samples run in duplicate were used for the generation of gene expression levels. The array data for two subjects were excluded because of sample annotation concerns, leaving a total of 265 samples. To minimize the potential con- founding effect of lung cancer, data from 19 subjects with a diagnosis of lung cancer as of January 2010 were excluded as were data from eight subjects who lacked lung function testing within 1 year of their study bronchoscopy, leaving a total of 238 samples.
####

clinical<-read.table(file.path(dirFormat, "clinical.txt"), stringsAsFactors = FALSE, header = T, sep="\t")
clinical<-clinical[ clinical$included=="yes",]
print(dim(clinical))

#Get list of CEL filenames for this path
CELfiles<-list.files(file.path(dirRaw), pattern  = "CEL.gz", full.names = TRUE)

install.packages(file.path(annotDir, "hugene10sthsentrezgprobe_19.0.0.tar.gz"), repos=NULL, type="source")

require(hugene10sthsentrezgprobe)


#Create wrapper to SCAN() function to run scan to normalize each dataset
SCANwrap<-function(GSM, CELfiles, outDir){
  celfilePath<-CELfiles[grepl(GSM, CELfiles)]
  normalized <- SCAN(celfilePath,probeSummaryPackage =  hugene10sthsentrezgprobe, outFilePath = file.path(outDir, GSM))}

GSMlist<-unique(clinical$sampGEO)

lapply(GSMlist[86:length(GSMlist)], SCANwrap, CELfiles, file.path(dirFormat, "SCAN_normalized"))

