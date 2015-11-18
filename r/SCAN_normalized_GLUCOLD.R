#Normalize GLUCOLD samples using SCAN

homeDir <- "C:/Users/bray/Dropbox/Boston University/Fall 2015/Challenge Project/GLUCOLD/"
dirRaw <- paste(homeDir, "GSE36221_RAW", sep="")
annotDir <- homeDir
dirFormat <- paste(homeDir, "FormattedData", sep="")

require(SCAN.UPC)

# Get list of sampGEO names
clinical <- read.table(file.path(dirFormat, "clinical.txt"), stringsAsFactors = FALSE, header = TRUE)

#Get list of CEL filenames for this path
CELfiles <- list.files(file.path(dirRaw), pattern  = "CEL.gz", full.names = TRUE)

install.packages(file.path(annotDir, "hugene10sthsentrezgprobe_19.0.0.tar.gz"), repos=NULL, type="source")

require(hugene10sthsentrezgprobe)


#Create wrapper to SCAN() function to run scan to normalize each dataset
SCANwrap <- function(GSM, CELfiles, outDir){
  celfilePath<-CELfiles[grepl(GSM, CELfiles)]
  normalized <- SCAN(celfilePath,probeSummaryPackage =  hugene10sthsentrezgprobe, outFilePath = file.path(outDir, GSM))}

GSMlist<-unique(clinical$sampGEO)

lapply(GSMlist[1:length(GSMlist)], SCANwrap, CELfiles, file.path(dirFormat, "SCAN_normalized"))
