####TCGA-Assembler


setwd("C:/Users/Administrator/Desktop/one database/TCGA-Assembler.2.0.6/TCGA-Assembler")##getwd()得到地址
# source("http://bioconductor.org/biocLite.R")
# biocLite("httr")
# biocLite("RCurl")
# biocLite("stringr")
# biocLite("HGNChelper")
# biocLite("rjson")a
library(httr)
library(RCurl)
library(stringr)
library(HGNChelper)
library(rjson)
#载入TCGA_assemble文件夹中的两个模块，其中A模块用来下载数据，B用来分析数据
source("Module_A.R")
source("Module_B.R")


##开始下数据
setwd("C:/Users/Administrator/Desktop/one database/TCGA-Assembler.2.0.6/TCGA-Assembler/COAD")
path=paste("C:/Users/Administrator/Desktop/one database/TCGA-Assembler.2.0.6/TCGA-Assembler/COAD/cases_table.csv",sep='')
Patient_ID<-read.csv(path,header=F)
vPatient_ID<-as.vector(as.array(Patient_ID[-1,2]))#####test
filename_READ_gene <- DownloadRNASeqData(cancerType = "COAD",
                                           assayPlatform = "gene_Array",
                                           tissueType = "TP",
                                           saveFolderName = "./gene",inputPatientIDs = vPatient_ID)


