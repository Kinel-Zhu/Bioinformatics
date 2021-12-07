# ##ClusterProfiler转换基因名ENSEMBL->SYMBOL
# source("http://bioconductor.org/biocLite.R")
# biocLite('clusterProfiler')
# biocLite('org.Hs.eg.db')

library(clusterProfiler)
library(org.Hs.eg.db)
ensg_test<-rownames(exp_d)
test <- bitr(ensg_test, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
head(test)