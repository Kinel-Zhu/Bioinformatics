 ##DESeq2
# source("http://bioconductor.org/biocLite.R")
# biocLite('clusterProfiler')
# biocLite('org.Hs.eg.db')
library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
setwd("C:/Users/Administrator/Desktop/one database")
#准备数据
P<-read.csv("TCGA-COAD-Counts.csv")
C<-read.csv("Colon.csv")
C<-C[,-2]
c_d<-data.frame(str_sub(as.character(C[,1]),1,15),C[,-1])
p_d<-data.frame(P)
colnames(p_d)<-c(rep("treat",length(p_d)))
colnames(c_d)[1]<-c("gene_id")
colnames(p_d)[1]<-c("gene_id")
exp_d<-merge(c_d,p_d,by="gene_id")
rownames(exp_d)<-exp_d[,1]
exp_d<-exp_d[,-1]
##设置因子
condition<-factor(c(rep("control",length(c_d)-1),rep("treat",length(p_d)-1)),levels = c("contol","treat"))
colData<-data.frame(row.names = colnames(exp_d),condition)
deg<-DESeqDataSetFromMatrix(exp_d,colData,design=~condition)
deg<-DESeq(deg)
res=results(deg,contrast = c("condition","control","treat"))
##筛选结果
diff_gene<-subset(res,padj<0.05 & abs(log2FoldChange)>1)
dim(diff_gene)
##转换基因名
ensembl_gene<-rownames(diff_gene)
diff_gene<-cbind(ensembl_gene,diff_gene)
colnames(diff_gene)[1]<-c("ENSEMBL")
test <- bitr(ensembl_gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
diff_name<-merge(diff_gene,test,by="ENSEMBL")
##输出
write.csv(diff_name,file="coad_DEG.csv")
