##source('getdata.R')

##01.连接数据库
library(RMySQL)
getdata<-function(tablename)
{
  conn <- dbConnect(MySQL(), dbname = "cancer", password="password",username="root", host="127.0.0.1", port=3306)
  b1<- dbGetQuery(conn,paste("SELECT distinct LEFT(LOWER(sample),12) FROM",tablename))
  c1<- dbGetQuery(conn,paste("SELECT distinct LEFT(gene_name,15) FROM",tablename))
  a1<- dbGetQuery(conn,paste("SELECT distinct sample FROM",tablename))
  sql<-paste("SELECT `value` FROM ",tablename," WHERE sample = '",a1[[1]][1],"' ORDER BY gene_name ASC",sep="")
  a<-dbGetQuery(conn,sql)
  for(i in c(2:lengths(a1)))
  {a2<-a1[[1]][i]
  csql<-paste("SELECT `value` FROM ",tablename," WHERE sample = '",a2,"' ORDER BY gene_name ASC",sep="")
  a3<- dbGetQuery(conn,csql)
  a<-data.frame(a,a3)
  print(paste(i,'/',lengths(a1),sep=""))}##06start
  rownames(a)<-c(c1[[1]])
  colnames(a)<-c(b1[[1]])
  dbDisconnect(conn)
  return(a)
}

filterdata<-function(stagename,stagevalue,cancername)
{
  conn <- dbConnect(MySQL(), dbname = "cancer", password="password",username="root", host="127.0.0.1", port=3306)
  #phase_type<-dbGetQuery(conn,"select distinct pathologic_m from cancer_sample_inf_stage")
  sub_sample<-dbGetQuery(conn,paste("select distinct sample from cancer_sample_inf_stage where ",stagename,"='",stagevalue,"' and cancer_name='",cancername,"'",sep=""))

  dbDisconnect(conn)
  return(sub_sample)
}

##02.取数据
a<-getdata("cancer_value_coad")##cancer_value_tgct
a_sample<-filterdata("pathologic_m","m0","COAD")
find_sub<-function(x){
  sub_index<-which(colnames(a)==x)
  return(sub_index)
}
sub_ind<-apply(a_sample,1,find_sub)
ind<-c(as.numeric(sub_ind))
ind<-ind[which(is.na(ind)!=TRUE)]
a_m0<-a[,ind]

b<-getdata("normal_value_colon")

#03.转矩阵
test1<-data.frame(rownames(a_m0),a_m0)
test2<-data.frame(rownames(b),b)
#colnames(test1)<-c(rep("treat",length(test1)))
#colnames(test2)<-c(rep("control",length(test1)))
colnames(test1)[1]<-"gene_id"
colnames(test2)[1]<-"gene_id"
exp_d<-merge(test2,test1,by="gene_id")
rownames(exp_d)<-exp_d[,1]
exp_d<-exp_d[,-1]


library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
#04.设因子、DE计算
condition<-factor(c(rep("control",length(test2)-1),rep("treat",length(test1)-1)),levels = c("control","treat"))
colData<-data.frame(row.names = colnames(exp_d),condition)
deg<-DESeqDataSetFromMatrix(exp_d,colData,design=~condition)
deg<-DESeq(deg)
res=results(deg,contrast = c("condition","control","treat"))

#05.筛选结果 
diff_gene<-subset(res,padj<0.05 & abs(log2FoldChange)>1)
dim(diff_gene)

#06.转换基因名加入分期信息
ensembl_gene<-rownames(diff_gene)
diff_gene<-cbind(ensembl_gene,diff_gene)
colnames(diff_gene)[1]<-c("ENSEMBL")
test <- bitr(ensembl_gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
diff_name<-merge(diff_gene,test,by="ENSEMBL")
phase_type<-rep("m0",dim(diff_name)[1])
diff<-data.frame(phase_type,diff_name[,1],diff_name[,8],diff_name[,2:6],sort)
colnames(diff)<-c("Phase_type","ENSEMBL","SYMBOL","baseMean","log2FoldChange" ,"lfcSE" ,"stat" , "pvalue" , "padj")
diff_m0<-diff[order(diff$pvalue),]


###相同方法计算m1，m1a,m1b,m1c，得diff_m1/m1a/m1b

###输出结果
setwd("C:/Users/Administrator/Desktop")
write.csv(diff_name,file="coad_DEG_m.csv")



