rm(list=ls()) 
Data <- read.csv("C:/Users/zjy/Desktop/coad_DEG_m.csv", header=T)###coad_DEG_m依旧是m0,m1,m1a,m1b的结果合并后的东西

down<-which(Data$log2FoldChange<0)
de = rep("Up", length(Data$log2FoldChange))
de[down] = "down"
subject<-1:length(Data$log2FoldChange)
subject<-factor(subject)
test<-data.frame(subject,Data$Phase_type,de,Data$padj)




aov2 <- aov(test$Data.padj ~ test$Data.Phase_type*test$de, data=test)
summary(aov2)

aov1 <- aov(test$Data.padj ~ test$de, data=test)
summary(aov1)

########结果截图一下吧，或者复制一下文本
aov11 <- aov(test$Data.padj ~ test$Data.Phase_type, data=test)
summary(aov11)

##ANOVA单双因素分析结果：相邻分期是影响p值的显著因素


install.packages('HH')
library(HH)
p_adj<-test$Data.padj
Up_down<-test$de
phase<-test$Data.Phase_type
png('test.png')
interaction2wt(p_adj~Up_down*phase,data = test)
dev.off()
