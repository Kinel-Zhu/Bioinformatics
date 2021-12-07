##00.完成目标：M分期差异性基因的可视化，故需要把m0,m1,m1a,m1b的结果合并起来

#01.清除工作区
rm(list=ls()) 

#02.准备数据
##读入已经合并好的数据
diff_gene <- read.csv("C:/Users/zjy/Desktop/coad_DEG_m.csv", header=T)
##得到m0,m1,m1a,m1b四个分期分别有多少基因，后面分组用
indicator<-duplicated(diff_gene$Phase_type)
Type_num<-table(diff_gene$Phase_type[indicator])+1
Type<-levels(diff_gene$Phase_type)
##目的：把coad_DEG_m.csv中的每个基因前面都加个因子(factor)
##做法：合并的时候按顺序合并的基因，所以因子的生成也是按顺序的
Group <-rep(Type,Type_num) #组别变量
Group <- factor(Group) #组别因子化
Attribute <- c(rep("Attribute",length(Group)))
Attribute <- factor(Attribute)
##把构建好的因子加上log2FoldChange一起放到数据框里
value <- diff_gene$log2FoldChange 
Data <- data.frame(Group=Group,Attribute=Attribute,value=value) #生成数据框


#03.计算标准误
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)##############
  datac$ci <- datac$se * ciMult
  return(datac)
}


#04.运行结果
Data_summary <- summarySE(Data, measurevar="value", groupvars=c("Group","Attribute"))

####################画图的前期数据准备已经完成################################################

#05.结果可视化
library(ggplot2)
setwd("C:/Users/zjy/Documents")
##款式1
P1<- ggplot(Data, aes(x=Group, y=value,fill=Group)) + 
  geom_violin(trim=FALSE,color="white") + #“color=”设置小提琴图的轮廓线的颜色(#trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线图
  scale_fill_manual(values = c("#FF6666", "#FFFF66","#99CC66", "#FFCC99"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=20), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Value")+xlab("") #设置x轴和y轴的标题

jpeg(file = "02Vioplot1.jpg",width =1600,height = 2000,units = "px",res =300) #结果保存
print(P1)
dev.off()

#款式2
polor2 = ggplot(Data, aes(x = Group,
                          y = value,
                          fill = Group))+
  geom_violin(alpha = 0.95, width = 1) + 
  theme_bw()+
  coord_polar()+
  scale_fill_discrete(c=100, l=100)+
  labs(title = "lint weight per ball")+
  xlab("cultivar")+
  ylab("weight(g)")
jpeg(file = "02Vioplot2.jpg",width =1600,height = 2000,units = "px",res =300) #结果保存
print(polor2)
dev.off()




###########这个要连接数据库，所以我没有测试，之前Attribute_a来自于这里，可能我把代码复制到前面的时候忘记把_a删除了= =。
##########接下来的几个都是半小提琴图


#款式3
rm(list=ls())
phase_gene<-diff_gene[which(diff_gene$Phase_type=="m0"),"ENSEMBL"]
library(RMySQL)

###使用下面定义的getdata()函数，目的：取出m0,m1,m1a,m1b分期对应的正常人的基因表达值
getdata<-function(tablename,Phase_type){
  conn <- dbConnect(MySQL(), dbname = "cancer", password="password",username="root", host="127.0.0.1", port=3306)
  a1<- diff_gene[which(diff_gene$Phase_type==Phase_type),"ENSEMBL"]
  sql<-paste("SELECT value FROM ",tablename," WHERE gene_name like  '",a1,"%' ORDER BY gene_name ASC",sep="")
  a<-dbGetQuery(conn,sql)
  value<-sum(a)/dim(a)[1]
  for(i in c(2:length(a1)))
  {a2<-a1[i]
  csql<-paste("SELECT value FROM ",tablename," WHERE gene_name like '",a2,"%' ORDER BY gene_name ASC",sep="")
  a3<- dbGetQuery(conn,csql)
  value<-c(value,sum(a3)/dim(a3)[1])
  print(paste(i,'/',length(a1),sep=""))}
  dbDisconnect(conn)
  return(value)
}
value_n_m0<-getdata("normal_value_colon","m0")
value_n_m1<-getdata("normal_value_colon","m1")
value_n_m1a<-getdata("normal_value_colon","m1a")
value_n_m1b<-getdata("normal_value_colon","m1b")
#Data <- data.frame(Group=Group,Attribute=Attribute_a,value=value) #生成数据框

###使用下面定义的getcancerdata()函数，目的：取出m0,m1,m1a,m1b分期对应的患者的基因表达值
getcancerdata<-function(Phase_type){
  conn <- dbConnect(MySQL(), dbname = "cancer", password="password",username="root", host="127.0.0.1", port=3306)
  a1<- diff_gene[which(diff_gene$Phase_type==Phase_type),"ENSEMBL"]
  sql<-paste("SELECT value FROM cancer_value_coad WHERE gene_name like  '",a1,"%' ORDER BY gene_name ASC",sep="")
  a<-dbGetQuery(conn,sql)
  value<-sum(a)/dim(a)[1]
  for(i in c(2:length(a1)))
  {a2<-a1[i]
  csql<-paste("SELECT value FROM cancer_value_coad WHERE gene_name like '",a2,"%' ORDER BY gene_name ASC",sep="")
  a3<- dbGetQuery(conn,csql)
  value<-c(value,sum(a3)/dim(a3)[1])
  print(paste(i,'/',length(a1),sep=""))}
  dbDisconnect(conn)
  return(value)
  
}
value_c_m0<-getcancerdata("m0")
value_c_m1<-getcancerdata("m1")
value_c_m1a<-getcancerdata("m1a")
value_c_m1b<-getcancerdata("m1b")

####按顺序组合正常人与患者的基因表达数据
value<-c(value_c_m0,value_n_m0,value_c_m1,value_n_m1,value_c_m1a,value_n_m1a,value_c_m1b,value_n_m1b)
#rm(list=ls())
#exp_value <- read.csv("C:/Users/Administrator/Desktop/coad_cn_m.csv", header=T)
num_m0<-length(value_c_m0)
num_m1<-length(value_c_m1)
num_m1a<-length(value_c_m1a)
num_m1b<-length(value_c_m1b)
Group <-rep(c('m0','m0','m1','m1','m1a','m1a','m1b','m1b'),c(num_m0,num_m0,num_m1,num_m1,num_m1a,num_m1a,num_m1b,num_m1b)) #组别变量
Group <- factor(Group) #组别因子化
Attribute <- rep(c('Attribute1','Attribute2','Attribute1','Attribute2','Attribute1','Attribute2','Attribute1','Attribute2'),c(num_m0,num_m0,num_m1,num_m1,num_m1a,num_m1a,num_m1b,num_m1b))
Attribute <- factor(Attribute)
exp_value <- data.frame(Group=Group,Attribute=Attribute,value=value) #生成数据框
#Data <- data.frame(Group=exp_value$Group,Attribute=exp_value$Attribute,value=exp_value$value) #生成数据框
Data_summary <- summarySE(exp_value, measurevar="value", groupvars=c("Group","Attribute"))

GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, 
                      xminv = x - violinwidth * (x - xmin), 
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}



P3 <- ggplot(data=Data, aes(x=Group, y=value,fill=Attribute)) + 
  geom_split_violin(trim=FALSE,color="white") + #绘制分半的小提琴图
  geom_point(data = Data_summary,aes(x=Group, y=value),pch=19,position=position_dodge(0.9),size=1.5)+ #绘制均值为点图
  geom_errorbar(data = Data_summary,aes(ymin = value-ci, ymax=value+ci,x=Group), #误差条表示95%的置信区间
                width=0.1,  #误差条末端短横线的宽度
                position=position_dodge(0.9), 
                color="black",
                alpha = 0.7,
                size=0.5) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=20), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Value")+xlab("") #设置x轴和y轴的标题

P3

jpeg(file = "results_Value_3.jpg",width =1600,height = 2000,units = "px",res =300) #结果保存
print(P3)
dev.off()

##款式4
exp_value <- read.csv("C:/Users/Administrator/Desktop/coad_cn_m.csv", header=T)
Data <- data.frame(Group=exp_value$Group,Attribute=exp_value$Attribute,value=exp_value$value) #生成数据框
Data_summary <- summarySE(Data, measurevar="value", groupvars=c("Group","Attribute"))

P4<-ggplot(Data_summary,aes(x=Group, y=value, fill=Attribute)) + #“fill=”设置填充颜色依据Attribute指定
  geom_point(aes(x=Group, y=value),pch=19,position=position_dodge(0.9),size=2.5)+ #绘制均值为点图
  geom_bar(stat = "identity",position = "dodge",alpha = 0.7) + #绘制条形图
  
  #如果误差条想表示标准差：请设置 ymin = value-sd, ymax=value+sd
  #如果误差条想表示标准误：请设置 ymin = value-se, ymax=value+se
  geom_errorbar(aes(ymin = value-ci, ymax=value+ci), #误差条表示95%的置信区间
                width=0.1, #误差条末端短横线的宽度
                position=position_dodge(0.9), 
                color="black",
                alpha = 0.7,
                size=0.5) +
  # scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ #设置填充颜色
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=20), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Value")+xlab("")+ #设置x轴和y轴的标题
  scale_fill_discrete(name="Attribute",breaks=c("Attribute_1","Attribute_2"),labels=c("cancer","normal"))

jpeg(file = "results_Value_4.jpg",width =1600,height = 2000,units = "px",res =300) #结果保存
print(P4)
dev.off()




##参考网址：https://blog.csdn.net/zhouhucheng00/article/details/86082760
##参考网址：https://www.jianshu.com/p/8235bb92ffa1


