path=paste("C:/Users/Administrator/Desktop/coad_mut_M.csv")
mut_data=read.csv(path,header=F)


##处理数据
mut_data=mut_data[which(mut_data[,5]!='MX'& mut_data[,3]!='Silent'& mut_data[,1]!='gene'),]

gene_list=as.vector(unique(mut_data[,1]))
mut_pic<-data.frame(gene_name=gene_list,m0=rep(0,length(gene_list)),m1=rep(0,length(gene_list)),sum=rep(0,length(gene_list)))
for(i in 1:length(gene_list)){
  mut_pic$m0[i]=length(which(mut_data[,1]==gene_list[i]&mut_data[,5]=='M0'))
  mut_pic$m1[i]=length(which(mut_data[,1]==gene_list[i]&mut_data[,5]=='M1'))
  mut_pic$sum[i]=length(which(mut_data[,1]==gene_list[i]))
}
mut_pic=mut_pic[order(mut_pic[,4],decreasing=T),]


##ggplot2 stack条形�?
png(file = "C:/Users/Administrator/Desktop/coad_m.png",width=574,height=336)

library(ggplot2)
mut_pic_gg<-data.frame(gene_name=c(gene_list[1:20],gene_list[1:20]),satge=c(rep('M0',20),rep('M1',20)),count=c(mut_pic$m0[1:20],mut_pic$m1[1:20]))
stage<-factor(mut_pic_gg$satge)##设定因子
ggplot(mut_pic_gg,aes(gene_name,count))+
  geom_bar(aes(fill=stage),stat = "identity",position="stack",width=0.8)+
  labs(title = "COAD_M_Mut",x=" ",y="counts")+
  theme(axis.text.x = element_text(angle = 60,hjust = 1,vjust = 0.9),
          axis.title.y =element_text(angle = 0,hjust = 1,vjust = 1,size = rel(1.1)) ,
          plot.title = element_text(hjust = 0.5, face = "plain"),
          panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"))

dev.off()
##上面为去掉灰色网�?
# ggplot(mut_pic_gg,aes(gene_name,count))+
#   geom_bar(aes(fill=stage),stat = "identity",position="stack",width=0.8)+
#   labs(title = "COAD_M_Mut")+
#   theme(axis.text.x = element_text(angle = 60,hjust = 1,vjust = 1))



