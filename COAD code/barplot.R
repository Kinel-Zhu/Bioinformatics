path=paste("C:/Users/Administrator/Desktop/coad_mut_M.csv")
mut_data=read.csv(path,header=F)
mut_data=mut_data[which(mut_data[,5]!='MX'& mut_data[,3]!='Silent'& mut_data[,1]!='gene'),]

gene_list=as.vector(unique(mut_data[,1]))
mut_pic<-data.frame(gene_name=gene_list,m0=rep(0,length(gene_list)),m1=rep(0,length(gene_list)),sum=rep(0,length(gene_list)))
for(i in 1:length(gene_list)){
  mut_pic$m0[i]=length(which(mut_data[,1]==gene_list[i]&mut_data[,5]=='M0'))
  mut_pic$m1[i]=length(which(mut_data[,1]==gene_list[i]&mut_data[,5]=='M1'))
  mut_pic$sum[i]=length(which(mut_data[,1]==gene_list[i]))
}


mut_pic=mut_pic[order(mut_pic[,4],decreasing=T),]

colors_list<- c("orange","brown")
stage_list<-c("M0",'M1')

png(file = "C:/Users/Administrator/Desktop/barchart_coad_m.png")
values<-t(mut_pic[1:20,c(2,3)])


par(cex.axis=0.8)
barplot(values,main = "coad_M_mut",names.arg =mut_pic[1:20,1],xlab = "gene",ylab = "rates", col = colors_list)
##用ggplot theme()解决字体斜置
legend("topright", stage_list, cex = 0.6, fill = colors_list)

dev.off()


