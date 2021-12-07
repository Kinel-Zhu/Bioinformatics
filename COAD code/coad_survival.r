#exp_mut有统计学差异有几个
#

#################################other_survival_analysis##########################################
library(survival)
clinical<-read.csv('C:/Users/Administrator/Desktop/COAD_project/coad/COAD.merged_only_clinical_clin_format.csv')
#ccc=c("M0","M1", "N0", "N1", "N2", "T1", "T2", "T3", "T4", "stage1", "stage2", "stage3", "stage4")
ccc=c("M1","N1", "N2", "T1", "T2","T4", "stage1", "stage2", "stage3", "stage4")
res<-list()
exp_mut1=exp_mut[c(2,4,5,6,7,9,10,11,12,13)]
pv1=c()
for(j in 1:length(ccc)){
pa=paste("C:/Users/Administrator/Desktop/COAD_project/coad/mutation/coad_mut_",ccc[j],".csv",sep='')
mut=read.csv(pa)
pdf(paste("C:/Users/Administrator/Desktop/COAD_project/coad/exp_mut/coad_",ccc[j],"_sur.pdf",sep=""),38,40)
par(mfrow=c(7,6))
gee=exp_mut1[[j]]

res1=c()
for(o in 1:length(gee)){
mut_g<-paste(gee[o],'_mut',sep='')
n_g<-paste('n',mut_g,sep='')
mut_p53<-mut[which(as.character(mut[,1])==gee[o] ),4]###TP3突变的sample
mut_samp<-as.character(mut[,4])####all sample
mut_p53<-unique(substr(mut_p53,1,12))####TP3突变的sample的bcr的 1-12位☆
mut_samp<-unique(substr(mut_samp,1,12))
np53=setdiff(mut_samp,mut_p53)##求向量x与向量y中不同的元素(只取x中不同的元素)###没有发生TP53突变了的sample
l=which(clinical[,1]=='patient.bcr_patient_barcode')
c_name=toupper(unlist(clinical[l,]))##大写☆clinical中所有的sample
f=dim(clinical)[2]-1##列数☆clinical中所有的sample数
sum_res=rep('not_clear',f)
m01=0
m01=intersect(mut_p53,c_name)###TP53发生突变的sample与clinical的sample做交集

for(k in 1:length(m01)){
p=0
p=which(c_name==m01[k])
sum_res[p]=mut_g
}
m02=0
m02=intersect(np53,c_name)
for(k in 1:length(m02)){
p=0
p=which(c_name==m02[k])
sum_res[p]=n_g
}
###↑，只有mut文件里的存在的sample才可以取消not_clear


if(length(m01)>3){##########m02要不要大于3
pp="Probability of disease free survival"

ind_keep <- grep(".*days_to_new_tumor_event_after_initial_treatment.*",clinical[,1])##grep查询函数☆
# Since there are numerous follow ups, let's collapse them together and keep the first value (the higher one) if more than one is available
new_tum <- as.matrix(clinical[ind_keep,])##########☆
new_tum=new_tum[,-1]#############
new_tum_collapsed <- c()
for (a in 1:dim(new_tum)[2]){##dim(new_tum)[2]样本数
  if(sum(is.na(new_tum[,a])) < dim(new_tum)[1]){##每个样本对于那6个days_to_new_tumor_event_after_initial_treatment至少存在一个数值
    m <- max(as.numeric(new_tum[,a]),na.rm=T)############☆接着上面的注释，取6个当中最大的数值，没有就NA
    new_tum_collapsed <- c(new_tum_collapsed,m)
  } else {
    new_tum_collapsed <- c(new_tum_collapsed,"NA")
  }
}
# do the same to death
ind_keep <- grep(".*days_to_death.*",clinical[,1])##########☆
death <- as.matrix(clinical[ind_keep,])##########☆
death=death[,-1]
death_collapsed <- c()
for (a in 1:dim(death)[2]){
  if(sum(is.na(death[,a])) < dim(death)[1]){
    m <- max(as.numeric(death[,a]),na.rm=T)##########☆sample存活最久的death值
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,"NA")
  }
}
# and days last follow up here we take the most recent which is the max number
ind_keep <- grep(".*days_to_last_followup.*",clinical[,1])
fl <- as.matrix(clinical[ind_keep,])
fl=fl[,-1]
fl_collapsed <- c()
for (a in 1:dim(fl)[2]){
  if(sum(is.na(fl[,a])) < dim(fl)[1]){
    m <- max(as.numeric(fl[,a]),na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,"NA")
  }
}
# and put everything together
all_clin <- data.frame(new_tum_collapsed,death_collapsed,fl_collapsed)
colnames(all_clin) <- c("new_tumor_days", "death_days", "followUp_days")
# create vector with time to new tumor containing data to censor for new_tumor
all_clin$new_time <- c()##你可以理解为要在all_clin里加一个这样的变量，$提取
for (i in 1:length(as.numeric(as.character(all_clin$new_tumor_days)))){
  all_clin$new_time[i] <- ifelse(is.na(as.numeric(as.character(all_clin$new_tumor_days))[i]),
                    as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$new_tumor_days))[i])###大概就是new_tumor_days为NA，就用followup_days代替，放在new_time里
}
# create vector time to death containing values to censor for death
all_clin$new_death <- c()
for (i in 1:length(as.numeric(as.character(all_clin$death_days)))){
  all_clin$new_death[i] <- ifelse(is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                 as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
}
# create vector for death censoring
#all_clin$death_event <- ifelse(clinical[which(clinical[,1]=="patient.vital_status"),] == "alive", 0,1)##活着0死去1，获取死亡事件
#finally add row.names to clinical
#rownames(all_clin) <- clinical[clinical[which(clinical[,1]=="patient.bcr_patient_barcode"),]###添加barcode
death_event1<-ifelse(clinical[which(clinical[,1]=="patient.vital_status"),] == "alive", 0,1)
death_event1=death_event1[-1]
all_clin$death_event<-death_event1
rownames(all_clin) <-c_name[-1] 

event_p53<-sum_res
ind_clin<-which(event_p53!="not_clear")
msam_num=length(which(event_p53==mut_g))
nsam_num=length(which(event_p53==n_g))


s <- survfit(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_p53[ind_clin])
s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_p53[ind_clin]), error = function(e) return(NA))
pv <- ifelse(is.na(s1),next,(round(1 - pchisq(s1$chisq, length(s1$n) - 1),3)))[[1]]##pv用列存储,卡方分布：qchisq(p, df=N,ncp=0),四舍五入round,
pv1=c(pv1,pv)#################pv存在pv1中，共609个，小于0.01的有11个
if(pv<0.05){
pv=paste('p=',pv)

if((o%%6)==1){####????????六个一组六个一组？

ylab="Probability of overall survival"

}else{
ylab=""

}
plot(s, lwd=1,main=gee[0],cex.lab=1.5, mark.time=T,xlab="Overall survival time (days)",ylab=ylab,col=c("red","blue"));legend('topright',pv,bty="n",cex=1);legend("bottomleft",c(paste('n=',nsam_num,sep=''),paste('n=',msam_num,sep='')),bty="n",cex=1.5,fill=c("blue","red"), horiz=F);
#legend("bottomleft",inset=.05,c("CDH1 wildtype","CDH1 mutated"),bty="n",cex=1,fill=c("red","blue"), horiz=F)

## create vector for new tumor censoring

line_w<-paste(gee[o],'-wildtype',sep='')
line_m<-paste(gee[o],'-mutated',sep='')
plot.new()
legend("left",cex=2,inset=0.05,c(line_w,line_m),bty="n",fill=c("blue","red"), horiz=F)

res1=c(res1,gee[o])
}
}
# plot.new()
# legend("left",inset=0.05,c("Immune genes wildtype","Immune genes mutated"),bty="n",cex=1.5,fill=c("blue","red"), horiz=F)
#dev.off()

}
res[[j]]=res1##有显著性差异的用res存储，总的用exp_mut存储

dev.off()
}

####################################M0N0T3_survival############################################
ccc=c("M0","N0","T3")
rres<-list()
exp_mut1=exp_mut[c(1,3,8)]
pv1=c()
for(j in 1:length(ccc)){
pa=paste("C:/Users/Administrator/Desktop/COAD_project/coad/mutation/coad_mut_",ccc[j],".csv",sep='')
mut=read.csv(pa)
pdf(paste("C:/Users/Administrator/Desktop/COAD_project/coad/exp_mut/",ccc[j],"_sur.pdf",sep=""),38,40)
par(mfrow=c(7,6))
gee=exp_mut1[[j]]
res1=c()
for(o in 1:length(gee)){
mut_g<-paste(gee[o],'_mut',sep='')
n_g<-paste('n',mut_g,sep='')
mut_p53<-mut[which(as.character(mut[,1])==gee[o] ),4]###TP3突变的sample
mut_samp<-as.character(mut[,4])####all sample
mut_p53<-unique(substr(mut_p53,1,12))####TP3突变的sample的bcr的 1-12位☆
mut_samp<-unique(substr(mut_samp,1,12))
np53=setdiff(mut_samp,mut_p53)##求向量x与向量y中不同的元素(只取x中不同的元素)###没有发生TP53突变了的sample
l=which(clinical[,1]=='patient.bcr_patient_barcode')
c_name=toupper(unlist(clinical[l,]))##大写☆
f=dim(clinical)[2]-1##列数☆
sum_res=rep('not_clear',f)
m01=0
m01=intersect(mut_p53,c_name)###TP53发生突变的sample与clinical的sample做交集
for(k in 1:length(m01)){
p=0
p=which(c_name==m01[k])
sum_res[p]=mut_g
}
m02=0
m02=intersect(np53,c_name)
for(k in 1:length(m02)){
p=0
p=which(c_name==m02[k])
sum_res[p]=n_g
}
if(length(m01)>3){##########m02要不要大于3
pp="Probability of disease free survival"
ind_keep <- grep(".*days_to_new_tumor_event_after_initial_treatment.*",clinical[,1])##grep查询函数☆
# Since there are numerous follow ups, let's collapse them together and keep the first value (the higher one) if more than one is available
new_tum <- as.matrix(clinical[ind_keep,])##########☆
new_tum=new_tum[,-1]#############
new_tum_collapsed <- c()
for (a in 1:dim(new_tum)[2]){
  if(sum(is.na(new_tum[,a])) < dim(new_tum)[1]){
    m <- max(as.numeric(new_tum[,a]),na.rm=T)############☆
    new_tum_collapsed <- c(new_tum_collapsed,m)
  } else {
    new_tum_collapsed <- c(new_tum_collapsed,"NA")
  }
}
# do the same to death
ind_keep <- grep(".*days_to_death.*",clinical[,1])##########☆
death <- as.matrix(clinical[ind_keep,])##########☆
death=death[,-1]
death_collapsed <- c()
for (a in 1:dim(death)[2]){
  if(sum(is.na(death[,a])) < dim(death)[1]){
    m <- max(as.numeric(death[,a]),na.rm=T)##########☆sample存活最久的death值
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,"NA")
  }
}
# and days last follow up here we take the most recent which is the max number
ind_keep <- grep(".*days_to_last_followup.*",clinical[,1])
fl <- as.matrix(clinical[ind_keep,])
fl=fl[,-1]
fl_collapsed <- c()
for (a in 1:dim(fl)[2]){
  if(sum(is.na(fl[,a])) < dim(fl)[1]){
    m <- max(as.numeric(fl[,a]),na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,"NA")
  }
}
# and put everything together
all_clin <- data.frame(new_tum_collapsed,death_collapsed,fl_collapsed)
colnames(all_clin) <- c("new_tumor_days", "death_days", "followUp_days")
# create vector with time to new tumor containing data to censor for new_tumor
all_clin$new_time <- c()##你可以理解为要在all_clin里加一个这样的变量，$提取
for (i in 1:length(as.numeric(as.character(all_clin$new_tumor_days)))){
  all_clin$new_time[i] <- ifelse(is.na(as.numeric(as.character(all_clin$new_tumor_days))[i]),
                    as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$new_tumor_days))[i])###大概就是new_tumor_days为NA，就用followup_days代替，放在new_time里
}
# create vector time to death containing values to censor for death
all_clin$new_death <- c()
for (i in 1:length(as.numeric(as.character(all_clin$death_days)))){
  all_clin$new_death[i] <- ifelse(is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                 as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
}
# create vector for death censoring
#all_clin$death_event <- ifelse(clinical[which(clinical[,1]=="patient.vital_status"),] == "alive", 0,1)##活着0死去1，获取死亡事件
#finally add row.names to clinical
#rownames(all_clin) <- clinical[clinical[which(clinical[,1]=="patient.bcr_patient_barcode"),]###添加barcode
death_event1<-ifelse(clinical[which(clinical[,1]=="patient.vital_status"),] == "alive", 0,1)
death_event1=death_event1[-1]
all_clin$death_event<-death_event1
rownames(all_clin) <-c_name[-1] 
event_p53<-sum_res
ind_clin<-which(event_p53!="not_clear")
msam_num=length(which(event_p53==mut_g))
nsam_num=length(which(event_p53==n_g))
s <- survfit(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_p53[ind_clin])
s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(all_clin$new_death))[ind_clin],all_clin$death_event[ind_clin])~event_p53[ind_clin]), error = function(e) return(NA))
pv <- ifelse(is.na(s1),next,(round(1 - pchisq(s1$chisq, length(s1$n) - 1),3)))[[1]]##pv用列存储,卡方分布：qchisq(p, df=N,ncp=0),四舍五入round,
pv1=c(pv1,pv)#################pv存在pv1中，共609个，小于0.01的有11个
if(pv<0.05){
pv=paste('p=',pv)
if((o%%6)==1){####????????六个一组六个一组？
ylab="Probability of overall survival"
}else{
ylab=""
}
plot(s, lwd=1,main=gee[0],cex.lab=1.5, mark.time=T,xlab="Overall survival time (days)",ylab=ylab,col=c("red","blue"));legend('topright',pv,bty="n",cex=1);legend("bottomleft",c(paste('n=',nsam_num,sep=''),paste('n=',msam_num,sep='')),bty="n",cex=1.5,fill=c("blue","red"), horiz=F);
#legend("bottomleft",inset=.05,c("CDH1 wildtype","CDH1 mutated"),bty="n",cex=1,fill=c("red","blue"), horiz=F)
## create vector for new tumor censoring
line_w<-paste(gee[o],'-wildtype',sep='')
line_m<-paste(gee[o],'-mutated',sep='')
plot.new()
legend("left",cex=2,inset=0.05,c(line_w,line_m),bty="n",fill=c("blue","red"), horiz=F)
res1=c(res1,gee[o])
}
}
# plot.new()
# legend("left",inset=0.05,c("Immune genes wildtype","Immune genes mutated"),bty="n",cex=1.5,fill=c("blue","red"), horiz=F)
#dev.off()
}
rres[[j]]=res1##有显著性差异的用res存储，总的用exp_mut存储
dev.off()
}

##############################################other_sur_high/low_exp########################################################################
bbb=c("M","N","T","stage")
kk=c(1,2,3,4)
c=read.csv("C:/Users/Administrator/Desktop/COAD_project/coad/coad_exp_nor_sym.csv",header=F)
num_low_exp=c()
for(i in 1:length(bbb)){
pa=paste("C:/Users/Administrator/Desktop/COAD_project/coad/exp_sym/",bbb[i],"_exp.csv",sep='')######stage来测试
b=read.csv(pa,header=F)
for(j in 1:kk[i]){
p=sum(kk[1:i-1])+j
##########nor/M/N/T/S的gene序列都是一一对应的
nor_exp=subset(c,c[,1]%in%res[[p]])
pat_exp=subset(b,b[,1]%in%res[[p]])
nor_exp=as.matrix(nor_exp)#dim 6 42
pat_exp=as.matrix(pat_exp)#dim 6 90
#######nor_exp/pat_exp[,1]gene顺序相同
if(dim(pat_exp)[1]!=0){
if(dim(pat_exp)[1]==1){
result=ifelse(mean(as.numeric(nor_exp[-1]))<mean(as.numeric(pat_exp[-1])),1,0)
num_low_exp=c(num_low_exp,sum(result==0))
ress=pat_exp[,1]
names(ress)<-NULL
names(ress)<-result
res[[p]]<-ress
}else{
nor_exp1=apply(nor_exp[,-1],2,as.numeric)#dim 6 41
pat_exp1=apply(pat_exp[,-1],2,as.numeric)
nor_exp1_mean=apply(nor_exp1,1,mean)
pat_exp1_mean=apply(pat_exp1,1,mean)
result=ifelse(nor_exp1_mean<pat_exp1_mean,1,0)
num_low_exp=c(num_low_exp,sum(result==0))
ress=pat_exp[,1]
names(ress)<-NULL
names(ress)<-result
res[[p]]<-ress
}
}
}
}

##########################################################M0N0T3_sur_high/low_exp########################################################
bbb=c("M","N","T")
c=read.csv("C:/Users/Administrator/Desktop/COAD_project/coad/coad_exp_nor_sym.csv",header=F)
num_low_exp=c()
for(i in 1:length(bbb)){
pa=paste("C:/Users/Administrator/Desktop/COAD_project/coad/exp_sym/",bbb[i],"_exp.csv",sep='')######stage来测试
b=read.csv(pa,header=F)
nor_exp=subset(c,c[,1]%in%rres[[i]])
pat_exp=subset(b,b[,1]%in%rres[[i]])
nor_exp=as.matrix(nor_exp)#dim 6 42
pat_exp=as.matrix(pat_exp)#dim 6 90
#######nor_exp/pat_exp[,1]gene顺序相同
if(dim(pat_exp)[1]==1){
result=ifelse(mean(as.numeric(nor_exp[-1]))<mean(as.numeric(pat_exp[-1])),1,0)
num_low_exp=c(num_low_exp,sum(result==0))
ress=pat_exp[,1]
names(ress)<-NULL
names(ress)<-result
res[[i]]<-ress
}else{
nor_exp1=apply(nor_exp[,-1],2,as.numeric)#dim 6 41
pat_exp1=apply(pat_exp[,-1],2,as.numeric)
nor_exp1_mean=apply(nor_exp1,1,mean)
pat_exp1_mean=apply(pat_exp1,1,mean)
result=ifelse(nor_exp1_mean<pat_exp1_mean,1,0)
num_low_exp=c(num_low_exp,sum(result==0))
ress=pat_exp[,1]
names(ress)<-NULL
names(ress)<-result
rres[[i]]<-ress
}
}

