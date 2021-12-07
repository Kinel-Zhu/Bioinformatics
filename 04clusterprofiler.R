rm(list=ls()) 
Data <- read.csv("C:/Users/zjy/Desktop/coad_DEG_m.csv", header=T)
##################################################################################################
library(UpSetR)
m0 <- subset(Data, Phase_type=="m0")
m1 <- subset(Data, Phase_type=="m1")
m1a <- subset(Data, Phase_type=="m1a")
m1b <- subset(Data, Phase_type=="m1b")

listinput <- list(m0 = m0$SYMBOL,
                  m1 = m1$SYMBOL,
                  m1a = m1a$SYMBOL,
                  m1b = m1b$SYMBOL)
p <- upset(fromList(listinput),nsets = 4, order.by = "freq", 
           queries = list(list(query=intersects, params=list("m1"), color="red", active=T), 
                          list(query=intersects, params=list("m1a"), color="green", active=T), 
                         # list(query=intersects, params=list("m1b"), color="yellow", active=T), 
                          list(query=intersects, params=list("m1b"), color="blue", active=T)))
          

jpeg(file = "03Upset1.jpg",width =1600,height = 2000,units = "px",res =300) #结果保存
print(p)
dev.off()

##################################################################################################
library(GOplot)
library(topGO)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)


entrez_id=mapIds(x=org.Hs.eg.db,keys=as.character(m0$SYMBOL),keytype="SYMBOL",colum="ENTREZID")
entrez_id=na.omit(entrez_id)


go <- enrichGO(entrez_id,,OrgDb=org.Hs.eg.db, ont='all',keyType ="ENTREZID",pvalueCutoff=0.01,qvalueCutoff=0.05,readable=T)
write.csv(go,file="coad_erich.go.CCBPMF_m.csv")
m0_ccbpmf<-dotplot(go, split='ONTOLOGY',showCategory = 10) + facet_grid(ONTOLOGY~., scale='free')
jpeg(file = "m0_ccbpmf3.jpg",width =2000,height = 5000,units = "px",res =300) #结果保存
print(m0_ccbpmf)

go.BP <- enrichGO(entrez_id, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01, qvalueCutoff = 0.05,keyType = 'ENTREZID')
bpgraph<-plotGOgraph(go.BP)
jpeg(file = "04Bpgraph2.jpg",width =2000,height = 1600,units = "px",res =300) #结果保存
print(bpgraph)


#oragnx <- setReadable(erich.go.BP, 'org.Hs.eg.db', 'ENTREZID')  ## 将 Gene ID 转换为 symbol
#test<-cnetplot(oragnx,showCategory = 2)
#jpeg(file = "test.jpg",width =2000,height = 1600,units = "px",res =300) #结果保存
#print(test)

#GOBubble(go, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)  
##################################################################################################
