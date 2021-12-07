#单细胞数据预处理
#### QC core code
data = Read10X(data.dir=datapath)
obj = CreateSeuratObject(counts=data)
obj[['percent.mt']] = PercentageFeatureSet(obj,pattern='MT-')	# QC
vln = VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave('VlnPlot.pdf',vln,width=60,height=30,units='cm')
plot1 = FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")	# 以线粒体比例和feature数目过滤细胞
plot2 = FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
cp = CombinePlots(plots = list(plot1, plot2))
ggsave('CombinePlots.pdf',cp,width=60,height=30,units='cm')
mt = 20
obj = subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < mt)	# 以线粒体比例和feature数目过滤细胞

#单细胞数据聚类
#### Preprocess core code
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)	# 规范化数据
obj = FindVariableFeatures(obj, selection.method = "vst")	# 识别高度可变的特征（特征选择）
obj = ScaleData(obj, features = rownames(obj))	# 缩放数据
obj = RunPCA(obj, features = VariableFeatures(object = obj))	# PCA
pca = DimPlot(obj, reduction = "pca")
ggsave('PCA.pdf',pca,width=30,height=30,units='cm')

obj = JackStraw(obj,num.replicate = 100)	# 找到具有统计学显著性的主成分
obj = ScoreJackStraw(obj,dims=1:20)
jackstraw = JackStrawPlot(obj, dims = 1:20)
ggsave('JackStraw.pdf',jackstraw,width=30,height=30,units='cm')

p = ElbowPlot(obj,ndims=50)
ggsave('ElbowPlot.pdf',p,width=30,height=20,units='cm')
obj1 = FindNeighbors(obj, dims = 1:25)
obj.tsne = FindClusters(obj1, resolution = 1)
obj.tsne = RunTSNE(obj.tsne,dims=1:25)
tnseplot = DimPlot(obj.tsne, reduction = "tsne",pt.size = 1.2,label = T,cols=cols)
ggsave('TSNE.pdf',tnseplot,width=25,height=20,units='cm')

#识别肿瘤微环境中的各细胞类型/某一细胞类型的亚型
#### FeaturePlot core code
FeaturePlot(dat, reduction = 'tsne', features = c('CD3E','CD3G','CD4','CD8A'), 
            min.cutoff = 0, coord.fixed=T,order=T)

#### core code
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se
## class: SummarizedExperiment 
## dim: 19363 713 
## metadata(0):
## assays(1): logcounts
## rownames(19363): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
## rowData names(0):
## colnames(713): GSM112490 GSM112491 ... GSM92233 GSM92234
## colData names(3): label.main label.fine label.ont
library(scRNAseq)
hESCs <- LaMannoBrainData('human-es')
hESCs <- hESCs[,1:100]
library(SingleR)
pred.hesc <- SingleR(test = hESCs, ref = hpca.se, assay.type.test=1,
                     labels = hpca.se$label.main)
#### SingleR core code
# Find markers for all clusters
all.markers <- FindAllMarkers(object = pbmc_small)
head(x = all.markers)

#计算各个细胞的细胞周期评分
#### CellCycle core code
ccs = list(cellcycleGene = c('ASPM', 'CENPE', 'CENPF', 'DLGAP5', 'MKI67', 'NUSAP1', 'PCLAF', 'STMN1', 'TOP2A', 'TUBB'))
ccs.dat = ccScore.analyze(dat.expr=mat,dat.gene.list=ccs,k=10,bin=25,controlSize=100)
dat$CellCycleScore = ccs.dat$cellcycleGene_raw

#估计各个细胞的细胞周期状态
#### CellCycleScoring core code
tmp = CellCycleScoring(dat, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
dat@meta.data[,c('S.Score','G2M.Score','Phase')] = tmp@meta.data[rownames(tmp@meta.data),c('S.Score','G2M.Score','Phase')]

#计算细胞的分化潜能/干性程度
#### CytoTRACE core code
step = 20000
mat.list = lapply(1:ceiling(ncol(dat)/step),function(x){
  min = (x-1) * step + 1
  max = ifelse(x*step>=ncol(dat), ncol(dat), x*step)
  return(as.matrix(dat@assays$RNA@data[,min:max])) })    
mat = do.call(cbind,mat.list)
res = CytoTRACE(mat,ncores = 4)
dat$stemness = res$CytoTRACE

#基于monocle2推断细胞变化的轨迹
#### monocle2 core code

#linux
#RNA速率分析
#### RNA velocyto core code
#linux velocyto run10x -m hg38_rmsk.gtf -@ 8 cellrangerRs genes.gtf

#基于slingshot推断细胞变化的轨迹
#### slingshot core code
#linux sim <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'PCA')

#基于PAGA推断亚群之间的关联性
#### PAGA core code
sc.tl.paga(adata, groups='leiden')

#基于Seurat-CCA整合单细胞数据
#### CCA core code
anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
combined <- IntegrateData(anchorset = anchors)

#基于Seurat-merge整合单细胞数据
#### Merge core code
sce.big <- merge(sceList[[1]],y = c(sceList[[2]],sceList[[3]],sceList[[4]],  sceList[[5]],sceList[[6]],sceList[[7]],sceList[[8]]),add.cell.ids = folders, project = "project")

#基于harmony整合单细胞数据
#### HARMONY core code
scRNA_harmony <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]], scRNAlist[[3]], scRNAlist[[4]], scRNAlist[[5]],scRNAlist[[6]], scRNAlist[[7]], scRNAlist[[8]], scRNAlist[[9]], scRNAlist[[10]]))
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")

#基于pySCENIC推断转录因子
#### pySCENIC core code
#linux pyscenic grn -o $grnboost2_output -m grnboost2 --seed 12345 --num_workers 6 $preprocessing_fname /public/workspace/caiyun/SCENIC/hs_hgnc_curated_tfs.txt
#linux pyscenic ctx step1.adjacencies.tsv /public/workspace/caiyun/SCENIC/hg19-500bp-upstream-7species.mc9nr.feather /public/workspace/caiyun/SCENIC/hg19-tss-centered-10kb-7species.mc9nr.feather --annotations_fname /public/workspace/caiyun/SCENIC/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname step0.csv --mode "custom_multiprocessing" --output step2.regulons.tsv --num_workers 8
#linux pyscenic aucell step0.csv step2.regulons.tsv -o step3.auc_mtx.csv --num_workers 6

#基于DoRothEA推测转录因子活性
#### DoRothEA core code
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
regulon <- dorothea_regulon_human %>%dplyr::filter(confidence %in% c("A","B","C"))
pbmc <- run_viper(pbmc, regulon,options = list(method = "scale", minsize = 4, eset.filter = FALSE, cores = 1, verbose = FALSE))

#基于CellPhoneDB进行细胞间通讯分析
#### CellPhoneDB core code
#linux cellphonedb method statistical_analysis yourmetafile.txt yourcountsfile.txt --iterations=10 --threads=2

#基于CellChat进行细胞间通讯分析
#### CellChat core code
cellchat <- createCellChat(data = data.input)