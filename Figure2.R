library(ggplot2)
library(Seurat)
library(dplyr)
library(DoubletFinder) 
library(future)
library(ArchR)
plan("multicore", workers = 20)
options(future.globals.maxSize = Inf)
setwd('/wd')
dir.create("outFigures")
dir.create("outData")
sample=c()
experi=read.csv('/wd/sample.csv',header=TRUE,sep=",",check.names=FALSE)
dim.usage=20
Find_doublet <- function(data){
  sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn[bcmvn$BCmetric==max(bcmvn$BCmetric),]$pK))
  DoubletRate = ncol(data)*8*1e-6
  nExp_poi <- round(DoubletRate*ncol(data))
  homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
 data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_low"
  data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_high"
  return(data)
}
allsample=NULL
for(i in 1:length(experi$sample)){
  data=Read10X_h5(experi$sample[i])
  data <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)
  data[["percent.mt"]] <- PercentageFeatureSet(data , pattern = "^MT-")
  for(j in 2:length(names(experi))){
    data[[names(experi)[j]]]=experi[i,j]
  }
  QC1=VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("/wd/outFigures/",sample[i],".QC1.pdf"),QC1,width=16,height=9)
  plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  QC2=CombinePlots(plots = list(plot1, plot2))
  ggsave(paste0("/wd/outFigures/",sample[i],".QC2.pdf"),QC2,width=16,height=9)
  #filter cells
  data <- subset(data, subset = nFeature_RNA > 500 & percent.mt < 10)
  data <- NormalizeData(object = data) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA()
  data <- FindNeighbors(data, dims = 1:20) %>%
    FindClusters(resolution = 0.5) %>%
    RunUMAP(dims = 1:20) %>% 
    RunTSNE(dims = 1:20)
  data <- Find_doublet(data)
  data@meta.data$DF_hi.lo <- data@meta.data$doublet_low
  data@meta.data$DF_hi.lo[which(data@meta.data$DF_hi.lo == "Doublet" & data@meta.data$doublet_high == "Singlet")] <- "Doublet-Low"
  data@meta.data$DF_hi.lo[which(data@meta.data$DF_hi.lo == "Doublet")] <- "Doublet-High"
  table(data@meta.data$DF_hi.lo)
  data=subset(data,DF_hi.lo=='Singlet')
  allsample=append(allsample,data)
}

list <- lapply(X = allsample, FUN = function(x) {
  x <- SCTransform(x, verbose = FALSE,do.scale = T)
})
save(list,file='list.RData')
list.features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000)
list <- PrepSCTIntegration(object.list = list, anchor.features = list.features, verbose = FALSE)
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT"
                                  ,anchor.features = list.features, verbose = FALSE)#,l2.norm=F,k.filter=20,k.score = 20)#dims = 1:10,k.score = 10
seurat<- IntegrateData(anchorset = anchors)#, k.weight = 20)
DefaultAssay(seurat) <- "integrated"     
seurat <- FindVariableFeatures(seurat,selection.method = "vst", nfeatures = 2000) 
seurat <-ScaleData(seurat)
seurat <- RunPCA(seurat,npcs = 100,features = VariableFeatures(object = seurat))
seurat <- FindNeighbors(seurat, dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:50) 
col=paletteContinuous(set = "paired",n=11)
p=DimPlot(seurat,reduction='umap',raster =F,shuffle=F,cols=col,label=T,repel = T)
ggsave('umap.pdf',p)
###############################
#####################
##############
luminal=c('MSMB','KLK2','KLK3','ACPP','NKX3-1')
basal=c('KRT5','KRT15')#12 13 7
SMC='MYH11'#6 14 
TC=c('CD3G','CD247','PTPRC')#10
Mac=c('CD74','C1QA')#18
EC=c('PECAM1')#15 21
Per=c('RGS5','IGFBP7','PDGFRB')#22

Neu=c('NCAM1','PMP22','SOX10','GPM6B')#19
Fib=c('PDGFRA','FN1','DCN')#5
pdf("feature_NE.pdf",width=10,height=10)
a=FeaturePlot(prostate_seurat, features = G2M.Score,cols=c("lightgrey", "red"))
dev.off()
genes_to_check=c('KRT15','KRT5',
                 'MSMB','NKX3-1',
                 'PECAM1','VWF',
                 'MYH11','ACTG2',
                 'DCN','PDGFRA',
                 'CD247','CD3G',
                 'CD74','C1QA',
                 'NCAM1','SOX10'
)
DefaultAssay(seurat)='RNA'
sobj=seurat

sobj@meta.data$ct=factor(sobj@meta.data$ct,levels=c("BE","LE","SMC","Fib",
                                                    "EC","TC",'Mac','Neu'))
p_all_markers <- DotPlot(sobj, features = genes_to_check,
                         assay='RNA' ,group.by = 'ct' )  + 
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
p_all_markers
data<-p_all_markers$data
colnames(data)<-c("AverageExpression_unscaled","Precent Expressed","Features","celltype","Average Expression")
unique(data$`Precent Expressed`)
ggplot(data,aes(celltype,Features,size = `Precent Expressed` ))+
  geom_point(shape=21,aes(fill= `Average Expression`),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,10))+theme_bw()+
  scale_fill_gradient(low = "grey", high = "#E54924")+
  theme(legend.position = "right",legend.box = "vertical", #图例位置
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=16,angle = 45, 
                                    vjust = 0.5, hjust=0.5),#x轴
        axis.text.y  = element_text(color="black",size=12),#y轴
        legend.text = element_text(size =12,color="black"),#图例
        legend.title = element_text(size =12,color="black"),#图例
        axis.title.y=element_text(vjust=1,  
                                  size=16)
  )+labs(x=" ",y = "Features")
######
####
#cell identity
library(tidyverse)
sobj_y <- sobj %>%
  subset(age == "Young")
Idents(sobj_y) <- sobj_y$ct
markers <- FindAllMarkers(sobj_y, only.pos = TRUE)
markers <- markers %>%
  as_tibble() %>%
  filter(p_val_adj < 0.05)
markers %>%
  filter(avg_log2FC > 1 & pct.1 > 0.1) %>%
  group_by(cluster) %>%
  top_n(50, wt = avg_log2FC)%>%
  summarise(num = n())
cell_identity_score <- function(data, feature) {
  data <- AddModuleScore(data,
                         features = feature,
                         name = "cell_identity_score__")
  
  df <- data@meta.data %>%
    as_tibble(rownames = "cellid") %>%
    select(class, age, sample, ct,treat,group, starts_with("cell_identity_score__")) %>%
    rename_at(vars(starts_with("cell_identity_score__")), ~{
      t <- names(feature)[as.numeric(gsub("cell_identity_score__","",.x))]
      str_c("cell_identity_score_", t)
    })
}

feature <- markers %>%
  mutate(
    feature_name = cluster
  ) %>%
  select(feature_name, gene) %>%
  group_nest(feature_name) %>%
  pmap(~ {
    t <- ..2[["gene"]]
    setNames(list(t), ..1)
  }) %>%
  unlist(recursive = F)
score_df <- cell_identity_score(sobj, feature)
score_long <- score_df %>%
  pivot_longer(
    starts_with("cell_identity_score"),
    names_prefix = "cell_identity_score_",
    names_to = "cell_identity_score_celltype",
    values_to = "cell_identity_score")
library(ggpubr)
score_long %>%
  filter(celltype == cell_identity_score_celltype & celltype %in% celltypes) %>%
  mutate(celltype = factor(celltype, levels = celltypes) ) %>%
  ggviolin(x = "group", y = "cell_identity_score", fill = "age", 
           palette = group_cols, 
           add = "boxplot", add.params = list(fill="white")) + 
  stat_compare_means(comparisons = age_comparisons) +
  theme(
    legend.position = "right",
    axis.title.x = element_blank()
  ) +
  facet_grid(.~celltype)
############################
###################
#####cv########
mat.df <- GetAssayData(object = sobj, slot = "data") %>% as.matrix
hvgs <- FindVariableFeatures(object = sobj, nfeatures = 0.1 * nrow(mat.df)) 
hvgs <-VariableFeatures(object = hvgs) 
mat.df %<>% .[hvgs, ]
sobj$ident=paste(sobj$ct,sobj$age,sep='_')
Idents(sobj)='ident'
CV.df <- NULL
celltype=c('LE','BE','SMC','Fib','Neu','TC','Mac','EC')
length_celltype <- length(celltype)
cv_wd='/wd'
for (i in 1:length(celltype)) {
  object_Y <- WhichCells(sub, ident = paste0(celltype[i],"_Young"))
  object_O <- WhichCells(sub, ident = paste0(celltype[i],"_Old"))
  df_Y <- mat.df[,object_Y]
  df_O <- mat.df[,object_O]
  
  Mat.abs <- NULL
  for (j in 1:ncol(df_Y)) {
    mat.abs <- abs(df_Y[,j] - df_O)
    Mat.abs %<>% cbind(mat.abs)
    rm(mat.abs)
    gc()
  }
  
  fwrite(Mat.abs %>% as.data.frame, paste(cv_wd,celltype[i], "_CV.csv", sep = ""), row.names=T, sep = "\t")
  gc()
  
  mean.df <- apply(Mat.abs, 1, mean)
  write.csv(mean.df,paste(cv_wd,celltype[i],"_mean_sub.csv",sep = ""))
  sd.df <- apply(Mat.abs, 1, sd)
  write.csv(sd.df,paste(cv_wd,celltype[i],"_sd_sub.csv",sep = ""))
  cv.df <-(sd.df / mean.df) * 100 %>%  as.matrix
  write.csv(cv.df,paste(cv_wd,celltype[i],"_cv_sub.csv",sep = ""))
  CV.df %<>% cbind(cv.df)
  colnames(CV.df)[i] <- celltype[i]
  gc()
  print(paste0(celltype[i]," has been calculated"))
}
write.csv(CV.df,paste(cv_wd,"_cv_sub.csv",sep = ""))
cv=CV.df
cv$gene=row.names(cv)
cv%<>%pivot_longer(1:8,names_to = 'celltype',values_to = 'cv')
library(CellChat)
library(ggplot2)
library(ggpubr)
group_cols=scPalette(8)
comparisons <- list(c("LE", "BE"),c("BE", "TC"),c("M_C", "O_C"), c("O_C", "O_CR"))
ggpubr::ggboxplot(cv,x = "celltype", y = "cv", fill = "celltype",
                  palette = group_cols, 
                  add.params = list(fill="white"))
stat_compare_means(comparisons = comparisons) +
  theme(
    legend.position = "right",
    axis.title.x = element_blank()
  )

##################
##############
####DEG####\
sobj@meta.data$ct_group=paste(sobj$age,sobj$ct,sep = '/')
Idents(sobj)=sobj$ct_group
DEG <- data.frame()
for (cell in levels(sobj$ct)){
  x.inv <- try(markers <- FindMarkers(sobj, ident.1 = paste('Old',cell,sep = '/'), 
                                      ident.2 = paste('Young',cell,sep = '/'), logfc.threshold = 0), silent=TRUE)
  if ('try-error' %in% class(x.inv))
    next
  markers$gene <- rownames(markers)
  markers$ct <- cell
  DEG <- rbind(DEG, markers)
  print(paste0(cell, ' is finished'))
}
DEG$change = as.factor(ifelse(DEG$p_val_adj < 0.05 & abs(DEG$avg_log2FC) >= 0.25
                                , ifelse(DEG$avg_log2FC> 0.25 ,'Upregulated','Downregulated'),'Notsignificant'))

###############
######
##GO enrichment
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tidyverse)
library(RColorBrewer)
up=subset(DEG,DEG$change=='Upregulated')
down=subset(DEG,DEG$change=='Downregulated')
GO_down <- compareCluster(
  gene ~ ct,
  data = down,
  fun = "enrichGO",
  OrgDb = "org.Hs.eg.db",
  ont = "ALL",keyType='SYMBOL'
)
GO_up <- compareCluster(
  gene ~ ct,
  data = up,
  fun = "enrichGO",
  OrgDb = "org.Hs.eg.db",
  ont = "ALL",keyType='SYMBOL'
)
View(summary(GO_up))
GO_up  %>% 
  pairwise_termsim() %>%
  emapplot(
    showCategory = 50,
    layout = 'kk',
    pie = "count",
    min_edge = 0.2,
    cex_label_category = 2,
    cex_line = 0.5,
    cex_category = 0.8,
    cex_pie2axis = 7,
    legend_n = 3) 
################
#######
####Geneset score
sobj <- AddModuleScore(
  object = sobj,
  features =list(inflammatory_response) ,
  ctrl = 5,
  name = 'inflammatory_response'
)
meta=sobj@meta.data
meta%<>%group_by(ct,age)%>%summarise(autophagy=median(autophagy1),DNA_repair=median(DNA_repair1),
                                    inflammatory_response=median(inflammatory_response1),
                                    cell_cycle=median(cell_cycle1),
                                    cellular_senescence=median(cellular_senescence1),
                                    Fibrosis=median(Fibrosis1))
df=meta[c(17:24),3:16]
df=as.data.frame(df)
rownames(df)=df$age
range(df)
bk <- c(seq(-0.05281067,0,by=0.001),seq(0,0.09435018,by=0.001))
bk
col = c(colorRampPalette(colors = c("#30706a","white"))(53),
        colorRampPalette(colors = c("white","#a26c3f"))(95))
pheatmap::pheatmap(df
                   , cluster_rows = F, cluster_cols = F, border_color = '#999999', 
                   color = col 
)
meta=sobj@meta.data
meta=meta[,c(5,40,50:62)]
meta=meta%>%pivot_longer(3:7,names_to = 'term',values_to = 'score')
ggpubr::ggviolin(meta, x='age', y='score', fill = 'age', width = 0.5)+
  geom_boxplot(width=0.3)+
 ggpubr::stat_compare_means(comparisons = list(c("Young", "Old")))+ #, label = "p.signif"
  theme(axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size=10,angle=90,color = 'black'),axis.text.y = element_text(size=10),
        legend.text=element_text(size=10),legend.title=element_text(size=0))+facet_wrap(~ct+term)


###############
######
#









