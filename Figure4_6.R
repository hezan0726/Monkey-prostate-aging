#####bulk RNA seq
#si-GRHL2 vs si-NC
df=read.table('rna_seq.featurecounts.full.tsv')
df=df[,-c(2:5)]
colnames(df)=df[1,]
df=df[-1,]
colnames(df)=c('GeneID','KD1','KD2','KD3','KD4','KD5','Ctrl1','Ctrl2','Ctrl3','Ctrl4','Ctrl5')
df=merge(ref38,df,by='GeneID')
df=df[,-1]
df=df[!duplicated(df$GeneName),]
rownames(df)=df$GeneName
df=df[,-1]
df1=apply(df,2, as.numeric)
rownames(df1)=rownames(df)
saveRDS(df1,'count.RDS')
library(tidyverse)
library(factoextra)
library(ggplot2)
library(genefilter)
df_count=df
df=edgeR::cpm(df)
group <- data.frame(sample =colnames(df),
                    group=c("KD","KD","KD","KD","KD","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl"))
rv <- genefilter::rowVars(df)
select <- order(rv, decreasing = TRUE)[seq_len(1000)]
pca_data <- cbind(t(log10(df[select,]+1)),group)
expr_pca <- prcomp(pca_data[,1:1000],scale = T,center = T)

fviz_screeplot(expr_pca, addlabels = TRUE, ylim = c(0, 40))

fviz_pca_ind(expr_pca,
             label = 'all',
             geom.ind = c('point','text'),
             habillage = group$group,    #分组变量
             addEllipses = T,
             ellipse.level = c(0.88),
             palette = c('#1F78B4','#33A02C','#FF7F00','#ff3f35')) +
  theme_bw((base_size=14))+  
  theme(text = element_text(size = 20),
        legend.margin = margin(-10),
        axis.text = element_text(size = 5, colour = 'black'),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm")
  )+
  ggtitle('')
####################################################################
group <- factor(c('KD','KD','KD','KD','KD','Ctrl','Ctrl','Ctrl','Ctrl','Ctrl'))
Data <- data.frame(row.names = colnames(df), 
                   group = group)
library(DESeq2)
df=df_count
dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = Data,
                              design = ~ group)
dds2 <- DESeq(dds)
res <- results(dds2,contrast = c("group",'KD','Ctrl'))
res=as.data.frame(res)
res$change = as.factor(ifelse(res$padj < 0.05 & abs(res$log2FoldChange) >= 0.5
                              , ifelse(res$log2FoldChange> 0.5 ,'UP','DOWN'),'NOT'))
res$GeneName<-rownames(res)
table(res$change)
res=subset(res,!res$change=='NOT')
write.csv(res,paste0('GRHL2_KD vs Ctrl',".DEG.csv"))
#go
require(DOSE)
require(clusterProfiler)#基于超几何分布检验的富集分析
require(enrichplot)
library(ggplot2)
library(ReactomePA)
library(meshes)
DEG_GO = compareCluster(data = DEG,GeneName ~ change, fun='enrichGO', OrgDb= 'org.Hs.eg.db',ont='ALL',keyType='SYMBOL')
dotplot(DEG_GO, showCategory=5) 
DEG_GO_summary=summary(DEG_GO)
write.csv(DEG_GO_summary,paste0('GRHL2_KD vs Ctrl',".GO.csv"))

#mouse
library(tidyverse)
library(factoextra)
library(ggplot2)
library(genefilter)
df2=edgeR::cpm(df)
group <- data.frame(sample =colnames(df),
                    group=c(rep("OC",11),rep("O_grhl2",10),rep("YC",10)))
rv <- genefilter::rowVars(df2)
select <- order(rv, decreasing = TRUE)[seq_len(1000)]
pca_data <- cbind(t(log10(df2[select,]+1)),group)
expr_pca <- prcomp(pca_data[,1:1000],scale = T,center = T)

fviz_pca_ind(expr_pca,
             label = 'all',
             geom.ind = c('point','text'),
             habillage = group$group,    #分组变量
             addEllipses = T,
             ellipse.level = c(0.88),
             palette = c('#1f78b4','#e2a461','#92b8c7','#ff3f35')) +
  theme_bw((base_size=14))+  
  theme(text = element_text(size = 20),
        legend.margin = margin(-10),
        axis.text = element_text(size = 5, colour = 'black'),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm")
  )+
  ggtitle('')
####################################################################
group <- factor(c(rep("OC",11),rep("O_grhl2",10),rep("YC",10)))
Data <- data.frame(row.names = colnames(df), 
                   group = group)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = Data,
                              design = ~ group)
dds2 <- DESeq(dds)
res <- results(dds2,contrast = c("group",'OC','YC'))
res=as.data.frame(res)
res$change = as.factor(ifelse(res$padj < 0.05 & abs(res$log2FoldChange) >= 0.5
                              , ifelse(res$log2FoldChange> 0.5 ,'UP','DOWN'),'NOT'))

res$GeneName<-rownames(res)

table(res$change)
OC_YC_DEG=res
Grhl2_OC_DEG=res
#EXPRESSION
library(magrittr)
library(tidyverse)
genes=c('Grhl2','Cdk19','Cdkn1a','Trp53')
df$GeneName=rownames(df)
df1<-as_tibble(df[genes,])%>%
  pivot_longer(-GeneName,names_to = 'Sample',values_to = 'Count')
Data$Sample=row.names(Data)
Data%<>%
  as_tibble()%>%
  left_join(df1)
Data$cpm=edgeR::cpm(Data$Count)
ggpubr::ggviolin(Data, x='group', y='cpm', fill = 'group', width = 0.8,facet.by = 'GeneName')+
  geom_boxplot(width=0.15)+
  ggpubr::stat_compare_means(comparisons = list(c("YC", "OC"),
                                                c("OC", "O_grhl2"),
                                                c("YC", "O_grhl2")))+ 
  theme(axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size=10,angle=90,color = 'black'),axis.text.y = element_text(size=10),legend.text=element_text(size=10),legend.title=element_text(size=0))
save(OC_YC_DEG,Grhl2_OC_DEG,file='DEG.RData')
OC_YC_DEG$Tag='OC_YC_DEG'
Grhl2_OC_DEG$Tag='Grhl2_OC_DEG'
DEG=rbind(Grhl2_OC_DEG,OC_YC_DEG)
DEG=subset(DEG,!DEG$change=='NOT')
rescue=subset(DEG,DEG$GeneName%in%rescue_gene)
range(rescue$foldchange)
bk=seq(-5.262734,8.003479, length.out = 100)
color2 = colorRampPalette(colors = c("#229187","#f8f8f8"))(40)
color1 = colorRampPalette(colors = c("#f8f8f8","#a96b36"))(length(41:100))
pheatmap(c[,c(2,3)],
         show_colnames =T,
         show_rownames = T, 
         cluster_rows =T,
         cluster_cols =F,
         breaks = bk,
         color = c(color2,color1))





