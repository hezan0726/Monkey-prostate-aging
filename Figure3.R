####scenic
library(SCENIC)
library(AUCell)
library(pheatmap)
library(SingleCellExperiment)
library(Seurat)
library(SCopeLoomR)
library(doParallel)
DefaultAssay(sobj)='RNA'
exprMat <- as(as.matrix(sobj@assays$RNA@data), 'sparseMatrix')
exprMat<- as.matrix(exprMat)
exprMat=exprMat[unique(DEG),]
cellInfo <- data.frame(seuratCluster=Idents(sobj))
org="hgnc" 
dbDir= 'cisTarget_databases/hg19'
myDatasetTitle="10x" 
data(defaultDbNames)
dbs <- list('500bp'= 'hg19-500bp-upstream-7species.mc9nr.feather', 
            '10kb' = 'hg19-tss-centered-10kb-7species.mc9nr.feather')
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs,
                                  datasetTitle=myDatasetTitle, 
                                  nCores=10)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions)
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) 
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
runSCENIC_4_aucell_binarize(scenicOptions)
export2loom(scenicOptions, exprMat_log)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
regulonTargetsInfo <- readRDS("/wd/2.5_regulonTargetsInfo.Rds")
regulonTargetsInfo=merge(regulonTargetsInfo,DEG,by='TF')#target

##################
####monocle2
library(monocle)
sub=subset(sobj,ct%in%c('BE','LE'))
matrix=sub@assays$RNA@data#or counts
matrix=matrix[which(rowSums(matrix) > 10),]
feature_ann<-data.frame(gene_id=rownames(matrix),gene_short_name=rownames(matrix))
rownames(feature_ann)<-rownames(matrix)
fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<- sub@meta.data
rownames(sample_ann)<-colnames(matrix)
pd<-new("AnnotatedDataFrame", data =sample_ann)
monocle_cds<-newCellDataSet(matrix,phenoData =pd,
                            featureData =fd,
                            expressionFamily=negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds) 
monocle_cds <- estimateDispersions(monocle_cds) 
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(monocle_cds),num_cells_expressed >= 10))
diff_celltype <- differentialGeneTest(monocle_cds[expressed_genes,],fullModelFormulaStr = "~ct")
diff_celltype<- diff_celltype[order(diff_celltype$qval),]
ordering_genes <- row.names(diff_celltype[1:150,]) 
monocle_cds <- setOrderingFilter(monocle_cds,ordering_genes = ordering_genes)
plot_ordering_genes(monocle_cds)
monocle_cds <- reduceDimension(monocle_cds, method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds)
plot_cell_trajectory(monocle_cds , color_by = "ct",show_branch_points=F)
plot_cell_trajectory(monocle_cds , color_by = "Pseudotime")
plot_cell_trajectory(monocle_cds , color_by = "age")
df <- pData(monocle_cds)
ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()+facet_grid(~age)
p=plot_pseudotime_heatmap(monocle_cds[ row.names(diff_celltype[1:600,]),],
                        num_clusters =2, 
                        cores = 1,return_heatmap=T,
                        show_rownames = T)
cds_subset = monocle_cds[row.names(diff_celltype[1:600,]),]
cluster_rows = TRUE
hclust_method = "ward.D2"
num_clusters = 2
hmcols = NULL
add_annotation_row = NULL 
add_annotation_col = NULL
show_rownames = FALSE
use_gene_short_name = TRUE
norm_method = c("vstExprs")
scale_max = 3
scale_min = -3 
trend_formula = "~sm.ns(Pseudotime, df=3)"
return_heatmap = TRUE
cores = 5
save.wd = "./"



num_clusters <- min(num_clusters, nrow(cds_subset))
pseudocount <- 1
newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), 
                                       max(pData(cds_subset)$Pseudotime), length.out = 100))
m <- genSmoothCurves(cds_subset, cores = 1, trend_formula = trend_formula, 
                     relative_expr = T, new_data = newdata)

m.copy <- m
m = m[!apply(m, 1, sum) == 0, ]
norm_method <- match.arg(norm_method)
if (norm_method == "vstExprs" && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) == 
    FALSE) {
  m = vstExprs(cds_subset, expr_matrix = m)
} else if (norm_method == "log") {
  m = log10(m + pseudocount)
}
m = m[!apply(m, 1, sd) == 0, ]
m = Matrix::t(scale(Matrix::t(m), center = TRUE))
m = m[is.na(row.names(m)) == FALSE, ]
m[is.nan(m)] = 0
m[m > scale_max] = scale_max
m[m < scale_min] = scale_min
heatmap_matrix <- m
row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1
if (is.null(hmcols)) {
  bks <- seq(-3.1, 3.1, by = 0.1)
  hmcols <- blue2green2red(length(bks) - 1)
} else {
  bks <- seq(-3.1, 3.1, length.out = length(hmcols))
}
ph <- pheatmap(heatmap_matrix, useRaster = T, cluster_cols = FALSE, 
               cluster_rows = cluster_rows, show_rownames = F, show_colnames = F, 
               clustering_distance_rows = row_dist, clustering_method = hclust_method, 
               cutree_rows = 2, silent = TRUE, filename = NA, 
               breaks = bks, border_color = NA, color = hmcols)
if (cluster_rows) {
  annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
                                                       num_clusters)))
} else {
  annotation_row <- NULL
}
if (!is.null(add_annotation_row)) {
  old_colnames_length <- ncol(annotation_row)
  annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), 
  ])
  colnames(annotation_row)[(old_colnames_length + 1):ncol(annotation_row)] <- colnames(add_annotation_row)
}
if (!is.null(add_annotation_col)) {
  if (nrow(add_annotation_col) != 100) {
    stop("add_annotation_col should have only 100 rows (check genSmoothCurves before you supply the annotation data)!")
  }
  annotation_col <- add_annotation_col
} else {
  annotation_col <- NA
}
if (use_gene_short_name == TRUE) {
  if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
    feature_label <- as.character(fData(cds_subset)[row.names(heatmap_matrix), 
                                                    "gene_short_name"])
    feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
    row_ann_labels <- as.character(fData(cds_subset)[row.names(annotation_row), 
                                                     "gene_short_name"])
    row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
  } else {
    feature_label <- row.names(heatmap_matrix)
    row_ann_labels <- row.names(annotation_row)
  }
} else {
  feature_label <- row.names(heatmap_matrix)
  if (!is.null(annotation_row)) 
    row_ann_labels <- row.names(annotation_row)
}
row.names(heatmap_matrix) <- feature_label
if (!is.null(annotation_row)) 
  row.names(annotation_row) <- row_ann_labels
colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
ph_res <- pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE, 
                   cluster_rows = cluster_rows, show_rownames = show_rownames, 
                   show_colnames = F, clustering_distance_rows = row_dist, 
                   clustering_method = hclust_method, cutree_rows = num_clusters, 
                   annotation_row = annotation_row, annotation_col = annotation_col, 
                   treeheight_row = 20, breaks = bks, fontsize = 6, color = hmcols, 
                   border_color = NA, silent = TRUE, filename = NA)

save(heatmap_matrix, file = paste0(save.wd, "plot_pseudotime_heatmap_matrix.rdata"))



ph_res$tree_row$labels[ph_res$tree_row$order]


clusters <- cutree(ph_res$tree_row, k = num_clusters)
clustering <- data.frame(clusters)
clustering[,1] <- as.numeric(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
clustering$Gene <- rownames(clustering)
clustering <- clustering[order(clustering$Gene_Clusters,decreasing = F),]
table(clustering$Gene_Clusters)


select.cluster <- 2

select.gene <- subset(clustering, Gene_Clusters %in% select.cluster)$Gene

plot.mat <- heatmap_matrix[select.gene, ]

plot.mat=as.data.frame(plot.mat)
plot.mat$gene=row.names(plot.mat)
plot.mat=plot.mat%>%pivot_longer(1:100,names_to = 'cell',values_to = 'exp')
plot.mat$cell %<>% factor(levels = seq(1, 100 ,1))
plot.mat$group <- "smooth"
ggplot(plot.mat, aes(x = cell, y = exp)) + 
  #geom_line(aes(group = gene)) +
  geom_smooth(aes(group = group))+
  theme_ArchR()
###########
pseudotime <- pseudotime(monocle_cds)
pseudotime <- pseudotime[rownames(sub@meta.data)]
sub$pseudotime <- pseudotime
meta=sub@meta.data
df=data.frame(ct=meta$ct,age=meta$age,pseu=meta$Pseudotime)
df <- df[order(df$pseu),]
df$ct=as.character(df$ct)
col=paletteContinuous(set = "paired",n=length(table(meta$ct)))
names(col)=levels(meta$ct)
df2=melt(df,id='pseu')
df2$pseu=factor(df2$pseu)
df_ct=subset(df2,df2$variable=='ct')
df_sub <- df_ct[(sample(rownames(df_ct), size =400, replace=F)),]
df_sub%>%ggplot(aes(x=pseu,y=variable))+
  geom_tile(aes(fill=value),size=0)+ #color和size分别指定方块边线的颜色和粗细
  scale_x_discrete("",expand = c(0,0))+ #不显示横纵轴的label文本；画板不延长
  scale_y_discrete("",expand = c(0,0))+
  scale_fill_manual(values = cols)+ #指定自定义的颜色
  theme(
    axis.text.x.bottom = element_text(size=0),axis.text.y.left = element_text(size = 12), #修改坐标轴文本大小
    axis.ticks = element_blank(), #不显示坐标轴刻度
    legend.title = element_blank() #不显示图例title
  )
##########################
cluster=cutree(p$ph$tree_row,k=2)
cluster=as.data.frame(cluster)
cluster$gene=row.names(cluster)
sub <- AddModuleScore(
  object = sub,
  features =list(cluster) ,
  ctrl = 5,
  name = 'cluster'
)
#################










