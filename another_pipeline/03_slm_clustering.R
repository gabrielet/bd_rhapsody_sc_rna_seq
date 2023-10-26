#' ---
#' title: "SLM clustering for zinbwaved data"
#' author: "gabrielet"
#' output: html_document
#' date: '`r format(Sys.Date(), "%d %B, %Y")`'
#' ---
#' 
#' ```{r setup, include=FALSE}
#' knitr::opts_chunk$set(
#'   tidy=TRUE,
#'   tidy.opts=list(width.cutoff=120),
#'   message=FALSE,
#'   warning=FALSE,
#'   eval=TRUE
#' )
#' ```

# Load libraries
library("bluster")
library("cluster")
library("compositions")
library("cowplot")
library("DropletUtils")
library("ggplot2")
library("gridExtra")
library("pheatmap")
library("plyr")
library("RColorBrewer")
library("scater")
library("scran")
library("Seurat")
library("SingleCellExperiment")

# set paths
# rawPath <- "/home/gabriele/cbmc/scrnaseq/raw_data_new/572_20/" ; expN <- "572" ; kZinb <- 20 ; perC <- 1000 ; # AD male
 rawPath <- "/home/gabriele/cbmc/scrnaseq/raw_data_new/566_20/" ; expN <- "566" ; kZinb <- 20 ; perC <- 1000 ; # AD female
 
saveIn <- paste0(rawPath, "results_agosto_nuovi/")
saveImg <- paste0(saveIn, "figures/")

print(paste0("analysing exp ", expN))

#check if saveIn directory exist, otherwise create it
ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
ifelse(dir.exists(saveImg), TRUE, dir.create(saveImg, recursive=F))

# set palette with 12 colours for colouring plots
colorsPal <- brewer.pal(12, "Paired")

# number of highly variable genes used for zinbwaving the data
zbWaved <- readRDS(paste0(saveIn, paste0("full_zbwaved_data_", expN, "_K_", kZinb, "_top_", perC, ".Rds")))

# transform SCE back to Seurat to use the clustering algorithms
zbwSeurat <- as.Seurat(zbWaved, data="counts")

# SNN graph
zbwSeurat <- FindNeighbors(object=zbwSeurat, reduction="zinbwave", dims=1:15)

# SLM clustering with different resolution
seurat_clusters_0.5 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=0.5, random.seed=1)$seurat_clusters
seurat_clusters_1 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=1, random.seed=1)$seurat_clusters
seurat_clusters_1.5 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=1.5, random.seed=1)$seurat_clusters

slmClusters <- as.matrix(data.frame(slm_0.5=seurat_clusters_0.5, slm_1=seurat_clusters_1, slm_1.5=seurat_clusters_1.5))

#store clustering results in zinbSce
zbWaved$slm_0.5 <- as.factor(slmClusters[, 1])
zbWaved$slm_1 <- as.factor(slmClusters[, 2])
zbWaved$slm_1.5 <- as.factor(slmClusters[, 3])

# a summary of clusters content
tbl <- table(zbWaved$slm_1, zbWaved$fine_to_main)

# plot tSNE
png(paste0(saveImg, "03_01_slm_0.5_clustering_", expN, "_K_", kZinb, "_top_", perC, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
grid.arrange(
	scater::plotTSNE(zbWaved, colour_by="slm_0.5", text_by="slm_0.5") + theme(legend.position="none") + labs(title="By cluster"),
	scater::plotReducedDim(zbWaved, "TSNE", colour_by="fine_to_main") + theme(legend.position="right") + labs(title="By cell type") + guides(color=guide_legend("Cell type")) ,
	ncol=2)
dev.off()

# plot tSNE
png(paste0(saveImg, "03_01_slm_1_clustering_", expN, "_K_", kZinb, "_top_", perC, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
grid.arrange(
	scater::plotTSNE(zbWaved, colour_by="slm_1", text_by="slm_1") + theme(legend.position="none") + labs(title="By cluster"),
	scater::plotReducedDim(zbWaved, "TSNE", colour_by="fine_to_main") + theme(legend.position="right") + labs(title="By cell type") + guides(color=guide_legend("Cell type")) ,
	ncol=2)
dev.off()

# plot tSNE
png(paste0(saveImg, "03_01_slm_1.5_clustering_", expN, "_K_", kZinb, "_top_", perC, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
grid.arrange(
	scater::plotTSNE(zbWaved, colour_by="slm_1.5", text_by="slm_1.5") + theme(legend.position="none") + labs(title="By cluster"),
	scater::plotReducedDim(zbWaved, "TSNE", colour_by="fine_to_main") + theme(legend.position="right") + labs(title="By cell type") + guides(color=guide_legend("Cell type")) ,
	ncol=2)
dev.off()

# finally, export the data
saveRDS(object=zbWaved, file=paste0(saveIn, "slm_clustered_data_", expN, "_K_", kZinb, "_top_", perC, ".Rds"))
