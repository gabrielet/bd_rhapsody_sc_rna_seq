#' ---
#' title: "SLM clustering for zinbwaved-batched AD and EAE experiments"
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
expsList <- list()
expsList[["AD"]] <- "/home/gabriele/cbmc/scrnaseq/raw_data_new/AD_batched/"
expsList[["EAE"]] <- "/home/gabriele/cbmc/scrnaseq/raw_data_new/EAE_batched/"

for (expN in names(expsList)) {

	rawPath <- expsList[[expN]]
	saveIn <- paste0(rawPath, "results/")
	saveImg <- paste0(saveIn, "figures/")

	print(paste0("analysing batched experiments from ", expN))

	#check if saveIn directory exist, otherwise create it
	ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
	ifelse(dir.exists(saveImg), TRUE, dir.create(saveImg, recursive=F))

	# set palette with 12 colours for colouring plots
	colorsPal <- brewer.pal(12, "Paired")

	# number of highly variable genes used for zinbwaving the data
	zinHV <- 1000
	zbWaved <- readRDS(paste0(rawPath, "zbwaved_data_batch_corrected_batched_K_20_top_1000.Rds"))

	# rename colnames since cells may have the same name
	colnames(zbWaved) <- seq(1, ncol(zbWaved), by=1)

	# obtain a seurat object
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

	# new labels that compress the fine name into something similar to main labels
	# zbWaved$fine_to_main <- unlist(lapply(strsplit(unlist(lapply(strsplit(zbWaved$pruned_fine, "\\."), `[[`, 1)), " "), `[[`, 1))

	zbWaved$fine_to_main <- gsub(" ", "", unlist(lapply(strsplit(zbWaved$pruned_fine, "\\("), `[[`, 1)))

	# a summary of clusters content
	tbl <- table(zbWaved$slm_1, zbWaved$fine_to_main)

	# plot tSNE
	png(paste0(saveImg, "03_01_slm_clustering_batched_K_1000_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	grid.arrange(
		scater::plotTSNE(zbWaved, colour_by="slm_1.5", text_by="slm_1.5") + theme(legend.position="none") + labs(title="By cluster"),
		scater::plotReducedDim(zbWaved, "TSNE", colour_by="fine_to_main") + theme(legend.position="right") + labs(title="By cell type") + guides(color=guide_legend("Cell type")) ,
		ncol=2)
	dev.off()

	# finally, export the data
	saveRDS(object=zbWaved, file=paste0(saveIn, "slm_clustered_data_batched_K_20_top_", zinHV, "_", expN, ".Rds"))
}
