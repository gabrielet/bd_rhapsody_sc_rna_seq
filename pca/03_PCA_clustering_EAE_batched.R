#' ---
#' title: "K-means clustering for EAE size-factored data"
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

#' 2) performing [K-means clustering](https://bioconductor.org/books/release/OSCA/clustering.html#k-means-clustering){target="_blank"}  

# Load libraries
library("bluster")
library("cluster")
library("cowplot")
library("ggplot2")
library("gridExtra")
library("limma")
library("pheatmap")
library("plyr")
library("RColorBrewer")
library("scater")
library("scran")
library("Seurat")
library("SingleCellExperiment")

# set seed
set.seed(113)

# set paths and set kS parameter, which is dependent on the experiment
rawPath <- "/home/gabriele/work/cbmc/scrnaseq/raw_data_new/" ; expN <- "EAE_batched"

saveIn <- paste0(rawPath, "batched_EAE/")
saveImg <- paste0(saveIn, "figures/")

print(paste0("analysing exp ", expN))

# load data
kmeansSce <- readRDS(paste0(saveIn, paste0("denoised_data_kmeans_", expN, ".Rds")))

#' 4) A more practical use of k-means is to deliberately set k to a large value to achieve overclustering. This will forcibly partition cells inside broad clusters that do not have well-defined internal structure. For example, we might be interested in the change in expression from one “side” of a cluster to the other, but the lack of any clear separation within the cluster makes it difficult to separate with graph-based methods, even at the highest resolution. k-means has no such problems and will readily split these broad clusters for greater resolution, though obviously one must be prepared for the additional work involved in interpreting a greater number of clusters.  

# set an arbitrary kS, 15 seems to be best for separating clusters outside of the neutrophils one, for 
kS <- 15

# and an arbitrary seed
set.seed(131)

# new labels that compress the fine name into something similar to main labels
kmeansSce$fine_to_main <- unlist(lapply(strsplit(unlist(lapply(strsplit(kmeansSce$pruned_fine, "\\."), `[[`, 1)), " "), `[[`, 1))

clustKM <- kmeans(reducedDim(kmeansSce, "PCA"), centers=kS)
colLabels(kmeansSce) <- factor(clustKM$cluster)
kmeansSce$clusterS <- factor(clustKM$cluster)

# a summary of clusters content
tbl <- table(clustKM$cluster, kmeansSce$fine_to_main)

# plot by cluster
# plotting tSNE using cluster colours and cell label colours, side by side
png(paste0(saveImg, "03_01_tSNE_kmeans_K_", kS, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
grid.arrange(
	scater::plotTSNE(kmeansSce, colour_by="clusterS", text_by="clusterS") + theme(legend.position="none") + labs(title="By cluster"),
	scater::plotReducedDim(kmeansSce, "TSNE", colour_by="fine_to_main") + theme(legend.position="right") + labs(title="By cell type") + guides(color=guide_legend("Cell type")) ,
	ncol=2)
dev.off()

#' evaluate clusters separation: #https://bioconductor.org/books/release/OSCA/clustering.html#assessing-cluster-separation-1  

#' The within-cluster sum of squares (WCSS) for each cluster is the most relevant diagnostic for k-means, given that the algorithm aims to find a clustering that minimizes the WCSS  

ncells <- tabulate(clustKM$cluster)
tab <- data.frame(wcss=clustKM$withinss, ncells=ncells)
tab$rms <- sqrt(tab$wcss/tab$ncells)
tab

cent.tree <- hclust(dist(clustKM$centers), "ward.D2")
png(paste0(saveImg, "03_02_dendrogram_kmeans_K_", kS, "_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	plot(cent.tree)
dev.off()

# print tSNE with all the cells
png(paste0(saveImg, "02_01_kmean_tSNE_ALL_CLUSTERS_conditionBYstage_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotReducedDim(kmeansSce, dimred="TSNE", colour_by="conditionBYstage") + theme(legend.position="right") + labs(title="By condition") + guides(color=guide_legend("Condition")))
dev.off()

png(paste0(saveImg, "02_02_kmean_tSNE_ALL_CLUSTERS_cell_types_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotReducedDim(kmeansSce, dimred="TSNE", colour_by="pruned_main") + theme(legend.position="right") + labs(title="By cell type") + guides(color=guide_legend("Cell type")))
dev.off()

# plot tSNE with all the cells
png(paste0(saveImg, "03_01_kmean_tSNE_ALL_CLUSTERS_stage_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotReducedDim(kmeansSce, dimred="TSNE", colour_by="stage") + theme(legend.position="right") + labs(title="By condition") + guides(color=guide_legend("Condition")))
dev.off()

# plot tSNE by condition
png(paste0(saveImg, "03_01_kmean_tSNE_ALL_CLUSTERS_CONDITION_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotTSNE(kmeansSce, colour_by="condition") + theme(legend.position="right") + labs(title="Conditions"))
dev.off()

# plot tSNE by cluster
png(paste0(saveImg, "03_01_kmean_tSNE_ALL_CLUSTERS_CLUSTER_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotTSNE(kmeansSce, colour_by="clusterS") + theme(legend.position="right") + labs(title="Clusters"))
dev.off()

# plot tSNE by cluster
png(paste0(saveImg, "03_01_kmean_tSNE_ALL_CLUSTERS_CLUSTER_con_etichette", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotTSNE(kmeansSce, colour_by="clusterS", text_by="clusterS") + theme(legend.position="none") + labs(title="Clusters"))
dev.off()

# finally, store the normalised objects
saveRDS(object=kmeansSce, file=paste0(saveIn, "kmeans_K_", kS, "_clustered_data_", expN, ".Rds"))
