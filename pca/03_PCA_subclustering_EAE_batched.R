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
rawPath <- "/home/gabriele/work/cbmc/scrnaseq/raw_data_new/" ; expN <- "EAE_batched"

saveIn <- paste0(rawPath, "batched_EAE/")
saveImg <- paste0(saveIn, "figures/")

print(paste0("analysing exp ", expN))

#check if saveIn directory exist, otherwise create it
ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
ifelse(dir.exists(saveImg), TRUE, dir.create(saveImg, recursive=F))

# set palette with 12 colours for colouring plots
colorsPal <- brewer.pal(12, "Paired")

# set cellTypes
cellType <- c("T.8|T.CD8", "Neutrophils", "\\(T.4|\\(T.CD4", "B cells", "DC|Macrophages|Microglia|Monocytes")
cellLab <-  c("CD8", "Neutrophils", "CD4", "B_cells", "Myeloid")

# load data, all cells
kmeansSce_original <- readRDS(paste0(saveIn, paste0("denoised_data_kmeans_", expN, ".Rds")))

# loop through cell type
for (cT in seq(1, length(cellType), by=1)){

	###################################################################### get a specific cell type to perform subclustering ######################################################################
	kmeansSce <- kmeansSce_original[, grep(cellType[cT], kmeansSce_original$pruned_fine)]
	###############################################################################################################################################################################################

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
	
	# plot condition only
	png(paste0(saveImg, "03_subcl_02_kmean_tSNE_SUBCLUSTERED_condition_only_", expN, "_", cellLab[cT], ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotTSNE(kmeansSce, colour_by="condition"))
	dev.off()
	
	# plot condition only
	png(paste0(saveImg, "03_subcl_02_kmean_tSNE_SUBCLUSTERED_by_clusters_", expN, "_", cellLab[cT], ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotTSNE(kmeansSce, colour_by="clusterS"))
	dev.off()

	# plot region and condition
	png(paste0(saveImg, "03_subcl_03_kmean_tSNE_SUBCLUSTERED_conditionBYstage_", expN, "_", cellLab[cT], ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotReducedDim(kmeansSce, dimred="TSNE", colour_by="conditionBYstage") + theme(legend.position="right") + labs(title="By condition") + guides(color=guide_legend("Condition")))
	dev.off()

	# plot together
	png(paste0(saveImg, "03_subcl_04_kmean_tSNE_SUBCLUSTERED_stage_only_", expN, "_", cellLab[cT], ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotReducedDim(kmeansSce, dimred="TSNE", colour_by="stage") + theme(legend.position="right") + labs(title="By condition") + guides(color=guide_legend("Condition")))
	dev.off()

	# plot
	png(paste0(saveImg, "03_subcl_05_kmean_tSNE_SUBCLUSTERED_", expN, "_", cellLab[cT], ".png"), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
		print(scater::plotReducedDim(kmeansSce, dimred="TSNE", colour_by="pruned_fine") + theme(legend.position="bottom") + labs(title="By cell type") + guides(color=guide_legend("Cell type")))
	dev.off()
	
	###################################################################################################################
	# ELEONORA PLOT BRAIN ONLY
	elePlotOnset <- kmeansSce[, which(kmeansSce$stage=="Onset")]

	# plot
	png(paste0(saveImg, "03_00_kmean_tSNE_as_SUBCLUSTERED_Eleonora_Said_CLUSTERS_Onset_ONLY_", expN, "_", cellLab[cT], ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotTSNE(elePlotOnset, colour_by="clusterS"))
	dev.off()

	# ELEONORA PLOT MENINGES ONLY
	elePlotChronic <- kmeansSce[, which(kmeansSce$stage=="Chronic")]

	# plot
	png(paste0(saveImg, "03_00_kmean_tSNE_as_SUBCLUSTERED_Eleonora_Said_CLUSTERS_Chronic_ONLY_", expN, "_", cellLab[cT], ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotTSNE(elePlotChronic, colour_by="clusterS"))
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

	# finally, store the normalised objects
	saveRDS(object=kmeansSce, file=paste0(saveIn, "kmeans_K_", kS, "_clustered_data_", expN, ".Rds"))
}
