#' ---
#' title: "SLM clustering for reduced data"
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
library("SingleCellExperiment")
library("Seurat")
library("pals")

# set paths
rawPath <- "/home/gabriele/work/cbmc/scrnaseq/raw_data_new/" ; expN <- "EAE_batched"

saveIn <- paste0(rawPath, "batched_EAE_PCA/")
saveImg <- paste0(saveIn, "figures/")

print(paste0("analysing exp ", expN))

#check if saveIn directory exist, otherwise create it
ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
ifelse(dir.exists(saveImg), TRUE, dir.create(saveImg, recursive=F))

# set palette with 12 colours for colouring plots
colorsPal <- brewer.pal(12, "Paired")

# set cellTypes
cellType <- c("T.8|T.CD8", "Tgd", "Neutrophils", "\\(T.4|\\(T.CD4", "B cells", "DC|Macrophages|Microglia|Monocytes")
cellLab <-  c("CD8", "Tgd", "Neutrophils", "CD4", "B_cells", "Myeloid")

# load data, all cells
pca_ed_original <- readRDS(paste0(saveIn, "denoised_data_kmeans_", expN, ".Rds"))

# loop through cell type
for (cT in seq(1, length(cellType), by=1)){

	print("")
	print(paste0("analysing all the genes ", cellLab[cT], " for exp ", expN, " clustered after PCA"))
	print("")

	###################################################################### get a specific cell type to perform subclustering ######################################################################
	pca_ed <- pca_ed_original[, grep(cellType[cT], pca_ed_original$pruned_fine)]
	###############################################################################################################################################################################################

	# transform SCE back to Seurat to use the clustering algorithms
	zbwSeurat <- as.Seurat(pca_ed, data="counts")

	# SNN graph
	#zbwSeurat <- FindNeighbors(object=zbwSeurat, reduction="zinbwave", dims=1:20) #QUESTO QUELLO CHE SI FACEVA CON ZINBWAVE
	zbwSeurat <- FindNeighbors(object=zbwSeurat, reduction="corrected", dims=1:20)

	# SLM clustering with different resolution
	seurat_clusters_0.1 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=0.1, random.seed=1)$seurat_clusters
	seurat_clusters_0.3 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=0.3, random.seed=1)$seurat_clusters
	seurat_clusters_0.5 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=0.5, random.seed=1)$seurat_clusters
	seurat_clusters_1 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=1, random.seed=1)$seurat_clusters

	slmClusters <- as.matrix(data.frame(slm_0.1=seurat_clusters_0.1, slm_0.3=seurat_clusters_0.3, slm_0.5=seurat_clusters_0.5, slm_1=seurat_clusters_1))

	#store clustering results
	pca_ed$slm_0.1 <- as.factor(slmClusters[, 1])
	pca_ed$slm_0.3 <- as.factor(slmClusters[, 2])
	pca_ed$slm_0.5 <- as.factor(slmClusters[, 3])
	pca_ed$slm_1 <- as.factor(slmClusters[, 4])
	
	# new labels that compress the fine name into something similar to main labels
	pca_ed$fine_to_main <- gsub(" ", "", unlist(lapply(strsplit(pca_ed$pruned_fine, "\\("), `[[`, 1)))

	# plot condition only
	png(paste0(saveImg, "03_subcl_02_PCA_tSNE_SUBCLUSTERED_condition_only_", expN, "_", cellLab[cT], ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotReducedDim(pca_ed, dimred="TSNE", colour_by="condition"))
	dev.off()
	
	# plot condition only
	png(paste0(saveImg, "03_subcl_02_PCA_tSNE_SUBCLUSTERED_by_clusters_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotReducedDim(pca_ed, dimred="TSNE", colour_by="slm_0.1"))
	dev.off()

	# plot region and condition
	png(paste0(saveImg, "03_subcl_03_PCA_tSNE_SUBCLUSTERED_conditionBYstage_", expN, "_", cellLab[cT], ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotReducedDim(pca_ed, dimred="TSNE", colour_by="conditionBYstage") + theme(legend.position="right") + labs(title="By condition") + guides(color=guide_legend("Condition")))
	dev.off()

	# plot together
	png(paste0(saveImg, "03_subcl_04_PCA_tSNE_SUBCLUSTERED_stage_only_", expN, "_", cellLab[cT], ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotReducedDim(pca_ed, dimred="TSNE", colour_by="stage") + theme(legend.position="right") + labs(title="By condition") + guides(color=guide_legend("Condition")))
	dev.off()

	# plot
	png(paste0(saveImg, "03_subcl_05_PCA_tSNE_SUBCLUSTERED_", expN, "_", cellLab[cT], ".png"), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
		print(scater::plotReducedDim(pca_ed, dimred="TSNE", colour_by="pruned_fine") + theme(legend.position="bottom") + labs(title="By cell type") + guides(color=guide_legend("Cell type")))
	dev.off()
	
	###################################################################################################################
	# ELEONORA PLOT BRAIN ONLY
	elePlotOnset <- pca_ed[, which(pca_ed$stage=="Onset")]

	# plot
	png(paste0(saveImg, "03_00_PCA_tSNE_as_SUBCLUSTERED_Eleonora_Said_CLUSTERS_Onset_ONLY_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotReducedDim(elePlotOnset, dimred="TSNE", colour_by="slm_0.1"))
	dev.off()

	# ELEONORA PLOT MENINGES ONLY
	elePlotChronic <- pca_ed[, which(pca_ed$stage=="Chronic")]

	# plot
	png(paste0(saveImg, "03_00_PCA_tSNE_as_SUBCLUSTERED_Eleonora_Said_CLUSTERS_Chronic_ONLY_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotReducedDim(elePlotChronic, dimred="TSNE", colour_by="slm_0.1"))
	dev.off()

	###################################################################################################################

	# finally, export the data
	saveRDS(object=pca_ed, file=paste0(saveIn, "slm_clustered_data_", expN, "_", cellLab[cT], "_ONLY.Rds"))

	# FINALLY, add clusters info to cleanSce
	cleanSce <- readRDS(paste0(saveIn, paste0("QCed_data_", expN, ".Rds")))

	###################################################################### get a specific cell type to perform subclustering ######################################################################
	cleanSce <- cleanSce[, grep(cellType[cT], cleanSce$pruned_fine)]
	###############################################################################################################################################################################################

	cleanSce$slm_0.1 <- pca_ed$slm_0.1
	cleanSce$slm_0.3 <- pca_ed$slm_0.3
	cleanSce$slm_0.5 <- pca_ed$slm_0.5
	cleanSce$slm_1 <- pca_ed$slm_1

	# then export cleanSce with clustering info
	saveRDS(object=cleanSce, file=paste0(saveIn, "QCed_data_with_clusters_", expN, "_", cellLab[cT], "_ONLY.Rds"))
}
