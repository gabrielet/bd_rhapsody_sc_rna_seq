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

# select experiment
expN <- "Onset"
#expN <- "Chronic"

print(paste0("analysing exp ", expN))

# set paths
raw_path <- "/home/gabriele/work/cbmc/scrnaseq/raw_data_new/"
saveIn <- paste0(raw_path, "EAE_", expN, "/")
saveImg <- paste0(saveIn, "figures/")

#check if saveIn directory exist, otherwise create it
ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
ifelse(dir.exists(saveImg), TRUE, dir.create(saveImg, recursive=F))

# set palette with 12 colours for colouring plots
colorsPal <- brewer.pal(12, "Paired")

# set cellTypes
cellType <- c("T.8|T.CD8", "Tgd", "Neutrophils", "\\(T.4|\\(T.CD4", "B cells", "DC|Macrophages|Microglia|Monocytes")
cellLab <-  c("CD8", "Tgd", "Neutrophils", "CD4", "B_cells", "Myeloid")

# load data, all cells
# best 1000 genes and k
#perC <- 1000 ; kZinb <- 20 ; zbWaved_original <- readRDS(paste0(saveIn, "full_zbwaved_data_", expN, "_K_", kZinb, "_top_", perC, ".Rds")) ; chosenRed <- "zinbwave" ; dimRedMet <- "zbwave" ; dimRedPar <- "TSNE" # for zimbwave

zbWaved_original <- readRDS(paste0(paste0(saveIn, "denoised_data_kmeans_", expN, ".Rds"))) ; chosenRed <- "corrected" ; dimRedMet <- "PCA" ; dimRedPar <- "PCA" # for PCA

# loop through cell type
for (cT in seq(1, length(cellType), by=1)){

	print("")
	print(paste0("analysing all the genes ", cellLab[cT], " for exp ", expN, " clustered after zbwave"))
	print("")

	###################################################################### get a specific cell type to perform subclustering ######################################################################
	zbWaved <- zbWaved_original[, grep(cellType[cT], zbWaved_original$pruned_fine)]
	###############################################################################################################################################################################################

	# transform SCE back to Seurat to use the clustering algorithms
	zbwSeurat <- as.Seurat(zbWaved, data="counts")

	# SNN graph
	zbwSeurat <- FindNeighbors(object=zbwSeurat, reduction=chosenRed, dims=1:20)

	# SLM clustering with different resolution
	seurat_clusters_0.1 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=0.1, random.seed=1)$seurat_clusters
	seurat_clusters_0.3 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=0.3, random.seed=1)$seurat_clusters

	slmClusters <- as.matrix(data.frame(slm_0.1=seurat_clusters_0.1, slm_0.3=seurat_clusters_0.3))

	#store clustering results
	zbWaved$slm_0.1 <- as.factor(slmClusters[, 1])
	zbWaved$slm_0.3 <- as.factor(slmClusters[, 2])
	
	# new labels that compress the fine name into something similar to main labels
	zbWaved$fine_to_main <- gsub(" ", "", unlist(lapply(strsplit(zbWaved$pruned_fine, "\\("), `[[`, 1)))

	# plot condition only
	png(paste0(saveImg, "03_subcl_02_", dimRedMet, "_", dimRedPar, "_SUBCLUSTERED_condition_only_", expN, "_", cellLab[cT], "_SLM_ZERO_TRE.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotReducedDim(zbWaved, dimred=dimRedPar, colour_by="condition"))
	dev.off()
	
	# plot condition only
	png(paste0(saveImg, "03_subcl_02_", dimRedMet, "_", dimRedPar, "_SUBCLUSTERED_by_clusters_", expN, "_", cellLab[cT], "_SLM_ZERO_TRE.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotReducedDim(zbWaved, dimred=dimRedPar, colour_by="slm_0.3"))
	dev.off()

	# plot
	png(paste0(saveImg, "03_subcl_05_", dimRedMet, "_", dimRedPar, "_SUBCLUSTERED_", expN, "_", cellLab[cT], "_SLM_ZERO_TRE.png"), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
		print(scater::plotReducedDim(zbWaved, dimred=dimRedPar, colour_by="pruned_fine") + theme(legend.position="bottom") + labs(title="By cell type") + guides(color=guide_legend("Cell type")))
	dev.off()

	###################################################################################################################

	# finally, export the data
	saveRDS(object=zbWaved, file=paste0(saveIn, "slm_clustered_data_", expN, "_", cellLab[cT], "_ONLY.Rds"))

	# FINALLY, add clusters info to cleanSce
	cleanSce <- readRDS(paste0(saveIn, paste0("QCed_data_", expN, ".Rds")))

	###################################################################### get a specific cell type to perform subclustering ######################################################################
	cleanSce <- cleanSce[, grep(cellType[cT], cleanSce$pruned_fine)]
	###############################################################################################################################################################################################

	cleanSce$slm_0.1 <- zbWaved$slm_0.1
	cleanSce$slm_0.3 <- zbWaved$slm_0.3

	# then export cleanSce with clustering info
	saveRDS(object=cleanSce, file=paste0(saveIn, "QCed_data_with_clusters_", dimRedMet, "_", expN, "_", cellLab[cT], "_ONLY.Rds"))
}
