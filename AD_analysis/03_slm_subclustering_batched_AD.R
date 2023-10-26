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
library("SingleCellExperiment")
library("Seurat")
library("pals")

# set paths
rawPath <- "/home/gabriele/work/cbmc/scrnaseq/raw_data_new/" ; expN <- "AD_batched"

saveIn <- paste0(rawPath, "batched_AD/")
saveImg <- paste0(saveIn, "figures/")

print(paste0("analysing exp ", expN))

#check if saveIn directory exist, otherwise create it
ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
ifelse(dir.exists(saveImg), TRUE, dir.create(saveImg, recursive=F))

# set palette with 12 colours for colouring plots
colorsPal <- brewer.pal(12, "Paired")

# number of highly variable genes used for zinbwaving the data
perC <- 1000
kZinb <- 20

# set cellTypes
cellType <- c("T.8|T.CD8", "Tgd", "Neutrophils", "\\(T.4|\\(T.CD4", "B cells", "DC|Macrophages|Microglia|Monocytes")
cellLab <-  c("CD8", "Tgd", "Neutrophils", "CD4", "B_cells", "Myeloid")

# load data, all cells
zbWaved_original <- readRDS(paste0(saveIn, "full_zbwaved_data_", expN, "_K_", kZinb, "_top_", perC, ".Rds"))

# recompute variable conditionBYgender
zbWaved_original$conditionBYgender <- paste(zbWaved_original$condition, zbWaved_original$gender, sep="_")

# loop through cell type
for (cT in seq(1, length(cellType), by=1)){

	print("")
	print(paste0("analysing all the genes ", cellLab[cT], " for exp ", expN, " clustered with zinbwave"))
	print("")

	###################################################################### get a specific cell type to perform subclustering ######################################################################
	zbWaved <- zbWaved_original[, grep(cellType[cT], zbWaved_original$pruned_fine)]
	###############################################################################################################################################################################################

	# transform SCE back to Seurat to use the clustering algorithms
	zbwSeurat <- as.Seurat(zbWaved, data="counts")

	# SNN graph
	zbwSeurat <- FindNeighbors(object=zbwSeurat, reduction="zinbwave", dims=1:20)

	# SLM clustering with different resolution
	seurat_clusters_0.1 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=0.1, random.seed=1)$seurat_clusters
	seurat_clusters_0.3 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=0.3, random.seed=1)$seurat_clusters
	seurat_clusters_0.5 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=0.5, random.seed=1)$seurat_clusters
	seurat_clusters_1 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=1, random.seed=1)$seurat_clusters

	slmClusters <- as.matrix(data.frame(slm_0.1=seurat_clusters_0.1, slm_0.3=seurat_clusters_0.3, slm_0.5=seurat_clusters_0.5, slm_1=seurat_clusters_1))

	#store clustering results in zinbSce
	zbWaved$slm_0.1 <- as.factor(slmClusters[, 1])
	zbWaved$slm_0.3 <- as.factor(slmClusters[, 2])
	zbWaved$slm_0.5 <- as.factor(slmClusters[, 3])
	zbWaved$slm_1 <- as.factor(slmClusters[, 4])

	# new labels that compress the fine name into something similar to main labels
	zbWaved$fine_to_main <- gsub(" ", "", unlist(lapply(strsplit(zbWaved$pruned_fine, "\\("), `[[`, 1)))

	# plot condition only
	png(paste0(saveImg, "03_subcl_02_zinbwave_tSNE_SUBCLUSTERED_condition_only_", expN, "_K_", kZinb, "_", cellLab[cT], "_SLM_ZERO_CINQUE.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotTSNE(zbWaved, colour_by="condition"))
	dev.off()
	
	# plot condition only
	png(paste0(saveImg, "03_subcl_02_zinbwave_tSNE_SUBCLUSTERED_by_clusters_", expN, "_K_", kZinb, "_", cellLab[cT], "_SLM_ZERO_CINQUE.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotTSNE(zbWaved, colour_by="slm_0.5"))
	dev.off()

	# plot together
	png(paste0(saveImg, "03_subcl_04_zinbwave_tSNE_SUBCLUSTERED_gender_only_", expN, "_K_", kZinb, "_", cellLab[cT], "_SLM_ZERO_CINQUE.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotReducedDim(zbWaved, dimred="TSNE", colour_by="gender") + theme(legend.position="right") + labs(title="By condition") + guides(color=guide_legend("Condition")))
	dev.off()

	# plot
	png(paste0(saveImg, "03_subcl_05_zinbwave_tSNE_SUBCLUSTERED_", expN, "_K_", kZinb, "_", cellLab[cT], "_SLM_ZERO_CINQUE.png"), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
		print(scater::plotReducedDim(zbWaved, dimred="TSNE", colour_by="pruned_fine") + theme(legend.position="bottom") + labs(title="By cell type") + guides(color=guide_legend("Cell type")))
	dev.off()
	
	# print tSNE with all the cells
	png(paste0(saveImg, "03_subcl_01_zinbwave_tSNE_SUBCLUSTERED_conditionBYgender_", expN, "_K_", kZinb, "_", cellLab[cT], "_SLM_ZERO_CINQUE.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(print(scater::plotReducedDim(zbWaved, dimred="TSNE", colour_by="conditionBYgender") + theme(legend.position="right") + labs(title="By condition") + guides(color=guide_legend("Condition"))))
	dev.off()
	
	###################################################################################################################
	# ELEONORA PLOT BRAIN ONLY
	elePlotMale <- zbWaved[, which(zbWaved$gender=="M")]

	# plot
	png(paste0(saveImg, "03_00_zinbwave_tSNE_as_SUBCLUSTERED_Eleonora_Said_CLUSTERS_MALES_ONLY_", expN, "_K_", kZinb, "_", cellLab[cT], "_SLM_ZERO_CINQUE.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotTSNE(elePlotMale, colour_by="slm_0.5"))
	dev.off()

	# ELEONORA PLOT MENINGES ONLY
	elePlotFemale <- zbWaved[, which(zbWaved$gender=="F")]

	# plot
	png(paste0(saveImg, "03_00_zinbwave_tSNE_as_SUBCLUSTERED_Eleonora_Said_CLUSTERS_FEMALES_ONLY_", expN, "_K_", kZinb, "_", cellLab[cT], "_SLM_ZERO_CINQUE.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotTSNE(elePlotFemale, colour_by="slm_0.5"))
	dev.off()

	###################################################################################################################
	
	# finally, export the data
	saveRDS(object=zbWaved, file=paste0(saveIn, "slm_clustered_data_", expN, "_K_", kZinb, "_top_", perC, "_", cellLab[cT], "_ONLY.Rds"))

	# FINALLY, add clusters info to cleanSce
	cleanSce <- readRDS(paste0(saveIn, paste0("QCed_data_", expN, ".Rds")))

	###################################################################### get a specific cell type to perform subclustering ######################################################################
	cleanSce <- cleanSce[, grep(cellType[cT], cleanSce$pruned_fine)]
	###############################################################################################################################################################################################

	cleanSce$slm_0.1 <- zbWaved$slm_0.1
	cleanSce$slm_0.3 <- zbWaved$slm_0.3
	cleanSce$slm_0.5 <- zbWaved$slm_0.5
	cleanSce$slm_1 <- zbWaved$slm_1

	# then export cleanSce with clustering info
	saveRDS(object=cleanSce, file=paste0(saveIn, "QCed_data_with_clusters_", expN, "_K_", kZinb, "_top_", perC, "_", cellLab[cT], "_ONLY.Rds"))
}
