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

print("")
print(paste0("analysing all the genes for exp ", expN, " clustered aftered PCA"))
print("")

# load data, all cells
pca_ed <- readRDS(paste0(saveIn, "denoised_data_kmeans_", expN, ".Rds"))

####################################################################################################################################
# set cellTypes
cellType <- c("T.8|T.CD8", "Tgd", "Neutrophils", "\\(T.4|\\(T.CD4", "B cells", "DC|Macrophages|Microglia|Monocytes")
cellLab <-  c("CD8", "Tgd", "Neutrophils", "CD4", "B_cells", "Myeloid")

# loop through cell type
for (cT in seq(1, length(cellType), by=1)){

	# ELEONORA PLOT BRAIN ONLY
	elePlotChronic <- pca_ed[, which(pca_ed$stage=="Chronic")]

	# find interesting cells
	poSBRN <- grep(cellType[cT], elePlotChronic$pruned_fine)

	# get labels ready
	elePlotChronic$onlyInt <- "allOther"
	elePlotChronic$onlyInt[poSBRN] <- cellLab[cT]

	# plot
	png(paste0(saveImg, "03_00_PCA_tSNE_ALL_CLUSTERS_as_Eleonora_Said_", cellLab[cT], "_Chronic_ONLY_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotReducedDim(elePlotChronic, dimred="TSNE", colour_by="onlyInt") + theme(legend.position="right") + labs(title="By condition") + guides(color=guide_legend("Condition")))
	dev.off()

	# ELEONORA PLOT MENINGES ONLY
	elePlotOnset <- pca_ed[, which(pca_ed$stage=="Onset")]

	# find interesting cells
	poSMNG <- grep(cellType[cT], elePlotOnset$pruned_fine)

	# get labels ready
	elePlotOnset$onlyInt <- "allOther"
	elePlotOnset$onlyInt[poSMNG] <- cellLab[cT]

	# plot
	png(paste0(saveImg, "03_00_PCA_tSNE_ALL_CLUSTERS_as_Eleonora_Said_", cellLab[cT], "_ONSET_ONLY_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotReducedDim(elePlotOnset, dimred="TSNE", colour_by="onlyInt") + theme(legend.position="right") + labs(title="By condition") + guides(color=guide_legend("Condition")))
	dev.off()

	# ELEONORA PLOT MENINGES + BRN
	elePlotALL <- pca_ed

	# find interesting cells
	poSALL <- grep(cellType[cT], elePlotALL$pruned_fine)

	# get labels ready
	elePlotALL$onlyInt <- "allOther"
	elePlotALL$onlyInt[poSALL] <- cellLab[cT]

	# plot
	png(paste0(saveImg, "03_00_PCA_tSNE_ALL_CLUSTERS_as_Eleonora_Said_", cellLab[cT], "_ONSET_and_CHRONIC_ONLY_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotReducedDim(elePlotALL, dimred="TSNE", colour_by="onlyInt") + theme(legend.position="right") + labs(title="By condition") + guides(color=guide_legend("Condition")))
	dev.off()

	############################################################## a cell type ##############################################################

	# set colours for plotting. remove RED and LIGHT-BLUE from the cols25() palette
	palettE <- unname(cols25()[-c(2, 8)][1:(length(unique(pca_ed$pruned_fine))-1)])
	# set a nice red for celltype
	names(palettE) <- unique(pca_ed$pruned_fine)[!is.na(unique(pca_ed$pruned_fine))]
	a_cell_type <- grep(cellType[cT], names(palettE))
	palettE[a_cell_type] <- "#E41A1C"

	png(paste0(saveImg, "03_08_PCA_tSNE_ALL_CLUSTERS_", cellLab[cT], "_ONLY.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotReducedDim(pca_ed, dimred="TSNE", colour_by="pruned_fine") + scale_color_manual(values=palettE) + guides(color=guide_legend("Cell type")) + theme(legend.position="none") + labs(title="Clustered by cell type"))
	dev.off()

}

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

# print tSNE with all the cells
png(paste0(saveImg, "02_01_PCA_tSNE_ALL_CLUSTERS_conditionBYstage_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotReducedDim(pca_ed, dimred="TSNE", colour_by="conditionBYstage") + theme(legend.position="right") + labs(title="By condition") + guides(color=guide_legend("Condition")))
dev.off()

png(paste0(saveImg, "02_02_PCA_tSNE_ALL_CLUSTERS_cell_types_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotReducedDim(pca_ed, dimred="TSNE", colour_by="pruned_main") + theme(legend.position="right") + labs(title="By cell type") + guides(color=guide_legend("Cell type")))
dev.off()

# plot tSNE with all the cells
png(paste0(saveImg, "02_03_PCA_tSNE_ALL_CLUSTERS_stage_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotReducedDim(pca_ed, dimred="TSNE", colour_by="stage") + theme(legend.position="right") + labs(title="By condition") + guides(color=guide_legend("Condition")))
dev.off()

# plot tSNE by condition
png(paste0(saveImg, "02_04_PCA_tSNE_ALL_CLUSTERS_CONDITION_", expN, "_K_1000.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotReducedDim(pca_ed, dimred="TSNE", colour_by="condition") + theme(legend.position="right") + labs(title="Conditions"))
dev.off()

# plot tSNE by cluster
png(paste0(saveImg, "02_05_PCA_tSNE_ALL_CLUSTERS_CLUSTER_", expN, "_K_1000_SLM_ZERO_UNO.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotReducedDim(pca_ed, dimred="TSNE", colour_by="slm_0.1") + theme(legend.position="right") + labs(title="Clusters"))
dev.off()

# plot tSNE by cluster
png(paste0(saveImg, "02_06_PCA_tSNE_ALL_CLUSTERS_CLUSTER_con_etichette", expN, "_K_1000_SLM_ZERO_UNO.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotReducedDim(pca_ed, dimred="TSNE", colour_by="slm_0.1", text_by="slm_0.1") + theme(legend.position="none") + labs(title="Clusters"))
dev.off()

# finally, export the data
saveRDS(object=pca_ed, file=paste0(saveIn, "slm_clustered_data_", expN, ".Rds"))

# FINALLY, add clusters info to cleanSce
cleanSce <- readRDS(paste0(saveIn, paste0("QCed_data_", expN, ".Rds")))

cleanSce$slm_0.1 <- pca_ed$slm_0.1
cleanSce$slm_0.3 <- pca_ed$slm_0.3
cleanSce$slm_0.5 <- pca_ed$slm_0.5
cleanSce$slm_1 <- pca_ed$slm_1

# then export cleanSce with clustering info
saveRDS(object=cleanSce, file=paste0(saveIn, "QCed_data_with_clusters_", expN, ".Rds"))
