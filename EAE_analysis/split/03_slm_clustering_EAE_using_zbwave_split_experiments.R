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

# set palette with 12 colours for colouring plots
colorsPal <- brewer.pal(12, "Paired")

print("")
print(paste0("analysing all the genes for exp ", expN, " clustered aftered zbwave"))
print("")

# load data, all cells
# best 1000 genes and k
#perC <- 1000 ; kZinb <- 20 ; zbWaved <- readRDS(paste0(saveIn, "full_zbwaved_data_", expN, "_K_", kZinb, "_top_", perC, ".Rds")) ; chosenRed <- "zinbwave" ; dimRedMet <- "zbwave" ; dimRedPar <- "TSNE" # for zimbwave

zbWaved <- readRDS(paste0(paste0(saveIn, "denoised_data_kmeans_", expN, ".Rds"))) ; chosenRed <- "corrected" ; dimRedMet <- "PCA" ; dimRedPar <- "PCA" # for PCA

####################################################################################################################################
# set cellTypes
cellType <- c("T.8|T.CD8", "Tgd", "Neutrophils", "\\(T.4|\\(T.CD4", "B cells", "DC|Macrophages|Microglia|Monocytes")
cellLab <-  c("CD8", "Tgd", "Neutrophils", "CD4", "B_cells", "Myeloid")

# loop through cell type
for (cT in seq(1, length(cellType), by=1)){

	############################################################## color a cell type ##############################################################
	
	zbWaved$celltype <- rep("Other Cells", length(zbWaved$pruned_fine))
	# set colours for plotting. remove RED and LIGHT-BLUE from the cols25() palette
	zbWaved$celltype[grep(cellLab[cT], zbWaved$pruned_fine)] <- cellLab[cT]
	zbWaved$celltype <- factor(zbWaved$celltype, levels=c("Other Cells", cellLab[cT]))
	
	png(paste0(saveImg, "03_08_", dimRedMet, "_", dimRedPar, "_ALL_CLUSTERS_", cellLab[cT], "_ONLY.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotReducedDim(zbWaved, dimred=dimRedPar, colour_by="celltype"))
	dev.off()
}

# transform SCE back to Seurat to use the clustering algorithms
zbwSeurat <- as.Seurat(zbWaved, data="counts")

# SNN graph
zbwSeurat <- FindNeighbors(object=zbwSeurat, reduction=chosenRed, dims=1:20)

# SLM clustering with different resolution
seurat_clusters_0.3 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=0.3, random.seed=1)$seurat_clusters

slmClusters <- as.matrix(data.frame(slm_0.3=seurat_clusters_0.3))

#store clustering results
zbWaved$slm_0.3 <- as.factor(slmClusters[, 1])

# new labels that compress the fine name into something similar to main labels
zbWaved$fine_to_main <- gsub(" ", "", unlist(lapply(strsplit(zbWaved$pruned_fine, "\\("), `[[`, 1)))

png(paste0(saveImg, "02_02_", dimRedMet, "_", dimRedPar, "_ALL_CLUSTERS_cell_types_", expN, "_SLM_ZERO_TRE.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotReducedDim(zbWaved, dimred=dimRedPar, colour_by="pruned_main") + theme(legend.position="right") + labs(title="By cell type") + guides(color=guide_legend("Cell type")))
dev.off()

# plot tSNE by condition
png(paste0(saveImg, "02_04_", dimRedMet, "_", dimRedPar, "_ALL_CLUSTERS_CONDITION_", expN, "_K_1000_SLM_ZERO_TRE.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotReducedDim(zbWaved, dimred=dimRedPar, colour_by="condition") + theme(legend.position="right") + labs(title="Conditions"))
dev.off()

# plot tSNE by cluster
png(paste0(saveImg, "02_05_", dimRedMet, "_", dimRedPar, "_ALL_CLUSTERS_CLUSTER_", expN, "_K_1000_SLM_ZERO_TRE.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotReducedDim(zbWaved, dimred=dimRedPar, colour_by="slm_0.3") + theme(legend.position="right") + labs(title="Clusters"))
dev.off()

# plot tSNE by cluster
png(paste0(saveImg, "02_06_", dimRedMet, "_", dimRedPar, "_ALL_CLUSTERS_CLUSTER_con_etichette", expN, "_K_1000_SLM_ZERO_TRE.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotReducedDim(zbWaved, dimred=dimRedPar, colour_by="slm_0.3", text_by="slm_0.3") + theme(legend.position="none") + labs(title="Clusters"))
dev.off()

# finally, export the data
saveRDS(object=zbWaved, file=paste0(saveIn, "slm_clustered_data_", expN, ".Rds"))

# FINALLY, add clusters info to cleanSce
cleanSce <- readRDS(paste0(saveIn, paste0("QCed_data_", expN, ".Rds")))

cleanSce$slm_0.3 <- zbWaved$slm_0.3

# then export cleanSce with clustering info
saveRDS(object=cleanSce, file=paste0(saveIn, "QCed_data_with_clusters_", dimRedMet, "_", expN, ".Rds"))
