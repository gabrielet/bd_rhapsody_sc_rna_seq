#' ---
#' title: "Analysing GammaDelta cells in AD exps"
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
library("limma")
library("clusterProfiler")
library("org.Mm.eg.db")
library("RColorBrewer")
library("scater")
library("scran")
library("Seurat")
library("SingleCellExperiment")
library("openxlsx")
library("fgsea")
library("pathview")
library("ReactomePA")
library("DESeq2")
library("enrichplot")
library("DOSE")
library("HelpersMG")
library("gg.gap")
library("ddpcr")
library("data.table")
library("ggpubr")
library("edgeR")

# import plotting functions
source("../00_functions.R")

# initialise label
lbl <- 1

# set path
rawPath <- "/home/gabriele/work/cbmc/scrnaseq/raw_data_new/" ; expN <- "EAE_batched"
saveIn <- paste0(rawPath, "batched_EAE_PCA/")

# NOTE: we are using QCed data file which logcounts were normalised using sizeFactors
# but this data file contains also the information from the SLM clustering
# performed using the pca_ed data, in order to have a good clustering but also all
# the genes

kS <- 15

# get normally clustered data
clusteredSce <- readRDS(paste0(saveIn, "QCed_data_with_clusters_", expN, ".Rds"))

# create fine_to_main labels as it was done for the zbWaved object in 03
clusteredSce$fine_to_main <- gsub(" ", "", unlist(lapply(strsplit(clusteredSce$pruned_fine, "\\("), `[[`, 1)))

##################### FULL DATASET

tab_alldata_labels <- table(clusteredSce$fine_to_main, clusteredSce$condition, clusteredSce$stage)

# get EAE full table ready to be exported
etabFEAE <- as.data.frame.matrix(tab_alldata_labels[, , "Chronic"])
eTabrnFEAE <- cbind.data.frame(cellType=rownames(etabFEAE), etabFEAE)
eTabcnFEAE <- rbind.data.frame(c("cluster", colnames(etabFEAE)), eTabrnFEAE)

# and CTRL
etabFCTRL <- as.data.frame.matrix(tab_alldata_labels[, , "Onset"])
eTabrnFCTRL <- cbind.data.frame(cellType=rownames(etabFCTRL), etabFCTRL)
eTabcnFCTRL <- rbind.data.frame(c("cluster", colnames(etabFCTRL)), eTabrnFCTRL)

# export the full table
write.xlsx(eTabcnFEAE, paste0(saveIn, "chronic_cluster_content_full_labels_ALL_DATASET_exp_", expN, "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=F)
write.xlsx(eTabcnFCTRL, paste0(saveIn, "onset_cluster_content_full_labels_ALL_DATASET_exp_", expN, "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=F)

# set cellTypes
cellType <- c("T.8|T.CD8", "Tgd", "Neutrophils", "\\(T.4|\\(T.CD4", "B cells", "DC|Macrophages|Microglia|Monocytes")
cellLab <-  c("CD8", "Tgd", "Neutrophils", "CD4", "B_cells", "Myeloid")

# loop through cell type
for (cT in seq(1, length(cellLab), by=1)){

	print("")
	print(paste0("analysing all the genes ", cellLab[cT], " for exp ", expN, " clustered with Zinbwave"))
	print("")
	
	# set paths
	saveUpDown <- paste0(saveIn, "enrichment_", cellLab[cT], "/UpDown/")
	saveImg <- paste0(saveIn, "figures/analysis_", cellLab[cT], "/")

	# check if saveIn directory exist, otherwise create it
	ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
	ifelse(dir.exists(saveUpDown), TRUE, dir.create(saveUpDown, recursive=T))
	ifelse(dir.exists(saveImg), TRUE, dir.create(saveImg, recursive=T))

	# set palette with 12 colours for colouring plots
	colorsPal <- brewer.pal(12, "Paired")
	
	# plotting heatmap with cell type
	tab_fine_to_main <- table(Assigned=clusteredSce$fine_to_main, Cluster=clusteredSce$slm_0.1, Condition=clusteredSce$condition)

	# get EAE table ready to be exported
	eTabEAE <- as.data.frame.matrix(tab_fine_to_main[, , "EAE"])
	eTabRNEAE <- cbind.data.frame(cellType=rownames(eTabEAE), eTabEAE)
	eTabCNEAE <- rbind.data.frame(c("cluster", colnames(eTabEAE)), eTabRNEAE)

	# same for CTRL
	eTabCTRL <- as.data.frame.matrix(tab_fine_to_main[, , "CTRL"])
	eTabRNCTRL <- cbind.data.frame(cellType=rownames(eTabCTRL), eTabCTRL)
	eTabCNCTRL <- rbind.data.frame(c("cluster", colnames(eTabCTRL)), eTabRNCTRL)

	# export both tables
	write.xlsx(eTabCNEAE, paste0(saveImg, "EAE_cluster_content_ALL_CLUSTERS_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=F)

	write.xlsx(eTabCNCTRL, paste0(saveImg, "CTRL_cluster_content_ALL_CLUSTERS_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=F)

	##################### FULL

	tab_full_labels <- table(Assigned=clusteredSce$pruned_fine, Cluster=clusteredSce$slm_0.1, Condition=clusteredSce$condition)

	# get EAE full table ready to be exported
	etabFEAE <- as.data.frame.matrix(tab_full_labels[, , "EAE"])
	eTabrnFEAE <- cbind.data.frame(cellType=rownames(etabFEAE), etabFEAE)
	eTabcnFEAE <- rbind.data.frame(c("cluster", colnames(etabFEAE)), eTabrnFEAE)

	# and CTRL
	etabFCTRL <- as.data.frame.matrix(tab_full_labels[, , "CTRL"])
	eTabrnFCTRL <- cbind.data.frame(cellType=rownames(etabFCTRL), etabFCTRL)
	eTabcnFCTRL <- rbind.data.frame(c("cluster", colnames(etabFCTRL)), eTabrnFCTRL)

	# export the full table
	write.xlsx(eTabcnFEAE, paste0(saveImg, "EAE_cluster_content_full_labels_ALL_CLUSTERS_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=F)

	write.xlsx(eTabcnFCTRL, paste0(saveImg, "CTRL_cluster_content_full_labels_ALL_CLUSTERS_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=F)

	##################### FULL STAGE

	tab_STAGE_labels <- table(clusteredSce$pruned_fine, clusteredSce$condition, clusteredSce$stage)

	# get EAE full table ready to be exported
	etabFEAE <- as.data.frame.matrix(tab_STAGE_labels[, , "Chronic"])
	eTabrnFEAE <- cbind.data.frame(cellType=rownames(etabFEAE), etabFEAE)
	eTabcnFEAE <- rbind.data.frame(c("cluster", colnames(etabFEAE)), eTabrnFEAE)

	# and CTRL
	etabFCTRL <- as.data.frame.matrix(tab_STAGE_labels[, , "Onset"])
	eTabrnFCTRL <- cbind.data.frame(cellType=rownames(etabFCTRL), etabFCTRL)
	eTabcnFCTRL <- rbind.data.frame(c("cluster", colnames(etabFCTRL)), eTabrnFCTRL)

	# export the full table
	write.xlsx(eTabcnFEAE, paste0(saveImg, "chronic_cluster_content_full_labels_ALL_CLUSTERS_BY_STAGE_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=F)

	write.xlsx(eTabcnFCTRL, paste0(saveImg, "onset_cluster_content_full_labels_ALL_CLUSTERS_BY_STAGE_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=F)

	########################################################################################################################################################################

	# get data set conditions, i.e. CTRL and EAE
	conD <- unique(clusteredSce$condition)

	# search for cellType
	cTPos <- grep(cellType[cT], clusteredSce$pruned_fine)

	print(paste0(cellLab[cT], " are found in ", length(unique(clusteredSce$slm_0.1[cTPos])), " different clusters"))

	fR <- table(clusteredSce$slm_0.1[cTPos])
	clustS <- names(fR[which(fR!=0)])

	################################################# SOME PLOTTING

	print(fR[which(fR!=0)])

	# Cluster composition plots using main labels

	df <- data.frame(table(clusteredSce$slm_0.1, clusteredSce$fine_to_main, clusteredSce$condition))

	colnames(df) <- c("Cluster", "CellType", "Condition", "Frequency")
	lims <- levels(factor(clusteredSce$fine_to_main))
	lims <- lims[order(table(clusteredSce$fine_to_main))]

	# plot!
	p <- plot_membership(df, cond=conD[2], col=colorsPal[2])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[2], "_cluster_composition_ALL_CLUSTERS_", expN, "_SLM_ZERO_UNO.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()

	p <- plot_membership(df, cond=conD[1], col=colorsPal[6])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[1], "_cluster_composition_ALL_CLUSTERS_", expN, "_SLM_ZERO_UNO.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()
	
	# Cluster composition plots SPLITTNG BY STAGE AND CONDITION

	conD <- unique(clusteredSce$conditionBYstage)

	dfReg <- data.frame(table(clusteredSce$slm_0.1, clusteredSce$fine_to_main, clusteredSce$conditionBYstage))
	colnames(dfReg) <- c("Cluster", "CellType", "Condition", "Frequency")
	
	write.xlsx(dfReg, paste0(saveIn, "composition_information_ALL_CLUSTERS_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)

	lims <- levels(factor(clusteredSce$fine_to_main))
	lims <- lims[order(table(clusteredSce$fine_to_main))]

	# plot!
	p <- plot_membership(dfReg, cond=conD[2], col=colorsPal[4])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[2], "_cluster_composition_ALL_CLUSTERS_", cellLab[cT], "_", expN, "_SLM_ZERO_UNO.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()

	p <- plot_membership(dfReg, cond=conD[1], col=colorsPal[8])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[1], "_cluster_composition_ALL_CLUSTERS_", cellLab[cT], "_", expN, "_SLM_ZERO_UNO.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()
	
	# plot!
	p <- plot_membership(dfReg, cond=conD[3], col=colorsPal[10])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[3], "_cluster_composition_ALL_CLUSTERS_", cellLab[cT], "_", expN, "_SLM_ZERO_UNO.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()

	p <- plot_membership(dfReg, cond=conD[4], col=colorsPal[12])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[4], "_cluster_composition_ALL_CLUSTERS_", cellLab[cT], "_", expN, "_SLM_ZERO_UNO.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()

	print(paste0("######## SEARCHING FOR ", cellType[cT], " IN THE WHOLE DATASET ########"))

	################################################# BACK TO THE EXPERIMENT

	# and subset the original single cell experiment using the name of the cellType cells.
	cTSubSet <- clusteredSce[, cTPos]
	
	# CONTE GREZZE ELEONORA
	raw_counts <- assay(cTSubSet, "counts")
	colnames(raw_counts) <- paste(colnames(cTSubSet), cTSubSet$condition, cTSubSet$stage, cTSubSet$pruned_fine, "cluster", cTSubSet$slm_0.1, sep="_")
	with_row_names <- cbind.data.frame(rownames(raw_counts), raw_counts)
	colnames(with_row_names) <- c("genes", colnames(raw_counts))
	
	write.table(with_row_names, paste0(saveIn, "counts_CONTE_GREZZE_CLUSTERING_GENERALE_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
	#######################

	################################################# EXPORT THE NEW TABLE ELEONORA ASKED FEBRUARY #################################################

	superfull_table <- assay(cTSubSet, "logcounts")
	colnames(superfull_table) <- paste(colnames(cTSubSet), cTSubSet$condition, cTSubSet$stage, cellLab[cT], "cluster", cTSubSet$slm_0.1, sep="_")
	with_row_names <- cbind.data.frame(rownames(superfull_table), superfull_table)
	colnames(with_row_names) <- c("genes", colnames(superfull_table))
	
	write.table(with_row_names, paste0(saveIn, "logcounts_all_the_selected_cells_version_barbara_elena_eleonora_only_clustering_generale_PRIMA_DI_LOGNORMCOUNTS_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")

	# remove genes with too many zeros
	cTSubSet <- cTSubSet[rowMeans(assays(cTSubSet)$counts > 0) > 0.05,]

	# alternative to the filter above
	#	cTSubSet <- cTSubSet[rowSums(counts(cTSubSet) > 0) > 10,]

	# set minmean equal to 0.5 since this is the value used in the previous analysis, i.e. the one in the manuscript
	# compute sizeFactors
	cTSubSet <- computeSumFactors(cTSubSet, min.mean=0.5)
	cTSubSet <- logNormCounts(cTSubSet)
	
	################################################# EXPORT THE NEW TABLE ELENA ASKED #################################################

	superfull_table <- assay(cTSubSet, "logcounts")
	colnames(superfull_table) <- paste(colnames(cTSubSet), cTSubSet$condition, cTSubSet$stage, cTSubSet$pruned_fine, "cluster", cTSubSet$slm_0.1, sep="_")
	with_row_names <- cbind.data.frame(rownames(superfull_table), superfull_table)
	colnames(with_row_names) <- c("genes", colnames(superfull_table))
	
	write.table(with_row_names, paste0(saveIn, "logcounts_clustering_generale_version_barbara_elena_eleonora_only_ALL_CLUSTERS_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
	
	################################################# EXPORT THE NEW TABLE ELENA ASKED #################################################

	
	# exporting logcounts for Elena's violins
	subCTRL <- cTSubSet[, which(cTSubSet$condition=="CTRL")]
	subEAE <- cTSubSet[, which(cTSubSet$condition=="EAE")]

	logCsCTRL <- cbind.data.frame(genes=rownames(subCTRL), as.data.frame(logcounts(subCTRL)))
	logCsEAE <- cbind.data.frame(genes=rownames(subEAE), as.data.frame(logcounts(subEAE)))
	
	write.table(logCsCTRL, paste0(saveIn, "logcounts_ALL_CLUSTERS_CTRL_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
	write.table(logCsEAE, paste0(saveIn, "logcounts_ALL_CLUSTERS_EAE_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")

#	write.xlsx(logCsCTRL, paste0(saveIn, "logcounts_ALL_CLUSTERS_CTRL_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)
#	write.xlsx(logCsEAE, paste0(saveIn, "logcounts_ALL_CLUSTERS_EAE_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)

	# exporting logcounts by condition for Elena's violins
	subCTRL <- cTSubSet[, which(cTSubSet$condition=="CTRL")]
	subEAE <- cTSubSet[, which(cTSubSet$condition=="EAE")]
	
	# Chronic ONLY
	subCTRL_M <- subCTRL[, which(subCTRL$stage=="Chronic")]
	subEAE_M <- subEAE[, which(subEAE$stage=="Chronic")]
	
	logCsCTRL_M <- cbind.data.frame(genes=rownames(subCTRL_M), as.data.frame(logcounts(subCTRL_M)))
	logCsEAE_M <- cbind.data.frame(genes=rownames(subEAE_M), as.data.frame(logcounts(subEAE_M)))

	write.table(logCsCTRL_M, paste0(saveIn, "logcounts_ALL_CLUSTERS_CTRL_M_ONLY_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
	write.table(logCsEAE_M, paste0(saveIn, "logcounts_ALL_CLUSTERS_EAE_M_ONLY_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
	
#	write.xlsx(logCsCTRL_M, paste0(saveIn, "logcounts_ALL_CLUSTERS_CTRL_M_ONLY_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)
#	write.xlsx(logCsEAE_M, paste0(saveIn, "logcounts_ALL_CLUSTERS_EAE_M_ONLY_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)

	# Onset ONLY
	subCTRL_F <- subCTRL[, which(subCTRL$stage=="Onset")]
	subEAE_F <- subEAE[, which(subEAE$stage=="Onset")]
	
	logCsCTRL_F <- cbind.data.frame(genes=rownames(subCTRL_F), as.data.frame(logcounts(subCTRL_F)))
	logCsEAE_F <- cbind.data.frame(genes=rownames(subEAE_F), as.data.frame(logcounts(subEAE_F)))

	write.table(logCsCTRL_F, paste0(saveIn, "logcounts_ALL_CLUSTERS_CTRL_F_ONLY_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
	write.table(logCsEAE_F, paste0(saveIn, "logcounts_ALL_CLUSTERS_EAE_F_ONLY_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
#	write.xlsx(logCsCTRL_F, paste0(saveIn, "logcounts_ALL_CLUSTERS_CTRL_F_ONLY_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)
#	write.xlsx(logCsEAE_F, paste0(saveIn, "logcounts_ALL_CLUSTERS_EAE_F_ONLY_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)
	
	###############################################################################	###############################################################################	###############################################################################	###############################################################################	###############################################################################	###############################################################################	
	
	# COMPARING CTRL VERSUS EAE

	cTSce <- list()
	conD <- unique(clusteredSce$condition)
	# then get CTRL cellType
	cTSce[[conD[1]]] <- cTSubSet[, which(cTSubSet$condition==conD[1])]
	# and EAE
	cTSce[[conD[2]]] <- cTSubSet[, which(cTSubSet$condition==conD[2])]
	
	# set CTRL as a reference, to compute fold change
	cTSubSet$condition <- relevel(as.factor(cTSubSet$condition), ref="CTRL")

	if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
	
		print("finding markers")
	
		# FINDMARKERS for DEGS
		markersG <- findMarkers(cTSubSet, groups=cTSubSet$condition, pval.type="any", direction="any")

		cC <- "EAE"
		chosen <- markersG[[cC]]

		# RELAXING THE THRESHOLD, is not done for meninges
		selfMade <- chosen[chosen$FDR < 0.05, ]
		#selfMade <- chosen[chosen$FDR < 0.1, ]

		# and export data
		finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

		# sort by logfoldchange
		sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

		# and export
#		write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ALL_CLUSTERS_", cellLab[cT], "_CELLS_", expN, "_SLM_ZERO_UNO.csv"), col.names=T, row.names=F, sep="\t")
		write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ALL_CLUSTERS_", cellLab[cT], "_CELLS_", expN, "_SLM_ZERO_UNO.xlsx"), overwrite=TRUE)
	}
	
	# GABRIELA Chronic VERSUS Onset
	
	cTSubSetGC <- cTSubSet
	
	cTSce <- list()
	conD <- unique(clusteredSce$stage)
	# then get CTRL cellType
	cTSce[[conD[1]]] <- cTSubSetGC[, which(cTSubSetGC$stage==conD[1])]
	# and EAE
	cTSce[[conD[2]]] <- cTSubSetGC[, which(cTSubSetGC$stage==conD[2])]
	
	# set onset as a reference, to compute fold change
	cTSubSetGC$stage <- relevel(as.factor(cTSubSetGC$stage), ref="Onset")
	
	if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
	
		print("finding markers")
	
		# FINDMARKERS for DEGS
		markersG <- findMarkers(cTSubSetGC, groups=cTSubSetGC$stage, pval.type="any", direction="any")

		cC <- "Chronic"
		chosen <- markersG[[cC]]

		# RELAXING THE THRESHOLD, is not done for meninges
		selfMade <- chosen[chosen$FDR < 0.05, ]
		#selfMade <- chosen[chosen$FDR < 0.1, ]

		# and export data
		finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

		# sort by logfoldchange
		sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

		# and export
#		write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ALL_CLUSTERS_chronic_VS_onset_", cellLab[cT], "_CELLS_", expN, "_SLM_ZERO_UNO.csv"), col.names=T, row.names=F, sep="\t")
		write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ALL_CLUSTERS_chronic_VS_onset_", cellLab[cT], "_CELLS_", expN, "_SLM_ZERO_UNO.xlsx"), overwrite=TRUE)
	}	
	
	# GABRIELA Chronic VERSUS Onset, EAE only
	
	# get EAE only
	cTSubSetGC_EAE_only <- cTSubSet[, which(cTSubSet$condition=="EAE")]
	
	# set CTRL as a reference, to compute fold change
	cTSubSetGC_EAE_only$stage <- relevel(as.factor(cTSubSetGC_EAE_only$stage), ref="Onset")
	
	cTSce <- list()
	conD <- unique(cTSubSetGC_EAE_only$stage)

	if (ncol(cTSubSetGC_EAE_only[, which(cTSubSetGC_EAE_only$stage==conD[1])]) >= 10 & ncol(cTSubSetGC_EAE_only[, which(cTSubSetGC_EAE_only$stage==conD[2])]) >= 10) {
	
		print("finding markers")
	
		# FINDMARKERS for DEGS
		markersG <- findMarkers(cTSubSetGC_EAE_only, groups=cTSubSetGC_EAE_only$stage, pval.type="any", direction="any")

		cC <- "Chronic"
		chosen <- markersG[[cC]]

		# RELAXING THE THRESHOLD, is not done for meninges
		selfMade <- chosen[chosen$FDR < 0.05, ]
		#selfMade <- chosen[chosen$FDR < 0.1, ]

		# and export data
		finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

		# sort by logfoldchange
		sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

		# and export
#		write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ALL_CLUSTERS_chronic_VS_onset_EAE_ONLY_", cellLab[cT], "_CELLS_", expN, "_SLM_ZERO_UNO.csv"), col.names=T, row.names=F, sep="\t")
		write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ALL_CLUSTERS_chronic_VS_onset_EAE_ONLY_", cellLab[cT], "_CELLS_", expN, "_SLM_ZERO_UNO.xlsx"), overwrite=TRUE)
	}
	
	# GABRIELA Chronic VERSUS Onset, CTRL only
	
	# get EAE only
	cTSubSetGC_CTRL_only <- cTSubSet[, which(cTSubSet$condition=="CTRL")]
	
	# set CTRL as a reference, to compute fold change
	cTSubSetGC_CTRL_only$stage <- relevel(as.factor(cTSubSetGC_CTRL_only$stage), ref="Onset")
	
	cTSce <- list()
	conD <- unique(cTSubSetGC_CTRL_only$stage)

	if (ncol(cTSubSetGC_CTRL_only[, which(cTSubSetGC_CTRL_only$stage==conD[1])]) >= 10 & ncol(cTSubSetGC_CTRL_only[, which(cTSubSetGC_CTRL_only$stage==conD[2])]) >= 10) {
	
		print("finding markers")
	
		# FINDMARKERS for DEGS
		markersG <- findMarkers(cTSubSetGC_CTRL_only, groups=cTSubSetGC_CTRL_only$stage, pval.type="any", direction="any")

		cC <- "Chronic"
		chosen <- markersG[[cC]]

		# RELAXING THE THRESHOLD, is not done for meninges
		selfMade <- chosen[chosen$FDR < 0.05, ]
		#selfMade <- chosen[chosen$FDR < 0.1, ]

		# and export data
		finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

		# sort by logfoldchange
		sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

		# and export
#		write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ALL_CLUSTERS_chronic_VS_onset_CTRL_ONLY_", cellLab[cT], "_CELLS_", expN, "_SLM_ZERO_UNO.csv"), col.names=T, row.names=F, sep="\t")
		write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ALL_CLUSTERS_chronic_VS_onset_CTRL_ONLY_", cellLab[cT], "_CELLS_", expN, "_SLM_ZERO_UNO.xlsx"), overwrite=TRUE)
	}
	
	# GABRIELA eae VERSUS ctrl, ONSET only
	
	# get EAE only
	cTSubSetGC_Onset_only <- cTSubSet[, which(cTSubSet$stage=="Onset")]
	
	# set CTRL as a reference, to compute fold change
	cTSubSetGC_Onset_only$condition <- relevel(as.factor(cTSubSetGC_Onset_only$condition), ref="CTRL")
	
	cTSce <- list()
	conD <- unique(cTSubSetGC_Onset_only$condition)

	if (ncol(cTSubSetGC_Onset_only[, which(cTSubSetGC_Onset_only$condition==conD[1])]) >= 10 & ncol(cTSubSetGC_Onset_only[, which(cTSubSetGC_Onset_only$condition==conD[2])]) >= 10) {
	
		print("finding markers")
	
		# FINDMARKERS for DEGS
		markersG <- findMarkers(cTSubSetGC_Onset_only, groups=cTSubSetGC_Onset_only$condition, pval.type="any", direction="any")

		cC <- "EAE"
		chosen <- markersG[[cC]]

		# RELAXING THE THRESHOLD, is not done for meninges
		selfMade <- chosen[chosen$FDR < 0.05, ]
		#selfMade <- chosen[chosen$FDR < 0.1, ]

		# and export data
		finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

		# sort by logfoldchange
		sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

		# and export
		write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ALL_CLUSTERS_eae_VS_ctrl_ONSET_ONLY_", cellLab[cT], "_CELLS_", expN, "_SLM_ZERO_UNO.xlsx"), overwrite=TRUE)
	}
	
	
	# GABRIELA eae VERSUS ctrl, CHRONIC only
	
	# get EAE only
	cTSubSetGC_Chronic_only <- cTSubSet[, which(cTSubSet$stage=="Chronic")]
	
	# set CTRL as a reference, to compute fold change
	cTSubSetGC_Chronic_only$condition <- relevel(as.factor(cTSubSetGC_Chronic_only$condition), ref="CTRL")
	
	cTSce <- list()
	conD <- unique(cTSubSetGC_Chronic_only$condition)

	if (ncol(cTSubSetGC_Chronic_only[, which(cTSubSetGC_Chronic_only$condition==conD[1])]) >= 10 & ncol(cTSubSetGC_Chronic_only[, which(cTSubSetGC_Chronic_only$condition==conD[2])]) >= 10) {
	
		print("finding markers")
	
		# FINDMARKERS for DEGS
		markersG <- findMarkers(cTSubSetGC_Chronic_only, groups=cTSubSetGC_Chronic_only$condition, pval.type="any", direction="any")

		cC <- "EAE"
		chosen <- markersG[[cC]]

		# RELAXING THE THRESHOLD, is not done for meninges
		selfMade <- chosen[chosen$FDR < 0.05, ]
		#selfMade <- chosen[chosen$FDR < 0.1, ]

		# and export data
		finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

		# sort by logfoldchange
		sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

		# and export
		write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ALL_CLUSTERS_eae_VS_ctro_CHRONIC_ONLY_", cellLab[cT], "_CELLS_", expN, "_SLM_ZERO_UNO.xlsx"), overwrite=TRUE)
	}
}
