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

# set cellTypes
cellType <- c("T.8|T.CD8", "Tgd", "Neutrophils", "\\(T.4|\\(T.CD4", "B cells", "DC|Macrophages|Microglia|Monocytes")
cellLab <-  c("CD8", "Tgd", "Neutrophils", "CD4", "B_cells", "Myeloid")

# loop through cell type
for (cT in seq(1, length(cellLab), by=1)){

	print("")
	print(paste0("analysing all the genes ", cellLab[cT], " for exp ", expN, " clustered with slm"))
	print("")

	# then do everything
	saveIn <- paste0(rawPath, "batched_EAE_PCA/")
	saveUpDown <- paste0(saveIn, "enrichment_", cellLab[cT], "/UpDown/")
	saveImg <- paste0(saveIn, "figures/analysis_", cellLab[cT], "/")

	# check if saveIn directory exist, otherwise create it
	ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
	ifelse(dir.exists(saveUpDown), TRUE, dir.create(saveUpDown, recursive=T))
	ifelse(dir.exists(saveImg), TRUE, dir.create(saveImg, recursive=T))

	# set palette with 12 colours for colouring plots
	colorsPal <- brewer.pal(12, "Paired")

	# NOTE: we are using QCed data file which logcounts were normalised using sizeFactors
	# but this data file contains also the information from the SLM clustering
	# performed using the pca_ed data, in order to have a good clustering but also all
	# the genes
	kS <- 15	
	
	# load data
	cTSubSet <- readRDS(paste0(saveIn, "QCed_data_with_clusters_", expN, "_", cellLab[cT], "_ONLY.Rds"))
	
	# CONTE GREZZE ELEONORA
	raw_counts <- assay(cTSubSet, "counts")
	colnames(raw_counts) <- paste(colnames(cTSubSet), cTSubSet$condition, cTSubSet$stage, cTSubSet$pruned_fine, "cluster", cTSubSet$slm_0.1, sep="_")
	with_row_names <- cbind.data.frame(rownames(raw_counts), raw_counts)
	colnames(with_row_names) <- c("genes", colnames(raw_counts))
	
	write.table(with_row_names, paste0(saveIn, "counts_CONTE_GREZZE_SUBCLUSTERED_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
	#######################

	################################################# EXPORT THE NEW TABLE ELEONORA ASKED FEBRUARY #################################################

	superfull_table <- assay(cTSubSet, "logcounts")
	colnames(superfull_table) <- paste(colnames(cTSubSet), cTSubSet$condition, cTSubSet$stage, cellLab[cT], "cluster", cTSubSet$slm_0.1, sep="_")
	with_row_names <- cbind.data.frame(rownames(superfull_table), superfull_table)
	colnames(with_row_names) <- c("genes", colnames(superfull_table))
	
	write.table(with_row_names, paste0(saveIn, "logcounts_all_the_selected_cells_version_barbara_elena_eleonora_only_SUBCLUSTERING_PRIMA_DI_LOGNORMCOUNTS_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
	
	# find the clusters
	fR <- table(cTSubSet$slm_0.1)
	clustS <- names(fR[which(fR!=0)])

	# get label ready
	cTSubSet$fine_to_main <- gsub(" ", "", unlist(lapply(strsplit(cTSubSet$pruned_fine, "\\("), `[[`, 1)))
	
	##################### MAIN labels	
	
	# plotting heatmap with cell type
	tab_fine_to_main <- table(Assigned=cTSubSet$fine_to_main, Cluster=cTSubSet$slm_0.1, Condition=cTSubSet$condition)

	eTabCNEAE <- NULL
	eTabCNCTRL <- NULL
	
	if (length(unique(cTSubSet$fine_to_main)) > 1 && length(unique(cTSubSet$slm_0.1)) > 1) {
		# get EAE full table ready to be exported
		eTabEAE <- as.data.frame(tab_fine_to_main[, , "EAE"])
		eTabRNEAE <- cbind.data.frame(cellType=rownames(eTabEAE), eTabEAE)
		eTabCNEAE <- eTabRNEAE[, c(2,3,4)]
		colnames(eTabCNEAE) <- c("CellType", "Cluster", "Freq")

		# and CTRL
		eTabCTRL <- as.data.frame(tab_fine_to_main[, , "CTRL"])
		eTabRNCTRL <- cbind.data.frame(cellType=rownames(eTabCTRL), eTabCTRL)
		eTabCNCTRL <- eTabRNCTRL[, c(2,3,4)]
		colnames(eTabCNCTRL) <- c("CellType", "Cluster", "Freq")
	} else if (length(unique(cTSubSet$fine_to_main)) == 1) {
		# get EAE full table ready to be exported
		eTabEAE <- as.data.frame(tab_fine_to_main[, , "EAE"])
		eTabCNEAE <- cbind.data.frame(rep(unique(cTSubSet$fine_to_main), length(eTabEAE)), rownames(eTabEAE), eTabEAE)
		colnames(eTabCNEAE) <- c("CellType", "Cluster", "Freq")

		# and CTRL
		eTabCTRL <- as.data.frame(tab_fine_to_main[, , "CTRL"])
		eTabCNCTRL <- cbind.data.frame(rep(unique(cTSubSet$fine_to_main), length(eTabCTRL)), rownames(eTabCTRL), eTabCTRL)
		colnames(eTabCNCTRL) <- c("CellType", "Cluster", "Freq")
	} else if (length(unique(cTSubSet$slm_0.1)) == 1) {
		# get EAE full table ready to be exported
		eTabEAE <- as.data.frame(tab_fine_to_main[, , "EAE"])
		eTabCNEAE <- cbind.data.frame(rownames(eTabEAE), rep(unique(cTSubSet$slm_0.1), length(eTabEAE)), eTabEAE)
		colnames(eTabCNEAE) <- c("CellType", "Cluster", "Freq")

		# and CTRL
		eTabCTRL <- as.data.frame(tab_fine_to_main[, , "CTRL"])
		eTabCNCTRL <- cbind.data.frame(rownames(eTabCTRL), rep(unique(cTSubSet$slm_0.1), length(eTabCTRL)), eTabCTRL)
		colnames(eTabCNCTRL) <- c("CellType", "Cluster", "Freq")
	}	
	
	# export both tables
	write.xlsx(eTabCNEAE, paste0(saveImg, "EAE_cluster_SUBCLUSTERED_content_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), colNames=F, overwrite=T)

	write.xlsx(eTabCNCTRL, paste0(saveImg, "CTRL_cluster_SUBCLUSTERED_content_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), colNames=F, overwrite=T)

	##################### FULL labels

	tab_full_labels <- table(Assigned=cTSubSet$pruned_fine, Cluster=cTSubSet$slm_0.1, Condition=cTSubSet$condition)

	# load objs
	eTabcnFEAE <- NULL
	eTabcnFCTRL <- NULL
	if (length(unique(cTSubSet$pruned_fine)) > 1 && length(unique(cTSubSet$slm_0.1)) > 1) {
		# get EAE full table ready to be exported
		etabFEAE <- as.data.frame(tab_full_labels[, , "EAE"])
		eTabrnFEAE <- cbind.data.frame(cellType=rownames(etabFEAE), etabFEAE)
		eTabcnFEAE <- eTabrnFEAE[, c(2,3,4)]
		colnames(eTabcnFEAE) <- c("CellType", "Cluster", "Freq")

		# and CTRL
		etabFCTRL <- as.data.frame(tab_full_labels[, , "CTRL"])
		eTabrnFCTRL <- cbind.data.frame(cellType=rownames(etabFCTRL), etabFCTRL)
		eTabcnFCTRL <- eTabrnFCTRL[, c(2,3,4)]
		colnames(eTabcnFCTRL) <- c("CellType", "Cluster", "Freq")
	} else if (length(unique(cTSubSet$pruned_fine)) == 1) {
		# get EAE full table ready to be exported
		etabFEAE <- as.data.frame(tab_full_labels[, , "EAE"])
		eTabcnFEAE <- cbind.data.frame(rep(unique(cTSubSet$pruned_fine), length(etabFEAE)), rownames(etabFEAE), etabFEAE)
		colnames(eTabcnFEAE) <- c("CellType", "Cluster", "Freq")

		# and CTRL
		etabFCTRL <- as.data.frame(tab_full_labels[, , "CTRL"])
		eTabcnFCTRL <- cbind.data.frame(rep(unique(cTSubSet$pruned_fine), length(etabFCTRL)), rownames(etabFCTRL), etabFCTRL)
		colnames(eTabcnFCTRL) <- c("CellType", "Cluster", "Freq")
	} else if (length(unique(cTSubSet$slm_0.1)) == 1) {
		# get EAE full table ready to be exported
		etabFEAE <- as.data.frame(tab_full_labels[, , "EAE"])
		eTabcnFEAE <- cbind.data.frame(rownames(etabFEAE), rep(unique(cTSubSet$slm_0.1), length(etabFEAE)), etabFEAE)
		colnames(eTabcnFEAE) <- c("CellType", "Cluster", "Freq")

		# and CTRL
		etabFCTRL <- as.data.frame(tab_full_labels[, , "CTRL"])
		eTabcnFCTRL <- cbind.data.frame(rownames(etabFCTRL), rep(unique(cTSubSet$slm_0.1), length(etabFCTRL)), etabFCTRL)
		colnames(eTabcnFCTRL) <- c("CellType", "Cluster", "Freq")
	}

	# export the full table
	write.xlsx(eTabcnFEAE, paste0(saveImg, "EAE_cluster_content_SUBCLUSTERED_full_labels_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), colNames=T, overwrite=T)

	write.xlsx(eTabcnFCTRL, paste0(saveImg, "CTRL_cluster_content_SUBCLUSTERED_full_labels_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), colNames=T, overwrite=T)

	##################### FULL REGION ###############################################

	tab_REGION_labels <- table(cTSubSet$pruned_fine, cTSubSet$condition, cTSubSet$stage)

	# load objs
	eTabcnFCHRONIC <- NULL
	eTabcnFONSET <- NULL

	if (length(unique(cTSubSet$pruned_fine)) > 1 && length(unique(cTSubSet$slm_0.1)) > 1) {
		# get CHRONIC
		etabFCHRONIC <- as.data.frame(tab_REGION_labels[, , "Chronic"])
		eTabrnFCHRONIC <- cbind.data.frame(cellType=rownames(etabFCHRONIC), etabFCHRONIC)
		eTabcnFCHRONIC <- eTabrnFCHRONIC[, c(2,3,4)]
		colnames(eTabcnFCHRONIC) <- c("cellType", "condition", "freq")

		# and ONSET
		etabFONSET <- as.data.frame(tab_REGION_labels[, , "Onset"])
		eTabrnFONSET <- cbind.data.frame(cellType=rownames(etabFONSET), etabFONSET)
		eTabcnFONSET <- eTabrnFONSET[, c(2,3,4)]
		colnames(eTabcnFONSET) <- c("cellType", "condition", "freq")
	} else if (length(unique(cTSubSet$pruned_fine)) == 1) {
		# get CHRONIC full table ready to be exported
		etabFCHRONIC <- as.data.frame(tab_REGION_labels[, , "Chronic"])
		eTabcnFCHRONIC <- cbind.data.frame(rep(unique(cTSubSet$pruned_fine), length(etabFCHRONIC)), rownames(etabFCHRONIC), etabFCHRONIC)
		colnames(eTabcnFCHRONIC) <- c("CellType", "Cluster", "Freq")

		# and ONSET
		etabFONSET <- as.data.frame(tab_REGION_labels[, , "Onset"])
		eTabcnFONSET <- cbind.data.frame(rep(unique(cTSubSet$pruned_fine), length(etabFONSET)), rownames(etabFONSET), etabFONSET)
		colnames(eTabcnFONSET) <- c("CellType", "Cluster", "Freq")
	} else if (length(unique(cTSubSet$slm_0.1)) == 1) {
		# get CHRONIC full table ready to be exported
		etabFCHRONIC <- as.data.frame(tab_REGION_labels[, , "Chronic"])
		eTabcnFCHRONIC <- cbind.data.frame(rownames(etabFCHRONIC), rep(unique(cTSubSet$slm_0.1), nrow(etabFCHRONIC)), etabFCHRONIC)
		colnames(eTabcnFCHRONIC) <- c("CellType", "Cluster", "Freq")

		# and ONSET
		etabFONSET <- as.data.frame(tab_REGION_labels[, , "Onset"])
		eTabcnFONSET <- cbind.data.frame(rownames(etabFONSET), rep(unique(cTSubSet$slm_0.1), nrow(etabFONSET)), etabFONSET)
		colnames(eTabcnFONSET) <- c("CellType", "Cluster", "Freq")
	}
	
	# export the full table
	write.xlsx(eTabcnFCHRONIC, paste0(saveImg, "CHRONIC_cluster_content_full_labels_SUBCLUSTERED_BY_REGION_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=F)

	write.xlsx(eTabcnFONSET, paste0(saveImg, "ONSET_cluster_content_full_labels_SUBCLUSTERED_BY_REGION_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=F)

	####################### CLEAN AND NORMALISE #######################
	
	# remove genes with too many zeros
	cTSubSet <- cTSubSet[rowMeans(assays(cTSubSet)$counts > 0) > 0.05, ]

	# set minmean equal to 0.5 since this is the value used in the previous analysis, i.e. the one in the manuscript
	# compute sizeFactors
	cTSubSet <- computeSumFactors(cTSubSet, min.mean=0.5)
	cTSubSet <- logNormCounts(cTSubSet)
	
	################################################# EXPORT THE NEW TABLE ELENA ASKED #################################################

	superfull_table <- assay(cTSubSet, "logcounts")
	colnames(superfull_table) <- paste(colnames(cTSubSet), cTSubSet$condition, cTSubSet$stage, cTSubSet$pruned_fine, "cluster", cTSubSet$slm_0.1, sep="_")
	with_row_names <- cbind.data.frame(rownames(superfull_table), superfull_table)
	colnames(with_row_names) <- c("genes", colnames(superfull_table))
	
	write.table(with_row_names, paste0(saveIn, "logcounts_all_the_selected_cells_version_barbara_elena_eleonora_only_SUBCLUSTERED_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
	
#	write.xlsx(with_row_names, paste0(saveIn, "logcounts_ALL_THE_SELECTED_CELLS_VERSION_BARBARA_ELENA_ELEONORA_ONLY_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T)
	
	################################################# EXPORT THE NEW TABLE ELENA ASKED #################################################
	
	# exporting logcounts for Elena's violins BY CONDITION
	subCTRL <- cTSubSet[, which(cTSubSet$condition=="CTRL")]
	subEAE <- cTSubSet[, which(cTSubSet$condition=="EAE")]
	
	# append info to colnames
	coln_CTRL <- paste(colnames(subCTRL), subCTRL$condition, subCTRL$gender, subCTRL$pruned_fine, "cluster", subCTRL$slm_0.1, sep="_")
	with_row_names_CTRL <- cbind.data.frame(rownames(subCTRL), logcounts(subCTRL))
	colnames(with_row_names_CTRL) <- c("genes", coln_CTRL)
	
	coln_EAE <- paste(colnames(subEAE), subEAE$condition, subEAE$gender, subEAE$pruned_fine, "cluster", subEAE$slm_0.1, sep="_")
	with_row_names_EAE <- cbind.data.frame(rownames(subEAE), logcounts(subEAE))
	colnames(with_row_names_EAE) <- c("genes", coln_EAE)

	write.table(with_row_names_CTRL, paste0(saveIn, "logcounts_SUBCLUSTERING_CTRL_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
	write.table(with_row_names_EAE, paste0(saveIn, "logcounts_SUBCLUSTERING_EAE_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")

#	write.xlsx(with_row_names_CTRL, paste0(saveIn, "logcounts_SUBCLUSTERING_CTRL_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)
#	write.xlsx(with_row_names_EAE, paste0(saveIn, "logcounts_SUBCLUSTERING_EAE_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)
	
	# exporting logcounts for Elena's violins BY REGION
	subChronic <- cTSubSet[, which(cTSubSet$stage=="Chronic")]
	subOnset <- cTSubSet[, which(cTSubSet$stage=="Onset")]
	
	# append info to colnames
	coln_Chronic <- paste(colnames(subChronic), subChronic$condition, subChronic$gender, subChronic$pruned_fine, "cluster", subChronic$slm_0.1, sep="_")
	with_row_names_Chronic <- cbind.data.frame(rownames(subChronic), logcounts(subChronic))
	colnames(with_row_names_Chronic) <- c("genes", coln_Chronic)
	
	coln_Onset <- paste(colnames(subOnset), subOnset$condition, subOnset$gender, subOnset$pruned_fine, "cluster", subOnset$slm_0.1, sep="_")
	with_row_names_Onset <- cbind.data.frame(rownames(subOnset), logcounts(subOnset))
	colnames(with_row_names_Onset) <- c("genes", coln_Onset)
	
	write.table(with_row_names_Chronic, paste0(saveIn, "logcounts_SUBCLUSTERING_CHRONIC_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
	write.table(with_row_names_Onset, paste0(saveIn, "logcounts_SUBCLUSTERING_ONSET_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")

#	write.xlsx(with_row_names_Chronic, paste0(saveIn, "logcounts_SUBCLUSTERING_CHRONIC_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)
#	write.xlsx(with_row_names_Onset, paste0(saveIn, "logcounts_SUBCLUSTERING_ONSET_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)
	
	# exporting logcounts for Elena's violins BY BOTH
	subCTRL_Chronic <- cTSubSet[, which(cTSubSet$conditionBYstage=="CTRL_Chronic")]
	subCTRL_Onset <- cTSubSet[, which(cTSubSet$conditionBYstage=="CTRL_Onset")]
	subEAE_Chronic <- cTSubSet[, which(cTSubSet$conditionBYstage=="EAE_Chronic")]
	subEAE_Onset <- cTSubSet[, which(cTSubSet$conditionBYstage=="EAE_Onset")]
	
	# append info to colnames
	if (dim(subCTRL_Chronic)[2] > 0) {
		coln_CTRL_Chronic <- paste(colnames(subCTRL_Chronic), subCTRL_Chronic$condition, subCTRL_Chronic$gender, subCTRL_Chronic$pruned_fine, "cluster", subCTRL_Chronic$slm_0.1, sep="_")
		with_row_names_CTRL_Chronic <- cbind.data.frame(rownames(subCTRL_Chronic), logcounts(subCTRL_Chronic))
		colnames(with_row_names_CTRL_Chronic) <- c("genes", coln_CTRL_Chronic)
		
		write.table(with_row_names_CTRL_Chronic, paste0(saveIn, "logcounts_SUBCLUSTERING_CTRL_CHRONIC_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
#		write.xlsx(with_row_names_CTRL_Chronic, paste0(saveIn, "logcounts_SUBCLUSTERING_CTRL_CHRONIC_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)	
	}
	
	if (dim(subCTRL_Onset)[2] > 0) {
		coln_CTRL_Onset <- paste(colnames(subCTRL_Onset), subCTRL_Onset$condition, subCTRL_Onset$gender, subCTRL_Onset$pruned_fine, "cluster", subCTRL_Onset$slm_0.1, sep="_")
		with_row_names_CTRL_Onset <- cbind.data.frame(rownames(subCTRL_Onset), logcounts(subCTRL_Onset))
		colnames(with_row_names_CTRL_Onset) <- c("genes", coln_CTRL_Onset)
		
		write.table(with_row_names_CTRL_Onset, paste0(saveIn, "logcounts_SUBCLUSTERING_CTRL_ONSET_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
#		write.xlsx(with_row_names_CTRL_Onset, paste0(saveIn, "logcounts_SUBCLUSTERING_CTRL_ONSET_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)
	}
			
	if (dim(subEAE_Chronic)[2] > 0) {	
		coln_EAE_Chronic <- paste(colnames(subEAE_Chronic), subEAE_Chronic$condition, subEAE_Chronic$gender, subEAE_Chronic$pruned_fine, "cluster", subEAE_Chronic$slm_0.1, sep="_")
		with_row_names_EAE_Chronic <- cbind.data.frame(rownames(subEAE_Chronic), logcounts(subEAE_Chronic))
		colnames(with_row_names_EAE_Chronic) <- c("genes", coln_EAE_Chronic)
		
		write.table(with_row_names_EAE_Chronic, paste0(saveIn, "logcounts_SUBCLUSTERING_EAE_CHRONIC_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
#		write.xlsx(with_row_names_EAE_Chronic, paste0(saveIn, "logcounts_SUBCLUSTERING_EAE_CHRONIC_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)
	}
		
	if (dim(subEAE_Onset)[2] > 0) {	
		coln_EAE_Onset <- paste(colnames(subEAE_Onset), subEAE_Onset$condition, subEAE_Onset$gender, subEAE_Onset$pruned_fine, "cluster", subEAE_Onset$slm_0.1, sep="_")
		with_row_names_EAE_Onset <- cbind.data.frame(rownames(subEAE_Onset), logcounts(subEAE_Onset))
		colnames(with_row_names_EAE_Onset) <- c("genes", coln_EAE_Onset)
		
		write.table(with_row_names_EAE_Onset, paste0(saveIn, "logcounts_SUBCLUSTERING_EAE_ONSET_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")	
#		write.xlsx(with_row_names_EAE_Onset, paste0(saveIn, "logcounts_SUBCLUSTERING_EAE_ONSET_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)
	}
		
	# Cluster composition plots using main labels
	conD <- unique(cTSubSet$conditionBYstage)

	dfReg <- data.frame(table(cTSubSet$slm_0.1, cTSubSet$pruned_fine, cTSubSet$conditionBYstage))

	colnames(dfReg) <- c("Cluster", "CellType", "Condition", "Frequency")
	lims <- levels(factor(cTSubSet$pruned_fine))
	lims <- lims[order(table(cTSubSet$pruned_fine))]

	# plot!
	p <- plot_membership(dfReg, cond=conD[2], col=colorsPal[2])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[2], "_cluster_composition_SUBCLUSTERED_", cellLab[cT], "_exp_", expN, "_SLM_ZERO_UNO.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()

	p <- plot_membership(dfReg, cond=conD[1], col=colorsPal[6])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[1], "_cluster_composition_SUBCLUSTERED_", cellLab[cT], "_exp_", expN, "_SLM_ZERO_UNO.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()
	
	# plot!
	p <- plot_membership(dfReg, cond=conD[3], col=colorsPal[3])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[3], "_cluster_composition_SUBCLUSTERED_", cellLab[cT], "_exp_", expN, "_SLM_ZERO_UNO.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()

	p <- plot_membership(dfReg, cond=conD[4], col=colorsPal[5])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[4], "_cluster_composition_SUBCLUSTERED_", cellLab[cT], "_exp_", expN, "_SLM_ZERO_UNO.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()
	
	####################### TEST 1 IN EACH CLUSTER COMPARE STUFF #######################
	
	# loop through all the clusters containing a specific cellType
	for (cS in seq(1, length(clustS), by=1)) {
	
		# PERFORM CHRONIC+ONSET, CTRL VS EAE
	
		print(paste("cluster position ", clustS[cS]))
	
		# subset SCE to get only cS cluster cells
		aClust <- cTSubSet[, cTSubSet$slm_0.1==clustS[cS]]
		
		print(paste0("######## COMPARING ", cellType[cT], " INSIDE CLUSTER ", clustS[cS], " ########"))

		if(unique(length(aClust$condition))==2) {

			# set reference condition
			aClust$condition <- relevel(as.factor(aClust$condition), ref="CTRL")

			cTSce <- list()
			conD <- unique(aClust$condition)
			# then get CTRL cellType
			cTSce[[conD[1]]] <- aClust[, which(aClust$condition==conD[1])]
			# and EAE
			cTSce[[conD[2]]] <- aClust[, which(aClust$condition==conD[2])]
			
			if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
				# export logcounts
	#			write.xlsx(cbind.data.frame(rownames(aClust), logcounts(aClust)), paste0(saveUpDown, "logcounts", cellType, "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T)
				
				cat("finding markers\n")
				
				# FINDMARKERS for DEGS
				markersG <- findMarkers(aClust, groups=aClust$condition, pval.type="any", direction="any")
		
				cC <- "EAE"
				chosen <- markersG[[cC]]

				# RELAXING THE THRESHOLD, is not done for Onset
				selfMade <- chosen[chosen$FDR < 0.05, ]
				#selfMade <- chosen[chosen$FDR < 0.1, ]

				# and export data
				finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

				# sort by logfoldchange
				if (length(finalGenes$logFC) != 0) {
					sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

					# and export
	#				write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_SUBCLUSTERED_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
					write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_SUBCLUSTERED_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.xlsx"), overwrite=T)
				}
			}
		}
		
		# PERFORM CHRONIC ONLY, CTRL VS EAE
		
		aClustByReg <- cTSubSet[, cTSubSet$slm_0.1==clustS[cS] & cTSubSet$stage=="Chronic"]
		
		print(paste0("######## COMPARING ", cellType[cT], " INSIDE CLUSTER CHRONIC ONLY ", clustS[cS], " ########"))

		cTSce <- list()
		conD <- unique(aClustByReg$condition)
		
		if (length(conD)>1) {

			# set reference condition
			aClustByReg$condition <- relevel(as.factor(aClustByReg$condition), ref="CTRL")
			
			# then get CTRL cellType
			
			cTSce[[conD[1]]] <- aClustByReg[, which(aClustByReg$condition==conD[1])]
			# and EAE
			cTSce[[conD[2]]] <- aClustByReg[, which(aClustByReg$condition==conD[2])]
			
			if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
				# export logcounts
#				write.xlsx(cbind.data.frame(rownames(aClustByReg), logcounts(aClustByReg)), paste0(saveUpDown, "logcounts", cellType, "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T)
				
				cat("finding markers\n")
				
				# FINDMARKERS for DEGS
				markersG <- findMarkers(aClustByReg, groups=aClustByReg$condition, pval.type="any", direction="any")

				cC <- "EAE"
				chosen <- markersG[[cC]]

				# RELAXING THE THRESHOLD, is not done for Onset
				selfMade <- chosen[chosen$FDR < 0.05, ]
				#selfMade <- chosen[chosen$FDR < 0.1, ]

				# and export data
				finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

				# sort by logfoldchange
				if (length(finalGenes$logFC) != 0) {
					# sort by logfoldchange
					sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

					# and export
	#				write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_SUBCLUSTERED_CHRONIC_ONLY_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
					write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_SUBCLUSTERED_CHRONIC_ONLY_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.xlsx"), overwrite=T)
				}
			}
		} else {
			print("not enough conditions to compare")
		}
		
		# PERFORM ONSET ONLY, CTRL VS EAE
	
		aClustByReg <- cTSubSet[, cTSubSet$slm_0.1==clustS[cS] & cTSubSet$stage=="Onset"]

		print(paste0("######## COMPARING ", cellType[cT], " INSIDE CLUSTER ONSET ONLY ", clustS[cS], " ########"))

		cTSce <- list()
		conD <- unique(aClustByReg$condition)

		if (length(conD)>1) {

			# set reference condition
			aClustByReg$condition <- relevel(as.factor(aClustByReg$condition), ref="CTRL")
			
			# then get CTRL cellType
			
			cTSce[[conD[1]]] <- aClustByReg[, which(aClustByReg$condition==conD[1])]
			# and EAE
			cTSce[[conD[2]]] <- aClustByReg[, which(aClustByReg$condition==conD[2])]
			
			if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
				
				cat("finding markers\n")
				
				# FINDMARKERS for DEGS
				markersG <- findMarkers(aClustByReg, groups=aClustByReg$condition, pval.type="any", direction="any")

				cC <- "EAE"
				chosen <- markersG[[cC]]

				# RELAXING THE THRESHOLD, is not done for Onset
				selfMade <- chosen[chosen$FDR < 0.05, ]
				#selfMade <- chosen[chosen$FDR < 0.1, ]

				# and export data
				finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

				# sort by logfoldchange

				# sort by logfoldchange
				if (length(finalGenes$logFC) != 0) {
					sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

					# and export
	#				write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_SUBCLUSTERED_ONSET_ONLY_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
					write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_SUBCLUSTERED_ONSET_ONLY_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.xlsx"), overwrite=T)
				}
			}
		} else {
			print("not enough conditions to compare")
		}
		
		################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################ AS ALESSANDRO ASKED, WE ARE ADDING FURTHER ANALYSIS
		
		# PERFORM CTRL+EAE, CHRONIC vs ONSET
	
		print(paste("cluster position ", clustS[cS]))
	
		# subset SCE to get only cS cluster cells
		aClust <- cTSubSet[, cTSubSet$slm_0.1==clustS[cS]]
		
		print(paste0("######## COMPARING ", cellType[cT], " INSIDE CLUSTER ", clustS[cS], " ########"))

		if(unique(length(aClust$stage))==2) {

			# set reference stage
			aClust$stage <- relevel(as.factor(aClust$stage), ref="Onset")

			cTSce <- list()
			conD <- unique(aClust$stage)
			# then get CTRL cellType
			cTSce[[conD[1]]] <- aClust[, which(aClust$stage==conD[1])]
			# and EAE
			cTSce[[conD[2]]] <- aClust[, which(aClust$stage==conD[2])]
			
			if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
				# export logcounts
	#			write.xlsx(cbind.data.frame(rownames(aClust), logcounts(aClust)), paste0(saveUpDown, "logcounts", cellType, "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T)
				
				cat("finding markers\n")
				
				# FINDMARKERS for DEGS
				markersG <- findMarkers(aClust, groups=aClust$stage, pval.type="any", direction="any")
		
				cC <- "Chronic"
				chosen <- markersG[[cC]]

				# RELAXING THE THRESHOLD, is not done for Onset
				selfMade <- chosen[chosen$FDR < 0.05, ]
				#selfMade <- chosen[chosen$FDR < 0.1, ]

				# and export data
				finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

				# sort by logfoldchange
				if (length(finalGenes$logFC) != 0) {

					# sort by logfoldchange
					sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

					# and export
	#				write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_SUBCLUSTERED_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
					write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_SUBCLUSTERED_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.xlsx"), overwrite=T)
				}
			}
		}
		
		# PERFORM CTRL ONLY, Chronic VS Onset
		
		aClustByReg <- cTSubSet[, cTSubSet$slm_0.1==clustS[cS] & cTSubSet$condition=="CTRL"]
		
		print(paste0("######## COMPARING ", cellType[cT], " INSIDE CLUSTER CONTROL ONLY ", clustS[cS], " ########"))

		cTSce <- list()
		conD <- unique(aClustByReg$stage)
		
		if (length(conD)>1) {

			# set reference stage
			aClustByReg$stage <- relevel(as.factor(aClustByReg$stage), ref="Onset")
			
			# then get CTRL cellType
			
			cTSce[[conD[1]]] <- aClustByReg[, which(aClustByReg$stage==conD[1])]
			# and EAE
			cTSce[[conD[2]]] <- aClustByReg[, which(aClustByReg$stage==conD[2])]
			
			if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
				# export logcounts
#				write.xlsx(cbind.data.frame(rownames(aClustByReg), logcounts(aClustByReg)), paste0(saveUpDown, "logcounts", cellType, "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T)
				
				cat("finding markers\n")
				
				# FINDMARKERS for DEGS
				markersG <- findMarkers(aClustByReg, groups=aClustByReg$stage, pval.type="any", direction="any")

				cC <- "Chronic"
				chosen <- markersG[[cC]]

				# RELAXING THE THRESHOLD, is not done for Onset
				selfMade <- chosen[chosen$FDR < 0.05, ]
				#selfMade <- chosen[chosen$FDR < 0.1, ]

				# and export data
				finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

				# sort by logfoldchange
				if (length(finalGenes$logFC) != 0) {

					# sort by logfoldchange
					sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

					# and export
	#				write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_SUBCLUSTERED_CTRL_ONLY_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
					write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_SUBCLUSTERED_CTRL_ONLY_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.xlsx"), overwrite=T)
				}
			}
		} else {
			print("not enough stages to compare")
		}
		
		# PERFORM EAE ONLY, Chronic vs Onset
	
		aClustByReg <- cTSubSet[, cTSubSet$slm_0.1==clustS[cS] & cTSubSet$condition=="EAE"]

		print(paste0("######## COMPARING ", cellType[cT], " INSIDE CLUSTER ONSET ONLY ", clustS[cS], " ########"))

		cTSce <- list()
		conD <- unique(aClustByReg$stage)

		if (length(conD)>1) {

			# set reference stage
			aClustByReg$stage <- relevel(as.factor(aClustByReg$stage), ref="Onset")
			
			# then get CTRL cellType
			
			cTSce[[conD[1]]] <- aClustByReg[, which(aClustByReg$stage==conD[1])]
			# and EAE
			cTSce[[conD[2]]] <- aClustByReg[, which(aClustByReg$stage==conD[2])]
			
			if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
				
				cat("finding markers\n")
				
				# FINDMARKERS for DEGS
				markersG <- findMarkers(aClustByReg, groups=aClustByReg$stage, pval.type="any", direction="any")

				cC <- "Chronic"
				chosen <- markersG[[cC]]

				# RELAXING THE THRESHOLD, is not done for Onset
				selfMade <- chosen[chosen$FDR < 0.05, ]
				#selfMade <- chosen[chosen$FDR < 0.1, ]

				# and export data
				finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

				# sort by logfoldchange
				if (length(finalGenes$logFC) != 0) {

					# sort by logfoldchange
					sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

					# and export
	#				write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_SUBCLUSTERED_EAE_ONLY_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
					write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_SUBCLUSTERED_EAE_ONLY_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.xlsx"), overwrite=T)
				}
			}
		} else {
			print("not enough stages to compare")
		}		
	}
		
	####################### TEST 1.5 COMPARE A CLUSTER AGAINST ALL OTHERS, CHRONIC ONLY #######################

	print("####################### TEST 1.5 CHRONIC #######################")
	
	# get chronic only
	a_stage_chronic <- cTSubSet[, which(cTSubSet$stage=="Chronic")]
	
	################################################################ BEGIN ELEONORA
	# get EAE full table ready to be exported
	etabFEAEChronic <- as.data.frame(table(a_stage_chronic$slm_0.1))
	colnames(etabFEAEChronic) <- c("cluster", "frequency")

	# export the full table
	write.xlsx(etabFEAEChronic, paste0(saveImg, "ELEONORA_CHRONIC_cluster_content_full_labels_SUBCLUSTERED_BY_REGION_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)
	################################################################ END ELEONORA
	
	# loop through all the clusters containing a specific cellType
	for (cS in seq(1, length(clustS), by=1)) {
		
		# we need more than one cluster to perform this comparison
		if (length(unique(clustS))!=1) {
	
			print(paste("cluster position ", clustS[cS]))
		
			# subset SCE to get only cS cluster cells
			# this is for comparing one cluster versus all the other clusters
			a_stage_chronic$whatCluster[a_stage_chronic$slm_0.1 %in% clustS[cS]] <- "thisOne"
			a_stage_chronic$whatCluster[a_stage_chronic$slm_0.1 %in% clustS[-cS]] <- "allOther"

			################################################################ BEGIN ELEONORA
			# get cells in a specific cluster
			eleonoraSubSet <- a_stage_chronic[, a_stage_chronic$slm_0.1 %in% clustS[cS]]
			# get EAE full table ready to be exported
			etabFEAE <- as.data.frame(table(eleonoraSubSet$pruned_fine, eleonoraSubSet$condition))
			colnames(etabFEAE) <- c("cellType", "condition", "frequency")

			# export the full table
			write.xlsx(etabFEAE, paste0(saveImg, "CHRONIC_cluster_content_full_labels_SUBCLUSTERED_BY_REGION_exp_", expN, "_", cellLab[cT], "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)
			################################################################ END ELEONORA
			
			# otherwise we may compare, inside a cluster, CTRL versus AD
			
			print(paste0("######## COMPARING ", cellType[cT], " INSIDE CHRONIC CLUSTER ", clustS[cS], " ########"))
				
			# set ALLOTHER as a reference, to compute fold change
			a_stage_chronic$whatCluster <- relevel(as.factor(a_stage_chronic$whatCluster), ref="allOther")

			cTSce <- list()
			conD <- unique(a_stage_chronic$whatCluster)
			
			if (length(conD)>1) {
			
				# then get CTRL cellType
				cTSce[[conD[1]]] <- a_stage_chronic[, which(a_stage_chronic$whatCluster==conD[1])]
				# and EAE
				cTSce[[conD[2]]] <- a_stage_chronic[, which(a_stage_chronic$whatCluster==conD[2])]
				
				if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
					# export logcounts
	#				write.xlsx(cbind.data.frame(rownames(a_stage_chronic), logcounts(a_stage_chronic)), paste0(saveUpDown, "logcounts", cellType, "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T)
					
					cat("finding markers\n")
					
					# FINDMARKERS for DEGS
					markersG <- findMarkers(a_stage_chronic, groups=a_stage_chronic$whatCluster, pval.type="any", direction="any")
			
					cC <- "thisOne"
					chosen <- markersG[[cC]]

					# RELAXING THE THRESHOLD, is not done for Onset
					selfMade <- chosen[chosen$FDR < 0.05, ]
					#selfMade <- chosen[chosen$FDR < 0.1, ]

					# and export data
					finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

					# sort by logfoldchange
					if (length(finalGenes$logFC) != 0) {

						# sort by logfoldchange
						sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

						# and export
	#					write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ONE_VERSUS_ALL_OTHER_CHRONIC_ONLY_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
						write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ONE_VERSUS_ALL_OTHER_CHRONIC_SUBCLUSTERED_ONLY_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.xlsx"), overwrite=T)
					}
				}
			} else {
				print("not enough conditions to compare")
			}
		} else {
			print("can't compare with only a cluster")
		}
	}
		
	# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY# PERFORM ONSET ONLY
		
	####################### TEST 1.5 COMPARE A CLUSTER AGAINST ALL OTHERS, ONSET ONLY #######################

	print("####################### TEST 1.5 ONSET #######################")
	
	# get Onset only	
	a_stage_onset <- cTSubSet[, which(cTSubSet$stage=="Onset")]
	
	################################################################ BEGIN ELEONORA
	# get EAE full table ready to be exported
	etabFEAEOnset <- as.data.frame(table(a_stage_onset$slm_0.1))
	colnames(etabFEAEOnset) <- c("cluster", "frequency")

	# export the full table
	write.xlsx(etabFEAEOnset, paste0(saveImg, "ELEONORA_ONSET_cluster_content_full_labels_SUBCLUSTERED_BY_REGION_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)
	################################################################ END ELEONORA	
	
	# loop through all the clusters containing a specific cellType
	for (cS in seq(1, length(clustS), by=1)) {
		
		# we need more than one cluster to perform this comparison
		if (length(unique(clustS))!=1) {
	
			print(paste("cluster position ", clustS[cS]))
		
			# subset SCE to get only cS cluster cells
			# this is for comparing one cluster versus all the other clusters
			a_stage_onset$whatCluster[a_stage_onset$slm_0.1 %in% clustS[cS]] <- "thisOne"
			a_stage_onset$whatCluster[a_stage_onset$slm_0.1 %in% clustS[-cS]] <- "allOther"
			
			################################################################ BEGIN ELEONORA
			# get cells in a specific cluster
			eleonoraSubSet <- a_stage_onset[, a_stage_onset$slm_0.1 %in% clustS[cS]]
			# get EAE full table ready to be exported
			etabFEAE <- as.data.frame(table(eleonoraSubSet$pruned_fine, eleonoraSubSet$condition))
			colnames(etabFEAE) <- c("cellType", "condition", "frequency")

			# export the full table
			write.xlsx(etabFEAE, paste0(saveImg, "ONSET_cluster_content_full_labels_SUBCLUSTERED_BY_REGION_exp_", expN, "_", cellLab[cT], "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.xlsx"), overwrite=T, colNames=T)
			################################################################ END ELEONORA

			# otherwise we may compare, inside a cluster, CTRL versus AD
			
			print(paste0("######## COMPARING ", cellType[cT], " INSIDE ONSET CLUSTER ", clustS[cS], " ########"))
				
			# set ALLOTHER as a reference, to compute fold change
			a_stage_onset$whatCluster <- relevel(as.factor(a_stage_onset$whatCluster), ref="allOther")

			cTSce <- list()
			conD <- unique(a_stage_onset$whatCluster)
			
			if (length(conD)>1) {
			
				# then get CTRL cellType
				cTSce[[conD[1]]] <- a_stage_onset[, which(a_stage_onset$whatCluster==conD[1])]
				# and EAE
				cTSce[[conD[2]]] <- a_stage_onset[, which(a_stage_onset$whatCluster==conD[2])]
				
				if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
					# export logcounts
	#				write.xlsx(cbind.data.frame(rownames(a_stage_onset), logcounts(a_stage_onset)), paste0(saveUpDown, "logcounts", cellType, "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T)
					
					cat("finding markers\n")
					
					# FINDMARKERS for DEGS
					markersG <- findMarkers(a_stage_onset, groups=a_stage_onset$whatCluster, pval.type="any", direction="any")
			
					cC <- "thisOne"
					chosen <- markersG[[cC]]

					# RELAXING THE THRESHOLD, is not done for Onset
					selfMade <- chosen[chosen$FDR < 0.05, ]
					#selfMade <- chosen[chosen$FDR < 0.1, ]

					# and export data
					finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

					# sort by logfoldchange
					if (length(finalGenes$logFC) != 0) {

						# sort by logfoldchange
						sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

						# and export
	#					write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ONE_VERSUS_ALL_OTHER_ONSET_ONLY_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
						write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ONE_VERSUS_ALL_OTHER_ONSET_SUBCLUSTERED_ONLY_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.xlsx"), overwrite=T)
					}
				}
			} else {
				print("not enough conditions to compare")
			}
		} else {
			print("can't compare with only a cluster")
		}
	}
	
	####################### TEST 2 COMPARE A CLUSTER AGAINST ALL OTHERS #######################
	
	print("####################### TEST 2 #######################")
	
	# get cTSubSet, i.e. all the dataset
	cTSubSetTestTwo <- cTSubSet
	
	# loop through all the clusters containing a specific cellType
	for (cS in seq(1, length(clustS), by=1)) {
	
		
		# we need more than one cluster to perform this comparison
		if (length(unique(clustS))!=1) {
		
			print(paste("cluster position ", clustS[cS]))
		
			# subset SCE to get only cS cluster cells
			# this is for comparing one cluster versus all the other clusters
			cTSubSetTestTwo$whatCluster[cTSubSet$slm_0.1 %in% clustS[cS]] <- "thisOne"
			cTSubSetTestTwo$whatCluster[cTSubSet$slm_0.1 %in% clustS[-cS]] <- "allOther"

			# otherwise we may compare, inside a cluster, CTRL versus AD
			
			print(paste0("######## COMPARING ", cellType[cT], " INSIDE CLUSTER ", clustS[cS], " ########"))
				
			# set ALLOTHER as a reference, to compute fold change
			cTSubSetTestTwo$whatCluster <- relevel(as.factor(cTSubSetTestTwo$whatCluster), ref="allOther")

			cTSce <- list()
			conD <- unique(cTSubSetTestTwo$whatCluster)
			
			if (length(conD)>1) {
			
				# then get CTRL cellType
				cTSce[[conD[1]]] <- cTSubSetTestTwo[, which(cTSubSetTestTwo$whatCluster==conD[1])]
				# and EAE
				cTSce[[conD[2]]] <- cTSubSetTestTwo[, which(cTSubSetTestTwo$whatCluster==conD[2])]
				
				if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
					# export logcounts
	#				write.xlsx(cbind.data.frame(rownames(cTSubSetTestTwo), logcounts(cTSubSetTestTwo)), paste0(saveUpDown, "logcounts", cellType, "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_", cellLab[cT], "_SLM_ZERO_UNO.xlsx"), overwrite=T)
					
					cat("finding markers\n")
					
					# FINDMARKERS for DEGS
					markersG <- findMarkers(cTSubSetTestTwo, groups=cTSubSetTestTwo$whatCluster, pval.type="any", direction="any")
			
					cC <- "thisOne"
					chosen <- markersG[[cC]]

					# RELAXING THE THRESHOLD, is not done for Onset
					selfMade <- chosen[chosen$FDR < 0.05, ]
					#selfMade <- chosen[chosen$FDR < 0.1, ]

					# and export data
					finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

					# sort by logfoldchange
					if (length(finalGenes$logFC) != 0) {

						# sort by logfoldchange
						sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

						# and export
	#					write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ONE_VERSUS_ALL_OTHER_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.csv"), row.names=F, quote=F, sep="\t")
						write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ONE_VERSUS_ALL_OTHER_SUBCLUSTERED_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_UNO.xlsx"), overwrite=T)
					}
				}
			} else {
				print("not enough conditions to compare")
			}
		} else {
			print("can't compare with only a cluster")
		}
	}
}
