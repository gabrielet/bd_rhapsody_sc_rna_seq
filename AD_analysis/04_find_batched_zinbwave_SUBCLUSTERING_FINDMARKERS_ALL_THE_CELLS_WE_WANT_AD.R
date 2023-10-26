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
rawPath <- "/home/gabriele/work/cbmc/scrnaseq/raw_data_new/" ; expN <- "AD_batched"

# set cellTypes
cellType <- c("T.8|T.CD8", "Tgd", "Neutrophils", "\\(T.4|\\(T.CD4", "B cells", "DC|Macrophages|Microglia|Monocytes")
cellLab <-  c("CD8", "Tgd", "Neutrophils", "CD4", "B_cells", "Myeloid")

# loop through cell type
for (cT in seq(1, length(cellLab), by=1)){

	print("")
	print(paste0("analysing all the genes ", cellLab[cT], " for exp ", expN, " clustered with Zinbwave"))
	print("")

	# then do everything
	saveIn <- paste0(rawPath, "batched_AD/")
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
	# performed using the zinbwaved data, in order to have a good clustering but also all
	# the genes, instead of the best 1000 used for zinbwave
	perC <- 1000
	kZinb <- 20	
	
	cTSubSet <- readRDS(paste0(saveIn, "QCed_data_with_clusters_", expN, "_K_", kZinb, "_top_", perC, "_", cellLab[cT], "_ONLY.Rds"))
	
	# CONTE GREZZE ELEONORA
	raw_counts <- assay(cTSubSet, "counts")
	colnames(raw_counts) <- paste(colnames(cTSubSet), cTSubSet$condition, cTSubSet$gender, cTSubSet$pruned_fine, "cluster", cTSubSet$slm_0.5, sep="_")
	with_row_names <- cbind.data.frame(rownames(raw_counts), raw_counts)
	colnames(with_row_names) <- c("genes", colnames(raw_counts))
	
	write.table(with_row_names, paste0(saveIn, "counts_CONTE_GREZZE_SUBCLUSTERED_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")
	#######################

	################################################# EXPORT THE NEW TABLE ELEONORA ASKED FEBRUARY #################################################

	superfull_table <- assay(cTSubSet, "logcounts")
	colnames(superfull_table) <- paste(colnames(cTSubSet), cTSubSet$condition, cTSubSet$gender, cellLab[cT], "cluster", cTSubSet$slm_0.5, sep="_")
	with_row_names <- cbind.data.frame(rownames(superfull_table), superfull_table)
	colnames(with_row_names) <- c("genes", colnames(superfull_table))
	
	write.table(with_row_names, paste0(saveIn, "logcounts_all_the_selected_cells_version_barbara_elena_eleonora_only_SUBCLUSTERING_PRIMA_DI_LOGNORMCOUNTS_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")

	# find the clusters
	fR <- table(cTSubSet$slm_0.5)
	clustS <- names(fR[which(fR!=0)])

	# get label ready
	cTSubSet$fine_to_main <- gsub(" ", "", unlist(lapply(strsplit(cTSubSet$pruned_fine, "\\("), `[[`, 1)))
	
	##################### MAIN labels	
	
	# plotting heatmap with cell type
	tab_fine_to_main <- table(Assigned=cTSubSet$fine_to_main, Cluster=cTSubSet$slm_0.5, Condition=cTSubSet$condition)

	# load objs
	eTabCNTX <- NULL
	eTabCNWT <- NULL
	if (length(unique(cTSubSet$fine_to_main)) > 1 && length(unique(cTSubSet$slm_0.5)) > 1) {
		# get 3xTG full table ready to be exported
		eTabTX <- as.data.frame(tab_fine_to_main[, , "3xTG"])
		eTabRNTX <- cbind.data.frame(cellType=rownames(eTabTX), eTabTX)
		eTabCNTX <- eTabRNTX[, c(2,3,4)]
		colnames(eTabCNTX) <- c("CellType", "Cluster", "Freq")

		# and WT
		eTabWT <- as.data.frame(tab_fine_to_main[, , "WT"])
		eTabRNWT <- cbind.data.frame(cellType=rownames(eTabWT), eTabWT)
		eTabCNWT <- eTabRNWT[, c(2,3,4)]
		colnames(eTabCNWT) <- c("CellType", "Cluster", "Freq")
	} else if (length(unique(cTSubSet$fine_to_main)) == 1) {
		# get 3xTG full table ready to be exported
		eTabTX <- as.data.frame(tab_fine_to_main[, , "3xTG"])
		eTabCNTX <- cbind.data.frame(rep(unique(cTSubSet$fine_to_main), length(eTabTX)), rownames(eTabTX), eTabTX)
		colnames(eTabCNTX) <- c("CellType", "Cluster", "Freq")

		# and WT
		eTabWT <- as.data.frame(tab_fine_to_main[, , "WT"])
		eTabCNWT <- cbind.data.frame(rep(unique(cTSubSet$fine_to_main), length(eTabWT)), rownames(eTabWT), eTabWT)
		colnames(eTabCNWT) <- c("CellType", "Cluster", "Freq")
	} else if (length(unique(cTSubSet$slm_0.5)) == 1) {
		# get 3xTG full table ready to be exported
		eTabTX <- as.data.frame(tab_fine_to_main[, , "3xTG"])
		eTabCNTX <- cbind.data.frame(rownames(eTabTX), rep(unique(cTSubSet$slm_0.5), length(eTabTX)), eTabTX)
		colnames(eTabCNTX) <- c("CellType", "Cluster", "Freq")

		# and WT
		eTabWT <- as.data.frame(tab_fine_to_main[, , "WT"])
		eTabCNWT <- cbind.data.frame(rownames(eTabWT), rep(unique(cTSubSet$slm_0.5), length(eTabWT)), eTabWT)
		colnames(eTabCNWT) <- c("CellType", "Cluster", "Freq")
	}

	# export both tables
	write.xlsx(eTabCNTX, paste0(saveImg, "TX_cluster_SUBCLUSTERED_content_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), colNames=F, overwrite=T)

	write.xlsx(eTabCNWT, paste0(saveImg, "WT_cluster_SUBCLUSTERED_content_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), colNames=F, overwrite=T)

	##################### FULL labels

	tab_full_labels <- table(Assigned=cTSubSet$pruned_fine, Cluster=cTSubSet$slm_0.5, Condition=cTSubSet$condition)

	# load objs
	eTabcnFTX <- NULL
	eTabcnFWT <- NULL
	if (length(unique(cTSubSet$pruned_fine)) > 1 && length(unique(cTSubSet$slm_0.5)) > 1) {
		# get 3xTG full table ready to be exported
		etabFTX <- as.data.frame(tab_full_labels[, , "3xTG"])
		eTabrnFTX <- cbind.data.frame(cellType=rownames(etabFTX), etabFTX)
		eTabcnFTX <- eTabrnFTX[, c(2,3,4)]
		colnames(eTabcnFTX) <- c("CellType", "Cluster", "Freq")

		# and WT
		etabFWT <- as.data.frame(tab_full_labels[, , "WT"])
		eTabrnFWT <- cbind.data.frame(cellType=rownames(etabFWT), etabFWT)
		eTabcnFWT <- eTabrnFWT[, c(2,3,4)]
		colnames(eTabcnFWT) <- c("CellType", "Cluster", "Freq")
	} else if (length(unique(cTSubSet$pruned_fine)) == 1) {
		# get 3xTG full table ready to be exported
		etabFTX <- as.data.frame(tab_full_labels[, , "3xTG"])
		eTabcnFTX <- cbind.data.frame(rep(unique(cTSubSet$pruned_fine), length(etabFTX)), rownames(etabFTX), etabFTX)
		colnames(eTabcnFTX) <- c("CellType", "Cluster", "Freq")

		# and WT
		etabFWT <- as.data.frame(tab_full_labels[, , "WT"])
		eTabcnFWT <- cbind.data.frame(rep(unique(cTSubSet$pruned_fine), length(etabFWT)), rownames(etabFWT), etabFWT)
		colnames(eTabcnFWT) <- c("CellType", "Cluster", "Freq")
	} else if (length(unique(cTSubSet$slm_0.5)) == 1) {
		# get 3xTG full table ready to be exported
		etabFTX <- as.data.frame(tab_full_labels[, , "3xTG"])
		eTabcnFTX <- cbind.data.frame(rownames(etabFTX), rep(unique(cTSubSet$slm_0.5), length(etabFTX)), etabFTX)
		colnames(eTabcnFTX) <- c("CellType", "Cluster", "Freq")

		# and WT
		etabFWT <- as.data.frame(tab_full_labels[, , "WT"])
		eTabcnFWT <- cbind.data.frame(rownames(etabFWT), rep(unique(cTSubSet$slm_0.5), length(etabFWT)), etabFWT)
		colnames(eTabcnFWT) <- c("CellType", "Cluster", "Freq")
	}

	# export the full table
	write.xlsx(eTabcnFTX, paste0(saveImg, "TX_cluster_content_SUBCLUSTERED_full_labels_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), colNames=T, overwrite=T)

	write.xlsx(eTabcnFWT, paste0(saveImg, "WT_cluster_content_SUBCLUSTERED_full_labels_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), colNames=T, overwrite=T)

	##################### FULL GENDER

	tab_GENDER_labels <- table(cTSubSet$pruned_fine, cTSubSet$condition, cTSubSet$gender)

	# load objs
	eTabcnFTX <- NULL
	eTabcnFWT <- NULL

	if (length(unique(cTSubSet$pruned_fine)) > 1 && length(unique(cTSubSet$slm_0.5)) > 1) {
		# get MALE
		etabFTX <- as.data.frame(tab_GENDER_labels[, , "M"])
		eTabrnFTX <- cbind.data.frame(cellType=rownames(etabFTX), etabFTX)
		eTabcnFTX <- eTabrnFTX[, c(2,3,4)]
		colnames(eTabcnFTX) <- c("cellType", "condition", "freq")

		# and FEMALE
		etabFWT <- as.data.frame(tab_GENDER_labels[, , "F"])
		eTabrnFWT <- cbind.data.frame(cellType=rownames(etabFWT), etabFWT)
		eTabcnFWT <- eTabrnFWT[, c(2,3,4)]
		colnames(eTabcnFWT) <- c("cellType", "condition", "freq")
	} else if (length(unique(cTSubSet$pruned_fine)) == 1) {
		# get 3xTG full table ready to be exported
		etabFTX <- as.data.frame(tab_GENDER_labels[, , "M"])
		eTabcnFTX <- cbind.data.frame(rep(unique(cTSubSet$pruned_fine), length(etabFTX)), rownames(etabFTX), etabFTX)
		colnames(eTabcnFTX) <- c("CellType", "Cluster", "Freq")

		# and WT
		etabFWT <- as.data.frame(tab_GENDER_labels[, , "F"])
		eTabcnFWT <- cbind.data.frame(rep(unique(cTSubSet$pruned_fine), length(etabFWT)), rownames(etabFWT), etabFWT)
		colnames(eTabcnFWT) <- c("CellType", "Cluster", "Freq")
	} else if (length(unique(cTSubSet$slm_0.5)) == 1) {
		# get 3xTG full table ready to be exported
		etabFTX <- as.data.frame(tab_GENDER_labels[, , "M"])
		eTabcnFTX <- cbind.data.frame(rownames(etabFTX), rep(unique(cTSubSet$slm_0.5), length(etabFTX)), etabFTX)
		colnames(eTabcnFTX) <- c("CellType", "Cluster", "Freq")

		# and WT
		etabFWT <- as.data.frame(tab_GENDER_labels[, , "F"])
		eTabcnFWT <- cbind.data.frame(rownames(etabFWT), rep(unique(cTSubSet$slm_0.5), length(etabFWT)), etabFWT)
		colnames(eTabcnFWT) <- c("CellType", "Cluster", "Freq")
	}

	# export the full table
	write.xlsx(eTabcnFTX, paste0(saveImg, "MALE_cluster_content_full_labels_SUBCLUSTERED_BY_GENDER_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T, colNames=F)

	write.xlsx(eTabcnFWT, paste0(saveImg, "FEMALE_cluster_content_full_labels_SUBCLUSTERED_BY_GENDER_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T, colNames=F)

	####################### CLEAN AND NORMALISE #######################
	
	# remove genes with too many zeros
	cTSubSet <- cTSubSet[rowMeans(assays(cTSubSet)$counts > 0) > 0.05,]

	# set minmean equal to 0.5 since this is the value used in the previous analysis, i.e. the one in the manuscript
	# compute sizeFactors
	cTSubSet <- computeSumFactors(cTSubSet, min.mean=0.5)
	cTSubSet <- logNormCounts(cTSubSet)
	
	################################################# EXPORT THE NEW TABLE ELENA ASKED #################################################

	superfull_table <- assay(cTSubSet, "logcounts")
	colnames(superfull_table) <- paste(colnames(cTSubSet), cTSubSet$condition, cTSubSet$gender, cTSubSet$pruned_fine, "cluster", cTSubSet$slm_0.5, sep="_")
	with_row_names <- cbind.data.frame(rownames(superfull_table), superfull_table)
	colnames(with_row_names) <- c("genes", colnames(superfull_table))
	
	write.table(with_row_names, paste0(saveIn, "logcounts_all_the_selected_cells_version_barbara_elena_eleonora_SUBCLUSTERED_ONLY_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")
#	write.xlsx(with_row_names, paste0(saveIn, "logcounts_ALL_THE_SELECTED_CELLS_VERSION_BARBARA_ELENA_ELEONORA_ONLY_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T)
	
	################################################# EXPORT THE NEW TABLE ELENA ASKED #################################################
	
	# exporting logcounts for Elena's violins BY CONDITION
	sub_WT <- cTSubSet[, which(cTSubSet$condition=="WT")]
	sub_TX <- cTSubSet[, which(cTSubSet$condition=="3xTG")]
	
	# append info to colnames
	coln_WT <- paste(colnames(sub_WT), sub_WT$condition, sub_WT$gender, sub_WT$pruned_fine, "cluster", sub_WT$slm_0.5, sep="_")
	with_row_names_WT <- cbind.data.frame(rownames(sub_WT), logcounts(sub_WT))
	colnames(with_row_names_WT) <- c("genes", coln_WT)
	
	coln_TX <- paste(colnames(sub_TX), sub_TX$condition, sub_TX$gender, sub_TX$pruned_fine, "cluster", sub_TX$slm_0.5, sep="_")
	with_row_names_TX <- cbind.data.frame(rownames(sub_TX), logcounts(sub_TX))
	colnames(with_row_names_TX) <- c("genes", coln_TX)

	write.table(with_row_names_WT, paste0(saveIn, "logcounts_SUBCLUSTERING_WT_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")
	write.table(with_row_names_TX, paste0(saveIn, "logcounts_SUBCLUSTERING_3xTG_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")

#	write.xlsx(with_row_names_WT, paste0(saveIn, "logcounts_SUBCLUSTERING_WT_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T, colNames=T)
#	write.xlsx(with_row_names_TX, paste0(saveIn, "logcounts_SUBCLUSTERING_3xTG_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T, colNames=T)
	
	# exporting logcounts for Elena's violins BY GENDER
	subMALE <- cTSubSet[, which(cTSubSet$gender=="M")]
	subFEMALE <- cTSubSet[, which(cTSubSet$gender=="F")]
	
	# append info to colnames
	coln_MALE <- paste(colnames(subMALE), subMALE$condition, subMALE$gender, subMALE$pruned_fine, "cluster", subMALE$slm_0.5, sep="_")
	with_row_names_MALE <- cbind.data.frame(rownames(subMALE), logcounts(subMALE))
	colnames(with_row_names_MALE) <- c("genes", coln_MALE)
	
	coln_FEMALE <- paste(colnames(subFEMALE), subFEMALE$condition, subFEMALE$gender, subFEMALE$pruned_fine, "cluster", subFEMALE$slm_0.5, sep="_")
	with_row_names_FEMALE <- cbind.data.frame(rownames(subFEMALE), logcounts(subFEMALE))
	colnames(with_row_names_FEMALE) <- c("genes", coln_FEMALE)

	write.table(with_row_names_MALE, paste0(saveIn, "logcounts_SUBCLUSTERING_MALE_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")
	write.table(with_row_names_FEMALE, paste0(saveIn, "logcounts_SUBCLUSTERING_FEMALE_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")
	
#	write.xlsx(with_row_names_MALE, paste0(saveIn, "logcounts_SUBCLUSTERING_MALE_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T, colNames=T)
#	write.xlsx(with_row_names_FEMALE, paste0(saveIn, "logcounts_SUBCLUSTERING_FEMALE_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T, colNames=T)
	
	# exporting logcounts for Elena's violins BY BOTH
	sub_WT_MALE <- cTSubSet[, which(cTSubSet$conditionBYgender=="WT_M")]
	sub_WT_FEMALE <- cTSubSet[, which(cTSubSet$conditionBYgender=="WT_F")]
	sub_TX_MALE <- cTSubSet[, which(cTSubSet$conditionBYgender=="3xTG_M")]
	sub_TX_FEMALE <- cTSubSet[, which(cTSubSet$conditionBYgender=="3xTG_F")]
	
	# append info to colnames
	if (dim(sub_WT_MALE)[2] > 0) {
		coln_WT_MALE <- paste(colnames(sub_WT_MALE), sub_WT_MALE$condition, sub_WT_MALE$gender, sub_WT_MALE$pruned_fine, "cluster", sub_WT_MALE$slm_0.5, sep="_")
		with_row_names_WT_MALE <- cbind.data.frame(rownames(sub_WT_MALE), logcounts(sub_WT_MALE))
		colnames(with_row_names_WT_MALE) <- c("genes", coln_WT_MALE)
		
		write.table(with_row_names_WT_MALE, paste0(saveIn, "logcounts_SUBCLUSTERING_WT_MALE_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")		
#		write.xlsx(with_row_names_WT_MALE, paste0(saveIn, "logcounts_SUBCLUSTERING_WT_MALE_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T, colNames=T)
	}
	if (dim(sub_WT_FEMALE)[2] > 0) {
		coln_WT_FEMALE <- paste(colnames(sub_WT_FEMALE), sub_WT_FEMALE$condition, sub_WT_FEMALE$gender, sub_WT_FEMALE$pruned_fine, "cluster", sub_WT_FEMALE$slm_0.5, sep="_")
		with_row_names_WT_FEMALE <- cbind.data.frame(rownames(sub_WT_FEMALE), logcounts(sub_WT_FEMALE))
		colnames(with_row_names_WT_FEMALE) <- c("genes", coln_WT_FEMALE)

		write.table(with_row_names_WT_FEMALE, paste0(saveIn, "logcounts_SUBCLUSTERING_WT_FEMALE_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")
#		write.xlsx(with_row_names_WT_FEMALE, paste0(saveIn, "logcounts_SUBCLUSTERING_WT_FEMALE_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T, colNames=T)
	}
	if (dim(sub_TX_MALE)[2] > 0) {		
		coln_TX_MALE <- paste(colnames(sub_TX_MALE), sub_TX_MALE$condition, sub_TX_MALE$gender, sub_TX_MALE$pruned_fine, "cluster", sub_TX_MALE$slm_0.5, sep="_")
		with_row_names_TX_MALE <- cbind.data.frame(rownames(sub_TX_MALE), logcounts(sub_TX_MALE))
		colnames(with_row_names_TX_MALE) <- c("genes", coln_TX_MALE)

		write.table(with_row_names_TX_MALE, paste0(saveIn, "logcounts_SUBCLUSTERING_3xTG_MALE_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")
#		write.xlsx(with_row_names_TX_MALE, paste0(saveIn, "logcounts_SUBCLUSTERING_3xTG_MALE_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T, colNames=T)
	}
	if (dim(sub_TX_FEMALE)[2] > 0) {		
		coln_TX_FEMALE <- paste(colnames(sub_TX_FEMALE), sub_TX_FEMALE$condition, sub_TX_FEMALE$gender, sub_TX_FEMALE$pruned_fine, "cluster", sub_TX_FEMALE$slm_0.5, sep="_")
		with_row_names_TX_FEMALE <- cbind.data.frame(rownames(sub_TX_FEMALE), logcounts(sub_TX_FEMALE))
		colnames(with_row_names_TX_FEMALE) <- c("genes", coln_TX_FEMALE)
		
		write.table(with_row_names_TX_FEMALE, paste0(saveIn, "logcounts_SUBCLUSTERING_3xTG_FEMALE_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")
#		write.xlsx(with_row_names_TX_FEMALE, paste0(saveIn, "logcounts_SUBCLUSTERING_3xTG_FEMALE_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T, colNames=T)
	}
	
	# Cluster composition plots using main labels
	conD <- unique(cTSubSet$conditionBYgender)

	dfReg <- data.frame(table(cTSubSet$slm_0.5, cTSubSet$pruned_fine, cTSubSet$conditionBYgender))

	colnames(dfReg) <- c("Cluster", "CellType", "Condition", "Frequency")
	lims <- levels(factor(cTSubSet$pruned_fine))
	lims <- lims[order(table(cTSubSet$pruned_fine))]

	# plot!
	p <- plot_membership(dfReg, cond=conD[2], col=colorsPal[2])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[2], "_cluster_composition_SUBCLUSTERED_", cellLab[cT], "_exp_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()

	p <- plot_membership(dfReg, cond=conD[1], col=colorsPal[6])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[1], "_cluster_composition_SUBCLUSTERED_", cellLab[cT], "_exp_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()
	
	# plot!
	p <- plot_membership(dfReg, cond=conD[3], col=colorsPal[3])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[3], "_cluster_composition_SUBCLUSTERED_", cellLab[cT], "_exp_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()

	p <- plot_membership(dfReg, cond=conD[4], col=colorsPal[5])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[4], "_cluster_composition_SUBCLUSTERED_", cellLab[cT], "_exp_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()
	
	####################### TEST 1 IN EACH CLUSTER COMPARE STUFF #######################
	
	# loop through all the clusters containing a specific cellType
	for (cS in seq(1, length(clustS), by=1)) {
	
		# PERFORM MALE+FEMALE, WT VS 3X
	
		print(paste("cluster position ", clustS[cS]))
	
		# subset SCE to get only cS cluster cells
		aClust <- cTSubSet[, cTSubSet$slm_0.5==clustS[cS]]
		
		print(paste0("######## COMPARING ", cellType[cT], " INSIDE CLUSTER ", clustS[cS], " ########"))

		if(unique(length(aClust$condition))==2) {

			# set reference condition
			aClust$condition <- relevel(as.factor(aClust$condition), ref="WT")
			
			cTSce <- list()
			conD <- unique(aClust$condition)
			# then get CTRL cellType
			cTSce[[conD[1]]] <- aClust[, which(aClust$condition==conD[1])]
			# and EAE
			cTSce[[conD[2]]] <- aClust[, which(aClust$condition==conD[2])]
			
			if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
				# export logcounts
	#			write.xlsx(cbind.data.frame(rownames(aClust), logcounts(aClust)), paste0(saveUpDown, "logcounts", cellType, "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T)
				
				cat("finding markers\n")
				
				# FINDMARKERS for DEGS
				markersG <- findMarkers(aClust, groups=aClust$condition, pval.type="any", direction="any")
		
				cC <- "3xTG"
				chosen <- markersG[[cC]]

				# RELAXING THE THRESHOLD, is not done for females
				selfMade <- chosen[chosen$FDR < 0.05, ]
				#selfMade <- chosen[chosen$FDR < 0.1, ]

				# and export data
				finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

				# sort by logfoldchange
				if (length(finalGenes$logFC) != 0) {
					sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

					# and export
					#write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_INSIDE_CLUSTER_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")
					write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_SUBCLUSTERED_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T)
				}
			}
		}
		
		# PERFORM MALE ONLY, WT VS 3X
		
		aClustByReg <- cTSubSet[, cTSubSet$slm_0.5==clustS[cS] & cTSubSet$gender=="M"]
		
		print(paste0("######## COMPARING ", cellType[cT], " INSIDE CLUSTER MALE ONLY ", clustS[cS], " ########"))

		cTSce <- list()
		conD <- unique(aClustByReg$condition)
		
		if (length(conD)>1) {

			# set reference condition
			aClustByReg$condition <- relevel(as.factor(aClustByReg$condition), ref="WT")
			
			# then get CTRL cellType
			
			cTSce[[conD[1]]] <- aClustByReg[, which(aClustByReg$condition==conD[1])]
			# and EAE
			cTSce[[conD[2]]] <- aClustByReg[, which(aClustByReg$condition==conD[2])]
			
			if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
				# export logcounts
#				write.xlsx(cbind.data.frame(rownames(aClustByReg), logcounts(aClustByReg)), paste0(saveUpDown, "logcounts", cellType, "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T)
				
				cat("finding markers\n")
				
				# FINDMARKERS for DEGS
				markersG <- findMarkers(aClustByReg, groups=aClustByReg$condition, pval.type="any", direction="any")

				cC <- "3xTG"
				chosen <- markersG[[cC]]

				# RELAXING THE THRESHOLD, is not done for females
				selfMade <- chosen[chosen$FDR < 0.05, ]
				#selfMade <- chosen[chosen$FDR < 0.1, ]

				# and export data
				finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

				# sort by logfoldchange
				if (length(finalGenes$logFC) != 0) {
					# sort by logfoldchange
					sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

					# and export
					#write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_INSIDE_CLUSTER_MALE_ONLY_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")
					write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_INSIDE_CLUSTER_MALE_SUBCLUSTERED_ONLY_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T)
				}
			}
		} else {
			print("not enough conditions to compare")
		}
		
		# PERFORM FEMALE ONLY, WT VS 3X
	
		aClustByReg <- cTSubSet[, cTSubSet$slm_0.5==clustS[cS] & cTSubSet$gender=="F"]

		print(paste0("######## COMPARING ", cellType[cT], " INSIDE CLUSTER FEMALE ONLY ", clustS[cS], " ########"))

		cTSce <- list()
		conD <- unique(aClustByReg$condition)

		if (length(conD)>1) {

			# set reference condition
			aClustByReg$condition <- relevel(as.factor(aClustByReg$condition), ref="WT")
			
			# then get CTRL cellType
			
			cTSce[[conD[1]]] <- aClustByReg[, which(aClustByReg$condition==conD[1])]
			# and EAE
			cTSce[[conD[2]]] <- aClustByReg[, which(aClustByReg$condition==conD[2])]
			
			if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
				
				cat("finding markers\n")
				
				# FINDMARKERS for DEGS
				markersG <- findMarkers(aClustByReg, groups=aClustByReg$condition, pval.type="any", direction="any")

				cC <- "3xTG"
				chosen <- markersG[[cC]]

				# RELAXING THE THRESHOLD, is not done for females
				selfMade <- chosen[chosen$FDR < 0.05, ]
				#selfMade <- chosen[chosen$FDR < 0.1, ]

				# and export data
				finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)
				
				if (length(finalGenes$logFC) != 0) {

					# sort by logfoldchange
					sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

					# and export
					#write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_INSIDE_CLUSTER_FEMALE_ONLY_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")
					write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_INSIDE_CLUSTER_FEMALE_SUBCLUSTERED_ONLY", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T)
				}
			}
		} else {
			print("not enough conditions to compare")
		}
	}
		
	####################### TEST 1.5 COMPARE A CLUSTER AGAINST ALL OTHERS, MALE ONLY #######################

	print("####################### TEST 1.5 MALE #######################")
	
	# get male only
	aGenderMALE <- cTSubSet[, which(cTSubSet$gender=="M")]
	
	################################################################ BEGIN ELEONORA
	# get 3xTG full table ready to be exported
	etabFTXMALE <- as.data.frame(table(aGenderMALE$slm_0.5))
	colnames(etabFTXMALE) <- c("cluster", "frequency")

	# export the full table
	write.xlsx(etabFTXMALE, paste0(saveImg, "ANCORA_PIU_NUOVO_MALE_cluster_content_full_labels_SUBCLUST_BY_GENDER_SUBCLUSTERED_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T, colNames=T)
	################################################################ END ELEONORA
	
	# loop through all the clusters containing a specific cellType
	for (cS in seq(1, length(clustS), by=1)) {
		
		# we need more than one cluster to perform this comparison
		if (length(unique(clustS))!=1) {
	
			print(paste("cluster position ", clustS[cS]))
		
			# subset SCE to get only cS cluster cells
			# this is for comparing one cluster versus all the other clusters
			aGenderMALE$whatCluster[aGenderMALE$slm_0.5 %in% clustS[cS]] <- "thisOne"
			aGenderMALE$whatCluster[aGenderMALE$slm_0.5 %in% clustS[-cS]] <- "allOther"

			################################################################ BEGIN ELEONORA
			# get cells in a specific cluster
			eleonoraSubSet <- aGenderMALE[, aGenderMALE$slm_0.5 %in% clustS[cS]]
			# get 3xTG full table ready to be exported
			etabFTX <- as.data.frame(table(eleonoraSubSet$pruned_fine, eleonoraSubSet$condition))
			colnames(etabFTX) <- c("cellType", "condition", "frequency")

			# export the full table
			write.xlsx(etabFTX, paste0(saveImg, "NUOVO_MALE_cluster_content_full_labels_SUBCLUST_BY_GENDER_SUBCLUSTERED_exp_", expN, "_", cellLab[cT], "by_cluster_", clustS[cS], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T, colNames=T)
			################################################################ END ELEONORA
			
			# otherwise we may compare, inside a cluster, WT versus AD
			
			print(paste0("######## COMPARING ", cellType[cT], " INSIDE MALE CLUSTER ", clustS[cS], " ########"))
				
			# set ALLOTHER as a reference, to compute fold change
			aGenderMALE$whatCluster <- relevel(as.factor(aGenderMALE$whatCluster), ref="allOther")

			cTSce <- list()
			conD <- unique(aGenderMALE$whatCluster)
			
			if (length(conD)>1) {
			
				# then get CTRL cellType
				cTSce[[conD[1]]] <- aGenderMALE[, which(aGenderMALE$whatCluster==conD[1])]
				# and EAE
				cTSce[[conD[2]]] <- aGenderMALE[, which(aGenderMALE$whatCluster==conD[2])]
				
				if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
					# export logcounts
	#				write.xlsx(cbind.data.frame(rownames(aGenderMALE), logcounts(aGenderMALE)), paste0(saveUpDown, "logcounts", cellType, "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T)
					
					cat("finding markers\n")
					
					# FINDMARKERS for DEGS
					markersG <- findMarkers(aGenderMALE, groups=aGenderMALE$whatCluster, pval.type="any", direction="any")
			
					cC <- "thisOne"
					chosen <- markersG[[cC]]

					# RELAXING THE THRESHOLD, is not done for females
					selfMade <- chosen[chosen$FDR < 0.05, ]
					#selfMade <- chosen[chosen$FDR < 0.1, ]

					# and export data
					finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

					# sort by logfoldchange
					sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

					# and export
					#write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ONE_VERSUS_ALL_OTHER_MALE_ONLY", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")
					write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ONE_VERSUS_ALL_OTHER_MALE_ONLY_SUBCLUSTERED_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T)
				}
			} else {
				print("not enough conditions to compare")
			}
		} else {
			print("can't compare with only a cluster")
		}
	}

	################################################################# ELENA, MALE ONLY WT VERSUS 3XTG

	# otherwise we may compare, inside a cluster, WT versus AD
	
	print(paste0("######## COMPARING ", cellType[cT], " INSIDE MALE ########"))
		
	# set ALLOTHER as a reference, to compute fold change
	aGenderMALE$condition <- relevel(as.factor(aGenderMALE$condition), ref="WT")

	cTSce <- list()
	conD <- unique(aGenderMALE$condition)
	
	if (length(conD)>1) {
	
		# then get WT cellType
		cTSce[[conD[1]]] <- aGenderMALE[, which(aGenderMALE$condition==conD[1])]
		# and 3xTG
		cTSce[[conD[2]]] <- aGenderMALE[, which(aGenderMALE$condition==conD[2])]
		
		if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
			# export logcounts
#				write.xlsx(cbind.data.frame(rownames(aGenderMALE), logcounts(aGenderMALE)), paste0(saveUpDown, "logcounts", cellType, "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T)
			
			cat("finding markers\n")
			
			# FINDMARKERS for DEGS
			markersG <- findMarkers(aGenderMALE, groups=aGenderMALE$condition, pval.type="any", direction="any")
	
			cC <- "3xTG"
			chosen <- markersG[[cC]]

			# RELAXING THE THRESHOLD, is not done for females
			selfMade <- chosen[chosen$FDR < 0.05, ]
			#selfMade <- chosen[chosen$FDR < 0.1, ]

			# and export data
			finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

			# sort by logfoldchange
			sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

			# and export
			#write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ONE_VERSUS_ALL_OTHER_MALE_ONLY", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")
			write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_3x_versus_WT_MALE_ONLY_SUBCLUSTERED_", cellLab[cT], "_CELLS_exp_", expN, "_SLM_ZERO_CINQUE.xlsx"), overwrite=T)
		}
	} else {
		print("not enough conditions to compare")
	}
		
	# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY# PERFORM FEMALE ONLY
		
	####################### TEST 1.5 COMPARE A CLUSTER AGAINST ALL OTHERS, FEMALE ONLY #######################

	print("####################### TEST 1.5 FEMALE #######################")
	
	# get females only	
	aGenderFEMALE <- cTSubSet[, which(cTSubSet$gender=="F")]
	
	################################################################ BEGIN ELEONORA
	# get 3xTG full table ready to be exported
	etabFTXFEMALE <- as.data.frame(table(aGenderFEMALE$slm_0.5))
	colnames(etabFTXFEMALE) <- c("cluster", "frequency")

	# export the full table
	write.xlsx(etabFTXFEMALE, paste0(saveImg, "ANCORA_PIU_NUOVO_FEMALE_cluster_content_full_labels_SUBCLUST_BY_GENDER_SUBCLUSTERED_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T, colNames=T)
	################################################################ END ELEONORA	
	
	# loop through all the clusters containing a specific cellType
	for (cS in seq(1, length(clustS), by=1)) {
		
		# we need more than one cluster to perform this comparison
		if (length(unique(clustS))!=1) {
	
			print(paste("cluster position ", clustS[cS]))
		
			# subset SCE to get only cS cluster cells
			# this is for comparing one cluster versus all the other clusters
			aGenderFEMALE$whatCluster[aGenderFEMALE$slm_0.5 %in% clustS[cS]] <- "thisOne"
			aGenderFEMALE$whatCluster[aGenderFEMALE$slm_0.5 %in% clustS[-cS]] <- "allOther"
			
			################################################################ BEGIN ELEONORA
			# get cells in a specific cluster
			eleonoraSubSet <- aGenderFEMALE[, aGenderFEMALE$slm_0.5 %in% clustS[cS]]
			# get 3xTG full table ready to be exported
			etabFTX <- as.data.frame(table(eleonoraSubSet$pruned_fine, eleonoraSubSet$condition))
			colnames(etabFTX) <- c("cellType", "condition", "frequency")

			# export the full table
			write.xlsx(etabFTX, paste0(saveImg, "NUOVO_FEMALE_cluster_content_full_labels_SUBCLUST_BY_GENDER_SUBCLUSTERED_exp_", expN, "_", cellLab[cT], "by_cluster_", clustS[cS], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T, colNames=T)
			################################################################ END ELEONORA

			# otherwise we may compare, inside a cluster, WT versus AD
			
			print(paste0("######## COMPARING ", cellType[cT], " INSIDE FEMALE CLUSTER ", clustS[cS], " ########"))
				
			# set ALLOTHER as a reference, to compute fold change
			aGenderFEMALE$whatCluster <- relevel(as.factor(aGenderFEMALE$whatCluster), ref="allOther")

			cTSce <- list()
			conD <- unique(aGenderFEMALE$whatCluster)
			
			if (length(conD)>1) {
			
				# then get CTRL cellType
				cTSce[[conD[1]]] <- aGenderFEMALE[, which(aGenderFEMALE$whatCluster==conD[1])]
				# and EAE
				cTSce[[conD[2]]] <- aGenderFEMALE[, which(aGenderFEMALE$whatCluster==conD[2])]
				
				if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
					# export logcounts
	#				write.xlsx(cbind.data.frame(rownames(aGenderFEMALE), logcounts(aGenderFEMALE)), paste0(saveUpDown, "logcounts", cellType, "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T)
					
					cat("finding markers\n")
					
					# FINDMARKERS for DEGS
					markersG <- findMarkers(aGenderFEMALE, groups=aGenderFEMALE$whatCluster, pval.type="any", direction="any")
			
					cC <- "thisOne"
					chosen <- markersG[[cC]]

					# RELAXING THE THRESHOLD, is not done for females
					selfMade <- chosen[chosen$FDR < 0.05, ]
					#selfMade <- chosen[chosen$FDR < 0.1, ]

					# and export data
					finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

					# sort by logfoldchange
					sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

					# and export
					#write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ONE_VERSUS_ALL_OTHER_FEMALE_ONLY_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")
					write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ONE_VERSUS_ALL_OTHER_FEMALE_ONLY_SUBCLUSTERED_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T)
				}
			} else {
				print("not enough conditions to compare")
			}
		} else {
			print("can't compare with only a cluster")
		}
	}
	
	################################################################# ELENA, FEMALE ONLY, WT VERSUS 3XTG
		
	# otherwise we may compare, inside a cluster, WT versus AD
	
	print(paste0("######## COMPARING ", cellType[cT], " in FEMALEs ########"))
		
	# set ALLOTHER as a reference, to compute fold change
	aGenderFEMALE$condition <- relevel(as.factor(aGenderFEMALE$condition), ref="WT")

	cTSce <- list()
	conD <- unique(aGenderFEMALE$condition)
	
	if (length(conD)>1) {
	
		# then get WT cellType
		cTSce[[conD[1]]] <- aGenderFEMALE[, which(aGenderFEMALE$condition==conD[1])]
		# and 3xTG
		cTSce[[conD[2]]] <- aGenderFEMALE[, which(aGenderFEMALE$condition==conD[2])]
		
		if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
			# export logcounts
#				write.xlsx(cbind.data.frame(rownames(aGenderFEMALE), logcounts(aGenderFEMALE)), paste0(saveUpDown, "logcounts", cellType, "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T)
			
			cat("finding markers\n")
			
			# FINDMARKERS for DEGS
			markersG <- findMarkers(aGenderFEMALE, groups=aGenderFEMALE$condition, pval.type="any", direction="any")
	
			cC <- "3xTG"
			chosen <- markersG[[cC]]

			# RELAXING THE THRESHOLD, is not done for females
			selfMade <- chosen[chosen$FDR < 0.05, ]
			#selfMade <- chosen[chosen$FDR < 0.1, ]

			# and export data
			finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

			# sort by logfoldchange
			sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

			# and export
			#write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ONE_VERSUS_ALL_OTHER_FEMALE_ONLY", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")
			write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_3x_versus_WT_FEMALE_ONLY_SUBCLUSTERED_", cellLab[cT], "_CELLS_exp_", expN, "_SLM_ZERO_CINQUE.xlsx"), overwrite=T)
		}
	} else {
		print("not enough conditions to compare")
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
			cTSubSetTestTwo$whatCluster[cTSubSet$slm_0.5 %in% clustS[cS]] <- "thisOne"
			cTSubSetTestTwo$whatCluster[cTSubSet$slm_0.5 %in% clustS[-cS]] <- "allOther"

			# otherwise we may compare, inside a cluster, WT versus AD
			
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
	#				write.xlsx(cbind.data.frame(rownames(cTSubSetTestTwo), logcounts(cTSubSetTestTwo)), paste0(saveUpDown, "logcounts", cellType, "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_", cellLab[cT], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T)
					
					cat("finding markers\n")
					
					# FINDMARKERS for DEGS
					markersG <- findMarkers(cTSubSetTestTwo, groups=cTSubSetTestTwo$whatCluster, pval.type="any", direction="any")
			
					cC <- "thisOne"
					chosen <- markersG[[cC]]

					# RELAXING THE THRESHOLD, is not done for females
					selfMade <- chosen[chosen$FDR < 0.05, ]
					#selfMade <- chosen[chosen$FDR < 0.1, ]

					# and export data
					finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

					# sort by logfoldchange
					sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

					# and export
					#write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ONE_VERSUS_ALL_OTHER_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_CINQUE.csv"), row.names=F, quote=F, sep="\t")
					write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ONE_VERSUS_ALL_OTHER_SUBCLUSTERED_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_CINQUE.xlsx"), overwrite=T)
				}
			} else {
				print("not enough conditions to compare")
			}
		} else {
			print("can't compare with only a cluster")
		}
	}
}
