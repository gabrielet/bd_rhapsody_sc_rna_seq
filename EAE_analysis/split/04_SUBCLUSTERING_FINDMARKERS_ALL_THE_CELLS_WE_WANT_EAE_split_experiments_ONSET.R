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
source("../../00_functions.R")

# initialise label
lbl <- 1

# select experiment
expN <- "Onset"

print(paste0("analysing exp ", expN))

# set paths
rawPath <- "/home/gabriele/work/cbmc/scrnaseq/raw_data_new/"
chosenRed <- "corrected" ; dimRedMet <- "PCA" ; dimRedPar <- "PCA" # for PCA

# set cellTypes
cellType <- c("T.8|T.CD8", "Tgd", "Neutrophils", "\\(T.4|\\(T.CD4", "B cells", "DC|Macrophages|Microglia|Monocytes")
cellLab <-  c("CD8", "Tgd", "Neutrophils", "CD4", "B_cells", "Myeloid")

# loop through cell type
for (cT in seq(1, length(cellLab), by=1)){

	print("")
	print(paste0("analysing all the genes ", cellLab[cT], " for exp ", expN, " clustered with slm"))
	print("")

	# then do everything
	saveIn <- paste0(rawPath, "EAE_", expN, "/")
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
	
	# load data
	cTSubSet <- readRDS(paste0(saveIn, "QCed_data_with_clusters_", dimRedMet, "_", expN, "_", cellLab[cT], "_ONLY.Rds"))
	
	# CONTE GREZZE ELEONORA
	raw_counts <- assay(cTSubSet, "counts")
	colnames(raw_counts) <- paste(colnames(cTSubSet), cTSubSet$condition, cTSubSet$pruned_fine, "cluster", cTSubSet$slm_0.3, sep="_")
	with_row_names <- cbind.data.frame(rownames(raw_counts), raw_counts)
	colnames(with_row_names) <- c("genes", colnames(raw_counts))
	
	write.table(with_row_names, paste0(saveIn, "counts_CONTE_GREZZE_SUBCLUSTERED_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_TRE.csv"), row.names=F, quote=F, sep="\t")
	#######################

	################################################# EXPORT THE NEW TABLE ELEONORA ASKED FEBRUARY #################################################

	superfull_table <- assay(cTSubSet, "logcounts")
	colnames(superfull_table) <- paste(colnames(cTSubSet), cTSubSet$condition, cTSubSet$pruned_fine, "cluster", cTSubSet$slm_0.3, sep="_")
	with_row_names <- cbind.data.frame(rownames(superfull_table), superfull_table)
	colnames(with_row_names) <- c("genes", colnames(superfull_table))
	
	write.table(with_row_names, paste0(saveIn, "logcounts_all_the_selected_cells_version_barbara_elena_eleonora_only_SUBCLUSTERING_PRIMA_DI_LOGNORMCOUNTS_", expN, "_", cellLab[cT], "_SLM_ZERO_TRE.csv"), row.names=F, quote=F, sep="\t")
	
	# find the clusters
	fR <- table(cTSubSet$slm_0.3)
	clustS <- names(fR[which(fR!=0)])

	# get label ready
	cTSubSet$fine_to_main <- gsub(" ", "", unlist(lapply(strsplit(cTSubSet$pruned_fine, "\\("), `[[`, 1)))
	
	##################### MAIN labels
	
	# plotting heatmap with cell type
	tab_fine_to_main <- table(Assigned=cTSubSet$fine_to_main, Cluster=cTSubSet$slm_0.3, Condition=cTSubSet$condition)

	eTabCNEAE <- NULL
	eTabCNCTRL <- NULL

	tab_full_labels <- table(Assigned=cTSubSet$pruned_fine, Cluster=cTSubSet$slm_0.3, Condition=cTSubSet$condition)

	# load objs
	eTabcnFEAE <- NULL
	eTabcnFCTRL <- NULL

	if (length(unique(cTSubSet$condition))==2) {	
		if (length(unique(cTSubSet$fine_to_main)) > 1 && length(unique(cTSubSet$slm_0.3)) > 1) {
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
		} else if (length(unique(cTSubSet$slm_0.3)) == 1) {
			# get EAE full table ready to be exported
			eTabEAE <- as.data.frame(tab_fine_to_main[, , "EAE"])
			eTabCNEAE <- cbind.data.frame(rownames(eTabEAE), rep(unique(cTSubSet$slm_0.3), length(eTabEAE)), eTabEAE)
			colnames(eTabCNEAE) <- c("CellType", "Cluster", "Freq")

			# and CTRL
			eTabCTRL <- as.data.frame(tab_fine_to_main[, , "CTRL"])
			eTabCNCTRL <- cbind.data.frame(rownames(eTabCTRL), rep(unique(cTSubSet$slm_0.3), length(eTabCTRL)), eTabCTRL)
			colnames(eTabCNCTRL) <- c("CellType", "Cluster", "Freq")
		}
	
		# export both tables
		write.xlsx(eTabCNEAE, paste0(saveImg, "EAE_cluster_SUBCLUSTERED_content_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_TRE.xlsx"), colNames=F, overwrite=T)

		write.xlsx(eTabCNCTRL, paste0(saveImg, "CTRL_cluster_SUBCLUSTERED_content_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_TRE.xlsx"), colNames=F, overwrite=T)
		
		##################### FULL labels
		
		if (length(unique(cTSubSet$pruned_fine)) > 1 && length(unique(cTSubSet$slm_0.3)) > 1) {
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
		} else if (length(unique(cTSubSet$slm_0.3)) == 1) {
			# get EAE full table ready to be exported
			etabFEAE <- as.data.frame(tab_full_labels[, , "EAE"])
			eTabcnFEAE <- cbind.data.frame(rownames(etabFEAE), rep(unique(cTSubSet$slm_0.3), length(etabFEAE)), etabFEAE)
			colnames(eTabcnFEAE) <- c("CellType", "Cluster", "Freq")

			# and CTRL
			etabFCTRL <- as.data.frame(tab_full_labels[, , "CTRL"])
			eTabcnFCTRL <- cbind.data.frame(rownames(etabFCTRL), rep(unique(cTSubSet$slm_0.3), length(etabFCTRL)), etabFCTRL)
			colnames(eTabcnFCTRL) <- c("CellType", "Cluster", "Freq")
		}

		# export the full table
		write.xlsx(eTabcnFEAE, paste0(saveImg, "EAE_cluster_content_SUBCLUSTERED_full_labels_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_TRE.xlsx"), colNames=T, overwrite=T)

		write.xlsx(eTabcnFCTRL, paste0(saveImg, "CTRL_cluster_content_SUBCLUSTERED_full_labels_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_TRE.xlsx"), colNames=T, overwrite=T)
	} else {
		# which is the only condition
		cond <- unique(cTSubSet$condition)
	
		if (length(unique(cTSubSet$fine_to_main)) > 1 && length(unique(cTSubSet$slm_0.3)) > 1) {
			# get EAE full table ready to be exported
			eTabEAE <- as.data.frame(tab_fine_to_main[, , cond])
			eTabRNEAE <- cbind.data.frame(cellType=rownames(eTabEAE), eTabEAE)
			eTabCNEAE <- eTabRNEAE[, c(2,3,4)]
			colnames(eTabCNEAE) <- c("CellType", "Cluster", "Freq")
		} else if (length(unique(cTSubSet$fine_to_main)) == 1) {
			# get EAE full table ready to be exported
			eTabEAE <- as.data.frame(tab_fine_to_main[, , cond])
			eTabCNEAE <- cbind.data.frame(rep(unique(cTSubSet$fine_to_main), length(eTabEAE)), rownames(eTabEAE), eTabEAE)
			colnames(eTabCNEAE) <- c("CellType", "Cluster", "Freq")
		} else if (length(unique(cTSubSet$slm_0.3)) == 1) {
			# get EAE full table ready to be exported
			eTabEAE <- as.data.frame(tab_fine_to_main[, , cond])
			eTabCNEAE <- cbind.data.frame(rownames(eTabEAE), rep(unique(cTSubSet$slm_0.3), length(eTabEAE)), eTabEAE)
			colnames(eTabCNEAE) <- c("CellType", "Cluster", "Freq")
		}
	
		# export both tables
		write.xlsx(eTabCNEAE, paste0(saveImg, cond, "_cluster_SUBCLUSTERED_content_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_TRE.xlsx"), colNames=F, overwrite=T)

		##################### FULL labels
		
		if (length(unique(cTSubSet$pruned_fine)) > 1 && length(unique(cTSubSet$slm_0.3)) > 1) {
			# get EAE full table ready to be exported
			etabFEAE <- as.data.frame(tab_full_labels[, , cond])
			eTabrnFEAE <- cbind.data.frame(cellType=rownames(etabFEAE), etabFEAE)
			eTabcnFEAE <- eTabrnFEAE[, c(2,3,4)]
			colnames(eTabcnFEAE) <- c("CellType", "Cluster", "Freq")
		} else if (length(unique(cTSubSet$pruned_fine)) == 1) {
			# get EAE full table ready to be exported
			etabFEAE <- as.data.frame(tab_full_labels[, , cond])
			eTabcnFEAE <- cbind.data.frame(rep(unique(cTSubSet$pruned_fine), length(etabFEAE)), rownames(etabFEAE), etabFEAE)
			colnames(eTabcnFEAE) <- c("CellType", "Cluster", "Freq")
		} else if (length(unique(cTSubSet$slm_0.3)) == 1) {
			# get EAE full table ready to be exported
			etabFEAE <- as.data.frame(tab_full_labels[, , cond])
			eTabcnFEAE <- cbind.data.frame(rownames(etabFEAE), rep(unique(cTSubSet$slm_0.3), length(etabFEAE)), etabFEAE)
			colnames(eTabcnFEAE) <- c("CellType", "Cluster", "Freq")
		}

		# export the full table
		write.xlsx(eTabcnFEAE, paste0(saveImg, cond, "_cluster_content_SUBCLUSTERED_full_labels_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_TRE.xlsx"), colNames=T, overwrite=T)
	
	}

	####################### CLEAN AND NORMALISE #######################
	
	# remove genes with too many zeros
	cTSubSet <- cTSubSet[rowMeans(assays(cTSubSet)$counts > 0) > 0.05, ]

	# set minmean equal to 0.5 since this is the value used in the previous analysis, i.e. the one in the manuscript
	# compute sizeFactors
	cTSubSet <- computeSumFactors(cTSubSet, min.mean=0.5)
	cTSubSet <- logNormCounts(cTSubSet)
	
	################################################# EXPORT THE NEW TABLE ELENA ASKED #################################################

	superfull_table <- assay(cTSubSet, "logcounts")
	colnames(superfull_table) <- paste(colnames(cTSubSet), cTSubSet$condition, cTSubSet$pruned_fine, "cluster", cTSubSet$slm_0.3, sep="_")
	with_row_names <- cbind.data.frame(rownames(superfull_table), superfull_table)
	colnames(with_row_names) <- c("genes", colnames(superfull_table))
	
	write.table(with_row_names, paste0(saveIn, "logcounts_all_the_selected_cells_version_barbara_elena_eleonora_only_SUBCLUSTERED_exp_DOPO_LOGNORMCOUNTS_", expN, "_", cellLab[cT], "_SLM_ZERO_TRE.csv"), row.names=F, quote=F, sep="\t")
		
	################################################# EXPORT THE NEW TABLE ELENA ASKED #################################################
	
####################	# exporting logcounts for Elena's violins BY CONDITION
####################	subCTRL <- cTSubSet[, which(cTSubSet$condition=="CTRL")]
####################	subEAE <- cTSubSet[, which(cTSubSet$condition=="EAE")]
####################	
####################	# append info to colnames
####################	coln_CTRL <- paste(colnames(subCTRL), subCTRL$condition, subCTRL$gender, subCTRL$pruned_fine, "cluster", subCTRL$slm_0.3, sep="_")
####################	with_row_names_CTRL <- cbind.data.frame(rownames(subCTRL), logcounts(subCTRL))
####################	colnames(with_row_names_CTRL) <- c("genes", coln_CTRL)
####################	
####################	coln_EAE <- paste(colnames(subEAE), subEAE$condition, subEAE$gender, subEAE$pruned_fine, "cluster", subEAE$slm_0.3, sep="_")
####################	with_row_names_EAE <- cbind.data.frame(rownames(subEAE), logcounts(subEAE))
####################	colnames(with_row_names_EAE) <- c("genes", coln_EAE)

####################	write.table(with_row_names_CTRL, paste0(saveIn, "logcounts_SUBCLUSTERING_CTRL_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_TRE.csv"), row.names=F, quote=F, sep="\t")
####################	write.table(with_row_names_EAE, paste0(saveIn, "logcounts_SUBCLUSTERING_EAE_exp_", expN, "_", cellLab[cT], "_SLM_ZERO_TRE.csv"), row.names=F, quote=F, sep="\t")
		
	# Cluster composition plots using main labels
	conD <- unique(cTSubSet$condition)

	dfReg <- data.frame(table(cTSubSet$slm_0.3, cTSubSet$pruned_fine, cTSubSet$condition))

	colnames(dfReg) <- c("Cluster", "CellType", "Condition", "Frequency")
	lims <- levels(factor(cTSubSet$pruned_fine))
	lims <- lims[order(table(cTSubSet$pruned_fine))]

	# plot!
	p <- plot_membership(dfReg, cond=conD[2], col=colorsPal[2])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[2], "_cluster_composition_SUBCLUSTERED_", cellLab[cT], "_exp_", expN, "_SLM_ZERO_TRE.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()

	p <- plot_membership(dfReg, cond=conD[1], col=colorsPal[6])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[1], "_cluster_composition_SUBCLUSTERED_", cellLab[cT], "_exp_", expN, "_SLM_ZERO_TRE.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()
	
	####################### TEST 1 IN EACH CLUSTER COMPARE STUFF #######################
	
	# loop through all the clusters containing a specific cellType
	for (cS in seq(1, length(clustS), by=1)) {
	
		# PERFORM CTRL VS EAE
	
		print(paste("cluster position ", clustS[cS]))
	
		# subset SCE to get only cS cluster cells
		aClust <- cTSubSet[, cTSubSet$slm_0.3==clustS[cS]]
		
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
				
				cat("finding markers\n")
				
				# FINDMARKERS for DEGS
				markersG <- findMarkers(aClust, groups=aClust$condition, pval.type="any", direction="any")
		
				cC <- "EAE"
				chosen <- markersG[[cC]]

				# RELAXING THE THRESHOLD
				selfMade <- chosen[chosen$FDR < 0.05, ]
				#selfMade <- chosen[chosen$FDR < 0.1, ]

				# and export data
				finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

				# sort by logfoldchange
				if (length(finalGenes$logFC) != 0) {
					sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

					# and export
					write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_SUBCLUSTERED_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_TRE.xlsx"), overwrite=T)
				}
			}
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
			cTSubSetTestTwo$whatCluster[cTSubSet$slm_0.3 %in% clustS[cS]] <- "thisOne"
			cTSubSetTestTwo$whatCluster[cTSubSet$slm_0.3 %in% clustS[-cS]] <- "allOther"

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
	#					write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ONE_VERSUS_ALL_OTHER_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_TRE.csv"), row.names=F, quote=F, sep="\t")
						write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ONE_VERSUS_ALL_OTHER_SUBCLUSTERED_", cellLab[cT], "_CELLS_exp_", expN, "by_cluster_", clustS[cS], "_SLM_ZERO_TRE.xlsx"), overwrite=T)
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
