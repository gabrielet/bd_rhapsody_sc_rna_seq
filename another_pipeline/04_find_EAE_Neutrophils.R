#' ---
#' title: "Analysing Neutrophils cells in EAE exps"
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

# import plotting functions
source("00_functions.R")

expS <- c("588", "589")

# initialise label
lbl <- 1

# loop throgh the experiments
for (exP in expS) {
	
	# load the correct one
	if (exP == "588") {
		# set paths
		rawPath <- "/home/gabriele/cbmc/scrnaseq/raw_data_new/588_20/" ; expN <- "588"; sex <- "onset" ; kS <- 10 ; # EAE onset
	} else if (exP == "589") {
		#set paths
		rawPath <- "/home/gabriele/cbmc/scrnaseq/raw_data_new/589_20/" ; expN <- "589"; sex <- "chronic" ; kS <- 15 ; # EAE chronic
	} else {
		print("exp not found")
	}

	cellType <- "Neutrophils"

	print("")
	print(paste0("analysing ", cellType, " for exp ", expN, " ", sex, " mice, reduced with PCA"))
	print("")
	
	# then do everything
	saveIn <- paste0(rawPath, "results/")
	saveImg <- paste0(saveIn, "figures/analysis_", cellType, "/")
	saveUpDown <- paste0(saveIn, "enrichment_", cellType, "/UpDown/")
	saveRanks <- paste0(saveIn, "enrichment_", cellType, "/ranks/")
	saveAB <- paste0(saveImg, "ABseq/")

	# check if saveIn directory exist, otherwise create it
	ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
	ifelse(dir.exists(saveImg), TRUE, dir.create(saveImg, recursive=T))
	ifelse(dir.exists(saveUpDown), TRUE, dir.create(saveUpDown, recursive=T))
	ifelse(dir.exists(saveRanks), TRUE, dir.create(saveRanks, recursive=T))
	ifelse(dir.exists(saveAB), TRUE, dir.create(saveAB, recursive=T))

	# get antibody information
	antiB <- read.csv("/home/gabriele/cbmc/scrnaseq/scrnaseq/working_pipe/sc_analysis/abseqlist.txt", sep="\t", header=T)
	antiSce <- readRDS(paste0(saveIn, paste0("abseq_data_", expN, ".Rds")))
	# merge the two condition
	mergedAB <- cbind(antiSce[["EAE"]], antiSce[["CTRL"]])
	
	# set palette with 12 colours for colouring plots
	colorsPal <- brewer.pal(12, "Paired")

	# load PCA data, with proper kS
	clusteredSce <- readRDS(paste0(saveIn, paste0("kmeans_K_", kSExp, "_clustered_data_", expN, ".Rds")))

	# plotting heatmap with cell type
	tab_fine_to_main <- table(Assigned=clusteredSce$fine_to_main, Cluster=clusteredSce$clusterS)

	# get the table ready to be exported
	eTab <- as.data.frame.matrix(tab_fine_to_main)
	eTabRN <- cbind.data.frame(cellType=rownames(eTab), eTab)
	eTabCN <- rbind.data.frame(c("cluster", colnames(eTab)), eTabRN)

	# export the table
	write.table(eTabCN, paste0(saveImg, "cluster_content_exp_", expN, "_", cellType, ".csv"), quote=F, row.names=F, col.names=F, sep="\t")
	write.xlsx(eTabCN, paste0(saveImg, "cluster_content_exp_", expN, "_", cellType, ".xlsx"), colNames=F)
	
	tab_full_labels <- table(Assigned=clusteredSce$pruned_fine, Cluster=clusteredSce$clusterS)

	# get the full table ready to be exported
	etabF <- as.data.frame.matrix(tab_full_labels)
	eTabrnF <- cbind.data.frame(cellType=rownames(etabF), etabF)
	eTabcnF <- rbind.data.frame(c("cluster", colnames(etabF)), eTabrnF)
	
	# export the full table
	write.table(eTabcnF, paste0(saveImg, "cluster_content_full_labels_exp_", expN, "_", cellType, ".csv"), quote=F, row.names=F, col.names=F, sep="\t")
	write.xlsx(eTabcnF, paste0(saveImg, "cluster_content_full_labels_exp_", expN, "_", cellType, ".xlsx"), colNames=F)

	# get data set conditions, i.e. CTRL and EAE
	conD <- unique(clusteredSce$condition)

	# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_clusters_matrix.png"), type="cairo", units="in", width=9, height=12, pointsize=12, res=300)
	pheatmap(log2(tab_fine_to_main+10), color=colorRampPalette(c(colorsPal[1], colorsPal[6]))(101))
	dev.off()

	# plotting tSNE using cluster colours and cell label colours, side by side
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_clusters_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	grid.arrange(
		scater::plotTSNE(clusteredSce, colour_by="clusterS", text_by="clusterS") + theme(legend.position="none") + labs(title="By cluster"),
		scater::plotReducedDim(clusteredSce, "TSNE", colour_by="fine_to_main", text_by="fine_to_main") + theme(legend.position="none") + labs(title="By cell type"),
		ncol=2)
	dev.off()

	# MARKER GENES DETECTION

	# three different approaches explained:
	# using ImmGen we are able to label raw, not normalised, cells to obtain their predicted cell type. Then we perform two different normalisations
	# that are Zinbwave and sizeFactors. At this point two possibilities are open. The actual clustering is performed using the data normalised by
	# Zinbwave while sizeFactors are not considered.
	# first option: using only cellType cells that fall into the clusters containing cellType cells
	# second option: using all cells that fall into the clusters containing cellType cells
	# third option: using all cells that were labelled by ImmGen

	# AB info is also available to validate the cells.

	# search for cellType
	cTPos <- grep(cellType, clusteredSce$fine_to_main)
	
	print(paste0(cellType, " are found in ", length(unique(clusteredSce$clusterS[cTPos])), " different clusters"))

	fR <- table(clusteredSce$clusterS[cTPos])
	clustS <- names(fR[which(fR!=0)])

	print(fR[which(fR!=0)])

	# assign "unknown" labels to all cells
	clusteredSce$plot_labs <- rep("-", length(clusteredSce$pruned_fine))
	# find cellTypes cells
	clusteredSce$plot_labs[cTPos] <- cellType

	# plotting tSNE using cluster colours and cell label colours, side by side
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_clusters_only_", cellType, ".png"), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
	grid.arrange(
		scater::plotTSNE(clusteredSce, colour_by="clusterS", text_by="clusterS") + theme(legend.position="none") + labs(title="By cluster"),
		scater::plotReducedDim(clusteredSce, "TSNE", colour_by="plot_labs") + theme(legend.position="bottom") + labs(title="By cell type"),
		ncol=2)
	dev.off()

	# Cluster composition plots using main labels

	df <- data.frame(table(clusteredSce$clusterS, clusteredSce$fine_to_main, clusteredSce$condition))

	colnames(df) <- c("Cluster", "CellType", "Condition", "Frequency")
	lims <- levels(factor(clusteredSce$fine_to_main))
	lims <- lims[order(table(clusteredSce$fine_to_main))]

	# plot!
	p <- plot_membership(df, cond=conD[2], col=colorsPal[2])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[2], "_cluster_composition.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()

	p <- plot_membership(df, cond=conD[1], col=colorsPal[6])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[1], "_cluster_composition.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()

	######################################################################################################################################################################################################################################################################################################################################################################################################################

	# declare some variables to store the markers of each cluster for each option
	optOneGenes <- list()
	optTwoGenes <- list()
	optThreeGenes <- list()

	# set threshold for FRD
	thresholD <- 0.05

	# set a personal threshold determining the minimum number of cells to compare
	personalT <- 20

	# loop through all the clusters containing cellType
	for (cS in clustS) {

		# OPTION ONE: findMarkers ON ALL THE IMMGEN LABELLED CELLS OF A SPECIFIC CELLTYPE IN EACH CLUSTER #

		# subset SCE to get only cS cluster cells
		optOne <- clusteredSce[, clusteredSce$clusterS==cS]

		print("#################################################################################")
		print(paste0("analysing cluster ", cS, " cells content: "))
		
		print(table(optOne$fine_to_main))
		
		print(paste0("######## SEARCHING FOR ", cellType, " IN CLUSTER ", cS, " ########"))

		# now find cellType only in cS
		cTPos <- grep(cellType, optOne$fine_to_main)

		# and subset the original single cell experiment using the name of the cellType cells.
		cTSubSet <- optOne[, cTPos]

		if (ncol(cTSubSet) > 0) {
			# run option one on cTSubSet
			optOneGenes <- optionOne(cTSubSet, expN, saveUpDown, personalT, conD, cellType, thresholD)
			# check if everything went fine
			if (length(optOneGenes)==0) {
				print("no genes, no enrichment!")
			} else {
				# now perform GSEA on cTSubSet rank file names
				lbl <- lbl + 1
				rN <- paste0(saveImg, "04_0", lbl, "_", conD[1])
				# run enrichment on cTSubset
				ranksVals <- gseaMyData(cTSubSet, conD, rN, saveGSEA, saveGes, savePages, expN, paste0(cellType, "_cluster_", cS))
				# export rankings to perform later enrichment
				rTab <- cbind.data.frame(gene=names(ranksVals), foldChange=ranksVals)
				write.table(rTab, paste0(saveRanks, "rankings_genes_for_", cellType, "_cluster_", cS, ".csv"), col.names=T, row.names=F, sep="\t")
			}
		} else {
			print(paste0("not enough cells to compare the two conditions for ", cellType, " in cluster ", cS))
		}

		#############################################################################################

		# OPTION TWO: ALL THE CELLS IN THE CLUSTERS CONTAINING cellType CELLS

		print(paste0("######## SEARCHING FOR ALL THE CELLS IN CLUSTER ", cS, " ########"))

		# get cells in cluster cS
		cTSubSet <- clusteredSce[, clusteredSce$clusterS==cS]
		cTSubSet$fine_to_main[is.na(cTSubSet$fine_to_main)] <- "unknown"

		if (ncol(cTSubSet) > 0) {
			# run option two on cTSubSet
			optTwoGenes <- optionTwo(cTSubSet, expN, saveAB, saveUpDown, personalT, conD, cellType, mergedAB, thresholD)
			# no GSEA here, as decided with Elena
		} else {
			print(paste0("not enough cells to compare the two conditions for cluster ", cS))
		}
	}

	#############################################################################################

	# OPTION THREE: ONLY cellType CELLS IN THE FULL DATASET

	print(paste0("######## SEARCHING FOR ", cellType, " IN THE WHOLE DATASET ########"))

	# search for cellType 
	cTPos <- grep(cellType, clusteredSce$fine_to_main)

	# and subset the original single cell experiment using the name of the cellType cells.
	cTSubSet <- clusteredSce[, cTPos]

	if (ncol(cTSubSet) > 0) {
		# run option three on cTSubSet		
		optThreeGenes <- optionThree(cTSubSet, expN, saveUpDown, personalT, conD, cellType, thresholD)

		# check if everything went fine
		if (length(optThreeGenes)==0) {
			print("no genes, no enrichment!")
		} else {
			# now perform GSEA on cTSubSet rank file names
			lbl <- lbl + 1
			rN <- paste0(saveImg, "04_0", lbl, "_", conD[1])
			# run enrichment on cTSubset
			ranksVals <- gseaMyData(cTSubSet, conD, rN, saveGSEA, saveGes, savePages, expN, paste0("ALL_", cellType, "_CELLS"))
			# export rankings to perform later enrichment
			rTab <- cbind.data.frame(gene=names(ranksVals), foldChange=ranksVals)
			write.table(rTab, paste0(saveRanks, "rankings_genes_for_ALL_", cellType, "_CELLS.csv"), col.names=T, row.names=F, sep="\t")
		}
	} else {
		print(paste0("not enough cells to compare the two conditions"))
	}

	print("                                                 ")
	print("                                                 ")
	print("_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_")
	print("                                                 ")
	print("@@@@@@@@@              END              @@@@@@@@@")
	print("                                                 ")
	print("_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_")
}
