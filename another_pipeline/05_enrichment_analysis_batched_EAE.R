#' ---
#' title: "Quality controls"
#' author: "gabrielet"
#' output: html_document
#' date: '`r format(Sys.Date(), "%d %B, %Y")`'
#' ---
#' 
#' ```{r setup, include=FALSE}
#' knitr::sharedGenes_chunk$set(
#'   tidy = TRUE,
#'   tidy.sharedGenes = list(width.cutoff = 120),
#'   message = FALSE,
#'   warning = FALSE,
#'   eval = TRUE
#' )
#' ```

library("org.Mm.eg.db")
library("clusterProfiler")
library("ggplot2")
library("RColorBrewer")
library("stringr")
library("ddpcr")
library("pathview")
library("DOSE")
library("ReactomePA")
library("enrichplot")
library("openxlsx")
library("UpSetR")

# import plotting functions
source("00_functions.R")

#set paths
rawPath <- "/home/gabriele/cbmc/scrnaseq/raw_data_new/EAE_batched/"

cellType <- "Neutrophils"

# select experiment name
expN <- "batched"
  
# then do everything
saveIn <- paste0(rawPath, "results/")
saveImg <- paste0(saveIn, "figures/analysis_", cellType, "/")
#saveKegg <- paste0(saveIn, "enrichment_", cellType, "/KEGG_plots/")
#savePView <- paste0(saveKegg, "/pathviews/")
saveGO <- paste0(saveIn, "enrichment_", cellType, "/GO_plots/")
saveReactome <- paste0(saveIn, "enrichment_", cellType, "/Reactome_plots/")
saveCProf <- paste0(saveIn, "enrichment_", cellType, "/clusterProfiler/")
saveCProfRanks <- paste0(saveCProf, "rankings/")
saveCProfGOBP <- paste0(saveCProf, "GO_bio_process/")
saveCProfGOMF <- paste0(saveCProf, "GO_mol_function/")
saveCProfRCTP <- paste0(saveCProf, "REACTOME_pathway/")
saveUpDown <- paste0(saveIn, "enrichment_", cellType, "/UpDown/")
getRanks <- paste0(saveIn, "enrichment_", cellType, "/ranks/")

# check if saveIn directory exist, otherwise create it
ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
ifelse(dir.exists(saveImg), TRUE, dir.create(saveImg, recursive=T))
#ifelse(dir.exists(saveKegg), TRUE, dir.create(saveKegg, recursive=F))
#ifelse(dir.exists(savePView), TRUE, dir.create(savePView, recursive=F))
ifelse(dir.exists(saveGO), TRUE, dir.create(saveGO, recursive=F))
ifelse(dir.exists(saveReactome), TRUE, dir.create(saveReactome, recursive=F))
ifelse(dir.exists(saveCProf), TRUE, dir.create(saveCProf, recursive=F))
ifelse(dir.exists(saveCProfRanks), TRUE, dir.create(saveCProfRanks, recursive=F))
ifelse(dir.exists(saveCProfGOBP), TRUE, dir.create(saveCProfGOBP, recursive=F))
ifelse(dir.exists(saveCProfGOMF), TRUE, dir.create(saveCProfGOMF, recursive=F))
ifelse(dir.exists(saveCProfRCTP), TRUE, dir.create(saveCProfRCTP, recursive=F))

# create list of pathways to pass to enrichmydata function, to properly store the generated csv files
# ORDER IS FUNDAMENTAL HERE
cProfPathsList <- c(saveCProfGOBP, saveCProfGOMF, saveCProfRCTP)

# set palette with 12 colours for colouring plots
colorsPal <- brewer.pal(12, "Paired")

# get list containing both up and down regulated genes, coming from findMarkers from one of the three optOne, optTwo, optThree functions
filesList <- list.files(saveUpDown)[grep("up_and_down.csv", list.files(saveUpDown))]

# create list of all the genes, to find those that are shared across all the clusters
sharedGenes <- list()

# loop through the gene lists
for (fName in filesList) {

	print(paste0("analysing ", fName))
	#load data
	originalData <- read.csv(paste0(saveUpDown, fName), header=T, sep="\t", stringsAsFactors=T)
	# get identifier for this experiment
	id <- strsplit(fName, "_up")[[1]][1]

	# check if originalData contains actual information
	if (nrow(originalData) != 0) {
		enrichMyData(originalData, cProfPathsList, id)
	} else {
		print("no genes to enrich")
	}
	# get cluster name for upset plotting
	fN <- strsplit(fName, "_exp")[[1]][1]
	# and store the genes
	sharedGenes[[fN]] <- originalData$gName
}

# upset plot with the intersection of the up_and_down genes

# are there shared genes?
intersection <- unique(unlist(sharedGenes))

# export the full list with shared genes
write.table(intersection, paste0(saveImg, "shared_up_and_down_genes.csv"), col.names=F, row.names=F, sep="\t")

# initialise label count
lbl <- 1

# to perform a comparison al least two lists are required
if (length(names(sharedGenes)) > 1) {

	# loop through the different sets of genes
	for (lsT in names(sharedGenes)) {

		# if at least an intersection is found
		if (length(intersection) != 0) {
		
			uP <- upset(fromList(sharedGenes), order.by="freq", nsets=length(names(sharedGenes)), decreasing=T, sets=names(sharedGenes), keep.order=T, text.scale=1.5, sets.x.label="Set size", mainbar.y.label="Shared proteins")	
			
			# update label
			lbl <- lbl + 1
			png(paste0(saveImg, "04_0", lbl, "_upset_", lsT), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
				print(uP)
			dev.off()

			# plotting heatmap
			geneS <- vector()
			clustS <- vector()
			belonG <- vector()
			# loop through the cols
			for (nC in seq(1, ncol(uP$New_data), by=1)) {
				# loop through the rows
				for (nR in seq(1, nrow(uP$New_data), by=1)) {
					# get cluster name
					clustS <- append(clustS, names(uP$New_data)[nC])
					# get gene name
					geneS <- append(geneS, rownames(uP$New_data)[nR])
					# get belonging 
					bG <- uP$New_data[nR, nC]
					if (bG == 0) {
						belonG <- append(belonG, "not found")
					} else {
						belonG <- append(belonG, "found")
					}
				}
			}
			# put everything together
			dataG <- cbind.data.frame(geneS, clustS, belonG)

			# and print 
			p <- ggplot(dataG, aes(clustS, geneS, fill=factor(belonG))) + 
				geom_tile(color="black") +
				xlab("cluster") +
				ylab("gene name") +
				guides(fill=guide_legend(title="Gene in cluster is: "), position="bottom") +
				scale_fill_manual(values=c("not found"=colorsPal[1], "found"=colorsPal[6]))

			# update label
			lbl <- lbl + 1
			png(paste0(saveImg, "04_0", lbl, "_heatmap_", lsT), type="cairo", units="in", width=6, height=32, pointsize=12, res=300)
				plot(p)
			dev.off()
		}
	}
} else {
	print("the number of list is not sufficient to perform comparisons!")
}

# get list containing both up and down regulated genes
filesList <- list.files(getRanks)

# set threshold for fold change
fcThresh <- 1.5

# loop through the gene lists
for (fName in filesList) {

	print(paste0("analysing ", fName))
	
	#load data
	originalData <- read.csv(paste0(getRanks, fName), header=T, sep="\t", stringsAsFactors=T)
	
	# check if we are working with a cluster
	if (length(grep("cluster", fName)) != 0) {
		# get cluster name
		cS <- strsplit(strsplit(fName, "cluster_")[[1]][2], "\\.")[[1]][1]
		# check if originalData contains actual information
		if (nrow(originalData) != 0) {
#			lbl <- keggPlotting(originalData, saveKegg, saveCProfRanks, savePView, lbl, fName, cS)
			lbl <- goPlotting(originalData, saveGO, saveCProfRanks, lbl, fName, cS)
			lbl <- reactomePlotting(originalData, saveReactome, saveCProfRanks, fcThresh, lbl, fName, cS)
		} else {
			print("no genes to enrich")
		}
	# or with all the cells
	} else if (length(grep("cluster", fName)) == 0) {
		# set cluster name
		cS <- "ALL_CELLS"
		# check if originalData contains actual information
		if (nrow(originalData) != 0) {
#			lbl <- keggPlotting(originalData, saveKegg, saveCProfRanks, savePView, lbl, fName, cS)
			lbl <- goPlotting(originalData, saveGO, saveCProfRanks, lbl, fName, cS)
			lbl <- reactomePlotting(originalData, saveReactome, saveCProfRanks, fcThresh, lbl, fName, cS)
		} else {
			print("no genes to enrich")
		}
	}
}

# shut all the residual open graphics down
graphics.off()

