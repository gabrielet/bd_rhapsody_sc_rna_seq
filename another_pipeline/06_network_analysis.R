#' ---
#' title: "Quality controls"
#' author: "gabrielet"
#' output: html_document
#' date: '`r format(Sys.Date(), "%d %B, %Y")`'
#' ---
#' 
#' ```{r setup, include=FALSE}
#' knitr::opts_chunk$set(
#'   tidy = TRUE,
#'   tidy.opts = list(width.cutoff = 120),
#'   message = FALSE,
#'   warning = FALSE,
#'   eval = TRUE
#' )
#' ```

library("igraph")
library("ggplot2")
library("RColorBrewer")
library("STRINGdb")

# import plotting functions
source("00_functions.R")

# cell types
cellType <- c("Neutrophils", "Tgd", "Bcells", "CD8", "CD4")  
  
# set AD available rawPaths and experiments
rPs <- c("/home/gabriele/cbmc/scrnaseq/raw_data_new/572_20/", "/home/gabriele/cbmc/scrnaseq/raw_data_new/566_20/")
exS <- c("572", "566")

# loop through AD experiments
for (exP in c(1, 2)) {
	# and through cellTypes
	for (cType in cellType){

		# set rawPath and expN
		rawPath <- rPs[exP]
		expN <- exS[exP]	
	
		# then do everything
		saveIn <- paste0(rawPath, "results/")
		saveNets <- paste0(saveIn, "enrichment_", cType, "/networks/")
		saveUpDown <- paste0(saveIn, "enrichment_", cType, "/UpDown/")

		# check if saveIn directory exist, otherwise create it
		ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
		ifelse(dir.exists(saveNets), TRUE, dir.create(saveNets, recursive=F))

		# set palette with 12 colours for colouring plots
		colorsPal <- brewer.pal(12, "Paired")

		# get list containing both up and down regulated genes but only on the full dataset, not by single cluster
		filesList <- list.files(saveUpDown)[grep("up_and_down.csv", list.files(saveUpDown))][grep("ONLY", list.files(saveUpDown)[grep("up_and_down.csv", list.files(saveUpDown))])]

		# load STRING to SYMBOL mapping file
		# this file is obtained by the info file downloaded from STRING. from this file extract the first two columns to avoid importing issues
		mappingNames <- read.table("/home/gabriele/cbmc/scrnaseq/scrnaseq/10090.protein.info_cols_reduced.v11.0.txt", header=T)

		# if at least a files was found
		if (length(filesList) > 0) {
			# loop through the gene lists
			for (fName in filesList) {

				fName <- filesList

				#load data
				originalData <- read.csv(paste0(saveUpDown, fName), header=T, sep="\t", stringsAsFactors=T)
				
				# 10090 refers to Mus musculus
				stringDB <- STRINGdb$new(version="11", species=10090, score_threshold=400, input_directory=getwd())

				# mapping the genes into the STRING database to build the networks
				mappedGenes <- stringDB$map(originalData, "gName", removeUnmappedRows=F)

				# found and not found genes
				notFound <- mappedGenes[is.na(mappedGenes$STRING_id), ]
				foundGenes <- mappedGenes[!is.na(mappedGenes$STRING_id), ]
				# export not found genes!
			#	write.table(notFound$gName, paste0(saveIn, "notFoundSTRINGGenes.csv"), row.names=FALSE, col.names=F, quote=FALSE)

				# get neighbours for each found gene
				fNeighs <- stringDB$get_neighbors(foundGenes$STRING_id)
				fDF <- as.data.frame(fNeighs)
				mappedFN <- stringDB$map(fDF, "fNeighs", removeUnmappedRows=F)

				# build the network, using the geneNames (using goGroup or keggGroup geneNames should be equivalent!)
				originalNet <- simplify(stringDB$get_subnetwork(foundGenes$STRING_id), remove.multiple=TRUE, remove.loops=TRUE)
				V(originalNet)$symbol <- mappingNames$preferred_name[mappingNames$protein_external_id %in% V(originalNet)$name]
					
				# build the networks using geneNames and their first neighbours
				firstNet <- simplify(stringDB$get_subnetwork(fNeighs), remove.multiple=TRUE, remove.loops=TRUE)
				V(firstNet)$symbol <- mappingNames$preferred_name[mappingNames$protein_external_id %in% V(firstNet)$name]

				# computing Degree on both networks
				deg <- list()
				deg[["original"]] <- degree(originalNet, v=V(originalNet), mode="all", normalized=F)
				deg[["firstNeig"]] <- degree(firstNet, v=V(firstNet), mode="all", normalized=F)
				
				# computing Betweenness on both networks
				betw <- list()
				betw[["original"]] <- betweenness(originalNet, v=V(originalNet), directed=F, weights=NULL, normalized=F)
				betw[["firstNeig"]] <- betweenness(firstNet, v=V(firstNet), directed=F, weights=NULL, normalized=F)

				# plot networks, without first neighbours

				png(paste0(saveNets , "04_01_degree_original.png"), type="cairo", units="in", width=9, height=12, pointsize=12, res=300)
#					plot(originalNet, vertex.size=deg[["original"]], vertex.label=V(originalNet)$symbol, vertex.frame.color="white", vertex.label.color="black", vertex.color=brewer.pal(12, "Paired")[5])
					plot(originalNet, vertex.size=10, vertex.label=V(originalNet)$symbol, vertex.frame.color="white", vertex.label.color="black", vertex.color=brewer.pal(12, "Paired")[5])
				dev.off()
				
				png(paste0(saveNets, "04_0_betweenness_original.png"), type="cairo", units="in", width=9, height=12, pointsize=12, res=300)
#					plot(originalNet, vertex.size=betw[["original"]], vertex.label=V(originalNet)$symbol, vertex.frame.color="white", vertex.label.color="black", vertex.color=brewer.pal(12, "Paired")[5])
					plot(originalNet, vertex.size=10, vertex.label=V(originalNet)$symbol, vertex.frame.color="white", vertex.label.color="black", vertex.color=brewer.pal(12, "Paired")[5])
					dev.off()
				
				# plot first neighbours networks

			#	png(paste0(saveNets, "06_03_degree_firstNeighs.png"), type="cairo", units="in", width=9, height=24, pointsize=24, res=300)
			#		plot(firstNet, vertex.size=deg[["firstNeigh"]], vertex.label=V(firstNet)$symbol, vertex.frame.color="white", vertex.label.color="black", vertex.color=brewer.pal(12, "Paired")[5])
			#	dev.off()
			#	
			#	png(paste0(saveNets, "06_04_betweenness_firstNeighs.png"), type="cairo", units="in", width=9, height=24, pointsize=24, res=300)
			#		plot(firstNet, vertex.size=betw[["firstNeigh"]], vertex.label=V(firstNet)$symbol, vertex.frame.color="white", vertex.label.color="black", vertex.color=brewer.pal(12, "Paired")[5])
			#	dev.off()

			}
		}
	}
}

#####################################################################
#########################  now analyse EAE  #########################
#####################################################################

# cell types
cellType <- c("Neutrophils")
  
# set EAE available rawPaths and experiments
rPs <- c("/home/gabriele/cbmc/scrnaseq/raw_data_new/588_20/", "/home/gabriele/cbmc/scrnaseq/raw_data_new/589_20/")
exS <- c("588", "589")

# loop through AD experiments
for (exP in c(1, 2)) {
	for (cType in cellType){

		# set rawPath and expN
		rawPath <- rPs[exP]
		expN <- exS[exP]	
	
		# then do everything
		saveIn <- paste0(rawPath, "results/")
		saveNets <- paste0(saveIn, "enrichment_", cType, "/networks/")
		saveUpDown <- paste0(saveIn, "enrichment_", cType, "/UpDown/")

		# check if saveIn directory exist, otherwise create it
		ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
		ifelse(dir.exists(saveNets), TRUE, dir.create(saveNets, recursive=F))

		# set palette with 12 colours for colouring plots
		colorsPal <- brewer.pal(12, "Paired")

		# get list containing both up and down regulated genes but only on the full dataset, not by single cluster
		filesList <- list.files(saveUpDown)[grep("up_and_down.csv", list.files(saveUpDown))][grep("ONLY", list.files(saveUpDown)[grep("up_and_down.csv", list.files(saveUpDown))])]

		# load STRING to SYMBOL mapping file
		# this file is obtained by the info file downloaded from STRING. from this file extract the first two columns to avoid importing issues
		mappingNames <- read.table("/home/gabriele/cbmc/scrnaseq/scrnaseq/10090.protein.info_cols_reduced.v11.0.txt", header=T)

		# if at least a files was found
		if (length(filesList) > 0) {

			# loop through the gene lists
			for (fName in filesList) {

				fName <- filesList

				#load data
				originalData <- read.csv(paste0(saveUpDown, fName), header=T, sep="\t", stringsAsFactors=T)
				
				# 10090 refers to Mus musculus
				stringDB <- STRINGdb$new(version="11", species=10090, score_threshold=400, input_directory=getwd())

				# mapping the genes into the STRING database to build the networks
				mappedGenes <- stringDB$map(originalData, "gName", removeUnmappedRows=F)

				# found and not found genes
				notFound <- mappedGenes[is.na(mappedGenes$STRING_id), ]
				foundGenes <- mappedGenes[!is.na(mappedGenes$STRING_id), ]
				# export not found genes!
			#	write.table(notFound$gName, paste0(saveIn, "notFoundSTRINGGenes.csv"), row.names=FALSE, col.names=F, quote=FALSE)

				# get neighbours for each found gene
				fNeighs <- stringDB$get_neighbors(foundGenes$STRING_id)
				fDF <- as.data.frame(fNeighs)
				mappedFN <- stringDB$map(fDF, "fNeighs", removeUnmappedRows=F)

				# build the network, using the geneNames (using goGroup or keggGroup geneNames should be equivalent!)
				originalNet <- simplify(stringDB$get_subnetwork(foundGenes$STRING_id), remove.multiple=TRUE, remove.loops=TRUE)
				V(originalNet)$symbol <- mappingNames$preferred_name[mappingNames$protein_external_id %in% V(originalNet)$name]
					
				# build the networks using geneNames and their first neighbours
				firstNet <- simplify(stringDB$get_subnetwork(fNeighs), remove.multiple=TRUE, remove.loops=TRUE)
				V(firstNet)$symbol <- mappingNames$preferred_name[mappingNames$protein_external_id %in% V(firstNet)$name]

				# computing Degree on both networks
				deg <- list()
				deg[["original"]] <- degree(originalNet, v=V(originalNet), mode="all", normalized=F)
				deg[["firstNeig"]] <- degree(firstNet, v=V(firstNet), mode="all", normalized=F)
				
				# computing Betweenness on both networks
				betw <- list()
				betw[["original"]] <- betweenness(originalNet, v=V(originalNet), directed=F, weights=NULL, normalized=F)
				betw[["firstNeig"]] <- betweenness(firstNet, v=V(firstNet), directed=F, weights=NULL, normalized=F)

				# plot networks, without first neighbours

				png(paste0(saveNets , "04_01_degree_original.png"), type="cairo", units="in", width=9, height=12, pointsize=12, res=300)
					plot(originalNet, vertex.size=deg[["original"]], vertex.label=V(originalNet)$symbol, vertex.frame.color="white", vertex.label.color="black", vertex.color=brewer.pal(12, "Paired")[5])
				dev.off()
				
				png(paste0(saveNets, "04_0_betweenness_original.png"), type="cairo", units="in", width=9, height=12, pointsize=12, res=300)
					plot(originalNet, vertex.size=betw[["original"]], vertex.label=V(originalNet)$symbol, vertex.frame.color="white", vertex.label.color="black", vertex.color=brewer.pal(12, "Paired")[5])
				dev.off()
				
				# plot first neighbours networks

			#	png(paste0(saveNets, "06_03_degree_firstNeighs.png"), type="cairo", units="in", width=9, height=24, pointsize=24, res=300)
			#		plot(firstNet, vertex.size=deg[["firstNeigh"]], vertex.label=V(firstNet)$symbol, vertex.frame.color="white", vertex.label.color="black", vertex.color=brewer.pal(12, "Paired")[5])
			#	dev.off()
			#	
			#	png(paste0(saveNets, "06_04_betweenness_firstNeighs.png"), type="cairo", units="in", width=9, height=24, pointsize=24, res=300)
			#		plot(firstNet, vertex.size=betw[["firstNeigh"]], vertex.label=V(firstNet)$symbol, vertex.frame.color="white", vertex.label.color="black", vertex.color=brewer.pal(12, "Paired")[5])
			#	dev.off()

			}
		}
	}
}


































	################source <- vector()
	################target <- vector()
	################for (l in seq(1, length(E(firstNet)), by=1)) {
	################	source <- append(source, get.edgelist(firstNet)[l, ][1])
	################	target <- append(target, get.edgelist(firstNet)[l, ][2])
	################}


















