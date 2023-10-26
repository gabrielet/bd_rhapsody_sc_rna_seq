#' ---
#' title: "A collection of functions that are used to plot and enrich"
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

# function that plots the histograms and heatmap together

plot_membership <- function(df, cond, col) {

#	df <- df; cond <- testingExp[2]; col <- colorsPal[2];
	
	df <- df[df$Condition == cond, c(1:2, ncol(df))]
	p21 <- ggplot(df,aes(x=Cluster, y=CellType)) +
		geom_tile(width=0.9, height=0.9, aes(alpha=log(Frequency+1)), fill=col) +
		# coord_fixed() +
		theme_minimal() +
		geom_text(aes(label=ifelse(Frequency>0, Frequency, ""), alpha=log(Frequency+1)), color="white") +
		theme(legend.position="none") +
		scale_y_discrete(limits=lims) +
		scale_fill_manual(values=col)
	df_tot_CellType <- ddply(df, .variables=~ CellType, .fun=function(x) sum(x$Frequency))
	myT <- max(df_tot_CellType$V1)
	colnames(df_tot_CellType) <- c("CellType", "Frequency")
	# computing values as percentages, to keep the histograms away from too big values (that may happen)
	df_tot_CellType$Frequency <- round(df_tot_CellType$Frequency/sum(df_tot_CellType$Frequency) * 100, digits = 1)
	p22 <- ggplot(df_tot_CellType,aes(x=CellType, y=Frequency)) +
		geom_col(aes(alpha=log(Frequency+1)), fill=col) +
		theme_minimal() +
		coord_flip() +
		geom_text(aes(y=Frequency + 30, label=ifelse(Frequency>0, Frequency, ""), alpha=log(Frequency+1)), hjust="left") +
		theme(legend.position="none",
			axis.title=element_blank(),
			axis.text=element_blank(), panel.grid=element_blank()) +
		scale_x_discrete(limits=lims) +
		scale_y_continuous(limits=c(0, myT)) +
		scale_fill_manual(values=col)
	df_tot_Cluster <- ddply(df, .variables=~ Cluster, .fun=function(x) sum(x$Frequency))
	myC <- max(df_tot_Cluster$V1)
	colnames(df_tot_Cluster) <- c("Cluster","Frequency")
	# computing values as percentages, to keep the histograms away from too big values (that may happen)
	df_tot_Cluster$Frequency <- round(df_tot_Cluster$Frequency/sum(df_tot_Cluster$Frequency) * 100, digits = 1)
	p11 <- ggplot(df_tot_Cluster,aes(x=Cluster, y=Frequency)) +
		geom_col(aes(alpha=log(Frequency+1)), fill=col) +
		theme_minimal() +
		geom_text(aes(y=Frequency + 150, label=ifelse(Frequency>0, Frequency, ""), alpha=log(Frequency+1))) +
		theme(legend.position="none",
			axis.title=element_blank(),
			axis.text=element_blank(), panel.grid=element_blank()) +
		scale_y_continuous(limits=c(0, myC)) +
		scale_fill_manual(values=col)
	empty <- ggplot() +
		geom_point(aes(1,1), colour="white") +
		theme(axis.ticks=element_blank(),
			panel.background=element_blank(),
			axis.text.x=element_blank(),
			axis.text.y=element_blank(),
			axis.title.x=element_blank(),
			axis.title.y=element_blank())
	plot_grid(p11, empty, p21, p22, rel_heights=c(0.2,1), rel_widths=c(1,0.4), ncol=2,nrow=2, align="hv",axis="lrtb")
}

# plotting upset plots

# using this version: https://github.com/hms-dbmi/UpSetR/issues/85 of fromList() function

fromList <- function (input) {
	# Same as original fromList()...
	elements <- unique(unlist(input))
	data <- unlist(lapply(input, function(x) {x <- as.vector(match(elements, x))}))
	data[is.na(data)] <- as.integer(0)
	data[data != 0] <- as.integer(1)
	data <- data.frame(matrix(data, ncol = length(input), byrow = F))
	data <- data[which(rowSums(data) != 0), ]
	names(data) <- names(input)
	# ...Except now it conserves its original value names!
	row.names(data) <- elements
	# finally, return results
	return(data)
}

# analysing clusters

################## OPTION ONE FOR NOT-BATCHED AD DATA

optionOne <- function(subSet, expN, saveGReg, personalT, conD, cTy, thresHD) {

#	subSet <- cTSubSet; saveGReg <- saveUpDown; conD <- testingExp

	cTSce <- list()

	# then get WT cellType
	cTSce[[conD[1]]] <- subSet[, which(subSet$condition==conD[1])]
	# and 3xTG
	cTSce[[conD[2]]] <- subSet[, which(subSet$condition==conD[2])]

	#print info about number of cells for each condition
	print(paste0("in cluster ", cS, " there are ", ncol(cTSce[[conD[1]]]), " ", cTy, " for ", conD[1]))
	print(paste0("and ", ncol(cTSce[[conD[2]]]), " ", cTy, " for ", conD[2]))

	# this if else check is useful here and not for the two other options since here all the cellType cells are considered
	if (ncol(cTSce[[conD[1]]]) > personalT && ncol(cTSce[[conD[2]]]) > personalT) {

		# then search for marker genes using pairwise by condition so that WT will be compared to 3x (and vice versa)
		# https://rdrr.io/bioc/scran/man/findMarkers.html
		# it is important to note that findMarkers uses logcounts as the default assay, to compute the differentially expressed genes.
		# however when we apply zinbwave, its normalised values are put into the logcounts assay
		# also, only 1000 highly variable genes are used
		
		markersG <- findMarkers(subSet, groups=subSet$condition, pval.type="any", direction="any")

		# http://bioconductor.org/books/release/OSCA/marker-detection.html#looking-for-any-differences
		# The default philosophy of findMarkers() is to identify a combination of marker genes that - together - uniquely define one cluster against the rest. 
		# To this end, we collect the top DE genes from each pairwise comparison involving a particular cluster to assemble a set of candidate markers for that cluster.
		# We will demonstrate on cluster 7; the relevant DataFrame contains log2-fold changes of expression in cluster 7 over each other cluster, along with several 
		# statistics obtained by combining p-values across the pairwise comparisons involving 7.
		# chosen <- "7"
		# interesting <- markers.pbmc[[chosen]]

		# hence, if cond[1] stands for "3xTG", we are obtaining the fold change for the comparison of 3xTG over WT

		# filter on false discovery rate considering conD[1] as the chosen condition
		cC <- conD[1]
		chosen <- markersG[[cC]]
		selfMade <- chosen[chosen$FDR < thresHD, ]

		print(paste0(nrow(selfMade), " statistical significant genes were found for all the ", cTy, " cells in cluster ", cS))
		print(paste0("for the comparison of ", cC, " over ", conD[2]))

		# and export data
		finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

		# order by fold change, which is the 6th column
		finalOrdered <- list()
		finalOrdered[["UpAndDown"]] <-  finalGenes[order(finalGenes[, 6], decreasing=T), ]

		# get up and down, separated
		finalOrdered[["Up"]] <- finalOrdered[["UpAndDown"]][which(finalOrdered[["UpAndDown"]][, 6]>=0), ]
		dTmp <- finalOrdered[["UpAndDown"]][which(finalOrdered[["UpAndDown"]][, 6]<0), ]
		finalOrdered[["Down"]] <- dTmp[order(dTmp[, 6], decreasing=F), ]

		write.table(finalOrdered[["Up"]], paste0(saveGReg, cTy, "_cluster_", cS, "_exp_", expN, "_UP_REG.csv"), quote=F, row.names=F, col.names=T, sep="\t")
		write.xlsx(finalOrdered[["Up"]], paste0(saveGReg, cTy, "_cluster_", cS, "_exp_", expN, "_UP_REG.xlsx"))

		write.table(finalOrdered[["Down"]], paste0(saveGReg, cTy, "_cluster_", cS, "_exp_", expN, "_DOWN_REG.csv"), quote=F, row.names=F, col.names=T, sep="\t")
		write.xlsx(finalOrdered[["Down"]], paste0(saveGReg, cTy, "_cluster_", cS, "_exp_", expN, "_DOWN_REG.xlsx"))

		write.table(finalOrdered[["UpAndDown"]], paste0(saveGReg, cTy, "_cluster_", cS, "_exp_", expN, "_up_and_down.csv"), quote=F, row.names=F, col.names=T, sep="\t")
		write.xlsx(finalOrdered[["UpAndDown"]], paste0(saveGReg, cTy, "_cluster_", cS, "_exp_", expN, "_up_and_down.xlsx"))

		# store the genes for this cluster in the list
		optOneGenes <- finalGenes$gName

	} else {
		# set return to NULL
		optOneGenes <- NULL
		print(paste0("not enough cells to compare the two conditions for ", cTy, " in cluster ", cS))
	}
	# finally, return the results
	return(optOneGenes)
}

################## OPTION TWO FOR NOT-BATCHED AD DATA

optionTwo <- function(subSet, expN, saveAbs=NULL, saveGReg, personalT, conD, cTy, mgdAB, thresHD) {

	cTSce <- list()

	# then get WT cellType
	cTSce[[conD[1]]] <- subSet[, which(subSet$condition==conD[1])]
	# and 3xTG
	cTSce[[conD[2]]] <- subSet[, which(subSet$condition==conD[2])]

	#print info about number of cells for each condition
	print(paste0("in cluster ", cS, " there are ", ncol(cTSce[[conD[1]]]), " cells for ", conD[1]))
	print(paste0("and ", ncol(cTSce[[conD[2]]]), " cells for ", conD[2]))

	if (ncol(cTSce[[conD[1]]]) > personalT && ncol(cTSce[[conD[2]]]) > personalT) {

		# then search for marker genes using pairwise by condition so that WT will be compared to 3x (and vice versa)
		# https://rdrr.io/bioc/scran/man/findMarkers.html
		# it is important to note that findMarkers uses logcounts as the default assay, to compute the differentially expressed genes.
		# however when we apply zinbwave, its normalised values are put into the logcounts assay
		# also, only 1000 highly variable genes are used
		
		markersG <- findMarkers(subSet, groups=subSet$condition, pval.type="any", direction="any")

		# http://bioconductor.org/books/release/OSCA/marker-detection.html#looking-for-any-differences
		# The default philosophy of findMarkers() is to identify a combination of marker genes that - together - uniquely define one cluster against the rest. 
		# To this end, we collect the top DE genes from each pairwise comparison involving a particular cluster to assemble a set of candidate markers for that cluster.
		# We will demonstrate on cluster 7; the relevant DataFrame contains log2-fold changes of expression in cluster 7 over each other cluster, along with several 
		# statistics obtained by combining p-values across the pairwise comparisons involving 7.
		# chosen <- "7"
		# interesting <- markers.pbmc[[chosen]]

		# hence, if cond[1] stands for "3xTG", we are obtaining the fold change for the comparison of 3xTG over WT

		# filter on false discovery rate considering conD[1] as the chosen condition
		cC <- conD[1]
		chosen <- markersG[[cC]]
		selfMade <- chosen[chosen$FDR < thresHD, ]

		print(paste0(nrow(selfMade), " statistical significant genes were found for cluster ", cS))
		print(paste0("for the comparison of ", cC, " over ", conD[2]))

		# and export data
		finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

		# order by fold change, which is the 6th column, which is the 6th column
		finalOrdered <- list()
		finalOrdered[["UpAndDown"]] <-  finalGenes[order(finalGenes[, 6], decreasing=T), ]

		# get up and down, separated
		finalOrdered[["Up"]] <- finalOrdered[["UpAndDown"]][which(finalOrdered[["UpAndDown"]][, 6]>=0), ]
		dTmp <- finalOrdered[["UpAndDown"]][which(finalOrdered[["UpAndDown"]][, 6]<0), ]
		finalOrdered[["Down"]] <- dTmp[order(dTmp[, 6], decreasing=F), ]
		
		# store the genes for this cluster in the list
		optTwoGenes <- finalGenes$gName
	} else {
		# set return to NULL
		optTwoGenes <- NULL
		print(paste0("not enough cells to compare the two conditions for cluster ", cS))
	}

	# if saveAbs is provided, then analyse ABSeq information
	if (!is.null(saveAbs)) {
	
		# since it is not possible to directly infer the cell type using the antibodies that tagged a cell, what it might be done is to retrieve all the cells in a specific cluster and create a plot, or more, that defines the kinds of antibodies that are tagging such cells. this can be done by each cell type in the cluster

		# hence get the cell composition from the ABseq perspective

		# get all cell types in the cluster
		cTTmp <- unique(subSet$fine_to_main)
		# remove NAs
		cTypes <- cTTmp[!is.na(cTTmp)]
		
		# summarise cell types
		suMM <- table(subSet$fine_to_main)

		# initialise a storage for the barplot information
		barTmp <- list()
		barTabs <- vector()
		# loop through the cell types
		for (cTp in cTypes) {

			# get the cell names of that specific kind of cells in the cluster under analysis
			cNames <- colnames(subSet[, which(subSet$fine_to_main==cTp)])

			# get the abseq info for those cells
			subAB <- mgdAB[, colnames(mgdAB) %in% cNames]

			#create the data storage
			barTmp[[cTp]] <- cbind.data.frame(gene=antiB$antibody, count=rep(0, length(antiB$antibody)))

			# loop through the retrieved cells to gather info for the barplot
			for (cll in seq(1, ncol(subAB), by=1)) {
				# reduce multi antibody and sum up their contribution
				tmpDf <- rbind.data.frame(barTmp[[cTp]], mgdAB$antiB[[cll]])
				barTmp[[cTp]] <- ddply(tmpDf, .(gene), function(x) c(count=sum(x$count)))
			}

			# create readable names
			tmpN <- sapply(strsplit(barTmp[[cTp]]$gene, "\\|"), "[[", 1)
			barTmp[[cTp]]$gene <- sapply(strsplit(tmpN, "\\:"), "[[", 1)
			# divide counts by the number of cells, to obain an averaged value
			barTmp[[cTp]]$avg <- barTmp[[cTp]]$count / ncol(subAB)
			
			# bind the new barTmp
			valCT <- suMM[which(names(suMM)==cTp)]
			barTabs <- rbind.data.frame(barTabs, cbind.data.frame(gene=barTmp[[cTp]]$gene, cellT=paste0(rep(cTp, nrow(barTmp[[cTp]])), "=", valCT), avgCount=barTmp[[cTp]]$avg))
		}

		# eventually print each barTmp
		p <- ggplot(barTabs, aes(fill=gene, y=avgCount, x=gene)) + 
			geom_bar(position="dodge", stat="identity") +
			scale_fill_manual(values=colorsPal) +
			ggtitle(paste0("abSeq average content for cluster ", cS)) +
			facet_wrap(~cellT) +
			theme(legend.position="none", axis.text.x=element_text(angle=45, vjust=0.5, hjust=1)) +
			xlab("")
			
		gg.gap(plot=p, segments=c(500, 550), ylim=c(0, 2000))

		# update label
		lbl <- lbl + 1
		png(paste0(saveAbs, "04_0", lbl, "_abSeq_average_content_for cluster_", cS, "_in_", cTy, ".png"), type="cairo", units="in", width=16, height=12, pointsize=12, res=300)
		plot(p)
		dev.off()
	} else {
		# if saveAbs is not provided, i.e. it is NULL, don't do anything
		print("abseq information is missing, not going to analyse it")
	}
	# finally, return the results
	return(optTwoGenes)
}

################## OPTION THREE FOR NOT-BATCHED AD DATA

optionThree <- function(subSet, expN, saveGReg, personalT, conD, cTy, thresHD) {

	cTSce <- list()

	# then get WT cellType
	cTSce[[conD[1]]] <- subSet[, which(subSet$condition==conD[1])]
	# and 3xTG
	cTSce[[conD[2]]] <- subSet[, which(subSet$condition==conD[2])]

	#print info about number of cells for each condition
	print(paste0("in the full dataset there are ", ncol(cTSce[[conD[1]]]), " ", cTy, " for ", conD[1]))
	print(paste0("and ", ncol(cTSce[[conD[2]]]), " ", cTy, " for ", conD[2]))

	if (ncol(cTSce[[conD[1]]]) > personalT && ncol(cTSce[[conD[2]]]) > personalT) {

		# then search for marker genes using pairwise by condition so that WT will be compared to 3x (and vice versa)
		# https://rdrr.io/bioc/scran/man/findMarkers.html
		# it is important to note that findMarkers uses logcounts as the default assay, to compute the differentially expressed genes.
		# however when we apply zinbwave, its normalised values are put into the logcounts assay
		# also, only 1000 highly variable genes are used
		
		markersG <- findMarkers(subSet, groups=subSet$condition, pval.type="any", direction="any")
		
		# http://bioconductor.org/books/release/OSCA/marker-detection.html#looking-for-any-differences
		# The default philosophy of findMarkers() is to identify a combination of marker genes that - together - uniquely define one cluster against the rest. 
		# To this end, we collect the top DE genes from each pairwise comparison involving a particular cluster to assemble a set of candidate markers for that cluster.
		# We will demonstrate on cluster 7; the relevant DataFrame contains log2-fold changes of expression in cluster 7 over each other cluster, along with several 
		# statistics obtained by combining p-values across the pairwise comparisons involving 7.
		# chosen <- "7"
		# interesting <- markers.pbmc[[chosen]]

		# hence, if cond[1] stands for "3xTG", we are obtaining the fold change for the comparison of 3xTG over WT

		# filter on false discovery rate considering conD[1] as the chosen condition
		cC <- conD[1]
		chosen <- markersG[[cC]]
		selfMade <- chosen[chosen$FDR < thresHD, ]

		print(paste0(nrow(selfMade), " statistical significant genes were found for all the ", cTy, " cells"))
		print(paste0("for the comparison of ", cC, " over ", conD[2]))

		# and export data
		finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

		# order by fold change, which is the 6th column
		finalOrdered <- list()
		finalOrdered[["UpAndDown"]] <-  finalGenes[order(finalGenes[, 6], decreasing=T), ]

		# get up and down, separated
		finalOrdered[["Up"]] <- finalOrdered[["UpAndDown"]][which(finalOrdered[["UpAndDown"]][, 6]>=0), ]
		dTmp <- finalOrdered[["UpAndDown"]][which(finalOrdered[["UpAndDown"]][, 6]<0), ]
		finalOrdered[["Down"]] <- dTmp[order(dTmp[, 6], decreasing=F), ]

		write.table(finalOrdered[["Up"]], paste0(saveGReg, cTy, "_ONLY_exp_", expN, "_UP_REG.csv"), quote=F, row.names=F, col.names=T, sep="\t")
		write.xlsx(finalOrdered[["Up"]], paste0(saveGReg, cTy, "_ONLY_exp_", expN, "_UP_REG.xlsx"))

		write.table(finalOrdered[["Down"]], paste0(saveGReg, cTy, "_ONLY_exp_", expN, "_DOWN_REG.csv"), quote=F, row.names=F, col.names=T, sep="\t")
		write.xlsx(finalOrdered[["Down"]], paste0(saveGReg, cTy, "_ONLY_exp_", expN, "_DOWN_REG.xlsx"))

		write.table(finalOrdered[["UpAndDown"]], paste0(saveGReg, cTy, "_ONLY_exp_", expN, "_up_and_down.csv"), quote=F, row.names=F, col.names=T, sep="\t")
		write.xlsx(finalOrdered[["UpAndDown"]], paste0(saveGReg, cTy, "_ONLY_exp_", expN, "_up_and_down.xlsx"))
		
		# store the genes for this cluster in the list
		optThreeGenes <- finalGenes$gName
	} else {
		# set return to NULL
		optThreeGenes <- NULL
		print(paste0("not enough ", cTy, " cells to compare the two conditions"))
	}
	return(optThreeGenes)
}

# perform enrichment using DESeq2
gseaMyData <- function(subSet, condS, saveRk, saveGTab, saveGES, savePages, expN, idGSEA, cTy) {

#	subSet <- cTSubSet; condS <- testingExp; saveRk <- rN; saveGTab <- saveGSEA; saveGES <- saveGes; savePages <- savePages; expN <- expN; idGSEA <- paste0(cellType, "_cluster_", cS);
#	subSet <- cTSubSet; condS <- conD; saveRk <- rN; saveGTab <- saveGSEA; saveGES <- saveGes; savePages <- savePages; expN <- expN; idGSEA <- paste0(cTy, "_cluster_", cS);
	# obtain a list of DEGs
	subSet$condition <- as.factor(subSet$condition)
	quiet(dds <- DESeqDataSetFromMatrix(as.matrix(assays(subSet)$counts), colData=colData(subSet), design=~condition), all=T)
	## Pre filter

	# We perform a minimal pre-filtering to keep only rows that have at least 10 reads total. Note that more strict filtering to increase power is automatically applied via independent filtering on the mean of normalized counts within the results function

	keep <- rowSums(counts(dds)) >= 10
	dds <- dds[keep, ]

	# After this first filter, we remove genes without an EntrezID
	rowData(dds)$Symbol <- rownames(dds)
	quiet(rowData(dds)$Entrezid <- mapIds(x=org.Mm.eg.db, keys=rowData(dds)$Symbol, column="ENTREZID", keytype="SYMBOL"), all=T)

	keep <- !is.na(rowData(dds)$Entrezid)
	dds <- dds[keep, ]

	# set "WT" or "CTRL" as the reference level
	# if AD
	if (expN == "566" | expN == "572") {
		dds$condition <- factor(dds$condition, levels=c("WT", "3xTG"))
	# else EAE
	} else if (expN == "588" | expN == "589") {
		dds$condition <- factor(dds$condition, levels=c("CTRL", "EAE"))
	} else if (expN == "batched") {
		dds$condition <- factor(dds$condition, levels=c("WT", "3xTG"))
	}
	
	## Differential Expression

	# We launch DESeq2
	quiet(ddsDesequed <- DESeq(dds), all=T)
	res <- results(ddsDesequed)

	# A positive value of Log-Fold Change indicates higher expression values in 3X_Alz compared to WT

	# The analysis is performed by:
	# - ranking all genes in the data set;
	# - identifying the rank positions of all members of the gene set in the ranked data set;
	# - calculating an enrichment score (ES) that represents the difference between the observed rankings and that which would be expected assuming a random rank distribution.

	## Create ranks
	
	if (all(rownames(res) == rowData(ddsDesequed)$Symbol)) {
		# Rank all genes based on their statistics
		rkVals <- res$stat
		# and name each gene
		names(rkVals) <- rowData(ddsDesequed)$Entrezid
	} else {
		print("problem with names")
	}
	
	# finally, return results
	return(rkVals)
}

enrichMyData <- function(oData, saveCProf, id) {

	#oData <- originalData; saveCProf <- saveCProf;

	geneList <- list()
	reactomeList <- list()
	goBP <- list()
	goMF <- list()
	
	lolliBP <- list()
	lolliMF <- list()

	validID <- list()
	comparisonMatGO <- list()

	#call the database from which we will get the info
	database <- org.Mm.eg.db

	# up and down together, only
	geneList[["UpAndDown"]] <- originalData
	
	# loop through the gene lists
	for (nm in names(geneList)) {
		#create a string array with the names of the genes
		geneNames <- as.character(geneList[[nm]]$gName)

		#getting ENTREZ IDs to perform the analysis https://davetang.org/muse/2013/12/16/bioconductor-annotation-packages/
		quiet(listOfGenes <- AnnotationDbi::select(database, keys=geneNames, columns=c("ENTREZID", "SYMBOL"), keytype="SYMBOL"), all=T)

		#found genes, with ENTREZ IDs
		validID[[nm]] <- listOfGenes[!is.na(listOfGenes$ENTREZID), ]

		# GO bio process
		quiet(goBP[[nm]] <- as.data.frame(enrichGO(validID[[nm]]$ENTREZID, database, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05)), all=T)
		#translate entrez in symbol before exporting the table
		symbl <- vector()
		for (gSet in goBP[[nm]]$geneID) {
			genes <- unlist(strsplit(gSet, "/"))
			symbl <- append(symbl, paste(validID[[nm]]$SYMBOL[validID[[nm]]$ENTREZ %in% genes], collapse="/"))
		}
		goBP[[nm]] <- cbind.data.frame(goBP[[nm]], symbol=symbl)
		# then export the table using the first slot of saveCProf paths list
		write.table(goBP[[nm]], paste0(saveCProf[1], "cProf_FindMarkers_GO_bio_process_", nm, "_", id, ".csv"), row.names=F, col.names=T, quote=F, sep="\t")
		write.xlsx(goBP[[nm]], paste0(saveCProf[1], "cProf_FindMarkers_GO_bio_process_", nm, "_", id, ".xlsx"))
		
		# GO mol function
		quiet(goMF[[nm]] <- as.data.frame(enrichGO(validID[[nm]]$ENTREZID, database, ont="MF", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05)), all=T)		
		
		#translate entrez in symbol before exporting the table
		symbl <- vector()
		for (gSet in goMF[[nm]]$geneID) {
			genes <- unlist(strsplit(gSet, "/"))
			symbl <- append(symbl, paste(validID[[nm]]$SYMBOL[validID[[nm]]$ENTREZ %in% genes], collapse="/"))
		}
		goMF[[nm]] <- cbind.data.frame(goMF[[nm]], symbol=symbl)
		# then export the table using the second slot of saveCProf paths list
		write.table(goMF[[nm]], paste0(saveCProf[2], "cProf_FindMarkers_GO_mol_function_", nm, "_", id, ".csv"), row.names=F, col.names=T, quote=F, sep="\t")
		write.xlsx(goMF[[nm]], paste0(saveCProf[2], "cProf_FindMarkers_GO_mol_function_", nm, "_", id, ".xlsx"))
		
		# REACTOME pathways
		quiet(reactomeList[[nm]] <- as.data.frame(enrichPathway(validID[[nm]]$ENTREZID, organism="mouse", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05)), all=T)
		#translate entrez into symbol befor exporting the table
		symbl <- vector()
		for (gSet in reactomeList[[nm]]$geneID) {
			genes <- unlist(strsplit(gSet, "/"))
			symbl <- append(symbl, paste(validID[[nm]]$SYMBOL[validID[[nm]]$ENTREZ %in% genes], collapse="/"))
		}
		reactomeList[[nm]] <- cbind.data.frame(reactomeList[[nm]], symbol=symbl)
		# then export the table using the third slot of saveCProf paths list
		write.table(reactomeList[[nm]], paste0(saveCProf[3], "cProf_FindMarkers_REACTOME_", nm, "_", id, ".csv"), row.names=F, col.names=T, quote=F, sep="\t")
		write.xlsx(reactomeList[[nm]], paste0(saveCProf[3], "cProf_FindMarkers_REACTOME_", nm, "_", id, ".xlsx"))		
	}
	# finally, return results
	return(TRUE)
}

# interestingly, it is always possibe to subset any EnrichmentObject using the information available here:

# https://github.com/YuLab-SMU/enrichplot/issues/17
# http://yulab-smu.top/clusterProfiler-book/chapter13.html

# subsetting is surely useful to show only those terms that are actually interesting to the process under investigation

# some KEGG
keggPlotting <- function(rkDataF, saveKggImg, saveCProf, savePV, lbl, fN, cS) {

	# rkDataF <- originalData; saveKggImg <- saveKegg; savePV <- savePView; lbl <- lbl; fN <- fName; cS <- cS

	# all the plots are found here: https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#enrichplot

	gene <- rkDataF$gene
	
	# generate Rankings data structure
	ranks <- rkDataF$foldChange
	names(ranks) <- rkDataF$gene

	#getting SYMBOLS IDs from ENTREZ
	quiet(gSymb <- AnnotationDbi::select(org.Mm.eg.db, keys=as.character(gene), columns=c("SYMBOL", "ENTREZID"), keytype="ENTREZID"), all=T)

	eKEGG <- enrichKEGG(gene=gene, organism ='mmu', pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05)
	
	# entrez to symbols
	symbs <- vector()
	expKEGG <- as.data.frame(eKEGG)
	
	if (nrow(expKEGG) > 0) {
		for (nr in seq(1, nrow(expKEGG), by=1)) {
			# get symbols
			symbs <- append(symbs, paste(gSymb$SYMBOL[gSymb$ENTREZID %in% strsplit(expKEGG[nr, 8], "/")[[1]]], collapse="/"))
		}
		# get final table ready
		expKEGG <- cbind.data.frame(as.data.frame(eKEGG), symbol=symbs)
		
		write.table(expKEGG, paste0(saveCProf, "cProf_DESEQ_on_rankings_Kegg_cluster_", cS, ".csv", sep=""), row.names=F, col.names=T, quote=F, sep="\t")
		write.xlsx(expKEGG, paste0(saveCProf, "cProf_DESEQ_on_rankings_Kegg_cluster_", cS, ".xlsx"))
	} else {
		print("no enrichment")
	}
	# get directory to get back here after pathview
	herE <- getwd()
	if (nrow(eKEGG) > 0) {
		
		# plot only first ten to make it faster
		setwd(savePV)
		for (idK in eKEGG$ID[1:10]) {
			pathview(gene.data=ranks, pathway.id=idK, species="mmu", limit=list(gene=max(abs(ranks)), cpd=1))
		}
		setwd(herE)

		lbl <- lbl + 1
		p <- dotplot(eKEGG)
		png(paste0(saveKggImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "dotplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
		plot(p)
		dev.off()

		# setReadable translate entrezID into symbols
		eKEGGx <- setReadable(eKEGG, 'org.Mm.eg.db', 'ENTREZID')
		
		lbl <- lbl + 1
		p <- cnetplot(eKEGGx, foldChange=ranks, categorySize="pvalue", circular = TRUE, colorEdge = TRUE)
		png(paste0(saveKggImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "cnetplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in", width=32, height=32, pointsize=12, res=300)
		plot(p)
		dev.off()
		
		lbl <- lbl + 1
		p <- heatplot(eKEGGx, foldChange=ranks)
		png(paste0(saveKggImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "heatplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in",width=40, height=32, pointsize=12, res=300)
		plot(p)
		dev.off()

		# to plot a pairwise plot, there should be at least a couple terms
		if (nrow(eKEGG) >= 2) {
			dev.new()
			eKEGG <- pairwise_termsim(eKEGG)
			lbl <- lbl + 1
			p <- emapplot(eKEGG, cex_category=1.5,layout="kk")
			png(paste0(saveKggImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "pairwiseplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in", width=20, height=18, pointsize=12, res=300)
			plot(p)
			dev.off()
		}
	}
	# return label
	return(lbl)
}

# some GO
goPlotting <- function(rkDataF, saveGOImg, saveCProf, lbl, fN, cS) {

	# rkDataF <- originalData; saveGOImg <- saveGO; lbl <- lbl;

	gene <- rkDataF$gene
	
	# generate Rankings data structure
	ranks <- rkDataF$foldChange
	names(ranks) <- rkDataF$gene

	#getting SYMBOLS IDs from ENTREZ
	quiet(gSymb <- AnnotationDbi::select(org.Mm.eg.db, keys=as.character(gene), columns=c("SYMBOL", "ENTREZID"), keytype="ENTREZID"), all=T)

	eGO <- enrichGO(gene=gene, OrgDb=org.Mm.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05)
	
	# entrez to symbols
	symbs <- vector()
	expGO <- as.data.frame(eGO)
	
	if (nrow(expGO) > 0) {
		for (nr in seq(1, nrow(expGO), by=1)) {
			# get symbols
			symbs <- append(symbs, paste(gSymb$SYMBOL[gSymb$ENTREZID %in% strsplit(expGO[nr, 8], "/")[[1]]], collapse="/"))
		}
		# get final table ready
		expGO <- cbind.data.frame(as.data.frame(eGO), symbol=symbs)
		write.table(expGO, paste0(saveCProf, "cProf_DESEQ_on_rankings_GO_Bio_Process_cluster_", cS, ".csv", sep=""), row.names=F, col.names=T, quote=F, sep="\t")
		write.xlsx(expGO, paste0(saveCProf, "cProf_DESEQ_on_rankings_GO_Bio_Process_cluster_", cS, ".xlsx"))
	} else {
		print("no enrichment")
	}

	if (nrow(eGO) > 0) {
		lbl <- lbl + 1
		p <- goplot(eGO)
		png(paste0(saveGOImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "goplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in", width=20, height=15, pointsize=12, res=300)
		plot(p)
		dev.off()

		lbl <- lbl + 1
		p <- dotplot(eGO)
		png(paste0(saveGOImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "dotplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
		plot(p)
		dev.off()

		# setReadable translate entrezID into symbols
		eGOx <- setReadable(eGO, 'org.Mm.eg.db', 'ENTREZID')
		lbl <- lbl + 1
		p <- cnetplot(eGOx, foldChange=ranks, categorySize="pvalue", circular = TRUE, colorEdge = TRUE)
		png(paste0(saveGOImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "cnetplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in", width=32, height=32, pointsize=12, res=300)
		plot(p)
		dev.off()

		lbl <- lbl + 1
		p <- heatplot(eGOx, foldChange=ranks)
		png(paste0(saveGOImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "heatplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in",width=40, height=32, pointsize=12, res=300)
		plot(p)
		dev.off()

		# to plot a pairwise plot, there should be at least a couple terms
		if (nrow(eGOx) >= 2) {
			eGO <- pairwise_termsim(eGOx)
			lbl <- lbl + 1
			p <- emapplot(eGO, cex_category=1.5,layout="kk")
			png(paste0(saveGOImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "pairwiseplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in", width=20, height=18, pointsize=12, res=300)
			plot(p)
			dev.off()
		}
	} else {
		print("no enrichment")
	}
	# return label
	return(lbl)
}

# some REACTOME
reactomePlotting <- function(rkDataF, saveRctImg, saveCProf, fcThresh, lbl, fN, cS) {
	
	# rkDataF <- originalData; saveRctImg <- saveReactome; fcThresh <- fcThresh; lbl <- lbl;

	de <- rkDataF$gene[abs(rkDataF$foldChange) > fcThresh]
	
	# generate Rankings data structure
	ranks <- rkDataF$foldChange
	names(ranks) <- rkDataF$gene

	#getting SYMBOLS IDs from ENTREZ
	quiet(gSymb <- AnnotationDbi::select(org.Mm.eg.db, keys=as.character(rkDataF$gene), columns=c("SYMBOL", "ENTREZID"), keytype="ENTREZID"), all=T)

	eRCTM <- enrichPathway(gene=de, organism="mouse", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05, readable=TRUE)
	
	# entrez to symbols
	symbs <- vector()
	expRCT <- as.data.frame(eRCTM)
	
	if (nrow(expRCT) > 0) {
		for (nr in seq(1, nrow(expRCT), by=1)) {
			# get symbols
			symbs <- append(symbs, paste(gSymb$SYMBOL[gSymb$ENTREZID %in% strsplit(expRCT[nr, 8], "/")[[1]]], collapse="/"))
		}
		# get final table ready
		expRCT <- cbind.data.frame(as.data.frame(eRCTM), symbol=symbs)
		write.table(expRCT, paste0(saveCProf, "cProf_DESEQ_on_rankings_Reactome_cluster_", cS, ".csv", sep=""), row.names=F, col.names=T, quote=F, sep="\t")
		write.xlsx(expRCT, paste0(saveCProf, "cProf_DESEQ_on_rankings_Reactome_cluster_", cS, ".xlsx"))
	} else {
		print("no enrichment")
	}

	if (nrow(eRCTM) > 0) {
		lbl <- lbl + 1
		p <- dotplot(eRCTM)
		png(paste0(saveRctImg, "04_0", strsplit(gsub("rankings_genes", "dotplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
		plot(p)
		dev.off()
		
		# setReadable translate entrezID into symbols
		eRCTx <- setReadable(eRCTM, 'org.Mm.eg.db', 'ENTREZID')
		lbl <- lbl + 1
		p <- cnetplot(eRCTx, foldChange=ranks, categorySize="pvalue", circular = TRUE, colorEdge = TRUE)
		png(paste0(saveRctImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "cnetplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in", width=32, height=32, pointsize=12, res=300)
		plot(p)
		dev.off()
		
		lbl <- lbl + 1
		p <- heatplot(eRCTx, foldChange=ranks)
		png(paste0(saveRctImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "heatplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in",width=40, height=32, pointsize=12, res=300)
		plot(p)
		dev.off()

		# to plot a pairwise plot, there should be at least a couple pathways
		if (nrow(eRCTx) >= 2) {
			eRCT <- pairwise_termsim(eRCTx)
			lbl <- lbl + 1
			p <- emapplot(eRCT, cex_category=1.5,layout="kk")
			png(paste0(saveRctImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "pairwiseplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in", width=20, height=18, pointsize=12, res=300)
			plot(p)
			dev.off()
		}
	} else {
		print("no enrichment")
	}
	# return label
	return(lbl)
}

################## OPTION ONE FOR BATCHED EAE DATA

optionOneBatched <- function(subSet, expN, saveGReg, personalT, batch, cTy, thresHD) {
	
#	subSet <- cTSubSet; saveGReg <- saveUpDown; batch <- testingExp;

	cTSce <- list()

	# then get Onset EAE celltype
	cTSce[[batch[1]]] <- subSet[, which(subSet$batch==batch[1])]
	# and Chronic EAE celltype
	cTSce[[batch[2]]] <- subSet[, which(subSet$batch==batch[2])]

	#print info about number of cells for each batch
	print(paste0("in cluster ", cS, " there are ", ncol(cTSce[[batch[1]]]), " ", cTy, " for ", batch[1]))
	print(paste0("and ", ncol(cTSce[[batch[2]]]), " ", cTy, " for ", batch[2]))

	# this if else check is useful here and not for the two other options since here all the cellType cells are considered
	if (ncol(cTSce[[batch[1]]]) > personalT && ncol(cTSce[[batch[2]]]) > personalT) {

		# then search for marker genes using pairwise by batch so that WT will be compared to 3x (and vice versa)
		# https://rdrr.io/bioc/scran/man/findMarkers.html
		# it is important to note that findMarkers uses logcounts as the default assay, to compute the differentially expressed genes.
		# however when we apply zinbwave, its normalised values are put into the logcounts assay
		# also, only 1000 highly variable genes are used
		
		markersG <- findMarkers(subSet, groups=subSet$batch, pval.type="any", direction="any")

		# http://bioconductor.org/books/release/OSCA/marker-detection.html#looking-for-any-differences
		# The default philosophy of findMarkers() is to identify a combination of marker genes that - together - uniquely define one cluster against the rest. 
		# To this end, we collect the top DE genes from each pairwise comparison involving a particular cluster to assemble a set of candidate markers for that cluster.
		# We will demonstrate on cluster 7; the relevant DataFrame contains log2-fold changes of expression in cluster 7 over each other cluster, along with several 
		# statistics obtained by combining p-values across the pairwise comparisons involving 7.
		# chosen <- "7"
		# interesting <- markers.pbmc[[chosen]]

		# hence, if cond[1] stands for "3xTG", we are obtaining the fold change for the comparison of 3xTG over WT

		# filter on false discovery rate considering batch[1] as the chosen batch
		cC <- batch[1]
		chosen <- markersG[[cC]]
		selfMade <- chosen[chosen$FDR < thresHD, ]

		print(paste0(nrow(selfMade), " statistical significant genes were found for all the ", cTy, " cells in cluster ", cS))
		print(paste0("for the comparison of ", cC, " over ", batch[2]))

		# and export data
		finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

		# order by fold change, which is the 6th column
		finalOrdered <- list()
		finalOrdered[["UpAndDown"]] <-  finalGenes[order(finalGenes[, 6], decreasing=T), ]

		# get up and down, separated
		finalOrdered[["Up"]] <- finalOrdered[["UpAndDown"]][which(finalOrdered[["UpAndDown"]][, 6]>=0), ]
		dTmp <- finalOrdered[["UpAndDown"]][which(finalOrdered[["UpAndDown"]][, 6]<0), ]
		finalOrdered[["Down"]] <- dTmp[order(dTmp[, 6], decreasing=F), ]

		write.table(finalOrdered[["Up"]], paste0(saveGReg, cTy, "_cluster_", cS, "_exp_", expN, "_UP_REG.csv"), quote=F, row.names=F, col.names=T, sep="\t")
		write.xlsx(finalOrdered[["Up"]], paste0(saveGReg, cTy, "_cluster_", cS, "_exp_", expN, "_UP_REG.xlsx"))

		write.table(finalOrdered[["Down"]], paste0(saveGReg, cTy, "_cluster_", cS, "_exp_", expN, "_DOWN_REG.csv"), quote=F, row.names=F, col.names=T, sep="\t")
		write.xlsx(finalOrdered[["Down"]], paste0(saveGReg, cTy, "_cluster_", cS, "_exp_", expN, "_DOWN_REG.xlsx"))

		write.table(finalOrdered[["UpAndDown"]], paste0(saveGReg, cTy, "_cluster_", cS, "_exp_", expN, "_up_and_down.csv"), quote=F, row.names=F, col.names=T, sep="\t")
		write.xlsx(finalOrdered[["UpAndDown"]], paste0(saveGReg, cTy, "_cluster_", cS, "_exp_", expN, "_up_and_down.xlsx"))

		# store the genes for this cluster in the list
		optOneGenes <- finalGenes$gName

	} else {
		# set return to NULL
		optOneGenes <- NULL
		print(paste0("not enough cells to compare the two batches for ", cTy, " in cluster ", cS))
	}
	# finally, return the results
	return(optOneGenes)
}

################## OPTION TWO FOR BATCHED EAE DATA

optionTwoBatched <- function(subSet, expN, saveAbs=NULL, saveGReg, personalT, batch, cTy, mgdAB, thresHD) {

	cTSce <- list()

	# then get Onset EAE celltype
	cTSce[[batch[1]]] <- subSet[, which(subSet$batch==batch[1])]
	# and Chronic EAE celltype
	cTSce[[batch[2]]] <- subSet[, which(subSet$batch==batch[2])]

	#print info about number of cells for each batch
	print(paste0("in cluster ", cS, " there are ", ncol(cTSce[[batch[1]]]), " cells for ", batch[1]))
	print(paste0("and ", ncol(cTSce[[batch[2]]]), " cells for ", batch[2]))

	if (ncol(cTSce[[batch[1]]]) > personalT && ncol(cTSce[[batch[2]]]) > personalT) {

		# then search for marker genes using pairwise by batch so that WT will be compared to 3x (and vice versa)
		# https://rdrr.io/bioc/scran/man/findMarkers.html
		# it is important to note that findMarkers uses logcounts as the default assay, to compute the differentially expressed genes.
		# however when we apply zinbwave, its normalised values are put into the logcounts assay
		# also, only 1000 highly variable genes are used
		
		markersG <- findMarkers(subSet, groups=subSet$batch, pval.type="any", direction="any")

		# http://bioconductor.org/books/release/OSCA/marker-detection.html#looking-for-any-differences
		# The default philosophy of findMarkers() is to identify a combination of marker genes that - together - uniquely define one cluster against the rest. 
		# To this end, we collect the top DE genes from each pairwise comparison involving a particular cluster to assemble a set of candidate markers for that cluster.
		# We will demonstrate on cluster 7; the relevant DataFrame contains log2-fold changes of expression in cluster 7 over each other cluster, along with several 
		# statistics obtained by combining p-values across the pairwise comparisons involving 7.
		# chosen <- "7"
		# interesting <- markers.pbmc[[chosen]]

		# hence, if cond[1] stands for "3xTG", we are obtaining the fold change for the comparison of 3xTG over WT

		# filter on false discovery rate considering batch[1] as the chosen batch
		cC <- batch[1]
		chosen <- markersG[[cC]]
		selfMade <- chosen[chosen$FDR < thresHD, ]

		print(paste0(nrow(selfMade), " statistical significant genes were found for cluster ", cS))
		print(paste0("for the comparison of ", cC, " over ", batch[2]))

		# and export data
		finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

		# order by fold change, which is the 6th column, which is the 6th column
		finalOrdered <- list()
		finalOrdered[["UpAndDown"]] <-  finalGenes[order(finalGenes[, 6], decreasing=T), ]

		# get up and down, separated
		finalOrdered[["Up"]] <- finalOrdered[["UpAndDown"]][which(finalOrdered[["UpAndDown"]][, 6]>=0), ]
		dTmp <- finalOrdered[["UpAndDown"]][which(finalOrdered[["UpAndDown"]][, 6]<0), ]
		finalOrdered[["Down"]] <- dTmp[order(dTmp[, 6], decreasing=F), ]
		
		# store the genes for this cluster in the list
		optTwoGenes <- finalGenes$gName
	} else {
		# set return to NULL
		optTwoGenes <- NULL
		print(paste0("not enough cells to compare the two batches for cluster ", cS))
	}

	# if saveAbs is provided, then analyse ABSeq information
	if (!is.null(saveAbs)) {
	
		# since it is not possible to directly infer the cell type using the antibodies that tagged a cell, what it might be done is to retrieve all the cells in a specific cluster and create a plot, or more, that defines the kinds of antibodies that are tagging such cells. this can be done by each cell type in the cluster

		# hence get the cell composition from the ABseq perspective

		# get all cell types in the cluster
		cTTmp <- unique(subSet$fine_to_main)
		# remove NAs
		cTypes <- cTTmp[!is.na(cTTmp)]
		
		# summarise cell types
		suMM <- table(subSet$fine_to_main)

		# initialise a storage for the barplot information
		barTmp <- list()
		barTabs <- vector()
		# loop through the cell types
		for (cTp in cTypes) {

			# get the cell names of that specific kind of cells in the cluster under analysis
			cNames <- colnames(subSet[, which(subSet$fine_to_main==cTp)])

			# get the abseq info for those cells
			subAB <- mgdAB[, colnames(mgdAB) %in% cNames]

			#create the data storage
			barTmp[[cTp]] <- cbind.data.frame(gene=antiB$antibody, count=rep(0, length(antiB$antibody)))

			# loop through the retrieved cells to gather info for the barplot
			for (cll in seq(1, ncol(subAB), by=1)) {
				# reduce multi antibody and sum up their contribution
				tmpDf <- rbind.data.frame(barTmp[[cTp]], mgdAB$antiB[[cll]])
				barTmp[[cTp]] <- ddply(tmpDf, .(gene), function(x) c(count=sum(x$count)))
			}

			# create readable names
			tmpN <- sapply(strsplit(barTmp[[cTp]]$gene, "\\|"), "[[", 1)
			barTmp[[cTp]]$gene <- sapply(strsplit(tmpN, "\\:"), "[[", 1)
			# divide counts by the number of cells, to obain an averaged value
			barTmp[[cTp]]$avg <- barTmp[[cTp]]$count / ncol(subAB)
			
			# bind the new barTmp
			valCT <- suMM[which(names(suMM)==cTp)]
			barTabs <- rbind.data.frame(barTabs, cbind.data.frame(gene=barTmp[[cTp]]$gene, cellT=paste0(rep(cTp, nrow(barTmp[[cTp]])), "=", valCT), avgCount=barTmp[[cTp]]$avg))
		}

		# eventually print each barTmp
		p <- ggplot(barTabs, aes(fill=gene, y=avgCount, x=gene)) + 
			geom_bar(position="dodge", stat="identity") +
			scale_fill_manual(values=colorsPal) +
			ggtitle(paste0("abSeq average content for cluster ", cS)) +
			facet_wrap(~cellT) +
			theme(legend.position="none", axis.text.x=element_text(angle=45, vjust=0.5, hjust=1)) +
			xlab("")
			
		gg.gap(plot=p, segments=c(500, 550), ylim=c(0, 2000))

		# update label
		lbl <- lbl + 1
		png(paste0(saveAbs, "04_0", lbl, "_abSeq_average_content_for cluster_", cS, "_in_", cTy, ".png"), type="cairo", units="in", width=16, height=12, pointsize=12, res=300)
		plot(p)
		dev.off()
	} else {
		# if saveAbs is not provided, i.e. it is NULL, don't do anything
		print("abseq information is missing, not going to analyse it")
	}
	# finally, return the results
	return(optTwoGenes)
}

################## OPTION THREE FOR BATCHED EAE DATA

optionThreeBatched <- function(subSet, expN, saveGReg, personalT, batch, cTy, thresHD) {

#	subSet <- cTSubSet; saveGReg <- saveUpDown; batch <- testingExp

	cTSce <- list()

	# then get Onset EAE celltype
	cTSce[[batch[1]]] <- subSet[, which(subSet$batch==batch[1])]
	# and Chronic EAE celltype
	cTSce[[batch[2]]] <- subSet[, which(subSet$batch==batch[2])]

	#print info about number of cells for each batch
	print(paste0("in the full dataset there are ", ncol(cTSce[[batch[1]]]), " ", cTy, " for ", batch[1]))
	print(paste0("and ", ncol(cTSce[[batch[2]]]), " ", cTy, " for ", batch[2]))

	if (ncol(cTSce[[batch[1]]]) > personalT && ncol(cTSce[[batch[2]]]) > personalT) {

		# then search for marker genes using pairwise by batch so that WT will be compared to 3x (and vice versa)
		# https://rdrr.io/bioc/scran/man/findMarkers.html
		# it is important to note that findMarkers uses logcounts as the default assay, to compute the differentially expressed genes.
		# however when we apply zinbwave, its normalised values are put into the logcounts assay
		# also, only 1000 highly variable genes are used
		
		markersG <- findMarkers(subSet, groups=subSet$batch, pval.type="any", direction="any")
		
		# http://bioconductor.org/books/release/OSCA/marker-detection.html#looking-for-any-differences
		# The default philosophy of findMarkers() is to identify a combination of marker genes that - together - uniquely define one cluster against the rest. 
		# To this end, we collect the top DE genes from each pairwise comparison involving a particular cluster to assemble a set of candidate markers for that cluster.
		# We will demonstrate on cluster 7; the relevant DataFrame contains log2-fold changes of expression in cluster 7 over each other cluster, along with several 
		# statistics obtained by combining p-values across the pairwise comparisons involving 7.
		# chosen <- "7"
		# interesting <- markers.pbmc[[chosen]]

		# hence, if cond[1] stands for "3xTG", we are obtaining the fold change for the comparison of 3xTG over WT

		# filter on false discovery rate considering batch[1] as the chosen batch
		cC <- batch[1]
		chosen <- markersG[[cC]]
		selfMade <- chosen[chosen$FDR < thresHD, ]

		print(paste0(nrow(selfMade), " statistical significant genes were found for all the ", cTy, " cells"))
		print(paste0("for the comparison of ", cC, " over ", batch[2]))

		# and export data
		finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

		# order by fold change, which is the 6th column
		finalOrdered <- list()
		finalOrdered[["UpAndDown"]] <-  finalGenes[order(finalGenes[, 6], decreasing=T), ]

		# datatable(extensions = "Buttons", finalOrdered,filter="top", options = list(dom = 'Blfrtip',buttons = c('copy','csv','excel','pdf','print'),lengthMenu = list(c(10,25,50,-1),c(10,25,50,"All")),pageLength = 10, scrollX=T))

		# get up and down, separated
		finalOrdered[["Up"]] <- finalOrdered[["UpAndDown"]][which(finalOrdered[["UpAndDown"]][, 6]>=0), ]
		dTmp <- finalOrdered[["UpAndDown"]][which(finalOrdered[["UpAndDown"]][, 6]<0), ]
		finalOrdered[["Down"]] <- dTmp[order(dTmp[, 6], decreasing=F), ]

		write.table(finalOrdered[["Up"]], paste0(saveGReg, cTy, "_ONLY_exp_", expN, "_UP_REG.csv"), quote=F, row.names=F, col.names=T, sep="\t")
		write.xlsx(finalOrdered[["Up"]], paste0(saveGReg, cTy, "_ONLY_exp_", expN, "_UP_REG.xlsx"))

		write.table(finalOrdered[["Down"]], paste0(saveGReg, cTy, "_ONLY_exp_", expN, "_DOWN_REG.csv"), quote=F, row.names=F, col.names=T, sep="\t")
		write.xlsx(finalOrdered[["Down"]], paste0(saveGReg, cTy, "_ONLY_exp_", expN, "_DOWN_REG.xlsx"))

		write.table(finalOrdered[["UpAndDown"]], paste0(saveGReg, cTy, "_ONLY_exp_", expN, "_up_and_down.csv"), quote=F, row.names=F, col.names=T, sep="\t")
		write.xlsx(finalOrdered[["UpAndDown"]], paste0(saveGReg, cTy, "_ONLY_exp_", expN, "_up_and_down.xlsx"))

		# store the genes for this cluster in the list
		optThreeGenes <- finalGenes$gName
	} else {
		# set return to NULL
		optThreeGenes <- NULL
		print(paste0("not enough", cTy, " cells to compare the two batches"))
	}
	return(optThreeGenes)
}

# perform enrichment using DESeq2
gseaMyDataBatched <- function(subSet, condS, saveRk, saveGTab, saveGES, savePages, expN, idGSEA, cTy) {

#	subSet <- cTSubSet; condS <- testingExp; saveRk <- rN; saveGTab <- saveGSEA; saveGES <- saveGes; savePages <- savePages; expN <- expN; idGSEA <- paste0(cellType, "_cluster_", cS);
	
	# obtain a list of DEGs
	subSet$batch <- as.factor(subSet$batch)
	quiet(dds <- DESeqDataSetFromMatrix(as.matrix(assays(subSet)$counts), colData=colData(subSet), design=~batch), all=T)
	## Pre filter

	# We perform a minimal pre-filtering to keep only rows that have at least 10 reads total. Note that more strict filtering to increase power is automatically applied via independent filtering on the mean of normalized counts within the results function

	keep <- rowSums(counts(dds)) >= 10
	dds <- dds[keep, ]

	# After this first filter, we remove genes without an EntrezID
	rowData(dds)$Symbol <- rownames(dds)
	quiet(rowData(dds)$Entrezid <- mapIds(x=org.Mm.eg.db, keys=rowData(dds)$Symbol, column="ENTREZID", keytype="SYMBOL"), all=T)

	keep <- !is.na(rowData(dds)$Entrezid)
	dds <- dds[keep, ]

	# set "Onset" as the reference level
	dds$batch <- factor(dds$batch, levels=c("Onset", "Chronic"))
	
	## Differential Expression

	# We launch DESeq2
	quiet(ddsDesequed <- DESeq(dds), all=T)
	res <- results(ddsDesequed)

	# A positive value of Log-Fold Change indicates higher expression values in 3X_Alz compared to WT

	# The analysis is performed by:
	# - ranking all genes in the data set;
	# - identifying the rank positions of all members of the gene set in the ranked data set;
	# - calculating an enrichment score (ES) that represents the difference between the observed rankings and that which would be expected assuming a random rank distribution.

	## Create ranks

	# Rank all genes based on their statistics
	rkVals <- res$stat
	names(rkVals) <- rowData(ddsDesequed)$Entrezid
	
	# finally, return results
	return(rkVals)
}

########################################### FUNCTIONAL SCORE ###########################################

#For all signatures except neutrophil aging, functional scores were defined as the 
#average normalized expression of corresponding genes.

# computing functional score
simpleFC <- function(assy, rfc) {

#	assy <- assay(cTSubSet, "logcounts"); rfc <- gSign

	# compute average expression for a gene and all the cells
	scoreS <- rowMeans(assy)
	res <- list()
	# loop through the signatures
	for (sg in seq(1, ncol(rfc), by=1)){
		# get gene names of those genes that belong to a signature
		foundG <- rownames(assy)[rownames(assy) %in% rfc[, sg]]
		# and get their average exp value
		valS <- scoreS[rownames(assy) %in% rfc[, sg]]
		res[[names(gSign)[sg]]] <- cbind.data.frame(geneNames=foundG, functScore=valS)
	}
	return(res)
}

########################################### AGING SCORE ###########################################
#Aging score was defined as the weighted average of Z scores of age-related genes, 
#where the Z scores were calculated by scaling the normalized expression of a gene
#across all cells. 


#Gene weights were set to either 1 or −1 to reflect positive or negative relationships.
#The neutrophil maturation signature was derived by identifying the top 50 DEGs (as listed
#in Supplementary Table 4) with the highest fold-changes and adjusted P values < 0.05 
#between the mature cluster (G4) and immature clusters (G0–G3). Granule signatures were
#from ref. 20. Other functional signatures were derived from the Gene Ontology database,
#with the full gene list provided in Supplementary Table 4. For instance, to access the 
#phagocytosis function at the transcript level, we determined a phagocytosis score by 
#calculating the average expression of genes in the Gene Ontology term:
#‘phagocytosis, engulfment’ (GO: 0006911).

zscoreAC <- function(assy, rfc) {

#	assy <- assay(cTSubSet, "logcounts"); rfc <- gAgin

	zS <- vector()

	for (rw in seq(1, nrow(assy), by=1)) {
		# compute average expression for a gene and all the cells
		sD <- sd(assy[rw, ])
		avG <- mean(assy[rw, ])
		zS <- rbind.data.frame(zS, (assy[rw, ] - avG) / sD)
	}
	colnames(zS) <- colnames(assy)
	# get gene names of those genes that belong to a signature
	foundG <- rownames(assy)[rownames(assy) %in% rfc[, 1]]
	# and get their average exp value
	valS <- zS[rownames(assy) %in% rfc[, 1], ]
	res <- cbind.data.frame(geneNames=foundG, functScore=valS)

	return(res)
}
