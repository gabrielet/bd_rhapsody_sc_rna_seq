# Load libraries
library("bluster")
library("cluster")
library("compositions")
library("cowplot")
library("ggplot2")
library("gridExtra")
library("pheatmap")
library("plyr")
library("RColorBrewer")
library("scater")
library("scran")
library("Seurat")
library("SingleCellExperiment")
library("BiocParallel")
library("zinbwave")
library("clusterExperiment")
library("batchelor")

# select experiment
expN <- "Chronic"

print(paste0("analysing exp ", expN))

# set paths
raw_path <- "/home/gabriele/work/cbmc/scrnaseq/raw_data_new/"
saveIn <- paste0(raw_path, "EAE_", expN, "/")
saveImg <- paste0(saveIn, "figures/")

#check if saveIn directory exist, otherwise create it
ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
ifelse(dir.exists(saveImg), TRUE, dir.create(saveImg, recursive=F))

# set palette with 12 colours for colouring plots
colorsPal <- brewer.pal(12, "Paired")

# load data
cleanSce <- readRDS(paste0(saveIn, paste0("QCed_data_", expN, ".Rds")))

#' MODELLING GENE VARIANCE

#' 1) [log-normalized expression values](https://bioconductor.org/books/release/OSCA/feature-selection.html#variance-of-the-log-counts){target="_blank"}  

#' The simplest approach to quantifying per-gene variation is to simply compute the variance of the log-normalized expression values (referred to as “log-counts” for simplicity) for each gene across all cells in the population. This has an advantage in that the feature selection is based on the same log-values that are used for later downstream steps. In particular, genes with the largest variances in log-values will contribute the most to the Euclidean distances between cells. By using log-values here, we ensure that our quantitative definition of heterogeneity is consistent throughout the entire analysis.

# computing variation
geneVar <- modelGeneVar(cleanSce)

# visualising the fit
fitGVar <- metadata(geneVar)
png(paste0(saveImg, "01_01_geneVar_kmeans.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	plot(fitGVar$mean, fitGVar$var, xlab="Mean of log-expression", ylab="Variance of log-expression")
	curve(fitGVar$trend(x), col=colorsPal[5], add=TRUE, lwd=2)
dev.off()

#' The trend fit has several useful parameters (see ?fitTrendVar) that can be tuned for a more appropriate fit. For example, the defaults can occasionally yield an overfitted trend when the few high abundance genes are also highly variable. In such cases, users can reduce the contribution of those high-abundance genes by turning off density weights  

noWeight <- modelGeneVar(cleanSce, density.weights=FALSE)
fitNoW <- metadata(noWeight)

# visualising the fit and the non-weighted fit
png(paste0(saveImg, "01_02_geneVar_no_weights_kmeans.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	plot(fitGVar$mean, fitGVar$var, xlab="Mean of log-expression", ylab="Variance of log-expression")
	curve(fitGVar$trend(x), col=colorsPal[5], add=TRUE, lwd=2)
	curve(fitNoW$trend(x), col=colorsPal[1], add=TRUE, lwd=2)
	legend("topleft", col=c(colorsPal[5], colorsPal[1]), legend=c("Default", "No weight"), lwd=2)
dev.off()

#' 2) [quantifying technical noise](https://bioconductor.org/books/release/OSCA/feature-selection.html#sec:spikeins){target="_blank"}  

#' TODO:  
#' how do we know if our data contains spike-ins or not?  

#' according to Edyta from BD: the PhiX is additional sequences that are spiked in Rhapsody library and sequenced on Illumina. Those additional sequences are there for the purpose to increase the complexity of Rhaposdy library during sequencing but should not appear afterwards in the raw data, thus they are removed by Illumina software. In Rhaposdy pipeline we are verifying if those reads were completely removed and if all went fine with PhiX spiked in. For that reason we include in reference PhiX sequence to map against it the reads.  

#' assuming they don't: In the absence of spike-in data, one can attempt to create a trend by making some distributional assumptions about the noise. For example, UMI counts typically exhibit near-Poisson variation if we only consider technical noise from library preparation and sequencing.  

set.seed(131)
mergedPois <- modelGeneVarByPoisson(cleanSce)
mergedPois <- mergedPois[order(mergedPois$bio, decreasing=TRUE),]
head(mergedPois)

#' Interestingly, trends based purely on technical noise tend to yield large biological components for highly-expressed genes. This often includes so-called “house-keeping” genes coding for essential cellular components such as ribosomal proteins, which are considered uninteresting for characterizing cellular heterogeneity. These observations suggest that a more accurate noise model does not necessarily yield a better ranking of HVGs, though one should keep an open mind - house-keeping genes are regularly DE in a variety of conditions, and the fact that they have large biological components indicates that there is strong variation across cells that may not be completely irrelevant.  

# plot poisson noise trend
png(paste0(saveImg, "01_03_techNoise_kmeans.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	plot(mergedPois$mean, mergedPois$total, pch=16, xlab="Mean of log-expression", ylab="Variance of log-expression")
	curve(metadata(mergedPois)$trend(x), col=colorsPal[5], add=TRUE)
dev.off()

# same but using blocking factors
mergedPois <- modelGeneVarByPoisson(cleanSce, block=cleanSce$condition)
mergedPois <- mergedPois[order(mergedPois$bio, decreasing=TRUE),]

png(paste0(saveImg, "01_04_techNoise_blocked_kmeans.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	par(mfrow=c(1,2))
	blocked.stats <- mergedPois$per.block
	for (i in colnames(blocked.stats)) {
		current <- blocked.stats[[i]]
		plot(current$mean, current$total, main=i, pch=16, cex=0.5, xlab="Mean of log-expression", ylab="Variance of log-expression")
		curfit <- metadata(current)
		points(curfit$mean, curfit$var, col=colorsPal[5], pch=16)
		curve(curfit$trend(x), col=colorsPal[2], add=TRUE, lwd=2) 
	}
dev.off()

#' [SELECTING HIGHLY VARIABLE GENES](https://bioconductor.org/books/release/OSCA/feature-selection.html#hvg-selection){target="_blank"}  

#' assuming we performed the modelling of gene variation using not weighted method  

# order by gene variation
orderedNoW <- noWeight[order(noWeight$bio, decreasing=TRUE),]

#' The simplest HVG selection strategy is to take the [top X](https://bioconductor.org/books/release/OSCA/feature-selection.html#based-on-the-largest-metrics){target="_blank"} genes with the largest values for the relevant variance metric. The main advantage of this approach is that the user can directly control the number of genes retained, which ensures that the computational complexity of downstream calculations is easily predicted  

# using first 1000 to speed-up the computation
highVar <- getTopHVGs(orderedNoW, n=1000)

#' Another approach to feature selection is to [set a fixed threshold](https://bioconductor.org/books/release/OSCA/feature-selection.html#based-on-significance){target="_blank"} of one of the metrics. This is most commonly done with the (adjusted) p-value reported by each of the above methods. The p-value for each gene is generated by testing against the null hypothesis that the variance is equal to the trend.  

# highVar <- getTopHVGs(orderedNoW, fdr.threshold=0.05)

#' Here, the aim is to only remove the obviously uninteresting genes with [variances below the trend](https://bioconductor.org/books/release/OSCA/feature-selection.html#feature-selection-positive){target="_blank"}. By doing so, we avoid the need to make any judgement calls regarding what level of variation is interesting enough to retain.  

# highVar <- getTopHVGs(orderedNoW, var.threshold=0)

#' [DIMENSIONALITY REDUCTION](https://bioconductor.org/books/release/OSCA/dimensionality-reduction.html#principal-components-analysis){target="_blank"}

#' 2) using [technical noise](https://bioconductor.org/books/release/OSCA/dimensionality-reduction.html#using-the-technical-noise){target="_blank"}  

#' Note that denoisePCA() imposes internal caps on the number of PCs that can be chosen in this manner. By default, the number is bounded within the “reasonable” limits of 5 and 50 to avoid selection of too few PCs (when technical noise is high relative to biological variation) or too many PCs (when technical noise is very low).  

set.seed(131)
denoisedSce <- denoisePCA(cleanSce, technical=noWeight, subset.row=highVar)

# apply batch correction using batchelor
set.seed(131)
corrected <- fastMNN(denoisedSce, batch=denoisedSce$condition) 

# get the corrected dimensionality
reducedDim(denoisedSce, "corrected") <- reducedDim(corrected)

set.seed(131)
# run TSNE
#zbWaved <- runTSNE(zbWaved, dimred="zinbwave", perplexity=30) # QUESTO QUELLO CHE SI FACEVA CON ZINBWAVE
denoisedSce <- runTSNE(denoisedSce, dimred="corrected", perplexity=30)

# and export a summary table
summ <- table(denoisedSce$pruned_main, denoisedSce$condition)
summarY <- data.frame(rownames(summ), summ[, 1], summ[, 2])

##################################################################### EAE
colnames(summarY) <- c("cType", "CTRL", "EAE")

write.table(summarY, paste0(saveImg, "cells_X_condition_kmeans.csv"), quote=F, sep="\t", row.names=F)

#' As a general rule, focusing on local neighborhoods provides the safest interpretation of t-SNE and UMAP plots. These methods spend considerable effort to ensure that each cell’s nearest neighbors in high-dimensional space are still its neighbors in the two-dimensional embedding. Thus, if we see multiple cell types or clusters in a single unbroken “island” in the embedding, we could infer that those populations were also close neighbors in higher-dimensional space. However, less can be said about the distances between non-neighboring cells; there is no guarantee that large distances are faithfully recapitulated in the embedding, given the distortions necessary for this type of dimensionality reduction. It would be courageous to use the distances between islands (measured, on occasion, with a ruler!) to make statements about the relative similarity of distinct cell types.  

#' [CLUSTERING](https://bioconductor.org/books/release/OSCA/clustering.html#clustering){target="_blank"}  

# use one of the three clustering files
# finally, store the normalised object with assigned cell type
saveRDS(object=denoisedSce, file=paste0(saveIn, "denoised_data_kmeans_", expN, ".Rds"))

#' SIDE NOTE
#' as mentioned here: http://bioconductor.org/books/release/OSCA/feature-selection.html#feature-selection-subsetting we decided to not subset the SCEexperiment using the TOP HGV, in order to be able to find other genes that may be of interested. In particular see the consideration n. 2: "We can keep the original SingleCellExperiment object and specify the genes to use for downstream functions via an extra argument like subset.row=. This is useful if the analysis uses multiple sets of HVGs at different steps, whereby one set of HVGs can be easily swapped for another in specific steps." 

# we can export the reducedSce using only HVGs
######reducedSce <- denoisedSce[rownames(denoisedSce) %in% highVar, ]
####### and export
######saveRDS(object=reducedSce, file=paste0(saveIn, "reduced_HVGs_data_", expN, ".Rds"))
