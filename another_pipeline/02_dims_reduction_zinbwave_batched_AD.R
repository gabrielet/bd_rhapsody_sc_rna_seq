#' ---
#' title: "Dimensionality reduction with zinbwave for batch corrected AD experiments"
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

# set AD paths
rawPathMale <- "/home/gabriele/cbmc/scrnaseq/raw_data_new/572_20/" ; expNMale <- "572" # AD male
rawPathFemale <- "/home/gabriele/cbmc/scrnaseq/raw_data_new/566_20/" ; expNFemale <- "566" # AD female
saveIn <- "/home/gabriele/cbmc/scrnaseq/raw_data_new/AD_batched/"

saveInMale <- paste0(rawPathMale, "results/")
saveInFemale <- paste0(rawPathFemale, "results/")

ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))

print(paste0("performing batch correction on exps ", expNMale, " and ", expNFemale))

# to perform batch correction it is possible to use the datasets that were generated after the dimensionality reduction and the clustering steps, as mentioned here:
# http://bioconductor.org/books/release/OSCA/integrating-datasets.html#setting-up-the-data
# Each dataset was obtained from the TENxPBMCData package and separately subjected to basic processing steps. Separate processing prior to the batch correction step is more convenient, scalable and (on occasion) more reliable. For example, outlier-based QC on the cells is more effective when performed within a batch. The same can also be said for trend fitting when modelling the mean-variance relationship

# In the example that is provided there are also PCA dimensionality reduction and clustering using graph-based knn

# declare useful variables
toZinb <- list()

# load clustered data
toZinb[["male"]] <- readRDS(paste0(saveInMale, file=paste0("QCed_data_", expNMale, ".Rds")))
toZinb[["female"]] <- readRDS(paste0(saveInFemale, file=paste0("QCed_data_", expNFemale, ".Rds")))

# labelling the batch	
toZinb[["male"]]$batch <- "male"
toZinb[["female"]]$batch <- "female"

# set palette with 12 colours for colouring plots
colorsPal <- brewer.pal(12, "Paired")
universe <- intersect(rownames(toZinb[["male"]]), rownames(toZinb[["female"]]))

# Subsetting the SingleCellExperiment object.
toZinb[["male"]] <- toZinb[["male"]][universe,]
toZinb[["female"]] <- toZinb[["female"]][universe,]

# prepare the data for the merging step, if all the rownames are the same across the two exps
if (all(rowData(toZinb[["male"]]) == rowData(toZinb[["female"]])) == T){
	rowData(toZinb[["male"]]) <- rowData(toZinb[["female"]])
} else {
	print("there's a problem with rownames")
}

# and merge
mergedZinb <- cbind(toZinb[["male"]], toZinb[["female"]])

#get best 1000 genes
perC <- 1000
topGMat <- as.matrix(assays(mergedZinb)$counts)
topGMat <- makeFilterStats(topGMat, filterStats="var", transFun=log1p)

# set perC
topGMat <- filterData(topGMat, percentile=perC, filterStats="var")

sceTop1000 <- mergedZinb[rownames(topGMat),]
topGMat$condition <- sceTop1000$condition
topGMat$pruned_main <- sceTop1000$pruned_main
topGMat$pruned_fine <- sceTop1000$pruned_fine
topGMat$batch <- sceTop1000$batch
assayNames(topGMat) <- "counts"

#and finally run zinbwave

# best parameter combination is perC = 1000 and kZinb = 20

# it may happen that a warning arise at the end of the computation stating that:
# Warning message: In simpute.als(x, J, thresh, lambda, maxit, trace.it, warm.start, : Convergence not achieved by 100 iterations
# this is not a problem, according to this page:
# https://github.com/drisso/zinbwave/issues/9

kZinb <- 20

set.seed(131)

# zinbwave automatically perform the batch correction, by passing as the two covariates both the condition AND the batch as
# X="~condition+batch"

#register(bpstart(MulticoreParam(workers=4)))
register(SerialParam())
print(system.time(zbWaved <- zinbwave(Y=topGMat, K=kZinb, X="~condition+batch", residuals=TRUE, normalizedValues=TRUE, verbose=TRUE)))

# from zinbwave manual:

#########################W <- reducedDim(zbWaved)

#########################tsne_data <- Rtsne(W, pca = FALSE, perplexity=30, max_iter=5000)

#########################data.frame(Dim1=tsne_data$Y[,1], Dim2=tsne_data$Y[,2], 
#########################   bio=colData(zbWaved)$condition) %>%
#########################ggplot(aes(Dim1, Dim2, colour=bio)) + geom_point() + 
#########################scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()

# this is the alternative
# assigning normalizedValues to logcounts may inflence other functions that use logcounts to perform any kind of operations
assays(zbWaved)$logcounts <- assays(zbWaved)$normalizedValues

# also, consider that using normalised values as logcounts will make sure these value will be used by findMarkers

# then run tsne
zbWaved <- runTSNE(zbWaved, use_dimred="zinbwave", perplexity=30)

# or by doing so, if logcounts of zbWaved are not set using normalizedValues
#runTSNE(sce, exprs_values = "normalizedValues", dimred = "zinbwave")

# and plot
p <- scater::plotTSNE(zbWaved, colour_by="condition")

png(paste0(saveIn, "02_01_zinbwave_tSNE_barched.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
plot(p)
dev.off()

# finally, export the data
saveRDS(object=zbWaved, file=paste0(saveIn, "zbwaved_data_batch_corrected_batched_K_", kZinb, "_top_", perC, ".Rds"))
