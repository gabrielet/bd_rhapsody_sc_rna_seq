#' ---
#' title: "K-means clustering for EAE size-factored data"
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

#' 2) performing [K-means clustering](https://bioconductor.org/books/release/OSCA/clustering.html#k-means-clustering){target="_blank"}  

# Load libraries
library("bluster")
library("cluster")
library("cowplot")
library("ggplot2")
library("gridExtra")
library("limma")
library("pheatmap")
library("plyr")
library("RColorBrewer")
library("scater")
library("scran")
library("Seurat")
library("SingleCellExperiment")

# set seed
set.seed(113)

# set paths and set kS parameter, which is dependent on the experiment

 rawPath <- "/home/gabriele/cbmc/scrnaseq/raw_data_new/588_20/" ; expN <- "588" ; kSExp <- 10 # EAE onset
# rawPath <- "/home/gabriele/cbmc/scrnaseq/raw_data_new/589_20/" ; expN <- "589" ; kSExp <- 15 # EAE chronic

saveIn <- paste0(rawPath, "results_bRanks/")
saveImg <- paste0(saveIn, "figures/")

#check if saveIn directory exist, otherwise create it
ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
ifelse(dir.exists(saveImg), TRUE, dir.create(saveImg, recursive=F))

print(paste0("analysing exp ", expN))

# load data
kmeansSce <- readRDS(paste0(saveIn, paste0("denoised_data_", expN, ".Rds")))

#' 4) A more practical use of k-means is to deliberately set k to a large value to achieve overclustering. This will forcibly partition cells inside broad clusters that do not have well-defined internal structure. For example, we might be interested in the change in expression from one “side” of a cluster to the other, but the lack of any clear separation within the cluster makes it difficult to separate with graph-based methods, even at the highest resolution. k-means has no such problems and will readily split these broad clusters for greater resolution, though obviously one must be prepared for the additional work involved in interpreting a greater number of clusters.  

# set an arbitrary kS, 15 seems to be best for separating clusters outside of the neutrophils one, for 
kS <- kSExp

# and an arbitrary seed
set.seed(131)

# new labels that compress the fine name into something similar to main labels
kmeansSce$fine_to_main <- unlist(lapply(strsplit(unlist(lapply(strsplit(kmeansSce$pruned_fine, "\\."), `[[`, 1)), " "), `[[`, 1))

clustKM <- kmeans(reducedDim(kmeansSce, "PCA"), centers=kS)
colLabels(kmeansSce) <- factor(clustKM$cluster)
kmeansSce$clusterS <- factor(clustKM$cluster)

# a summary of clusters content
tbl <- table(clustKM$cluster, kmeansSce$fine_to_main)

# plot by cluster
# plotting tSNE using cluster colours and cell label colours, side by side
png(paste0(saveImg, "03_01_tSNE_kmeans_K_", kS, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
grid.arrange(
	scater::plotTSNE(kmeansSce, colour_by="clusterS", text_by="clusterS") + theme(legend.position="none") + labs(title="By cluster"),
	scater::plotReducedDim(kmeansSce, "TSNE", colour_by="fine_to_main") + theme(legend.position="right") + labs(title="By cell type") + guides(color=guide_legend("Cell type")) ,
	ncol=2)
dev.off()

#' evaluate clusters separation: #https://bioconductor.org/books/release/OSCA/clustering.html#assessing-cluster-separation-1  

#' The within-cluster sum of squares (WCSS) for each cluster is the most relevant diagnostic for k-means, given that the algorithm aims to find a clustering that minimizes the WCSS  

ncells <- tabulate(clustKM$cluster)
tab <- data.frame(wcss=clustKM$withinss, ncells=ncells)
tab$rms <- sqrt(tab$wcss/tab$ncells)
tab

cent.tree <- hclust(dist(clustKM$centers), "ward.D2")
png(paste0(saveImg, "03_02_dendrogram_K_", kS, "_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	plot(cent.tree)
dev.off()

# finally, store the normalised objects
saveRDS(object=kmeansSce, file=paste0(saveIn, "kmeans_K_", kS, "_clustered_data_", expN, ".Rds"))
