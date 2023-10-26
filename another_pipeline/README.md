# SINGLE CELL RNA SEQUENCING USING BD RHAPSODY TECHNOLOGY AND ABSeq

### Step 1 - Demultiplexing with blc2fastq (on CPT High Performance Computing platform)

To run the demultiplexing, several file were generated. In particular there are four csv files summarising each RUN's characteristics.

Run_364.csv:

```
[Header],,,,
Experiment Name,RUN364,,,
Date,2020-11-09,,,
Module,GenerateFASTQ - 2.0.1,,,
Workflow,GenerateFASTQ,,,
Library Prep Kit,BD_Rhapsody,,,
Description,Constantin_BD_Rhapsody,,,
Chemistry,Default,,,
[Reads],,,,
76,,,,
76,,,,
[Settings],,,,
adapter,CTGTCTCTTATACACATCT,,,
[Data],,,,
Sample_ID,Description,I7_Index_ID,index,Sample_Project
588_20,WTA_SampleTag_AbSeq,I3,AAGAGGCA,
```

Run_363.csv:

```
[Header],,,,
Experiment Name,RUN363,,,
Date,2020-11-09,,,
Module,GenerateFASTQ - 2.0.1,,,
Workflow,GenerateFASTQ,,,
Library Prep Kit,BD_Rhapsody,,,
Description,Constantin_BD_Rhapsody,,,
Chemistry,Default,,,
[Reads],,,,
76,,,,
76,,,,
[Settings],,,,
adapter,CTGTCTCTTATACACATCT,,,
[Data],,,,
Sample_ID,Description,I7_Index_ID,index,Sample_Project
566_20,WTA_SampleTag_AbSeq,I1,GCTACGCT,
```

Run_366.csv:

```
[Header],,,,
Experiment Name,RUN366,,,
Date,2020-11-09,,,
Module,GenerateFASTQ - 2.0.1,,,
Workflow,GenerateFASTQ,,,
Library Prep Kit,BD_Rhapsody,,,
Description,Constantin_BD_Rhapsody,,,
Chemistry,Default,,,
[Reads],,,,
76,,,,
76,,,,
[Settings],,,,
adapter,CTGTCTCTTATACACATCT,,,
[Data],,,,
Sample_ID,Description,I7_Index_ID,index,Sample_Project
572_20,WTA_SampleTag_AbSeq,I2,CGAGGCTG,
```

Run_368.csv:

```
[Header],,,,
Experiment Name,RUN368,,,
Date,2020-11-09,,,
Module,GenerateFASTQ - 2.0.1,,,
Workflow,GenerateFASTQ,,,
Library Prep Kit,BD_Rhapsody,,,
Description,Constantin_BD_Rhapsody,,,
Chemistry,Default,,,
[Reads],,,,
76,,,,
76,,,,
[Settings],,,,
adapter,CTGTCTCTTATACACATCT,,,
[Data],,,,
Sample_ID,Description,I7_Index_ID,index,Sample_Project
589_20,WTA_SampleTag_AbSeq,I4,GTAGAGGA,
```

Also, for each experiments a file was generated and called **demulti_*.sh** such as:

```sh
#!/bin/bash
bcl2fastq --runfolder-dir /home/gabriele.tosadori/sequenziamento/200701_NB500897_0363_AH5GCNBGXG/ --output-dir /home/gabriele.tosadori/sequenziamento/200701_NB500897_0363_AH5GCNBGXG_base_calls/ --sample-sheet RUN363.csv --no-lane-splitting
```

```
#!/bin/bash
bcl2fastq --runfolder-dir /home/gabriele.tosadori/sequenziamento/200702_NB500897_0364_AH5GGLBGXG/ --output-dir /home/gabriele.tosadori/sequenziamento/200702_NB500897_0364_AH5GGLBGXG_base_calls/ --sample-sheet RUN364.csv --no-lane-splitting
```

```
#!/bin/bash
bcl2fastq --runfolder-dir /home/gabriele.tosadori/sequenziamento/200701_NB500897_0363_AH5GCNBGXG/ --output-dir /home/gabriele.tosadori/sequenziamento/200701_NB500897_0363_AH5GCNBGXG_base_calls/ --sample-sheet RUN363.csv --no-lane-splitting
```

```
#!/bin/bash
bcl2fastq --runfolder-dir /home/gabriele.tosadori/sequenziamento/200706_NB500897_0366_AH5GKMBGXG/ --output-dir /home/gabriele.tosadori/sequenziamento/200706_NB500897_0366_AH5GKMBGXG_base_calls/ --sample-sheet RUN366.csv --no-lane-splitting
```

All of them were run using SLURM, by these script that were run from the command line of the HPC:

```
srun -c 8 --mem-per-cpu=4096 --output=demulti_200701.log --priority high demulti_200701.sh

srun -c 8 --mem-per-cpu=4096 --output=demulti_200702.log --priority high demulti_200702.sh

srun -c 8 --mem-per-cpu=4096 --output=demulti_200706.log --priority high demulti_200706.sh

srun -c 8 --mem-per-cpu=4096 --output=demulti_200713.log --priority high demulti_200713.sh
```

### Step 2 - Alignment using SevenBridges platform by BD

Once the demultiplexing step is performed, for each run 2 FASTQ files were generated. These FASTQ files contain information about the reads and their quality. All the reads, at this point, must be aligned to a reference genome which, in our case, is the GRCm39 downloaded from [this link at NCBI website](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27). From this same link, the transcriptome annotations were downloaded too, in form of **gtf** file.

```
Description:     Genome Reference Consortium Mouse Build 39
Organism name:    Mus musculus (house mouse)
Infraspecific name:    Strain: C57BL/6J
BioProject:    PRJNA20689
Submitter:    Genome Reference Consortium
Date:    2020/06/24
Synonyms:    mm39
Assembly level:    Chromosome
Genome representation:    full
RefSeq category:    reference genome
GenBank assembly accession:    GCA_000001635.9 (latest)
RefSeq assembly accession:    GCF_000001635.27 (latest)
RefSeq assembly and GenBank assembly identical:    yes
```

This genome, since it is not provided by BD, went through a preparation step performed by BD staff.

Using the SevenBridges platform, a Whole Transcriptome Analysis was run. Several files were used:

 - FASTQ files from the bcl2fastq step;
 - Reference genome;
 - PhiX phage sequences;
 - ABSeq references;
 - Transcriptome annotations for GRCm39.
 
SevenBridges provided several file to evaluate the goodness of the alignment, to understand how many reads and how many cells were meet the quality criteria and other information.

Also, the count matrices representing the number of transcripts found for each gene and each cell were generated and used for the Step 3.

### Step 3 - Cleaning and merging the samples using R

The first step of the R pipeline concerns loading the count matrices using Seurat, and converting them to be used as SingleCellExperiments.

Then some quality-control (QC) checks are performed and some low quality cells may be removed.

It is important to note that cells containing mitochondrial genes were removed. To achieve this goal, mitochondrial genes were retrieved using [a list available from the NCBI website](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27).

Also ABSeq data were extracted:

 - as they come;
 - after normalising them using centered-log-ratio.
 
To generate the figures that are mentioned in Step 5, the not normalised ABSeq data were used.

###### OUTPUTS:
 
 - raw ABseq data: **abseq_data.Rds**;
 - centered log ratio ABseq data: **centeredAB_data.Rds**;
 - raw data to be used for zinbwave: **QCed_to_zinbwave.Rds**;
 - normalised data to be used for PCA: **QCed_data.Rds**.

### Step 4 - Dimensionality reduction and clustering using R

Step 4 is performed by using Principal Component Analysis, when EAE data are analysed, or Zinbwave, when AD data were analysed. This difference is due to the fact that in AD, retrieving Gamma Delta T Cells required to use a more complex dimensionality reduction algorithm.

To take advantage of Zinbwave, only the 1000 most variable genes were considered. This was done in order to actually compute the dimensionality reduction since Zinbwave is computationally expensive and may take a very long time to provide an output.

Then clustering was performed. For EAE k-means was used, with K set at 20. This parameter was set after multiple tests and trials and it looks like it is the best solution. For AD, a **slm** clustering from the Seurat package was used. Here, different resolutions were tested, i.e. res=0, 1, 1.5. Eventually, 1.5 was chosen for the results.

###### OUTPUTS for dimensionality reduction:

 - from size factors analysis:
	 + exporting the reduced data with all the genes: **denoised_data.Rds**;
	 + exporting the reduced data with highly variable genes only: **reduced_HVGs_data.Rds**.
  
 - from zinbwave analyis:
	 + zinbwave reduced data with K=20 and 1000 most variable genes: **zinbwaved_data.Rds**.

###### OUTPUTS for clustering:

 - **kmeans_clustered_data.Rds** for size factors normalised data;
 - **kmeans_HVGs_clustered_data.Rds** for size factors with Highly Variable Genes only;
 - **slm_clustered_data_K_20_top_1000.Rds** for zinbwaved data and slm clustering;
 - **slm_clustered_data_batched_K_20_top_1000.Rds** for zinbwaved and batched experiments, and slm clustering.


### Step 5 - Finding cell types and their marker genes using R

Each cell type may require different data processing, e.g. different dimensionality reduction i.e. PCA, Zinbwave, and clustering, in order to be properly clustered.

Generally speaking, this step is divided in three parts that allowed analysing a specific cell type in three different settings:
 
 1. analyse each cluster that contain a specific cellType. For each cluster only cells of that specific cell type are considered;
 2. analyse each cluster that contain a specific cellType. For each cluster all the cells in the cluster are considered, no matter what their cellType is;
 3. analyse all and only cellType cells, without considering the clusters they may belong to.
 
For each one of these setting, a list of differentially expressed genes (DEGs) is obtained, if there are any. These lists are called **UP**, **DOWN**, and **up\_and\_Down** and are stored into the corresponding cellType **UpDown** directory.

It is possible that some cellType cells subset is found for a specific condition only, e.g. only for WT, making it impossible to obtain any DEG. Also, it may happen that there are too few cells hence the comparison between 3xTG vs WT or EAE vs CTRL cannot be performed. The threshold to determine if there are not enough cells is set to 25 for each of the two condition. This means that, to perform a comparison between EAE and CTRL there must be at least 50 cells, 25 belonging to each of the two conditions.

In order to find DEGs, all the three functions needs a **cTSubSet**, which is also used for a later GSEA analysis, and two parameters were set. One is the False Discovery Rate (FDR) that allows to obtain stronger candidates and is preferred to the p-value. **FDR is set at 0.05**. The other parameter is used to determine the minimum amount of cells that are allowed for each condition, e.g. WT or EAE, in order to perform the comparison. This value, called **personalT, is set to 25**.

For each set of found DEGs three tables are returned are stored in **expCode/results/enrichment\_cellType/UpDown/**. These tables represent Up, Down and Up and Down together DEGs, ordered by their fold-change:

 - Up from most positive to 0;
 - Down from more negative to 0;
 - Up and Down from most positive to most negative.

<!-- Then, a GSEA analysis, performed using the fgsea() algorithm is performed. GSEA is performed only on cells of a specific cellType for all the clusters in which they are found and for all the cells of a specific cellType found in all the dataset.

The output parameter pAdj, i.e. p-value adjusted, is a BH-adjusted p-value according to (the manual)[http://www.bioconductor.org/packages/release/bioc/manuals/fgsea/man/fgsea.pdf]. This parameter is used to order the resulting tables.

To perform GSEA, all the genes in the cTSubSet are previously ranked using DESeq2 statistics then the ranked DEGs are fed to the enrichment algorithm.

Also, two parameters are set (see link above):

 - minSizeMinimal size of a gene set to test.  All pathways below the threshold are excluded;
 - maxSizeMaximal size of a gene set to test.  All pathways above the threshold are excluded.

Two different GSEA analysis are performed: Reactome and C7 from the MSig database, compiled for *Mus musculus*.

The C7 dataset was downloaded from http://bioinf.wehi.edu.au/MSigDB/ and used by the enrichment function. It is necessary to download **Mm.c7.all.v7.1.entrez.rds**, to run the analysis.

Since C7 enrichment names refer to published research articles, the abstract of each found term is downloaded, to enable a complete review of all the retrieved information. Abstracts are found in the ***expCode/results/enrichment\_cellType/abstract*** directory.

For each set of DEGs **the ten best**, or less if there are less, and **all the terms** are returned, as two separate tables called **enrich\_ten\_best** and **enrich\_all** , respectively. These files are found in the ***expCode/results/enrichment_cellType/GSEA/*** directory. It contains also and *enrichmentScore* showing the ES obtained by the GSEA for each term. -->

The ranking of the genes, performed with DESeq2 as mentioned above and supposed to be used for GSEA, is returned to perform further analyses. These rankings are stored into the **ranks directory** and are used in the next step.

Some plots and a table are also provided and stored into ***expCode/results/figures/analysis\_cellType***:

- Cluster Matrix representing the content in cell types of the whole dataset;
- tSNE with clusters and cell types;
- tSNE showing only a specific cell type (in orange);
- two cluster composition plots representing the cell types, the number of that cell type in a specific cluster and two histograms one up and one right. The up histogram represents the percentage of a cell type for that cluster. The right histogram represent the percentage of a cell type for the dataset;
- a table called **_cluster_content_** representing the actual numbers of cells for each cluster.

**ABSeq analysis** was performed and ABSeq plots, one for each cluster, were generated (and located in the same **cellType figures directory** as above).

These plots represent the contribution of each antibody for each cell type found in a specific cluster. Ideally, the histograms for a certain cellType should be similar (if ImmGen misclassified some cells). If the histograms vary, then it might be that the clustering put together cells, in the same cluster, that are actually different.

To analyse ABseq data this goal a file listing the name of the antibodies was used and is called **abseqlist.txt**.

###### OUTPUTS:
 
 - rankings found in the **rank directory**;
 - up and down regulated genes for cellType cells, stored into the **UpDown directory**;
 - some plots.

### Step 6 - Analysis and enrichment using R

In this step it is fundamental to remember from where the genes lists are coming from. In particular there are two analysis:
 - FindMarkers() used to generate **UP**, **DOWN**, and **up\_and\_down** files that are later enriched using EnrichMyData()
 - DeSeq2 that is used inside the gseaMyData() function, that generates the **ranking\_genes\_for** files, that are used to GO, KEGG, and REACTOMEPlotting

The enrichment analysis is performed in two steps. The first step takes advantage of ClusterProfiler to perform three different analysis: a GO molecular function, GO biological process, <!-- and a KEGG pathway analysis --> and a Reactome pathway analysis. All of them are performed on all the sets of UpAndDown DEGs that were retrieved in Step 5. For each analysis a table called _clusterProfiler*_ is generated and stored in **_expCode/results/enrichment_cellType/cProfiler/_**.

The tables for Up and Down, obtained in the previous step, are enriched using clusterProfiler and the results are listed as a set of tables which name is **clusterProfiler_something**. The results are stored into the **cProfiler directory** that is divided in several subdirectories depending on the enrichment that was performed:

 - GO biological processes;
 - GO molecular functions;
 - Reactome Pathways.
 
Each directory may contain several tables depending on the files that were genereated and stored in the **UpDown directory**.

In this step, all the **up\_and\_down** lists are collected to determine the shared genes. A list of shared genes is returned. Also:

 - an **upset plot**;
 - a **heatmap**.
 
describing how exactly the intersections occur, are returned.

After this first step, another step is performed involving the plotting of many figures, one for each of these enrichment: <!-- KEGG, --> GO and Reactome.

To perform this further analysis, the ranking files are used. The rankings were generated using DESeq2 during the GSEA analysis in Step 5. Finally, for each enrichment, a table is exported showing all the results. They are stored in **\_expCode/results/enrichment_cellType/cProfiler/\_**. These files have the label **cProf\_DESEQ** and are obtained using the **ranking files** obtained in the previous step.

Each ranking is analysed and (several plots)[https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html] are produced, for each enrichment, i.e. <!-- KEGG -->GO and Reactome.

The plots, described in the link above, are slightly modified to make them more readable. This is the full list of the generated plots:

- **Dot plot** depicts the enrichment scores, p-adjusted, gene count and ratio as dot color, dimension, and position on the x axis;
- **Cnetplot** depicts the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) as a network;
- **Heatplot** is similar to cnetplot, while displaying the relationships as a heatmap. The gene-concept network may become too complicated if user want to show a large number significant terms. The heatplot can simplify the result and more easy to identify expression patterns;
- **Enrichment map** organizes enriched terms into a network with edges connecting overlapping gene sets. In this way, mutually overlapping gene sets are tend to cluster together, making it easy to identify functional module.

All these figures are stored in two directories that are **GO\_plots** and **Reactome\_plots**, found in each cellType's analysis directory.

###### OUTPUTS:
 
 - the figures described above.

### Summarising the full R pipeline

To perform the actual R analysis, there are many files that must be run.

  - **00\_functions.R**: this file contains the functions that are used throughout the pipeline, in particular in 04 and 05 steps;
 - **01\_cleaning\_AD.R**: clean AD data and get them ready to be processed by Zinbwave;
 - **01\_cleaning\_EAE.R**: clean EAE data and get them ready to be processed by PCA;
 - **02\_dims\_reduction\_sizeFactors.R**: dimensionality reduction using PCA, mainly used for EAE data but can be used for AD too, if they are processed using size factors;
 - **02\_dims\_reduction\_zinbwave.R**: dimensionality reduction using Zinbwave, mainly used for AD data but can be used for EAE too, if they are not processed using size factors;
 - **02\_dims\_reduction\_zinbwave\_batched\_AD.R**: dimensionality reduction using Zinbwave, for AD data;
 - **02\_dims\_reduction\_zinbwave\_batched\_EAE.R**: dimensionality reduction using Zinbwave, for EAE data;
 - **03\_kmeans\_clustering.R**: clustering for EAE using kmeans, but works for AD too after PCA, using the full dataset;
 - **03\_kmeans\_reducedHVG\_clustering.R**: clustering for EAE using kmeans, but works for AD too after PCA, using the 1000 most variable genes;
 - **03\_slm\_clustering.R**: clustering for AD using slm algorithm, but works for EAE too after Zinbwave;
 - **03\_slm\_clustering\_batched.R**: clustering for AD and EAE batched data using slm algorithm;
 - **04\_find\_AD\_Bcells\_zinbwave.R**: analysing Bcells for AD data;
 - **04\_find\_AD\_CD4\_zinbwave.R**;
 - **04\_find\_AD\_CD8\_zinbwave.R**;
 - **04\_find\_AD\_Neutrophils\_zinbwave.R**;
 - **04\_find\_AD\_Tgd\_zinbwave.R**;
 - **04\_find\_batched\_AD.R**: analysing all cell types for batched AD;
 - **04\_find\_batched\_EAE.R**: analysing Neutrophils cells for batched EAE;
 - **04\_find\_EAE\_Neutrophils\_PCA.R** : analysing Neutrophils cells for EAE data;
 - **05\_enrichment\_analysis\_AD.R**: performing enrichment with ClusterProfiler and generating several, different plots;
 - **05\_enrichment\_analysis\_batched\_AD.R**;
 - **05\_enrichment\_analysis\_batched\_EAE.R**;
 - **05\_enrichment\_analysis\_EAE.R**;
 - **06\_network\_analysis.R**: performing network analysis using DEGs.
 
**NOTE**

To run Step 4, run all the scripts through the command line. The fastest way is in parallel, but serial is ok too (it takes way more time). These commands also remove unwanted characters.

```
# AD

Rscript 04_find_AD_CD4_zinbwave.R | sed 's#\[1\] ##g' | sed 's/"//g' > text_res/results_AD_CD4.txt
Rscript 04_find_AD_Bcells_zinbwave.R | sed 's#\[1\] ##g' | sed 's/"//g' > text_res/results_AD_Bcells.txt
Rscript 04_find_AD_CD8_zinbwave.R | sed 's#\[1\] ##g' | sed 's/"//g' > text_res/results_AD_CD8.txt
Rscript 04_find_AD_Neutrophils_zinbwave.R | sed 's#\[1\] ##g' | sed 's/"//g' > text_res/results_AD_Neutrophils.txt
Rscript 04_find_AD_Tgd_zinbwave.R | sed 's#\[1\] ##g' | sed 's/"//g' > text_res/results_AD_Tgd.txt

# EAE
Rscript 04_find_EAE_Neutrophils_PCA.R | sed 's#\[1\] ##g' | sed 's/"//g' > text_res/results_Neutrophils_EAE.txt

# BATCHED
Rscript 04_find_batched_AD.R | sed 's#\[1\] ##g' | sed 's/"//g' > text_res/results_batched_AD.txt
Rscript 04_find_batched_EAE.R | sed 's#\[1\] ##g' | sed 's/"//g' > text_res/results_batched_EAE.txt
Rscript 04_find_batched_AD_CD4_CD8.R | sed 's#\[1\] ##g' | sed 's/"//g' > text_res/results_batched_AD_CD4_CD8.txt
```

Then, run 05 in the same way:

```
Rscript 05_enrichment_analysis_AD.R
Rscript 05_enrichment_analysis_batched_AD.R
Rscript 05_enrichment_analysis_EAE.R
Rscript 05_enrichment_analysis_batched_EAE.R
```

### The full pipeline should be:

```
#################################		STEP ONE

01_cleaning_AD.R	exp 566	+
01_cleaning_AD.R	exp 572	+

01_cleaning_EAE.R	exp 588	+
01_cleaning_EAE.R	exp 589	+

01_cleaning_AD_10x.R

#################################		STEP TWO

02_dims_reduction_sizeFactors.R	exp 588
02_dims_reduction_sizeFactors.R	exp 589

02_dims_reduction_zinbwave.R	exp 566 +
02_dims_reduction_zinbwave.R	exp 572 +

02_dims_reduction_zinbwave_batched_AD.R

02_dims_reduction_zinbwave_batched_EAE.R

02_dims_reduction_zinbwave_10x.R

#################################		STEP THREE

03_kmeans_clustering.R

03_slm_clustering.R

03_slm_clustering_batched.R

03_slm_clustering_10x.R

#################################		STEP FOUR

04_find_10x_zinbwave.R

04_find_AD_Bcells_zinbwave.R

04_find_AD_CD4_zinbwave.R

04_find_AD_CD8_zinbwave.R

04_find_AD_Neutrophils_zinbwave.R

04_find_AD_Tgd_zinbwave.R

04_find_EAE_Neutrophils_PCA.R

04_find_batched_AD.R

04_find_batched_AD_CD4_CD8.R

04_find_batched_EAE.R

#################################		STEP FIVE

05_enrichment_analysis_10x.R

05_enrichment_analysis_AD.R

05_enrichment_analysis_EAE.R

05_enrichment_analysis_batched_AD.R

05_enrichment_analysis_batched_EAE.R

#################################		STEP SIX

06_network_analysis. R
```


### R Version and Libraries

```
R version 4.0.4 (2021-02-15) -- "Lost Library Book"
Copyright (C) 2021 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)
```

 - zinbwave_1.12.0
 - UpSetR_1.4.0               
 - stringr_1.4.0
 - SingleR_1.4.0              
 - Seurat_3.2.3
 - scran_1.18.3               
 - scater_1.18.3
 - ReactomePA_1.34.0          
 - RColorBrewer_1.1-2
 - plyr_1.8.6                 
 - pheatmap_1.0.12
 - pathview_1.30.1            
 - org.Mm.eg.db_3.12.0
 - openxlsx_4.2.3
 - limma_3.46.0               
 - HelpersMG_4.4
 - gridExtra_2.3              
 - ggplot2_3.3.3
 - gg.gap_1.3                 
 - fgsea_1.16.0
 - enrichplot_1.10.1          
 - DropletUtils_1.10.2
 - DOSE_3.16.0                
 - DESeq2_1.30.0
 - ddpcr_1.15                 
 - cowplot_1.1.1
 - compositions_2.0-1         
 - clusterProfiler_3.18.0
 - clusterExperiment_2.10.0   
 - SingleCellExperiment_1.12.0
 - cluster_2.1.1              
 - celldex_1.0.0
 - bluster_1.0.0
 - BiocParallel_1.24.1