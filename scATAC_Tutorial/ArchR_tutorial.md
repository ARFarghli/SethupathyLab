scATAC Tutorial
================

``` r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
ArchR::installExtraPackages()
```

``` r
#set threads specific to your machine
addArchRThreads(threads = 35) 
```

    ## Setting default number of Parallel threads to 35.

We need a reference genome for downstream analyses. ArchR natively
supports hg19, hg38, mm9, and mm10.

``` r
inputFiles <- getTutorialData(tutorial = "Hematopoiesis")
inputFiles
```

    ##                                       scATAC_BMMC_R1 
    ##      "HemeFragments/scATAC_BMMC_R1.fragments.tsv.gz" 
    ##                                  scATAC_CD34_BMMC_R1 
    ## "HemeFragments/scATAC_CD34_BMMC_R1.fragments.tsv.gz" 
    ##                                       scATAC_PBMC_R1 
    ##      "HemeFragments/scATAC_PBMC_R1.fragments.tsv.gz"

``` r
addArchRGenome(genome = "hg19", install = FALSE)
```

    ## Setting default genome to Hg19.

### Creating Arrow Files

ArchR uses arrow files, the base unit of an ArchR analytical project.
Every arrow file stores all the data associated with an individual
sample (i.e. a single replicate of a particular condition). During
creation and as additional analyses are performed, ArchR updates and
edits each Arrow file to contain additional layers of information. IAn
Arrow file is actually just a path to an external file stored on disk.
More explicitly, an Arrow file is not an R-language object that is
stored in memory but rather an HDF5-format file stored on disk.

``` r
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
```

    ## Using GeneAnnotation set by addArchRGenome(Hg19)!
    ## Using GeneAnnotation set by addArchRGenome(Hg19)!

    ## ArchR logging to : ArchRLogs/ArchR-createArrows-25d241aa0385-Date-2021-05-05_Time-15-22-12.log
    ## If there is an issue, please report to github with logFile!

    ## Cleaning Temporary Files

    ## 2021-05-05 15:22:12 : Batch Execution w/ safelapply!, 0 mins elapsed.

    ## Warning in mclapply(..., mc.cores = threads, mc.preschedule = preschedule): 1
    ## function calls resulted in an error

    ## createArrowFiles has encountered an error, checking if any ArrowFiles completed..

    ## 2021-05-05 15:22:18 :

    ## ArchR logging successful to : ArchRLogs/ArchR-createArrows-25d241aa0385-Date-2021-05-05_Time-15-22-12.log

``` r
#We can inspect the ArrowFiles object to see that it is actually just a character vector of Arrow file paths.
ArrowFiles
```

    ## [1] "scATAC_PBMC_R1.arrow"      "scATAC_CD34_BMMC_R1.arrow"

Quality Control of scATAC-seq data is essential to remove cells that
contribute to low-quality data. There are three characteristics that
scATAC considers 1. The number of unique nuclear fragments (i.e. not
mapping to mitochondrial DNA). 2. The signal-to-background ratio. Low
signal-to-background ratio is often attributed to dead or dying cells
which have de-chromatinzed DNA which allows for random transposition
genome-wide. 3. The fragment size distribution. Due to nucleosomal
periodicity, we expect to see depletion of fragments that are the length
of DNA wrapped around a nucleosome (approximately 147 bp).

ArchR also infers doublets, a single droplet that contains multiple
cells.

``` r
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)
```

    ## ArchR logging to : ArchRLogs/ArchR-addDoubletScores-25d22c37ea28-Date-2021-05-05_Time-15-22-18.log
    ## If there is an issue, please report to github with logFile!

    ## 2021-05-05 15:22:18 : Batch Execution w/ safelapply!, 0 mins elapsed.

    ## 2021-05-05 15:22:18 : scATAC_PBMC_R1 (1 of 2) :  Computing Doublet Statistics, 0 mins elapsed.

    ## scATAC_PBMC_R1 (1 of 2) : UMAP Projection R^2 = 0.99262

    ## scATAC_PBMC_R1 (1 of 2) : UMAP Projection R^2 = 0.99262

    ## 2021-05-05 15:23:45 : scATAC_CD34_BMMC_R1 (2 of 2) :  Computing Doublet Statistics, 1.442 mins elapsed.

    ## scATAC_CD34_BMMC_R1 (2 of 2) : UMAP Projection R^2 = 0.98225

    ## scATAC_CD34_BMMC_R1 (2 of 2) : UMAP Projection R^2 = 0.98225

    ## ArchR logging successful to : ArchRLogs/ArchR-addDoubletScores-25d22c37ea28-Date-2021-05-05_Time-15-22-18.log

\#\#\#Creating an ArchRProject

``` r
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)
```

    ## Using GeneAnnotation set by addArchRGenome(Hg19)!
    ## Using GeneAnnotation set by addArchRGenome(Hg19)!

    ## Validating Arrows...

    ## Getting SampleNames...

    ## 

    ## Copying ArrowFiles to Ouptut Directory! If you want to save disk space set copyArrows = FALSE

    ## 1 2 
    ## Getting Cell Metadata...
    ## 
    ## Merging Cell Metadata...
    ## Initializing ArchRProject...
    ## 
    ##                                                    / |
    ##                                                  /    \
    ##             .                                  /      |.
    ##             \\\                              /        |.
    ##               \\\                          /           `|.
    ##                 \\\                      /              |.
    ##                   \                    /                |\
    ##                   \\#####\           /                  ||
    ##                 ==###########>      /                   ||
    ##                  \\##==......\    /                     ||
    ##             ______ =       =|__ /__                     ||      \\\
    ##         ,--' ,----`-,__ ___/'  --,-`-===================##========>
    ##        \               '        ##_______ _____ ,--,__,=##,__   ///
    ##         ,    __==    ___,-,__,--'#'  ==='      `-'    | ##,-/
    ##         -,____,---'       \\####\\________________,--\\_##,/
    ##            ___      .______        ______  __    __  .______      
    ##           /   \     |   _  \      /      ||  |  |  | |   _  \     
    ##          /  ^  \    |  |_)  |    |  ,----'|  |__|  | |  |_)  |    
    ##         /  /_\  \   |      /     |  |     |   __   | |      /     
    ##        /  _____  \  |  |\  \\___ |  `----.|  |  |  | |  |\  \\___.
    ##       /__/     \__\ | _| `._____| \______||__|  |__| | _| `._____|
    ## 

``` r
#We can query which data matrices are available in the ArchRProject. At this point in time, we should have “GeneScoreMatrix” and “TileMatrix”. As we continue to work and add to the ArchRProject, we can use the following function to query which matricies are added to the project.
getAvailableMatrices(proj)
```

    ## [1] "GeneScoreMatrix" "TileMatrix"

``` r
#Next we can filter out putative doublets based on the scores established in the `infer doublets` chunk. Importantly, this does not delete the data from the Arrow files, but rather forces ArchRProject to ignore these cells. 
proj <- filterDoublets(ArchRProj = proj)
```

    ## Filtering 139 cells from ArchRProject!
    ##  scATAC_PBMC_R1 : 60 of 2453 (2.4%)
    ##  scATAC_CD34_BMMC_R1 : 79 of 2818 (2.8%)

\#\#\#Dimensionality Reduction and Clustering

``` r
#ArchR implements an iterative LSI dimensionality reduction via the addIterativeLSI() function.
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
```

    ## Checking Inputs...

    ## ArchR logging to : ArchRLogs/ArchR-addIterativeLSI-25d27f1b49be-Date-2021-05-05_Time-15-24-57.log
    ## If there is an issue, please report to github with logFile!

    ## 2021-05-05 15:24:57 : Computing Total Across All Features, 0.002 mins elapsed.

    ## 2021-05-05 15:25:00 : Computing Top Features, 0.04 mins elapsed.

    ## ###########
    ## 2021-05-05 15:25:01 : Running LSI (1 of 2) on Top Features, 0.061 mins elapsed.
    ## ###########

    ## 2021-05-05 15:25:01 : Creating Partial Matrix, 0.061 mins elapsed.

    ## 2021-05-05 15:25:05 : Computing LSI, 0.129 mins elapsed.

    ## 2021-05-05 15:25:22 : Identifying Clusters, 0.417 mins elapsed.

    ## 2021-05-05 15:25:32 : Identified 6 Clusters, 0.581 mins elapsed.

    ## 2021-05-05 15:25:32 : Saving LSI Iteration, 0.581 mins elapsed.

    ## 2021-05-05 15:25:42 : Creating Cluster Matrix on the total Group Features, 0.745 mins elapsed.

    ## 2021-05-05 15:25:48 : Computing Variable Features, 0.849 mins elapsed.

    ## ###########
    ## 2021-05-05 15:25:48 : Running LSI (2 of 2) on Variable Features, 0.853 mins elapsed.
    ## ###########

    ## 2021-05-05 15:25:48 : Creating Partial Matrix, 0.853 mins elapsed.

    ## 2021-05-05 15:25:54 : Computing LSI, 0.94 mins elapsed.

    ## 2021-05-05 15:26:10 : Finished Running IterativeLSI, 1.221 mins elapsed.

``` r
#To call clusters in this reduced dimension sub-space, we use the addClusters() function which uses Seurat’s graph clustering as the default clustering method.
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
```

    ## ArchR logging to : ArchRLogs/ArchR-addClusters-25d22a5b6bc4-Date-2021-05-05_Time-15-26-10.log
    ## If there is an issue, please report to github with logFile!

    ## 2021-05-05 15:26:11 : Running Seurats FindClusters (Stuart et al. Cell 2019), 0.002 mins elapsed.

    ## Computing nearest neighbor graph

    ## Computing SNN

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 5132
    ## Number of edges: 314076
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8183
    ## Number of communities: 9
    ## Elapsed time: 0 seconds

    ## 2021-05-05 15:26:20 : Testing Outlier Clusters, 0.158 mins elapsed.

    ## 2021-05-05 15:26:20 : Assigning Cluster Names to 9 Clusters, 0.158 mins elapsed.

    ## 2021-05-05 15:26:20 : Finished addClusters, 0.159 mins elapsed.

``` r
#We can visualize our scATAC-seq data using a 2-dimensional representation such as Uniform Manifold Approximation and Projection (UMAP). To do this, we add a UMAP embedding to our ArchRProject object with the addUMAP() function. This function uses the uwot package to perform UMAP.
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")
```

    ## 15:26:20 UMAP embedding parameters a = 0.7669 b = 1.223

    ## 15:26:20 Read 5132 rows and found 30 numeric columns

    ## 15:26:20 Using Annoy for neighbor search, n_neighbors = 40

    ## 15:26:20 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 15:26:21 Writing NN index file to temp file /tmp/RtmpblA0Am/file25d2148a8c03
    ## 15:26:21 Searching Annoy index using 20 threads, search_k = 4000
    ## 15:26:21 Annoy recall = 100%
    ## 15:26:22 Commencing smooth kNN distance calibration using 20 threads
    ## 15:26:23 Found 2 connected components, falling back to 'spca' initialization with init_sdev = 1
    ## 15:26:23 Initializing from PCA
    ## 15:26:23 PCA: 2 components explained 66.23% variance
    ## 15:26:23 Commencing optimization for 500 epochs, with 310616 positive edges
    ## 15:26:40 Optimization finished
    ## 15:26:40 Creating temp model dir /tmp/RtmpblA0Am/dir25d269f4337e
    ## 15:26:40 Creating dir /tmp/RtmpblA0Am/dir25d269f4337e
    ## 15:26:40 Changing to /tmp/RtmpblA0Am/dir25d269f4337e
    ## 15:26:40 Creating /home/af547/Git/SethupathyLab/scATAC_Tutorial/HemeTutorial/Embeddings/Save-Uwot-UMAP-Params-IterativeLSI-25d236ee0a2a-Date-2021-05-05_Time-15-26-40.tar

``` r
#We can visualize the UMAP in a number of ways by calling various attributes of the cells stored in the `cellColData` matrix. Here, we can visualize the UMAP by sample, or clusters.
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
```

    ## ArchR logging to : ArchRLogs/ArchR-plotEmbedding-25d2224770e4-Date-2021-05-05_Time-15-26-40.log
    ## If there is an issue, please report to github with logFile!
    ## Getting UMAP Embedding
    ## ColorBy = cellColData
    ## Plotting Embedding
    ## 1 
    ## ArchR logging successful to : ArchRLogs/ArchR-plotEmbedding-25d2224770e4-Date-2021-05-05_Time-15-26-40.log

``` r
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
```

    ## ArchR logging to : ArchRLogs/ArchR-plotEmbedding-25d2a341300-Date-2021-05-05_Time-15-26-41.log
    ## If there is an issue, please report to github with logFile!
    ## Getting UMAP Embedding
    ## ColorBy = cellColData
    ## Plotting Embedding
    ## 1 
    ## ArchR logging successful to : ArchRLogs/ArchR-plotEmbedding-25d2a341300-Date-2021-05-05_Time-15-26-41.log

``` r
ggAlignPlots(p1, p2, type = "h")
```

![](ArchR_tutorial_files/figure-gfm/Dimensionality%20Reduction,%20UMAP%20Generation%20and%20Visualization-1.png)<!-- -->

``` r
#To save an editable vectorized version of this plot, we use the plotPDF() function.
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
```

    ## Plotting Ggplot!
    ## Plotting Ggplot!

\#\#\#Assigning Clusters with Gene Scores

The novelty of single cell approaches is to be able to resolve cellular
heterogeneity in complex tissues. We can identify cells population by
assigning cell-type specific markers to them.

``` r
#First, we add imputation weights using MAGIC **read up on MAGIC** to help smooth the dropout noise in our gene scores
proj <- addImputeWeights(proj)
```

    ## ArchR logging to : ArchRLogs/ArchR-addImputeWeights-25d229d72965-Date-2021-05-05_Time-15-26-48.log
    ## If there is an issue, please report to github with logFile!

    ## 2021-05-05 15:26:48 : Computing Impute Weights Using Magic (Cell 2018), 0 mins elapsed.

``` r
#Now we can overlay our marker gene scores on our 2D UMAP embedding.
markerGenes  <- c(
    "CD34",  #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", "MME", #B-Cell Trajectory
    "CD14", "MPO", #Monocytes
    "CD3D", "CD8A"#TCells
  )

p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)
```

    ## Getting ImputeWeights

    ## ArchR logging to : ArchRLogs/ArchR-plotEmbedding-25d22a1110f7-Date-2021-05-05_Time-15-26-53.log
    ## If there is an issue, please report to github with logFile!

    ## Getting UMAP Embedding

    ## ColorBy = GeneScoreMatrix

    ## Getting Matrix Values...

    ## 2021-05-05 15:26:54 :

    ## 

    ## Imputing Matrix

    ## Using weights on disk

    ## 1 of 1

    ## Using weights on disk

    ## 1 of 1

    ## Plotting Embedding

    ## 1 2 3 4 5 6 7 8 9 
    ## ArchR logging successful to : ArchRLogs/ArchR-plotEmbedding-25d22a1110f7-Date-2021-05-05_Time-15-26-53.log

``` r
#To plot a specific gene we can subset this plot list using the gene name.
p$CD14
```

![](ArchR_tutorial_files/figure-gfm/Assign%20genes%20to%20clusters-1.png)<!-- -->

``` r
#Plot all genes defined in markerGenes
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
```

![](ArchR_tutorial_files/figure-gfm/Assign%20genes%20to%20clusters-2.png)<!-- -->

``` r
#Save an editable PDF version
plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)
```

    ## Plotting Ggplot!
    ## Plotting Ggplot!
    ## Plotting Ggplot!
    ## Plotting Ggplot!
    ## Plotting Ggplot!
    ## Plotting Ggplot!
    ## Plotting Ggplot!
    ## Plotting Ggplot!
    ## Plotting Ggplot!
