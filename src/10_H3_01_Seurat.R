require("project.init")
require(Seurat)
require(Matrix)
require(methods)
project.init2("artemether")


sample.labels <- fread("metadata/Aggregate_Human3.csv")
data.path <- dirout("07_04_CleanHuman3_CleanSpikeIns/CorrectedData.RData")
outS <- "10_H3_01_Seurat/"
dir.create(dirout(outS))

seurat.file <- dirout(outS, "SeuratObject",".RData")
if(!file.exists(seurat.file)){
  source("src/10_SCRIPT_01_Seurat.R")
} else {
  load(seurat.file)
}


# SEURAT ANALYSIS ---------------------------------------------------------
source("src/FUNC_Enrichr.R")
str(extra.genes.to.plot <- unique(fread("metadata//PancreasMarkers.tsv")$Human_GeneName))
seurat.diff.test <- "roc"
figExt <- ".jpg"
# max.nr.of.cores <- 5
clusterings.exclude <- c("sample")
enrichrDBs <- c("NCI-Nature_2016", "WikiPathways_2016", "Human_Gene_Atlas", "Chromosome_Location", "NCI-60_Cancer_Cell_Lines")
# source("src/FUNC_Seurat2_fast_simple.R",echo=TRUE)
max.nr.of.cores <- 1
source("src/FUNC_Seurat2_fast_simple.R",echo=TRUE)




# SUBCLUSTER --------------------------------------------------------------
outSOrig <- outS
outS <- paste0(outS, "Subcluster/")
cellsToKeep <- row.names(subset(pbmc@meta.data, res.0.5 %in% as.character(c(0:2,4,5,7,15))))
cluster.precision <- 0.5
sample.x <- "subcluster"
source("src/FUNC_Seurat_subcluster.R",echo=TRUE)

# SEURAT ANALYSIS ---------------------------------------------------------
max.nr.of.cores <- 5
source("src/FUNC_Seurat2_fast_simple.R",echo=TRUE)
max.nr.of.cores <- 1
source("src/FUNC_Seurat2_fast_simple.R",echo=TRUE)




