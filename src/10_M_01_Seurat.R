require("project.init")
require(Seurat)
require(Matrix)
require(methods)
project.init2("artemether")


sample.labels <- fread("metadata/Aggregate_Mouse_withLTI.csv")
data.path <- dirout("07_11_CleanMouse_CleanSpikeIns/CorrectedData.RData")
outS <- "10_M_01_Seurat/"

seurat.file <- dirout(outS, "SeuratObject",".RData")
if(!file.exists(seurat.file)){
  source("src/10_SCRIPT_01_Seurat.R")
} else {
  load(seurat.file)
}

# SEURAT ANALYSIS ---------------------------------------------------------
source("src/FUNC_Enrichr.R")
str(extra.genes.to.plot <- unique(fread("metadata//PancreasMarkers.tsv")$Mouse_GeneName))
seurat.diff.test <- "roc"
figExt <- ".jpg"
max.nr.of.cores <- 5
clusterings.exclude <- c("sample")
enrichrDBs <- c("NCI-Nature_2016", "WikiPathways_2016", "Human_Gene_Atlas", "Chromosome_Location", "NCI-60_Cancer_Cell_Lines")
source("src/FUNC_Seurat2_fast_simple.R",echo=TRUE)
max.nr.of.cores <- 1
source("src/FUNC_Seurat2_fast_simple.R",echo=TRUE)




# SUBCLUSTER --------------------------------------------------------------
outSOrig <- outS
outS <- paste0(outS, "Subcluster/")
cellsToKeep <- row.names(subset(pbmc@meta.data, res.0.5 %in% as.character(c(0:2,5:6))))
cluster.precision <- 0.5
sample.x <- "subcluster"
seurat.file.sub <- dirout(outS, "data",".RData")
source("src/FUNC_Seurat_subcluster.R",echo=TRUE)


# SEURAT 
max.nr.of.cores <- 5
source("src/FUNC_Seurat2_fast_simple.R",echo=TRUE)
max.nr.of.cores <- 1
source("src/FUNC_Seurat2_fast_simple.R",echo=TRUE)


