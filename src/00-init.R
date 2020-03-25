require(data.table)
require(pheatmap)
require(Seurat)
require(Matrix)
require(methods)


ENRICHR.DBS <- c("NCI-Nature_2016", "WikiPathways_2016", "Human_Gene_Atlas", "Chromosome_Location", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "TRANSFAC_and_JASPAR_PWMs")

SCALE.FACTOR <- 1e6
MIN.GENES <- 200
MIN.UMIS <- 500

source("src/FUNCTIONS_HELPERS.R")
source("src/FUNCTIONS.R")

options(stringsAsFactors=F)

annMeta <- function(meta){
  meta[,replicate := gsub(".Islets_(.+)_(.+)", "\\1", sample)]
  meta[,treatment := gsub(".Islets_(.+)_(.+)", "\\2", sample)]
  meta
}

# human
loadHMeta <- function(){annMeta(fread(dirout("10_H_01_Seurat/", "/MetaData_AnnotatedCells.tsv")))}
loadHData <- function(){load(dirout("10_H_01_Seurat/SeuratObject.RData")); return(pbmc)}


# human 3
loadH3Meta <- function(){
  meta <- fread(dirout("10_H3_01_Seurat/", "/MetaData_AnnotatedCells.tsv"))
  meta[,replicate := gsub("hIslets_III_(.+)_(.+)_(\\d)", "III\\3_\\2", sample)]
  meta[,treatment := gsub("hIslets_III_(.+)_(.+)_(\\d)", "\\1", sample)]
  return(meta)
}
loadH3Data <- function(){load(dirout("10_H3_01_Seurat/SeuratObject.RData")); return(pbmc)}

# Mouse 
loadMMeta <- function(){annMeta(fread(dirout("10_M_01_Seurat/", "/MetaData_AnnotatedCells.tsv")))}
loadMData <- function(){load(dirout("10_M_01_Seurat/SeuratObject.RData")); return(pbmc)}

# Reassigned
loadHMeta2 <- function(){fread(dirout("15_Counts/", "/Reassigned_Human.tsv"))}
loadMMeta2 <- function(){fread(dirout("15_Counts/", "/Reassigned_Mouse.tsv"))}
loadH3Meta2 <- function(){fread(dirout("15_Counts/", "/Reassigned_Human3.tsv"))}
# loadH3Meta2.filtered <- function(){fread(dirout("18_01_Filtered_H3/", "Meta_Filtered.tsv"))}

# Reassigned unfiltered
loadHMetaUnfiltered <- function(){fread(dirout("15_Counts/", "/Reassigned_Human_unfiltered.tsv"))}
loadMMetaUnfiltered <- function(){fread(dirout("15_Counts/", "/Reassigned_Mouse_unfiltered.tsv"))}
loadH3MetaUnfiltered <- function(){fread(dirout("15_Counts/", "/Reassigned_Human3_unfiltered.tsv"))}

# Filtered Data
loadHData2 <- function(){load(dirout("15_Counts/SeuratObject_human.RData")); return(pbmc)}
loadH3Data2 <- function(){load(dirout("15_Counts/SeuratObject_human3.RData")); return(pbmc)}
loadMData2 <- function(){load(dirout("15_Counts/SeuratObject_mouse.RData")); return(pbmc)}


# Colors for plots --------------------------------------------------------
COLORS.TREATMENTS <- c("#33a02c", "#a6cee3", "#1f78b4", "#ff7f00", "#000080", "#000000")
names(COLORS.TREATMENTS) <- c("GABA", "A1", "A10", "FoxO", "Aold", "DMSO")
