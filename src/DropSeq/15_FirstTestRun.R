require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)

project.init2("artemether")

sample.x <- "IsletTestRun"
data.path <- paste0(Sys.getenv("PROCESSED"), "10x_datasets/results_pipeline/cellranger_count/Islets_10x_test/outs/raw_gene_bc_matrices/GRCh38/")
outS <- "15_TestRun/"
dir.create(dirout(outS))
source(paste0(Sys.getenv("CODEBASE"), "10x_datasets/src/FUNC_Seurat_preprocess.R"),echo=TRUE)
# extra.genes.to.plot <- unique(fread("metadata//PancreasMarkers.tsv")$Human_GeneName)
# x <- fread("metadata//PancreasMarkers.tsv")
# x[Mouse_GeneName == "", Mouse_GeneName := paste0(substr(Human_GeneName,0,1), tolower(substr(Human_GeneName,2,100)))]
# write.table(x, "metadata//PancreasMarkers.tsv", sep="\t", quote=F, row.names=F)
source(paste0(Sys.getenv("CODEBASE"), "10x_datasets/src/FUNC_Enrichr.R"))
source(paste0(Sys.getenv("CODEBASE"), "10x_datasets/src/FUNC_Seurat2.R"),echo=TRUE)




