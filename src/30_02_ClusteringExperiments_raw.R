require("project.init")
require(Seurat)
require(Matrix)
require(methods)
project.init2("artemether")

source("src/FUNC_Enrichr.R")

data.raw <- Read10X(paste0(getOption("PROCESSED.PROJECT"), "results_pipeline/cellranger_count/","Human1_2","/outs/raw_gene_bc_matrices_mex/GRCh38/"))

out <- "30_02_ClusteringExperiments_raw/"
dir.create(dirout(out))

# SUBCLUSTER --------------------------------------------------------------
meta <- loadHMeta2()
meta <- meta[!sample %in% c("hIslets_I_pla2g16", "hIslets_II_Aold")]
cellsToKeep.list <- with(meta, split(rn, factor(paste(treatment, replicate, sep="_"))))
str(cellsToKeep.list)

xnam <- names(cellsToKeep.list)[1]
for(xnam in names(cellsToKeep.list)){
  # do this for each
  outS <- paste0(out, xnam, "/")
  if(file.exists(dirout(outS, "SeuratObject",".RData"))) next
  cluster.precision <- c()
  sample.x <- xnam
  pbmc.data <- data.raw[,cellsToKeep.list[[xnam]]]
  tryCatch({
    source("src/30_02_ClusteringExperiments_raw_SCRIPT.R",echo=TRUE)
    ggplot(data.table(pbmc@dr$tsne@cell.embeddings), aes(x=tSNE_1, y=tSNE_2)) + geom_point(alpha=0.1)
    ggsave(dirout(outS, "TSNE.jpg"), w=4,h=4)
  }, error=function(e) message(e))
}



