require("project.init")
require(Seurat)
require(Matrix)
require(methods)
project.init2("artemether")

data.raw <- Read10X(paste0(getOption("PROCESSED.PROJECT"), "results_pipeline/cellranger_count/","Human1_2","/outs/raw_gene_bc_matrices_mex/GRCh38/"))

out <- "30_03_ClusteringExperiments_raw_subcluster/"
dir.create(dirout(out))

# SUBCLUSTER --------------------------------------------------------------
meta <- loadHMeta2()
meta <- meta[!sample %in% c("hIslets_I_pla2g16", "hIslets_II_Aold")]
cellsToKeep.list <- with(meta, split(rn, factor(paste(treatment, replicate, sep="_"))))
str(cellsToKeep.list)


# Create one base seurat object -------------------------------------------
pbmcOLD <- CreateSeuratObject(raw.data = data.raw[, meta$rn], min.cells = 3, min.genes = MIN.GENES, scale.factor=SCALE.FACTOR, project = "X")

xnam <- names(cellsToKeep.list)[1]
for(xnam in names(cellsToKeep.list)){
  # do this for each
  outS <- paste0(out, xnam, "/")
  cellsToKeep <- cellsToKeep.list[[xnam]]
  cluster.precision <- c()
  
  pbmc <- pbmcOLD
  source("src/FUNC_Seurat_subcluster.R",echo=TRUE)
  ggplot(data.table(pbmc@dr$tsne@cell.embeddings), aes(x=tSNE_1, y=tSNE_2)) + geom_point(alpha=0.1)
  ggsave(dirout(outS, "TSNE.jpg"), w=4,h=4)
}



