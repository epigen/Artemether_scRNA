require("project.init")
require(Seurat)
require(Matrix)
require(methods)
project.init2("artemether")

out <- "30_04_ClusteringExperiments/"
dir.create(dirout(out))


# LOAD META DATA  --------------------------------------------------------------
meta <- loadHMeta2()
meta <- meta[!sample %in% c("hIslets_I_pla2g16", "hIslets_II_Aold")]
cellsToKeep.list <- with(meta, split(rn, factor(paste(treatment, replicate, sep="_"))))
str(cellsToKeep.list)


# LOAD COUNT DATA -----------------------------------------------------------
(load(dirout("07_03_CleanHuman_CleanSpikeIns/","Uncorrected.Data.RData")))
(load(dirout("07_03_CleanHuman_CleanSpikeIns/", "CorrectedData.RData")))

# All meta barcodes are in both matrices
stopifnot(all(meta$rn %in% colnames(data.raw.cells)))
stopifnot(all(meta$rn %in% colnames(corr.export)))
# they have the same format
stopifnot(class(corr.export) == "dgCMatrix")
data.raw.cells <- as(data.raw.cells, "dgCMatrix")
# The have the same column names
stopifnot(length(intersect(colnames(corr.export),colnames(data.raw.cells))) == ncol(corr.export))
data.raw.cells <- data.raw.cells[,colnames(corr.export)]
stopifnot(all(colnames(corr.export) == colnames(data.raw.cells)))
# The have the same column sums
stopifnot(all(round(Matrix::colSums(data.raw.cells[,1:10])) == round(Matrix::colSums(corr.export[,1:10]))))


# Genes to use --------------------------------------------------------------
genes.to.use.file <- dirout(out, "Genes.to.use.tsv")
if(!file.exists(genes.to.use.file)){
  data.final.h <- loadHData()
  genes.to.use <- row.names(data.final.h@data)
  write.tsv(data.table(genes=genes.to.use), genes.to.use.file)
} else {
  genes.to.use <- fread(genes.to.use.file)$genes
}
stopifnot(all(genes.to.use %in% row.names(corr.export)))
stopifnot(all(genes.to.use %in% row.names(data.raw.cells)))


# raw or corrected data ---------------------------------------------------
typex <- "raw"
for(typex in c("raw", "corrected")){
  out2 <- paste0(out, "/", typex)
  dir.create(dirout(out2))
  
  xnam <- names(cellsToKeep.list)[1]
  for(xnam in names(cellsToKeep.list)){
    # do this for each
    outS <- paste0(out2, "/", xnam, "/")
    if(file.exists(dirout(outS, "SeuratObject",".RData"))) next
    
    pbmc.data <- if(typex == "raw") data.raw.cells else corr.export
    pbmc.data <- pbmc.data[genes.to.use,cellsToKeep.list[[xnam]]]
    
    cluster.precision <- c()
    sample.x <- xnam
    tryCatch({
      source("src/30_ClusteringExperiments_SCRIPT.R",echo=TRUE)
      
      ggplot(data.table(pbmc@dr$tsne@cell.embeddings), aes(x=tSNE_1, y=tSNE_2)) + geom_point(alpha=0.1)
      ggsave(dirout(outS, "TSNE.jpg"), w=4,h=4)
      
      ggplot(data.table(INS=pbmc@data["INS",]), aes(x=INS)) + geom_density()
      ggsave(dirout(outS, "INS.jpg"), w=4,h=4)
    }, error=function(e) message(e))
  }
}



