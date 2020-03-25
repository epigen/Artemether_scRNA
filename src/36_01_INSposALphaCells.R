require("project.init")
require(Seurat)
require(Matrix)
require(methods)
require(Rtsne)
project.init2("artemether")

out <- "36_01_AlphaInsHighVsLow/"
dir.create(dirout(out))



metaH <- loadHMeta2()
metaM <- loadMMeta2()
pbmcH <- loadHData2()
pbmcM <- loadMData2()

metaH <- assignGenesMeta(c("INS", "GCG"), metaH, pbmcH)
metaM <- assignGenesMeta(c("Ins2", "Gcg"), metaM, pbmcM)
metaM[,INS := Ins2]


ct <- "Alpha"
orgx <- "human"

# GET DIFF GENES ----------------------------------------------------------
cluster.markers <- list()
for(orgx in c("mouse", "human")){
  metaX <- metaH
  pbmc <- pbmcH
  if(orgx == "mouse"){
    metaX <- metaM
    pbmc <- pbmcM
  }
  metaX <- metaX[celltype2 == "Alpha"]
  metaX$group <- c("INShigh", "INSlow")[(metaX$INS == 0) + 1]
  pbmc <- SubsetData(pbmc, cells.use=metaX$rn)
  pbmc@meta.data[["cluster"]] <- metaX[match(row.names(pbmc@meta.data), metaX$rn)]$group
  pbmc@ident <- factor(pbmc@meta.data[["cluster"]])
  names(pbmc@ident) <- pbmc@cell.names # it needs those names apparently
  cluster.markers[[orgx]] <- FindMarkers(pbmc, test.use="negbinom",ident.1 = "INShigh", ident.2 = "INSlow")
  write.table(cluster.markers[[orgx]], dirout(out, orgx, "Markers.tsv"), sep="\t", quote=F, row.names=TRUE)
}
