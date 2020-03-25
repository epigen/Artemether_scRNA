require("project.init")
require(Seurat)
require(Matrix)
require(methods)
project.init2("artemether")

out <- "05_SpikeIns_Clean/"
dir.create(dirout(out))


path <- paste0(getOption("PROCESSED.PROJECT"),'results_pipeline/cellranger_count/MF179')


# LOAD DATA INTO SEURAT OBJECT --------------------------------------------
scDat.file <- dirout(out, "scData.RData")
if(!file.exists(scDat.file)){
  datH <- Read10X(paste0(path, '_hmG/outs/raw_gene_bc_matrices/hg19/'))
  datM <- Read10X(paste0(path, '_hmG/outs/raw_gene_bc_matrices/mm10/'))
  str(datH); str(datM)
  stopifnot(all(sort(colnames(datH)) == sort(colnames(datM))))
  pbmc.data <- rbind(datH, datM)
  # Apply UMI threshold
  cell.umis <- Matrix::colSums(pbmc.data)
  pbmc.data <- pbmc.data[,cell.umis >= MIN.UMIS]
  cell.genes <- Matrix::colSums(pbmc.data != 0)
  pbmc.data <- pbmc.data[,cell.genes >= MIN.GENES]
  # Set up Seurat Object
  pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = MIN.GENES, scale.factor=SCALE.FACTOR, project = "MH")
  #pbmc@data@x <- log(((exp(pbmc@data@x) - 1) / 100) + 1)
  save(pbmc, file=scDat.file)
} else {
  load(scDat.file)
}


# IDENTIFY ORGANISM -------------------------------------------------------
str(pbmc@data)
org.idx <- substr(row.names(pbmc@data), 0,4) == "hg19"
org.reads.H <- Matrix::colSums(pbmc@raw.data[org.idx, colnames(pbmc@data)])
org.reads.M <- Matrix::colSums(pbmc@raw.data[!org.idx, colnames(pbmc@data)])
pDat <- data.table(rh = org.reads.H, rm = org.reads.M, barcode=colnames(pbmc@data))
pDat$nGene <- Matrix::colSums(pbmc@raw.data[, colnames(pbmc@data)] != 0)
pDat[,org := "NA"]
pDat[,ratio := log2(rh/rm)]
pDat[ratio > 2, org := "human"]
pDat[ratio < -2, org := "mouse"]

ggplot(pDat, aes(x=rh, y=rm, color=org, size=nGene)) + geom_point(alpha=0.3) + geom_abline()+
  scale_y_log10() + scale_x_log10() + xlab("Human reads") + ylab("Mouse reads")
ggsave(dirout(out, "IdentifyOrganism.pdf"), height=5, width=6)



# EXPORT ANNOTATION AS FILE -----------------------------------------------
pDat <- pDat[org != "NA"]
write.tsv(pDat[,c("barcode", "org", "ratio")], file=dirout(out, "CrossSpecies_Alignment.tsv"))



# SIGNATURE OF SPIKE INS IN cross alignment --------------------------------------------------
human.sig <- data.table(gene=row.names(pbmc@data)[org.idx],
                        value=Matrix::rowMeans(pbmc@data[org.idx, pDat[org == "human"]$barcode]))
write.tsv(human.sig, dirout(out, "Human_hmG_signature.tsv"))

mouse.sig <- data.table(gene=row.names(pbmc@data)[!org.idx],
                        value=Matrix::rowMeans(pbmc@data[!org.idx, pDat[org == "mouse"]$barcode]))
write.tsv(mouse.sig, dirout(out, "Mouse_hmG_signature.tsv"))



# EXPORT SPIKE IN DATA ALIGNED TO TWO ORGANISMS ----------------------------------------
org.to.export <- "human"
for(org.to.export in c("mouse", "human")){
  export.file <- dirout(out, "SpikeData_", org.to.export, ".RData")
  if(file.exists(export.file)) next
  dat <- Read10X(paste0(path, "_", org.to.export, '/outs/raw_gene_bc_matrices/', ifelse(org.to.export == "mouse", "mm10", "GRCh38") ,"/"))
  dat.spike <- dat[, colnames(dat) %in% pDat$barcode]
  str(dat.spike)
  cell.umis <- Matrix::colSums(dat.spike)
  dat.spike <- dat.spike[,cell.umis >= MIN.UMIS]
  cell.genes <- Matrix::colSums(dat.spike != 0)
  dat.spike <- dat.spike[,cell.genes >= MIN.GENES]
  tpm.dat <- t(t(dat.spike) * SCALE.FACTOR/Matrix::colSums(dat.spike))
  dat.log.spike <- toLogTPM(tpm.dat)
  metaDat <- data.table(barcode=colnames(tpm.dat), nUMI=Matrix::colSums(dat.spike), nGene=Matrix::colSums(dat.spike != 0))
  metaDat <- merge(metaDat, pDat[,c("org", "barcode"), with=F], by="barcode")
  write.tsv(metaDat[order(org)], dirout(out, "SpikeMeta_", org.to.export,  ".tsv"))
  save(dat.log.spike, file=export.file)
}