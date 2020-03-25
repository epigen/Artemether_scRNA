require("project.init")
require(Seurat)
require(Matrix)
require(methods)
project.init2("artemether")

out <- "06_SpikeIns_Islets/"
dir.create(dirout(out))

# args = commandArgs(trailingOnly=TRUE)
# orgX <- args[1]
# if(is.na(args[1])) orgX <- "mouse"

path <- paste0(getOption("PROCESSED.PROJECT"),'results_pipeline/cellranger_count/')
f <- list.files(path, pattern = paste0("^.Islets", ".+_hmG$"))
f <- list.files(path, pattern = paste0("^hIslets_III_", ".+_hmG$"))

require(doMC); registerDoMC(cores=10)
fX <- f[8]
foreach(fX = f) %dopar% {
  orgX <- if(grepl("hIslets", fX)) "human" else "mouse"
  
  message(fX)
  outF <- paste0(out, fX, "/")
  dir.create(dirout(outF))
  
  # LOAD DATA INTO SEURAT OBJECT --------------------------------------------
  scDat.file <- dirout(outF, "scData.RData")
  if(!file.exists(scDat.file)){
    datH <- Read10X(paste0(path, fX, '/outs/raw_gene_bc_matrices/hg19/'))
    datM <- Read10X(paste0(path, fX, '/outs/raw_gene_bc_matrices/mm10/'))
    str(datH); str(datM)
    stopifnot(all(sort(colnames(datH)) == sort(colnames(datM))))
    pbmc.data <- rbind(datH, datM)
    # Apply UMI threshold
    cell.umis <- Matrix::colSums(pbmc.data)
    pbmc.data <- pbmc.data[,cell.umis >= 500]
    cell.genes <- Matrix::colSums(pbmc.data != 0)
    pbmc.data <- pbmc.data[,cell.genes >= 200]
    # Set up Seurat Object
    pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, scale.factor=1e6, project = "MH")
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
  pDat[,nUMI := rm + rh]
  pDat[org == "human",contamination := rm/nUMI]
  pDat[org == "mouse",contamination := rh/nUMI]
  pDat[org == orgX, contamination := NA]
  
  (p <- ggplot(pDat, aes(x=rh, y=rm, color=org, size=nGene)) + geom_point(alpha=0.3) + geom_abline()+
    scale_y_log10() + scale_x_log10() + xlab("Human reads") + ylab("Mouse reads") + theme_bw(12))
  ggsave(dirout(outF, "IdentifyOrganism.pdf"), height=4, width=5)
  ggsave(dirout(outF, "IdentifyOrganism.jpg"), height=4, width=5)
  ggsave(dirout(outF, "IdentifyOrganism2.pdf"), height=4, width=4, plot=p+guides(color=F, fill=F, size=F))
  ggsave(dirout(outF, "IdentifyOrganism2.jpg"), height=4, width=4, plot=p+guides(color=F, fill=F, size=F))
  
  # Annotate spike ins via correlation
  spike.sig <- fread(dirout("05_SpikeIns_Clean/",ifelse(orgX == "human", "Human", "Mouse"), "_hmG_signature.tsv"))
  spike.sig <- spike.sig[gene %in% row.names(pbmc@data)]
  corr <- cor(spike.sig$value, as.matrix(pbmc@data[spike.sig$gene,]), use="pairwise.complete.obs")
  pDat$corrToSpike <- corr[1,pDat$barcode]
  pDat[(corrToSpike > 0.9 & org == orgX) | ! org %in% c("NA", orgX), spikeIn := org]
  table(pDat$spikeIn)
  write.tsv(pDat[,c("barcode", "org", "ratio", "contamination","corrToSpike", "spikeIn")], file=dirout(outF, "IdentifyOrganism.tsv"))
  write.tsv(pDat, file=dirout(outF, "IdentifyOrganism_full.tsv"))
  
  
  # LOOK AT CONTAMINATION ---------------------------------------------------
  genes.use <- if(orgX == "human") row.names(pbmc@data)[org.idx] else row.names(pbmc@data)[!org.idx]
  cont <- data.table(inHuman = Matrix::rowMeans(pbmc@data[genes.use,pDat[org == "human"]$barcode]))
  cont$inMouse <- Matrix::rowMeans(pbmc@data[genes.use,pDat[org == "mouse"]$barcode])
  cont$gene <- if(orgX == "human") gsub("hg19_", "", genes.use) else gsub("mm10_", "", genes.use)
  write.tsv(cont, file=dirout(outF, "Contamination.tsv"))
  
  xAxis <- "inMouse";yAxis <- "inHuman"
  if(orgX == "human"){
    xAxis <- "inHuman";yAxis <- "inMouse"
  }
  ggplot(cont, aes_string(x=xAxis, y=yAxis)) + geom_point(alpha=0.5, color="grey") #+ geom_text(data=cont[inHuman > 130,], aes(label=gene))
  ggsave(dirout(outF, "Contamination_vs_Expression.pdf"), height=5, width=6)
  
  ggplot(cont, aes_string(x=xAxis, y=yAxis)) + geom_point() +
    scale_y_log10() + scale_x_log10() 
  ggsave(dirout(outF, "Contamination_vs_Expression_log.pdf"), height=5, width=6)
}