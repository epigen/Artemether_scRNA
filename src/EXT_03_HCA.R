require(project.init)
require(Matrix)
require("rhdf5")

require("project.init")
project.init2("artemether")

out <- "EXT_HCA/"
dir.create(dirout(out))

#require(cellrangerRkit)

read_h5 <- function(file, genome){
  dset <- h5read(file, name = genome)
  sparse_mat <- sparseMatrix(i = dset$indices + 1, p = dset$indptr, x = as.numeric(dset$data), dims = dset$shape, giveCsparse = FALSE)
  row.names(sparse_mat) <- as.character(dset$gene_names)
  colnames(sparse_mat) <- as.character(dset$barcodes)
  H5close()
  return(sparse_mat)
}

# READ AND MERGE TWO SAMPLES ----------------------------------------------
# Datasets downloaded from: https://preview.data.humancellatlas.org/
datBM <- read_h5(file=paste0("~/projects_shared/pathway_learning/external_data/HCA/ica_bone_marrow_h5.h5"), genome="GRCh38")
datCB <- read_h5(file=paste0("~/projects_shared/pathway_learning/external_data/HCA/ica_cord_blood_h5.h5"), genome="GRCh38")
str(datCB)
str(datBM)
quantile(as(datCB, "dgCMatrix")@i) # 0 - 33691
dim(datCB)
dim(datBM)
stopifnot(length(union(colnames(datCB), colnames(datBM))) == sum(ncol(datCB), ncol(datBM)))
head(colnames(datCB))
tail(colnames(datCB))
tail(colnames(datBM))
fullMat <- cbind(datCB, datBM)
#colnames(fullMat) <- gsub("Manton", "", gsub("_HiSeq_", "", colnames(fullMat)))
fullMat.counts <- fullMat
dim(fullMat)



# IDENTIFY CELLS -----------------------------------------------------
cell.umis <- Matrix::colSums(fullMat)
fullMat <- fullMat[,cell.umis >= 500]
cell.genes <- Matrix::colSums(fullMat != 0)
fullMat <- fullMat[,cell.genes >= 200]
mito.genes <- grep(pattern = "^MT-", x = rownames(fullMat), value = TRUE, ignore.case=TRUE)
percent.mito <- Matrix::colSums(fullMat[mito.genes, ])/Matrix::colSums(fullMat)
fullMat <- fullMat[,percent.mito < 0.15]
gc()
dim(fullMat)

fullLogTPM <- toTPM_fromRaw(fullMat, 1e6)
fullLogTPM <- toLogTPM(fullLogTPM)
hormoneDT <- data.table(
  barcode=colnames(fullLogTPM),
  INS=fullLogTPM["INS",],
  GCG=fullLogTPM["GCG",]
  )
write.tsv(hormoneDT, file=dirout(out, "Hormones.tsv"))
quantile(fullLogTPM["CST3",])
quantile(fullLogTPM["CD19",])
quantile(fullLogTPM["IL32",])
