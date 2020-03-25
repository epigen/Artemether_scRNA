
toTPM_fromRaw <- function(m, scale.factor){t(t(m)/Matrix::colSums(m)) * scale.factor}
toTPM <- function(m){m@x <- exp(m@x) - 1; return(m)}
toLogTPM <- function(m){m@x <- log(m@x + 1); return(m)}

corS <- function(...){cor(..., method="spearman")}

add.genes.sparseMatrix <- function(m, genes){
  gX <- genes[which(!genes %in% row.names(m))]
  m2 <- matrix(0, nrow=length(gX), ncol=ncol(m))
  m2 <- as(m2, "dgTMatrix")
  row.names(m2) <- gX
  rbind(m, m2)
}

toSparseMT <- function(mt){
  ret <- as(mt, "dgTMatrix")
  row.names(ret) <- row.names(mt)
  colnames(ret) <- colnames(mt)
  ret
}

ggS <- function(..., w=4, h=4){ggsave(..., width=4, height=4)}


assignGenesMeta <- function(genes, metaX, seuratObj, strict=T){ 
  for(gg in genes){
    if(!gg %in% row.names(seuratObj@data)){
      if(strict){
        stop(gg, " not found in sc data")
      }
    }
    metaX[[gg]] <- seuratObj@data[gg,metaX$rn]
  }
  metaX
}