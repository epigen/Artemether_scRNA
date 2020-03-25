require("project.init")
require(Seurat)
require(Matrix)
require(methods)
project.init2("artemether")

out <- "21_01_CellCorrelation/"
dir.create(dirout(out))

metaH <- loadHMeta2()
metaM <- loadMMeta2()
metaH3 <- loadH3Meta2()

(load(dirout("20_01_AlphaCells/", "Markers.RData")))
markers <- markers[power > 0.5]
write.tsv(markers, dirout(out, "Markers.tsv"))

xx <- "Markers"
for(xx in c("Markers")){
  outX <- paste0(out, xx, "/")
  dir.create(dirout(outX))

  org <- "human3"
  for(org in c("human", "mouse", "human3")){

    
    org.file <- dirout(outX, org, ".RData")
    if(!file.exists(org.file)){
      pbmc <- NULL
      meta <- NULL
      if(org == "mouse"){
        pbmc <- loadMData2()
        meta <- metaM
      }
      if(org == "human"){
        pbmc <- loadHData2()
        meta <- metaH
      }
      if(org == "human3"){
        pbmc <- loadH3Data2()
        meta <- metaH3
      }
      stopifnot(!is.null(pbmc) & !is.null(meta))
      
      genes <- row.names(pbmc@data)
      if(xx == "Markers") genes <- markers[V2 == gsub("\\d+$", "", org)]$V1
      
      meta <- meta[celltype %in% c("Alpha", "Beta", "Delta", "Gamma", "Endocrine")]
      meta2 <- data.table()
      for(rep in unique(meta$replicate)){
        metaR <- meta[replicate == rep]
        beta.mean <- Matrix::rowMeans(pbmc@data[genes, metaR[treatment == "DMSO" & celltype == "Beta"]$rn])
        alpha.mean <- Matrix::rowMeans(pbmc@data[genes, metaR[treatment == "DMSO" & celltype == "Alpha"]$rn])
        metaR$betaCor <- cor(as.matrix(pbmc@data[genes, metaR$rn]), beta.mean, method="spearman")
        metaR$alphaCor <- cor(as.matrix(pbmc@data[genes, metaR$rn]), alpha.mean, method="spearman")
        meta2 <- rbind(meta2, metaR)
      }
      save(meta2, genes, file=org.file)
    } else {
      load(org.file)
    }
    meta2[,corDiff := betaCor - alphaCor]
    write.tsv(meta2, dirout(outX, "MetaData_", org, ".tsv"))
    
    treatments <- as.character(unique(meta2$treatment))
    treatments <- treatments[treatments != "DMSO"]
    meta2$treatment <- factor(meta2$treatment, levels=c("DMSO", treatments))
    treat <- "A10"
    nrR <- length(unique(meta2$replicate))
    for(treat in treatments){
      # Density plots
      for(cc in c("betaCor", "alphaCor")){
        (p <- ggplot(meta2[treatment %in% c("DMSO", treat)], aes_string(x=cc, color="treatment")) + 
           theme_bw(12) + geom_density() + scale_color_manual(values=c("black", "red")))
        ggsave(dirout(outX, paste(org, treat, cc, sep="_"),"comb.pdf"), width=5, height=4)
        (p + facet_grid(. ~ replicate))
        ggsave(dirout(outX, paste(org, treat, cc, sep="_"), ".pdf"), width=5*nrR, height=4)
      }
      (p <- ggplot(meta2[treatment %in% c("DMSO", treat)], aes_string(x="corDiff", color="treatment")) + 
         theme_bw(12) + geom_density() + scale_color_manual(values=c("black", "red")))
      ggsave(dirout(outX, paste(org, treat, sep="_"),"_diff_comb.pdf"), width=5, height=4)
      (p + facet_grid(. ~ replicate))
      ggsave(dirout(outX, paste(org, treat, sep="_"), "_diff.pdf"), width=5*nrR, height=4)
      
      (p <- ggplot(meta2[treatment %in% c("DMSO", treat)], aes(x=alphaCor, y=betaCor, color=treatment)) + geom_point(alpha=0.05) +
        scale_color_manual(values=c("black", "red")))
      ggsave(dirout(outX, paste(org, treat, sep="_"), "_points.jpg"), width=5, height=4)
      (p + facet_grid(. ~ replicate))
      ggsave(dirout(outX, paste(org, treat, sep="_"), "_points_repl.jpg"), width=5*nrR, height=4)
    }
    for(cc in c("betaCor", "alphaCor")){
      ggplot(meta2, aes_string(x="treatment", y=cc)) + geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) +
        facet_grid(. ~ replicate) + theme_bw(12) + xRot()
      ggsave(dirout(outX, paste(org, "violin", cc, sep="_"), ".pdf"), width=5*nrR, height=4)
    }
  }
}