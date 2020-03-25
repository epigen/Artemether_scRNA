# This script does the same as 20_01 but just slightly different


require("project.init")
require(Seurat)
require(Matrix)
require(methods)
project.init2("artemether")

out <- "20_01_AlphaCells/"
dir.create(dirout(out))

metaH <- loadHMeta2()
metaM <- loadMMeta2()

pbmcM <- loadMData2()
pbmcH <- loadHData2()

# Identify Cell type markers ----------------------------------------------------------------
markers.file <- dirout(out, "Markers.RData")
if(!file.exists(markers.file)){
  require(doMC); registerDoMC(cores=2)
  
  org <- "mouse"
  for(org in c("human", "mouse")){
    pbmc <- pbmcM
    meta <- metaM
    if(org == "human"){
      pbmc <- pbmcH
      meta <- metaH
    }
    
    meta <- meta[treatment == "DMSO"]
    meta <- meta[rn %in% unique(do.call(c, lapply(split(meta$rn, factor(paste(meta$replicate, meta$celltype2))), function(x) sample(x, 200, replace=T))))]
    table(meta$celltype2)
    pbmc <- SubsetData(pbmc, cells.use=meta$rn)
    meta <- meta[match(row.names(pbmc@meta.data), meta$rn)]
    stopifnot(all(meta$rn == row.names(pbmc@meta.data)))
    pbmc@ident <- factor(meta$celltype2)
    names(pbmc@ident) <- pbmc@cell.names
    (celltypes <- unique(as.character(pbmc@ident)))
    cells <- "Alpha"
    foreach(cells = c("Alpha", "Beta")) %dopar% {
      if(!file.exists(dirout(out, cells, "_", org, ".tsv"))){
        cluster.markers <- FindMarkers(pbmc,
                                       thresh.use=0.7,
                                       test.use="roc",
                                       ident.1 = cells, 
                                       ident.2 = celltypes[celltypes != cells])
        write.table(cluster.markers, dirout(out, cells, "_", org, ".tsv"), sep="\t", quote=F, row.names=TRUE)
      }
    }
  }
  
  markers <- do.call(rbind, lapply(list.files(dirout(out), pattern="\\.tsv"), function(x) data.table(fread(dirout(out, x)), name=x)))
  markers[, name := gsub(".tsv", "", name)]
  markers <- cbind(markers, do.call(rbind, strsplit(markers$name, "_")))
  colnames(markers) <- make.unique(colnames(markers))
  save(markers, file=markers.file)
} else {
  load(markers.file)
}

# Calculate marker based scores ----------------------------------------------------------------
scores.file <- dirout(out, "Scores.RData")
if(!file.exists(scores.file)){
  org <- "mouse"
  for(org in c("mouse", "human")){
    pbmc <- pbmcM
    meta <- metaM
    if(org == "human"){
      pbmc <- pbmcH
      meta <- metaH
    }
    mark <- markers[V2 == org & myAUC > 0.75]
    ct <- "Alpha"
    for(ct in unique(mark$V1.1)){
      x <- as.matrix(pbmc@data[mark[V1.1 == ct]$V1,meta$rn])
      x <- x - rowMin(x)
      x <- x / rowMax(x)
      meta[[paste0(ct, "_score")]] <- colMeans(x)
      h <- sapply(split(meta$rn, factor(with(meta, paste(celltype2, treatment, replicate, sep="_")))),function(rns) rowMeans(x[,rns, drop=F]))
      pdf(dirout(out, "GeneAverages_", org, "_", ct, ".pdf"), width=25, height=8, onefile=F)
      pheatmap(h, cluster_cols=F)
      dev.off()
    }
    if(org == "mouse") metaM <- meta
    if(org == "human") metaH <- meta
  }
  save(metaH, metaM, file=scores.file)
} else { load(scores.file)}


# Plot scores ----------------------------------------------------------------
for(org in c("mouse", "human")){
  meta <- if(org == "human") metaH  else metaM
  meta$treatment <- factor(meta$treatment, levels = c("DMSO", unique(meta[treatment != "DMSO"]$treatment)))
  for(score in c("Beta_score", "Alpha_score")){
    for(treat in unique(meta$treatment)){
      ggplot(meta[treatment %in% c(treat, "DMSO") & celltype2 %in% c("Alpha", "Beta")], aes_string(x=score, color="treatment")) + theme_bw(12) + 
        stat_ecdf() + facet_grid(replicate ~ celltype2)+ scale_color_manual(values=c("black", "red"))
      ggsave(dirout(out, "Density_", org, "_", treat, "_", score, ".pdf"), width=8, height=12)
    }
    ggplot(meta[celltype2 %in% c("Alpha", "Beta")], aes_string(y=score, x="treatment")) + theme_bw(12) + 
      geom_violin() + facet_grid(replicate ~ celltype2)
    ggsave(dirout(out, "Violin_", org, "_", "_", score, ".pdf"), width=8, height=12)
  }
}