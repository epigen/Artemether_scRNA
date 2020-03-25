# THIS SCRIPT REQUIRES:
# pbmc object (seurat)
# outS (output dir)


if(!exists("enrichrDBs")) enrichrDBs <- c("NCI-Nature_2016", "WikiPathways_2016", "Human_Gene_Atlas", "Chromosome_Location")
if(!exists("extra.genes.to.plot")) extra.genes.to.plot <- c()
if(!exists("seurat.diff.test")) seurat.diff.test <- "roc"
if(!exists("seurat.thresh.use")) seurat.thresh.use <- 0.25
if(!exists("max.nr.of.cores")) max.nr.of.cores <- 3
if(!exists("seurat.min.pct")) seurat.min.pct <- 0.25
if(!exists("cellMarker.file")) cellMarker.file <- "metadata/CellMarkers.csv"
if(!exists("enrichGeneList.oddsRatio")) tryCatch(source("src/FUNC_Enrichr.R"), error=function(e) stop("Cannot load Enrichr Functions, need to preload (FUNC_SEURAT2.R)"))
if(!exists("figExt")) figExt <- ".pdf"
if(!exists("clusterings.exclude")) clusterings.exclude <- NULL

print(dim(pbmc@data))

require(data.table)
require(pheatmap)
require(ggplot2)
require(doMC)
require(pryr)


# manage memory and number of tasks used
mem_u <- as.numeric(mem_used())/10**6
cores_u = max(1, min(12, floor(180000/mem_u)-1))
if(is.na(mem_u)) cores_u <- 3
cores_u <- min(cores_u, max.nr.of.cores)
message("Memory used: ", mem_u, " cores: ",cores_u)
registerDoMC(cores=cores_u)


# use colnames from pbmc@meta.data as clusterings
# only those with >= 2 groups with > 1 cell
clusterings <- colnames(pbmc@meta.data)[apply(pbmc@meta.data, 2, function(x) sum(table(x[x!="IGNORED"])>2)>1)]
clusterings <- clusterings[!clusterings %in% c("nUMI", "nGene", "orig.ident", "percent.mito", "var.ratio.pca")]
# reorder to do the kmeans at the end
clusterings <- clusterings[!grepl("ClusterNames", clusterings)]
# reorder to do the kmeans at the end
(clusterings <- c(clusterings[!grepl("res", clusterings)], clusterings[grepl("res", clusterings)]))
clusterings <- clusterings[!clusterings %in% clusterings.exclude]
message("Clusterings: ", paste(clusterings, collapse=" "))

for(cl.x in clusterings){
  print(table(pbmc@meta.data[[cl.x]]))
}


# PLOT MARKERS 2
message("Plotting Known marker genes")
outMarkers <- paste0(outS, "Markers/")
dir.create(dirout(outMarkers))
if(file.exists(cellMarker.file)){
  markers <- fread(cellMarker.file)[Marker != ""]
  for(i in 1:nrow(markers)){
    if(markers[i]$GeneSymbol %in% rownames(pbmc@data)){
      if(file.exists(dirout(outMarkers, markers[i]$GeneSymbol,figExt))) next
      marDat <- data.table(pbmc@dr$tsne@cell.embeddings, Expression=pbmc@data[markers[i]$GeneSymbol,])
      ggplot(marDat, aes(x=tSNE_1, y=tSNE_2, color=Expression)) + geom_point(alpha=0.5) + 
        scale_color_gradientn(colors=c("grey", "blue")) + theme(legend.position = 'none') +
        ggtitle(paste0(markers[i]$GeneSymbol, "/", markers[i]$Marker, "\n", markers[i]$CellType))
      ggsave(dirout(outMarkers, markers[i]$GeneSymbol,figExt), height=7, width=7)
    }
  }
}
for(geneSyn in extra.genes.to.plot){
  if(file.exists(dirout(outMarkers, geneSyn,figExt))) next
  if(geneSyn %in% row.names(pbmc@data)){
    marDat <- data.table(pbmc@dr$tsne@cell.embeddings, Expression=pbmc@data[geneSyn,])
    ggplot(marDat, aes(x=tSNE_1, y=tSNE_2, color=Expression)) + geom_point(alpha=0.5) +
      scale_color_gradientn(colors=c("grey", "blue")) + theme(legend.position = 'none') +
      ggtitle(paste0(geneSyn))
    ggsave(dirout(outMarkers, geneSyn,figExt), height=7, width=7)
  }
}




# PLOT UMIS ---------------------------------------------------------------
message("Plotting UMIs")
try({
  umip <- ggplot(data.table(pbmc@dr$tsne@cell.embeddings, UMIs=pbmc@meta.data$nUMI), aes(x=tSNE_1,y=tSNE_2, color=log10(UMIs))) +
    scale_color_gradient(low="blue", high="red") +
    geom_point() + ggtitle(paste(nrow(pbmc@data), "genes\n", ncol(pbmc@data), "cells")) + theme_bw(24)
  ggsave(dirout(outS, "UMI", figExt),plot=umip, height=7, width=7)
}, silent = TRUE)

# PLOT CLUSTERS
message("Plotting Clusters")
pDat <- data.table(pbmc@dr$tsne@cell.embeddings)
cl.x <- "sample"
for(cl.x in clusterings){
  print(cl.x)
  pDat[[cl.x]] <-pbmc@meta.data[[cl.x]]
  labelCoordinates <- pDat[,.(tSNE_1=median(tSNE_1), tSNE_2=median(tSNE_2)),by=cl.x]
  
  ggplot(pDat[get(cl.x) != "IGNORED"], aes_string(x=cl.x)) + geom_bar() + coord_flip()
  ggsave(dirout(outS, "Cluster_counts_", cl.x, ".pdf"), height=7, width=7)
  
  ggplot(pDat, aes_string(x="tSNE_1",y="tSNE_2", color=cl.x)) + geom_point(alpha=0.5) + # ggtitle(sample.x) +
    geom_label(data=labelCoordinates, aes_string(x="tSNE_1", y="tSNE_2", label=cl.x), color="black", alpha=0.5)
  ggsave(dirout(outS, "Cluster_tSNE_", cl.x, figExt), height=7, width=7)
}
write.table(pDat, dirout(outS,"Cluster.tsv"), sep="\t", quote=F, row.names=F)
clusterCounts <- pDat[,-c("tSNE_1", "tSNE_2"), with=TRUE]
clusterCounts <- do.call(rbind, lapply(names(clusterCounts), function(nam) data.table(clusterCounts[[nam]], nam)))
write.table(clusterCounts[, .N, by=c("V1", "nam")], dirout(outS, "ClusterCounts.tsv"), sep="\t", quote=F,row.names=F)


# Markers for each cluster ------------------------------------------------
message("Plotting cluster Markers")
(cl.x <- clusterings[1])
x <- foreach(cl.x = clusterings) %do% { 
  if(!is.null(pbmc@meta.data[[cl.x]])){
    print(cl.x)
    x <- gsub("ClusterNames_","", cl.x)
    out.cl <- paste0(outS, "Cluster_",x, "/")
    dir.create(dirout(out.cl))
    pbmc@ident <- factor(pbmc@meta.data[[cl.x]])
    names(pbmc@ident) <- pbmc@cell.names # it needs those names apparently
    clusters <- names(table(pbmc@ident))[table(pbmc@ident)>2]
    clusters <- clusters[clusters != "IGNORED"]
    (cl.i <- clusters[10])
    foreach(cl.i = clusters) %dopar% {
      # message(cl.i)
      if(!file.exists(dirout(out.cl, "Markers_Cluster",cl.i, ".tsv"))){
          #           print(tt <- Sys.time())
          cluster.markers <- FindMarkers(pbmc,  ident.1 = cl.i, ident.2 = clusters[clusters != cl.i],
                                         test.use=seurat.diff.test, 
                                         min.pct = seurat.min.pct,
                                         thresh.use=seurat.thresh.use)
          #           print(tt - Sys.time())
          write.table(cluster.markers, dirout(out.cl, "Markers_Cluster",cl.i, ".tsv"), sep="\t", quote=F, row.names=TRUE)
        try({
          pdf(dirout(out.cl,"Markers_Cluster",cl.i,".pdf"), height=15, width=15)
          FeaturePlot(pbmc, row.names(cluster.markers)[1:min(nrow(cluster.markers),9)],cols.use = c("grey","blue"))
          dev.off()
        }, silent=T)
        # Plot for all clusters heatmap
      }
    }
  }
}



# Enrichr on first Markers -----------------------------------------------------------------
message("EnrichR on first markers")
cl.x <- "patient"
cl.x <- clusterings[1]
for(cl.x in clusterings){
  if(!is.null(pbmc@meta.data[[cl.x]])){
    print(cl.x)
    x <- gsub("ClusterNames_","", cl.x)
    #     if(!file.exists(dirout(outS, "Enrichr_",x,".tsv"))){
    out.cl <- paste0(outS, "Cluster_",x, "/")
    f <- list.files(dirout(out.cl), pattern="^Markers_.+?[^2]\\.tsv$")
    genes <- lapply(f, function(fx) fread(dirout(out.cl, fx)))
    names(genes) <- gsub("Markers_", "", gsub(".tsv","",f))
    if(seurat.diff.test == "roc"){
      genes <- genes[sapply(genes, ncol) == 6]
      genes <- lapply(genes, function(fx) fx[V4 > 0.7]$V1)
    } else {
      genes <- genes[sapply(genes, ncol) == 5]
      genes <- lapply(genes, function(fx) fx[V2 < 0.05 & V3 > 0.3]$V1)
    }
    genes <- genes[sapply(genes, length) > 4]
    
    # HEATMAP
    hm.file <- dirout(outS, "Cluster_simple_HM_",cl.x, ".pdf")
    #if(!file.exists(hm.file)){
    try({
      cells.to.plot <- do.call(c, lapply(split(row.names(pbmc@meta.data), factor(pbmc@meta.data[[cl.x]])), function(x) sample(x, min(100, length(x)))))
      genes.to.plot <- do.call(c, lapply(genes, function(x) x[1:20]))
      genes.to.plot <- names(sort(table(genes.to.plot)))[1:min(200, length(genes.to.plot))]
      data.to.plot <- pbmc@data[genes.to.plot, cells.to.plot]
      data.to.plot <- data.to.plot/apply(data.to.plot, 1, max)
      pdf(hm.file, width=15, height=min(29, 1+length(genes.to.plot)*0.1),onefile=F)
      pheatmap(data.to.plot, show_colnames=F, cluster_cols=F, fontsize_row=6,
               annotation_col=pbmc@meta.data[,cl.x, drop=F], color=colorRampPalette(c("lightgrey", "blue", "black"))(50))
      dev.off()
    }, silent=TRUE)
    #}
    
    # ENRICHER RESULTS
    enrichRes <- data.table()
    for(grp.x in names(genes)){
      ret=try(as.data.table(enrichGeneList.oddsRatio(genes[[grp.x]],databases = enrichrDBs)),silent = FALSE)
      if(!any(grepl("Error",ret)) && nrow(ret) > 0){
        enrichRes <- rbind(enrichRes, data.table(ret, grp = grp.x))
      }
    }
    
    if(nrow(enrichRes) > 0){
      enrichRes <- enrichRes[hitLength >= 3]
      write.table(enrichRes, file=dirout(outS, "Enrich_simpleMarkers_",x,".tsv"), sep="\t", quote=F, row.names=F)
      # plot
      enrClass <- enrichRes$database[1]
      for(enrClass in unique(enrichRes$database)){
        enrichResX <- enrichRes[database == enrClass]
        if(nrow(enrichResX) > 0){
          enrichr.plot(enrichResX)
          ggsave(dirout(outS, "Enrich_simpleMarkers_", enrClass, "_",x,".pdf"), width=min(29, 6+ length(unique(enrichResX$grp))*0.3), height=min(29, length(unique(enrichResX$category))*0.3 + 4))
        }
      }
    } else {
      message("No enrichments for ", cl.x)
    }
  }
}


message("Analysis pipeline completed successfully!")

