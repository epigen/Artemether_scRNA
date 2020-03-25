require("project.init")
project.init2("artemether")

out <- "FIG_05_Human_3/"
dir.create(dirout(out))


# meta.filtered <- loadH3Meta2.filtered()
meta <- loadH3Meta2()
sobj <- loadH3Data2()


# Cell Numbers H3 ---------------------------------------------------------
write.tsv(data.table(meta[,.N, by=c("replicate", "treatment")], organism = "human"),dirout(out, "Supp_CellNumbers.tsv"))


# GET AVERAGES -----------------------------------------
for(typex in c("")){
  avg.file <- dirout(out, "Avgerages",typex,".csv")
  if(!file.exists(avg.file)){
    #     metaX <- if(typex == "") meta else meta.filtered
    metaX <- meta
    metaX <- metaX[celltype2 %in% c("Beta", "Alpha")]
    avg <- sapply(split(metaX$rn, factor(with(metaX, paste(celltype2, treatment, replicate)))), function(rns){
      log(Matrix::rowMeans(toTPM(sobj@data[,rns])) + 1)
    })
  
    # write files
    write.csv(avg, avg.file)
    if(typex == "") save(avg, file=dirout(out, "Averages", typex, ".RData"))
  }
}
(load(dirout(out, "Averages.RData")))


# CORRELATE AVERAGES ------------------------------------------------------
ctx <- "Alpha"
for(ctx in c("Alpha", "Beta")){
  idx <- grepl(ctx, colnames(avg))
  
  cMT <- corS(avg[,idx])
  dMT <- as.dist(1-cMT)
  diag(cMT) <- NA
  stopifnot(all(colnames(cMT) == colnames(avg[,idx])))
  stopifnot(all(attr(dMT, "Labels") == colnames(avg[,idx])))
  stopifnot(all(colnames(cMT) == row.names(cMT)))
  
  cleanDev(); pdf(dirout(out, "Clustering_Replicates_",ctx,".pdf"), w=5,h=5)
  pheatmap(cMT, clustering_distance_rows=dMT, clustering_distance_cols=dMT)
  dev.off()
}



# ALPHA + cells -----------------------------------------------------------
res <- data.table()
metaO <- copy(meta)
metaO[, time := gsub(".+_", "", metaO$replicate)]
metaO[, replicate := gsub("_\\d+", "", replicate)]
for(timex in unique(metaO$time)){
  metaO[INS > 0, insulin := "INS+"]
  metaO[!INS > 0, insulin := "INS-"]
  for(treatx in c("A10")){
#     for(repx in unique(metaO$replicate)){
#       if(!timex %in% unique(metaO[treatment == treatx & replicate == repx]$time)) next
#       m <- with(metaO[treatment %in% c("DMSO", treatx)][replicate == repx][celltype == "Alpha"][time == timex], table(insulin, treatment))
#       m <- m[c("INS+", "INS-"),c(treatx, "DMSO")]
#       f <- fisher.test(m)
#       res <- rbind(res, data.table(time=timex, treatment = treatx, replicate=repx, pvalue=f$p.value, oddsRatio=f$estimate, insulin="INS"))
#     }
    m <- with(metaO[treatment %in% c("DMSO", treatx)][celltype == "Alpha"][time == timex], table(insulin, treatment))
    m <- m[c("INS+", "INS-"),c(treatx, "DMSO")]
    f <- fisher.test(m)
    res <- rbind(res, data.table(time=timex, treatment = treatx, replicate="All", pvalue=f$p.value, oddsRatio=f$estimate, insulin="INS"))
  }
}
res
write.tsv(res, dirout(out, "INS_Pos_FisherTest_Human3.tsv"))
