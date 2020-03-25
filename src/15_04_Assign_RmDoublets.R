require("project.init")
require(Seurat)
require(Matrix)
require(methods)

project.init2("artemether")

out <- "15_Counts/"
dir.create(dirout(out))

celltypes.of.interest <- c("Alpha", "Beta", "Gamma", "Delta", "Endocrine")


# ORIGINAL COUNTS ---------------------------------------------------------
metaH <- loadHMeta()
metaM <- loadMMeta()
metaH3 <- loadH3Meta()

ggplot(metaH, aes(x=celltype, fill=treatment)) + geom_bar(position="dodge") + facet_grid(replicate ~ .)
ggsave(dirout(out, "Raw_Human.pdf"), height=8,width=15)

ggplot(metaH3, aes(x=celltype, fill=treatment)) + geom_bar(position="dodge") + facet_grid(replicate ~ .)
ggsave(dirout(out, "Raw_Human3.pdf"), height=8,width=15)

ggplot(metaM, aes(x=celltype, fill=treatment)) + geom_bar(position="dodge") + facet_grid(replicate ~ .)
ggsave(dirout(out, "Raw_Mouse.pdf"), height=8,width=15)



# PREDICTED GROUPS --------------------------------------------------------
inDir <- "14_03_Classifier_moreCelltypes_noEndocrine/"
list.files(dirout(inDir))

# HUMAN
predMetaH <- fread(dirout(inDir, "human/Prediction_allCells.tsv"), header=T)
cts <- unique(intersect(colnames(predMetaH), metaH$celltype))
predMetaH$PredictedCell <- apply(as.matrix(predMetaH[,cts, with=F]),1,function(row) names(row)[which(row == max(row))])
predMetaH$PredictionMax <- apply(as.matrix(predMetaH[,cts, with=F]),1, max)
predMetaH$PredictionConf <- apply(as.matrix(predMetaH[,cts, with=F]),1,function(row) max(row)/sort(row, decreasing=T)[2])
metaH <- merge(metaH, predMetaH[,c("rn", "PredictedCell", "PredictionConf", "PredictionMax")], by="rn", all.x=T)
stopifnot(nrow(metaH[is.na(PredictionMax)]) == 0)
# metaH[,celltype2 := celltype]
# metaH[celltype == "Endocrine", celltype2 := PredictedCell]

# HUMAN 3
predMetaH3 <- fread(dirout(inDir, "human3/Prediction_allCells.tsv"), header=T)
cts <- unique(intersect(colnames(predMetaH3), metaH3$celltype))
predMetaH3$PredictedCell <- apply(as.matrix(predMetaH3[,cts, with=F]),1,function(row) names(row)[which(row == max(row))])
predMetaH3$PredictionMax <- apply(as.matrix(predMetaH3[,cts, with=F]),1, max)
predMetaH3$PredictionConf <- apply(as.matrix(predMetaH3[,cts, with=F]),1,function(row) max(row)/sort(row, decreasing=T)[2])
metaH3 <- merge(metaH3, predMetaH3[,c("rn", "PredictedCell", "PredictionConf", "PredictionMax")], by="rn", all.x=T)
stopifnot(nrow(metaH3[is.na(PredictionMax)]) == 0)
# metaH3[,celltype2 := celltype]
# metaH3[celltype == "Endocrine", celltype2 := PredictedCell]

# MOUSE
predMetaM <- fread(dirout(inDir, "mouse/Prediction_allCells.tsv"), header=T)
cts <- unique(intersect(colnames(predMetaM), metaH3$celltype))
predMetaM$PredictedCell <- apply(as.matrix(predMetaM[,cts, with=F]),1,function(row) names(row)[which(row == max(row))])
predMetaM$PredictionMax <- apply(as.matrix(predMetaM[,cts, with=F]),1, max)
predMetaM$PredictionConf <- apply(as.matrix(predMetaM[,cts, with=F]),1,function(row) max(row)/sort(row, decreasing=T)[2])
metaM <- merge(metaM, predMetaM[,c("rn", "PredictedCell", "PredictionConf", "PredictionMax")], by="rn", all.x=T)
stopifnot(nrow(metaM[is.na(PredictionMax)]) == 0)
# metaM[,celltype2 := celltype]
# metaM[celltype == "Endocrine", celltype2 := PredictedCell]


# Identify cutoff for doublets --------------------------------------------
metas <- list(human=metaH, mouse=metaM, human3=metaH3)
metas$human <- metas$human[!grepl("II_GABA", sample)]
for(mnam in names(metas)){
  metas[[mnam]][,celltype2 := celltype]
  metas[[mnam]][celltype == "Endocrine", celltype2 := PredictedCell]
  
  for(sx in unique(metas[[mnam]]$sample)){
    for(ctx in unique(metas[[mnam]][sample == sx]$celltype2)){
      x <- metas[[mnam]][sample == sx][celltype2 == ctx]$nGene
      nGene.peak.x <- if(length(x) > 10){
        d <- density(x)
        d$x[which(d$y == max(d$y))]
      } else {
        mean(x)
      }
      metas[[mnam]][sample == sx & celltype2 == ctx, nGene.peak := nGene.peak.x]
    }
  }
  
  metas[[mnam]][, nGene.doublet := "Keep"]
  metas[[mnam]][(nGene > (nGene.peak + (nGene.peak))), nGene.doublet := "nGene.too.high"]
  metas[[mnam]][(nGene < (nGene.peak - (nGene.peak / 2))), nGene.doublet := "nGene.too.low"]
}

for(mnam in names(metas)){
  ggplot(metas[[mnam]][celltype2 %in% celltypes.of.interest], aes(x=nGene, fill=nGene.doublet)) + 
    geom_histogram() + facet_grid(treatment + replicate ~  celltype2)
  ggsave(dirout(out, "Doublets_nGene_Cutoff_",mnam,".pdf"), w=25,h=25)
  
  ggplot(metas[[mnam]][celltype != "Endocrine"], aes(x=PredictionConf, color=nGene.doublet)) + 
    geom_density() + facet_grid(replicate ~ treatment) + scale_x_log10()
  ggsave(dirout(out, "Doublets_PredictionConf_",mnam,".pdf"), w=15,h=15)
}

for(mnam in names(metas)){
  ggplot(metas[[mnam]][celltype2 %in% celltypes.of.interest], aes(x=treatment, fill=nGene.doublet)) + geom_bar() + facet_grid(celltype2 ~ replicate, scales="free_y")
  ggsave(dirout(out, "Doublets_nGene_Numbers_",mnam,".pdf"), w=15,h=15)
  
  metas[[mnam]][is.na(PredictionConf), PredictionConf := 10]
  ggplot(metas[[mnam]][celltype2 %in% celltypes.of.interest], aes(x=(PredictionConf > 3), fill=nGene.doublet)) + geom_bar(position="stack") + 
    facet_grid(celltype2 ~ replicate + treatment, scales="free_y") + xRot()
  ggsave(dirout(out, "Doublets_nGene_PredictionConf_Numbers_",mnam,".pdf"), w=15,h=15)
}



# # PLOT ALPHA PROBABILITY --------------------------------------------------
# ggplot(predMetaM, aes(x=Alpha, color=treatment)) + stat_ecdf() + facet_grid(replicate ~ .)
# ggsave(dirout(out, "Alpha_cell_Probablity_Mouse.pdf"), width=4, height=12)
# 
# ggplot(predMetaH, aes(x=Alpha, color=treatment)) + stat_ecdf() + facet_grid(replicate ~ .)
# ggsave(dirout(out, "Alpha_cell_Probablity_Human.pdf"), width=4, height=8)
# 
# ggplot(predMetaH3, aes(x=Alpha, color=treatment)) + stat_ecdf() + facet_grid(replicate ~ .)
# ggsave(dirout(out, "Alpha_cell_Probablity_Human3.pdf"), width=4, height=8)
# 



# Compare to ambiguous cells ----------------------------------------------
amb.file <- dirout(out, "AmbientCellNumbers.tsv")
if(!file.exists(amb.file)){  
  ff <- list.files(dirout("06_SpikeIns_Islets/"), pattern="IdentifyOrganism_full.tsv", recursive=T)
  names(ff) <- gsub("(\\_\\d)?(\\_S\\d+)?\\_hmG", "", dirname(ff))
  ambDT <- do.call(rbind, lapply(names(ff), function(fnam) data.table(fread(dirout("06_SpikeIns_Islets/",ff[[fnam]]))[,.N, by="org"], sample=fnam)))
  write.tsv(ambDT, amb.file)
} else {
  ambDT <- fread(amb.file)
}
ambDT[, orgSample := ifelse(grepl("hIslets", sample), "human", "mouse")]
ambDT <- merge(ambDT[is.na(org)], ambDT[org != orgSample | is.na(org)][,.(TotalSum = sum(N)), by="sample"])
ambDT[, percentAmbiguous := N/TotalSum * 100]
ambDT[, sample := gsub("hIslets_III_A(\\d+)_(.+)hr", "hIslets_III_A10_\\2_\\1", sample)]
ambDT[, sample := gsub("hIslets_III_DMSO(\\d+)_(.+)hr", "hIslets_III_DMSO_\\2_\\1", sample)]

# Count removed cells from cutoff
mDT <- do.call(rbind, lapply(names(metas), function(xnam) data.table(metas[[xnam]], ds=xnam)))
mDT2 <- merge(mDT[PredictionConf > 3 & nGene.doublet != "nGene.too.high"][,.N, by="sample"], mDT[,.N, by="sample"], by="sample")
mDT2[, doubletsRmPercent := 100*(1-N.x/N.y)]

# plot, test, export
pDT <- merge(mDT2, ambDT, by="sample")[sample != "hIslets_II_FoxO"]
pDT <- merge(pDT, unique(mDT[,c("sample", "ds"),with=F]))
ggplot(pDT, aes(x=doubletsRmPercent, y=percentAmbiguous)) + geom_point() + facet_wrap(~ds, scales="free") + 
  theme_bw(12)
ggsave(dirout(out, "Removed_vs_Ambiguous.pdf"), w=10,h=4)
write.tsv(pDT, dirout(out, "Removed_vs_Ambiguous.tsv"))
write.tsv(pDT[, cor(doubletsRmPercent, percentAmbiguous), by="ds"], dirout(out, "Removed_vs_Ambiguous_Correlation.tsv"))



# Plot predicted and Numbers ----------------------------------------------
metaH <- metas$human
metaH3 <- metas$human3
metaM <- metas$mouse

pDat <- metaH[,.N, by=c("celltype2", "replicate", "treatment", "celltype")]
pDat[,sum := sum(N), by=c("replicate", "treatment")]
pDat[,percent := N / sum * 100]
ggplot(pDat, aes(x=treatment,y=percent, fill=celltype=="Endocrine")) + geom_bar(stat="identity", position="stack") + facet_grid(replicate ~ celltype2)  +
  scale_fill_manual(values=c("grey", "blue")) + theme_bw(12)+ xRot()
ggsave(dirout(out, "Predicted_Human.pdf"), height=8,width=15)

pDat <- metaH3[,.N, by=c("celltype2", "replicate", "treatment", "celltype")]
pDat[,sum := sum(N), by=c("replicate", "treatment")]
pDat[,percent := N / sum * 100]
ggplot(pDat, aes(x=treatment,y=percent, fill=celltype=="Endocrine")) + geom_bar(stat="identity", position="stack") + facet_grid(replicate ~ celltype2)  +
  scale_fill_manual(values=c("grey", "blue")) + theme_bw(12)+ xRot()
ggsave(dirout(out, "Predicted_Human3.pdf"), height=8,width=15)

pDat <- metaM[,.N, by=c("celltype2", "replicate", "treatment", "celltype")]
pDat[,sum := sum(N), by=c("replicate", "treatment")]
pDat[,percent := N / sum * 100]
ggplot(pDat, aes(x=treatment,y=percent, fill=celltype=="Endocrine")) + geom_bar(stat="identity", position="stack") + facet_grid(replicate ~ celltype2)  +
  scale_fill_manual(values=c("grey", "blue")) + theme_bw(12)+ xRot()
ggsave(dirout(out, "Predicted_Mouse.pdf"), height=8,width=15)


# ANNOTATE WITH INS / GCG EXPRESSION DATA -------------------------------------------------
pbmcH <- loadHData()
pbmcM <- loadMData()
pbmcH3 <- loadH3Data()

assignGenesMeta <- function(genes, metaX, seuratObj){ 
  for(gg in genes){
    if(!gg %in% row.names(seuratObj@data)) stop(gg, " not found in sc data")
    metaX[[gg]] <- seuratObj@data[gg,metaX$rn]
  }
  metaX
}
metaH <- assignGenesMeta(c("INS", "GCG"), metaH, pbmcH)
metaH3 <- assignGenesMeta(c("INS", "GCG"), metaH3, pbmcH3)
metaM <- assignGenesMeta(c("Ins2", "Gcg"), metaM, pbmcM)

# EXPORT UNFILTERED -------------------------------------------------
write.tsv(metaH, dirout(out,"Reassigned_Human_unfiltered.tsv"))
write.tsv(metaM, dirout(out,"Reassigned_Mouse_unfiltered.tsv"))
write.tsv(metaH3, dirout(out,"Reassigned_Human3_unfiltered.tsv"))


# FILTER -------------------------------------------------
metaH <- metaH[PredictionConf > 3 & nGene.doublet != "nGene.too.high"]
metaH3 <- metaH3[PredictionConf > 3 & nGene.doublet != "nGene.too.high"]
metaM <- metaM[PredictionConf > 3 & nGene.doublet != "nGene.too.high"]

write.tsv(metaH, dirout(out,"Reassigned_Human.tsv"))
write.tsv(metaM, dirout(out,"Reassigned_Mouse.tsv"))
write.tsv(metaH3, dirout(out,"Reassigned_Human3.tsv"))


# PLOT DISTRIBUTIONS -------------------------------------------------
celltypes.of.interest <- c("Alpha", "Beta", "Gamma", "Delta", "Endocrine")
treat <- "A10"
gg <- "INS"
ct <- "Beta"
org <- "human3"
outGG <- paste0(out, "Details/")
dir.create(dirout(outGG))
for(org in c("human", "mouse", "human3")){
  meta <- metaH3
  if(org == "mouse") meta <- metaM
  if(org == "human") meta <- metaH
  for(treat in unique(meta[treatment != "DMSO"]$treatment)){
    for(gg in c("INS", "GCG", "Ins2", "Gcg")){
      if(!gg %in% colnames(meta)) next
      pDat <- meta[treatment %in% c("DMSO", treat) & celltype2 %in% celltypes.of.interest]
      pDat$treatment <- factor(pDat$treatment, levels = c("DMSO", treat))
      ggplot(pDat, aes_string(x=gg, color="treatment")) + geom_density() +
        theme_bw(12) + facet_grid(replicate ~ celltype2, scale="free") +
        scale_color_manual(values=c("black", "red"))
      ggsave(dirout(out, org, "_", treat, "_", gg, ".pdf"),height=9,width=16)
      
      ggplot(pDat, aes_string(x=gg, color="treatment")) + stat_ecdf() +
        theme_bw(12) + facet_grid(replicate ~ celltype2, scale="free") +
        scale_color_manual(values=c("black", "red"))
      ggsave(dirout(out, org, "_", treat, "_", gg, "_ECDF.pdf"),height=9,width=16)
      
      pDat$Expression <- pDat[[gg]]
      ggplot(pDat, aes(x=Expression, color=treatment, alpha=(celltype != "Endocrine"))) + stat_ecdf() +
        theme_bw(12) + facet_grid(replicate ~ celltype2, scale="free") +
        scale_color_manual(values=c("black", "red")) + xlab(gg)
      ggsave(dirout(out, org, "_", treat, "_", gg, "_ECDF_Endocrine.pdf"),height=9,width=16)
      
      for(ct in celltypes.of.interest){
        if(!gg %in% colnames(meta)) next
        pDat <- meta[treatment %in% c("DMSO", treat) & celltype2 %in% ct]
        if(nrow(pDat) == 0) next
        pDat$treatment <- factor(pDat$treatment, levels = c("DMSO", treat))
        pDat$celltype <- factor(pDat$celltype, levels = c(ct, "Endocrine"))
        ggplot(pDat, aes_string(x=gg, color="treatment", linetype = "celltype")) + geom_density() +
          theme_bw(12) + facet_grid(replicate ~ celltype2, scale="free") +
          scale_color_manual(values=c("black", "red"))
        ggsave(dirout(outGG, org, "_", treat, "_", ct, "_", gg, ".pdf"),height=6,width=6)
        
        ggplot(pDat, aes_string(x=gg, color="treatment", linetype = "celltype")) + stat_ecdf() +
          theme_bw(12) + facet_grid(replicate ~ celltype2, scale="free") +
          scale_color_manual(values=c("black", "red"))
        ggsave(dirout(outGG, org, "_", treat, "_", ct, "_", gg, ".pdf"),height=6,width=6)
        
      }
    }
  }
}


pbmc <- SubsetData(object=pbmcH, cells.use=metaH$rn)
save(pbmc, file=dirout(out, "SeuratObject_human.RData"))

pbmc <- SubsetData(object=pbmcM, cells.use=metaM$rn)
save(pbmc, file=dirout(out, "SeuratObject_mouse.RData"))

pbmc <- SubsetData(object=pbmcH3, cells.use=metaH3$rn)
save(pbmc, file=dirout(out, "SeuratObject_human3.RData"))

# metaH$org <- "human"
# metaM$org <- "mouse"
# 
# metaX <- rbind(metaH, metaM, fill=T)
# metaX[,Insulin :=  sum(INS, Ins2, na.rm=T), by=c("rn", "org")]
# metaX[,Glucagon :=  sum(GCG, Gcg, na.rm=T), by=c("rn", "org")]
# ggplot(metaX[treatment %in% c("DMSO", "A10") & celltype2 %in% c("Beta", "Alpha")], aes(x=replicate, y=exp(Insulin) - 1, fill=treatment)) + geom_violin() +
#   facet_grid(celltype2 ~ org, scales="free") + theme_bw(12)
# ggsave(dirout(out, "Org_Compare_Insulin.pdf"), width=8, height=4)
# 
# ggplot(metaX[treatment %in% c("DMSO", "A10") & celltype2 %in% c("Beta", "Alpha")], aes(x=replicate, y=exp(Glucagon) - 1, fill=treatment)) + geom_violin() +
#   facet_grid(celltype2 ~ org, scales="free") + theme_bw(12)
# ggsave(dirout(out, "Org_Compare_Glucagon.pdf"), width=8, height=4)
# 
# write.tsv(metaX, dirout(out,"Reassigned_Combined.tsv"))
