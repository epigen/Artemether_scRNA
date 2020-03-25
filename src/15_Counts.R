require("project.init")
require(Seurat)
require(Matrix)
require(methods)

project.init2("artemether")

out <- "15_Counts/"
dir.create(dirout(out))



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
predMetaH <- fread(dirout(inDir, "human/Prediction.tsv"), header=T)
predMetaH$PredictedCell <- apply(as.matrix(predMetaH[,unique(intersect(colnames(predMetaH), metaH$celltype)), with=F]),1,function(row) names(row)[which(row == max(row))])
metaH <- merge(metaH, predMetaH[,c("rn", "PredictedCell")], by="rn", all.x=T)
metaH[,celltype2 := celltype]
metaH[celltype == "Endocrine", celltype2 := PredictedCell]

predMetaH3 <- fread(dirout(inDir, "human3/Prediction.tsv"), header=T)
predMetaH3$PredictedCell <- apply(as.matrix(predMetaH3[,unique(intersect(colnames(predMetaH3), metaH3$celltype)), with=F]),1,function(row) names(row)[which(row == max(row))])
metaH3 <- merge(metaH3, predMetaH3[,c("rn", "PredictedCell")], by="rn", all.x=T)
metaH3[,celltype2 := celltype]
metaH3[celltype == "Endocrine", celltype2 := PredictedCell]

predMetaM <- fread(dirout(inDir, "mouse/Prediction.tsv"), header=T)
x <- as.matrix(predMetaM[,unique(intersect(colnames(predMetaM), metaM$celltype)), with=F])
head(x)
head(apply(x,1,max))
x2 <- apply(x,1,sort)
x2 <- t(x2)
head(x2)
quantile(x2[,8]-x2[,7])
predMetaM$PredictedCell <- apply(as.matrix(predMetaM[,unique(intersect(colnames(predMetaM), metaM$celltype)), with=F]),1,function(row) names(row)[which(row == max(row))])
metaM <- merge(metaM, predMetaM[,c("rn", "PredictedCell")], by="rn", all.x=T)
metaM[,celltype2 := celltype]
metaM[celltype == "Endocrine", celltype2 := PredictedCell]


ggplot(predMetaM, aes(x=Alpha, color=treatment)) + stat_ecdf() + facet_grid(replicate ~ .)
ggsave(dirout(out, "Alpha_cell_Probablity_Mouse.pdf"), width=4, height=12)

ggplot(predMetaH, aes(x=Alpha, color=treatment)) + stat_ecdf() + facet_grid(replicate ~ .)
ggsave(dirout(out, "Alpha_cell_Probablity_Human.pdf"), width=4, height=8)

ggplot(predMetaH3, aes(x=Alpha, color=treatment)) + stat_ecdf() + facet_grid(replicate ~ .)
ggsave(dirout(out, "Alpha_cell_Probablity_Human3.pdf"), width=4, height=8)


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




# CORRECTED DISTRIBUTIONS -------------------------------------------------
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

write.tsv(metaH, dirout(out,"Reassigned_Human.tsv"))
write.tsv(metaM, dirout(out,"Reassigned_Mouse.tsv"))
write.tsv(metaH3, dirout(out,"Reassigned_Human3.tsv"))

metaH$org <- "human"
metaM$org <- "mouse"

metaX <- rbind(metaH, metaM, fill=T)
metaX[,Insulin :=  sum(INS, Ins2, na.rm=T), by=c("rn", "org")]
metaX[,Glucagon :=  sum(GCG, Gcg, na.rm=T), by=c("rn", "org")]
ggplot(metaX[treatment %in% c("DMSO", "A10") & celltype2 %in% c("Beta", "Alpha")], aes(x=replicate, y=exp(Insulin) - 1, fill=treatment)) + geom_violin() +
  facet_grid(celltype2 ~ org, scales="free") + theme_bw(12)
ggsave(dirout(out, "Org_Compare_Insulin.pdf"), width=8, height=4)

ggplot(metaX[treatment %in% c("DMSO", "A10") & celltype2 %in% c("Beta", "Alpha")], aes(x=replicate, y=exp(Glucagon) - 1, fill=treatment)) + geom_violin() +
  facet_grid(celltype2 ~ org, scales="free") + theme_bw(12)
ggsave(dirout(out, "Org_Compare_Glucagon.pdf"), width=8, height=4)

write.tsv(metaX, dirout(out,"Reassigned_Combined.tsv"))
