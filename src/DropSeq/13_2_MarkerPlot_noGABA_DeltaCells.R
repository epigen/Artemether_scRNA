require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(gridExtra)
require(glmnet)
require(ROCR)
project.init2("artemether")


out <- "13_2_CellTypes_noGABA_DeltaCells/"
dir.create(dirout(out))


load(file=dirout("10_Seurat/", "DropSeq/","DropSeq.RData"))

(load(dirout("12_2_INS_GLC_noGABA_DeltaCells/Artemether_groups.RData")))

table(pDat$res.0.5)
pDat[,Celltype := "NA"]
pDat[res.0.5 == 2, Celltype := "Endothelial"]
pDat[res.0.5 %in% c(3), Celltype := "Ductal"]
pDat[res.0.5 %in% c(4), Celltype := "Acinar"]
pDat[res.0.5 == 5, Celltype := "Unclear1"]
pDat[res.0.5 == 6 & tSNE_1 < 0, Celltype := "Tcell"]
pDat[res.0.5 == 6 & tSNE_1 > 0, Celltype := "Dendrites"]
pDat[res.0.5 == 7, Celltype := "Unclear2"]
pDat[group == "GCG+", Celltype := "GCG+"]
pDat[group == "INS+", Celltype := "INS+"]
pDat[group == "INS+ GCG+", Celltype := "INS+ GCG+"]
pDat[Celltype == "NA", Celltype := NA]
pDat[SST > 6, Celltype := "Delta"]
table(pDat$Celltype)
write.table(pDat, file=dirout(out, "MetaData.tsv"), sep="\t", quote=F, row.names=F)


# T-SNE

cl.x <- "Celltype"
labelCoordinates <- pDat[,.(tSNE_1=median(tSNE_1), tSNE_2=median(tSNE_2)),by=cl.x]
ggplot(pDat, aes_string(x="tSNE_1",y="tSNE_2", color=cl.x)) + geom_point(alpha=0.5) + theme_bw(24) + 
  geom_label(data=labelCoordinates, aes_string(x="tSNE_1", y="tSNE_2", label=cl.x), color="black", alpha=0.5)
ggsave(dirout(out, "Cluster_tSNE_", cl.x, ".pdf"), width=8, height=7)


# Percentages

xx <- pDat[,.N, by=c("Celltype", "sample")]
xx[, sum := sum(N), by="sample"]
xx[,Percentage := N/sum * 100]
ggplot(xx, aes(x=Celltype, fill=sample, y=Percentage)) + geom_bar(position="dodge", stat="identity") +
  theme_bw(24) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) 
ggsave(dirout(out, "Celltype_Percentages.pdf"), height=7, width=7)

pDat <- pDat[!is.na(Celltype)]# & !Celltype %in% c("Dendrites", "Unclear2", "Tcell")]


# Plot unbiased markers

markers <- c("CD3D", "ARX","CRYBA2", "DLK1","GCG","IAPP","INS","KRT8", "KRT19", "LOXL4", "SST", "NPTX2", "PRSS1", "RBP4","REG1A","SPARC", "SPP1","TM4SF4","TPSAB1","TTR")
# markers <- c(fread("metadata/PancreasMarkers.tsv")$Human_GeneName,
# fread("metadata/CellMarkers.csv")$GeneSymbol)
pList <- list()
for(mm in unique(markers)){
  x <- data.table(pDat,Expression=pbmc@data[mm, pDat$cell])
  pList[[mm]] <- ggplot(x, aes(x=tSNE_1, y=tSNE_2, color=Expression)) + 
    geom_point(alpha=0.5) + ggtitle(mm) + theme_bw(12) + guides(color=F) +
    scale_color_gradient(low="lightgrey", high="blue")
}
plt <- grid.arrange(grobs=pList, ncol=4)
ggsave(dirout(out, "Marker_Tsne.jpg"),width=16, height=20, plot=plt)
ggsave(dirout(out, "Marker_Tsne.pdf"),width=16, height=20, plot=plt)

inPath <- dirout("10_Seurat/", "DropSeq/", "Cluster_res.0.5/")
  for(f in list.files(inPath, pattern="_version2.tsv")){
  markers <- c(markers, fread(paste0(inPath, f))[V4 > 0.1]$V3)
}
markers <- unique(markers[markers %in% row.names(pbmc@data)])
meanMar <- foreach(cl = unique(pDat$Celltype)) %do% {
  apply(pbmc@data[markers,pDat[Celltype == cl]$cell], 1, mean)
}
names(meanMar) <- unique(pDat$Celltype)
meanMar <- do.call(cbind, meanMar)
pdf(dirout(out, "MarkerMatrix.pdf"), height=15, width=7, onefile=F)
pheatmap(meanMar[,colnames(meanMar) != "NA"], scale="row")
dev.off()


# DELTA MARKERS
xx <- pDat[is.na(Celltype) | grepl("\\+", Celltype) | Celltype == "Delta"]
xx[is.na(Celltype),Celltype := "NA"]
pdf(dirout(out, "DeltaMarkers.pdf"), height=5, width=10, onefile=F)
pheatmap(pbmc@data[c("RBP4", "SST", "PCSK1", "HHEX", "INS", "GCG", "TTR", "IAPP","DLK1"), xx[order(Celltype)]$cell], 
         border_color=NA,cluster_rows=F,cluster_cols=F, show_colnames=F,color=colorRampPalette(c("black", "yellow"))(50),
         annotation_col = data.frame(Celltype = xx$Celltype, row.names=xx$cell))
dev.off()
for(g in c("RBP4", "SST", "INS")){
  xx[[g]] <- pbmc@data[g, xx$cell]
}
ggplot(xx, aes(x=RBP4, y=SST,color=Celltype, shape=sample, size=nGene)) + geom_point(alpha=0.8) + theme_bw(24)
ggsave(dirout(out, "DeltaMarkers_SST_RBP4.pdf"), height=7, width=7)
ggplot(xx, aes(x=INS, y=SST,color=Celltype, shape=sample, size=nGene)) + geom_point(alpha=0.8) + theme_bw(24)
ggsave(dirout(out, "DeltaMarkers_SST_INS.pdf"), height=7, width=7)
ggplot(xx, aes(x=GCG, y=SST,color=Celltype, shape=sample, size=nGene)) + geom_point(alpha=0.8) + theme_bw(24)
ggsave(dirout(out, "DeltaMarkers_SST_GCG.pdf"), height=7, width=7)
ggplot(xx, aes(x=SST)) + geom_density(fill="lightgrey") + theme_bw(16)+ ylab("Density")
ggsave(dirout(out, "DeltaMarkers_SST.pdf"), height=7, width=7) 




# Multinomial model -------------------------------------------------------
predGeneExclude <- c("INS", "GCG", "SST")
predMeta <- pDat[sample == "DMSO" & !is.na(Celltype) & Celltype != "INS+ GCG+"]
predMeta <- predMeta[Celltype %in% predMeta[,.N, by="Celltype"][N > 20]$Celltype]
table(predMeta$Celltype)
dat <- t(as.matrix(pbmc@data[-which(rownames(pbmc@data) %in% predGeneExclude),predMeta$cell]))
set.seed(5001)
glmFit <- glmnet(y=factor(predMeta$Celltype), x=dat,family="multinomial",alpha=1,lambda=0.05)
str(coef(glmFit))
save(glmFit, file=dirout(out, "ModelFit.RData"))
load(dirout(out, "ModelFit.RData"))

treat <- "DMSO"
for(treat in unique(pDat$sample)){
  message(treat)
  predMeta <- pDat[sample == treat & Celltype %in% names(coef(glmFit))]
  dat <- t(as.matrix(pbmc@data[-which(rownames(pbmc@data) %in% predGeneExclude),predMeta$cell]))
  predMeta$class <- predict(glmFit, dat, type="class")
  jac <- matrix(NA, length(unique(predMeta$class)), length(unique(predMeta$Celltype)))
  row.names(jac) <- unique(predMeta$class)
  colnames(jac) <- unique(predMeta$Celltype)
  for(cl1 in unique(predMeta$class)){ for(cl2 in unique(predMeta$Celltype)){
    x1 <- predMeta[class == cl1]$cell 
    x2 <- predMeta[Celltype == cl2]$cell
    jac[cl1, cl2] <- length(intersect(x1, x2))/length(union(x1,x2))
  }}
  pdf(dirout(out, "Prediction_Accuracy_", treat,".pdf"), width=5, height=4.5, onefile=F)
  pheatmap(jac[sort(row.names(jac)), sort(colnames(jac))],main= paste(treat, "\nAccuracy = ", round(sum(predMeta$Celltype == predMeta$class)/nrow(predMeta),3)),
           color=colorRampPalette(c("lightgrey", "blue"))(50), cluster_rows=F, cluster_cols=F)
  dev.off()
}

cf <- do.call(cbind, lapply(coef(glmFit), function(m) as.matrix(m)))
colnames(cf) <- names(coef(glmFit))
pdf(dirout(out, "Prediction_Coefficient.pdf"), width=5, height=10, onefile=F)
pheatmap(as.matrix(cf[apply(cf, 1, sum) != 0 & row.names(cf) != "",]), color=colorRampPalette(c("lightgrey", "blue"))(50))
dev.off()

treat <- "Artemether"
for(treat in unique(pDat$sample)){
  predMeta <- pDat[sample == treat & !is.na(Celltype)]
  dat <- t(as.matrix(pbmc@data[-which(rownames(pbmc@data) %in% predGeneExclude),predMeta$cell]))
  resp <- predict(glmFit, dat, type="response")[,,1]
  rowAnnot <- data.frame(predMeta[,"Celltype", with=F], row.names=predMeta$cell)
  pdf(dirout(out, "Prediction_",treat,"_Probabilities.pdf"), width=4, height=7, onefile=F)
  pheatmap(resp[order(predMeta$Celltype),], annotation_row=rowAnnot, cluster_rows=F, show_rownames=F,
           color=colorRampPalette(c("lightgrey", "blue"))(50))
  dev.off()

  stopifnot(all(predMeta$cell == row.names(resp)))
  predMeta$Beta_Probability <- resp[,"INS+"]
  predMeta$Alpha_Probability <- resp[,"GCG+"]
  
  pList <- lapply(unique(predMeta$Celltype), function(gg){
    ggplot(predMeta, aes(x=Beta_Probability, y=Alpha_Probability)) + 
    geom_point(alpha=0.7, color="lightgrey") + theme_bw(16) +
    geom_point(data=predMeta[Celltype == gg], alpha=0.7, color="black") +
    xlab("Beta cell class probability") + ylab("Alpha cell class probability") + ggtitle(gg)
    })
  names(pList) <- unique(predMeta$Celltype)
  plt <- gridExtra::grid.arrange(grobs=c(pList[c("INS+", "GCG+", "INS+ GCG+", "Delta")],deltaPlot), ncol=5)
  ggsave(dirout(out, "Prediction_",treat,"Probabilities_Delta.pdf"), width=16, height=4, plot=plt)
  
  plt <- gridExtra::grid.arrange(grobs=pList[c("INS+", "GCG+", "INS+ GCG+", "Endothelial")], ncol=4)
  ggsave(dirout(out, "Prediction_",treat,"Probabilities_Endo.pdf"), width=16, height=4, plot=plt)
  
  plt <- gridExtra::grid.arrange(grobs=pList, ncol=3)
  ggsave(dirout(out, "Prediction_",treat,"Probabilities_All.pdf"), width=13, height=13, plot=plt)
  
  pD <- data.table(resp, Celltype=predMeta$Celltype)
  colnames(pD) <- gsub("\\+", "", colnames(pD))
  pList2 <- lapply(gsub("\\+", "", colnames(resp)), function(grp){
    print(grp)
    ggplot(pD[grepl("\\+", Celltype) | Celltype == "Delta"], aes_string(x=grp, color="Celltype")) + 
    stat_ecdf(size=2) + ylab("Fraction") + theme_bw(16) + xlab("Delta cell class probability") + ggtitle(grp)
  })
  names(pList2) <- colnames(resp)
  ggsave(dirout(out, "Prediction_", treat, "_DeltaProb.pdf"), height=7, width=8, plot=pList2$Delta)
  plt <- gridExtra::grid.arrange(grobs=lapply(pList[c("INS+", "GCG+", "Delta")], function(p) p + guides(color=F)), ncol=3)
  
  ggplot(data.table(melt(pD,id.vars="Celltype",measure.vars=c("GCG", "INS")))[Celltype %in% c("GCG+", "INS+", "INS+ GCG+")], 
         aes(x=Celltype, y=value)) + geom_violin() + facet_grid(. ~ variable)
  
  pList[["DeltaProb"]] <- pList2[["Delta"]] + guides(color=F)
  plt <- gridExtra::grid.arrange(grobs=pList[c("INS+", "GCG+", "INS+ GCG+", "Endothelial", "DeltaProb")], ncol=5)
  ggsave(dirout(out, "Prediction_",treat,"Probabilities_Endo2.pdf"), width=16, height=4, plot=plt)
}