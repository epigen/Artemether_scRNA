require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(gridExtra)
require(glmnet)
require(ROCR)
project.init2("artemether")


out <- "13_CellTypes/"
dir.create(dirout(out))


load(file=dirout("10_Seurat/", "DropSeq/","DropSeq.RData"))

pDat <- data.table(pbmc@meta.data, data.table(pbmc@dr$tsne@cell.embeddings, keep.rownames=TRUE))

table(pDat$res.0.5)
pDat[res.0.5 == 0, Celltype := "Alpha"]
pDat[res.0.5 == 1, Celltype := "Beta"]
pDat[res.0.5 == 2, Celltype := "Endothelial"]
pDat[res.0.5 %in% c(3), Celltype := "Ductal"]
pDat[res.0.5 %in% c(4), Celltype := "Acinar"]
pDat[res.0.5 == 5, Celltype := "Unclear1"]
pDat[res.0.5 == 6 & tSNE_1 < 0, Celltype := "Tcell"]
pDat[res.0.5 == 6 & tSNE_1 > 0, Celltype := "Dendrites"]
pDat[res.0.5 == 7, Celltype := "Unclear2"]
pDat
write.table(pDat, file=dirout(out, "MetaData.tsv"), sep="\t", quote=F, row.names=F)


markers <- c(fread("metadata/PancreasMarkers.tsv")$Human_GeneName,
fread("metadata/CellMarkers.csv")$GeneSymbol)

inPath <- dirout("10_Seurat/", "DropSeq/", "Cluster_res.0.5/")
for(f in list.files(inPath, pattern="_version2.tsv")){
  markers <- c(markers, fread(paste0(inPath, f))[V4 > 0.1]$V3)
}
markers <- unique(markers[markers %in% row.names(pbmc@data)])
meanMar <- foreach(cl = unique(pDat$Celltype)) %do% {
  apply(pbmc@data[markers,pDat[Celltype == cl]$rn], 1, mean)
}
names(meanMar) <- unique(pDat$Celltype)
meanMar <- do.call(cbind, meanMar)
# meanMar <- meanMar - apply(meanMar, 1, min)
# meanMar <- meanMar / apply(meanMar, 1, max)
pdf(dirout(out, "MarkerMatrix.pdf"), height=20, width=7, onefile=F)
# pheatmap(meanMar, color=colorRampPalette(c("lightgrey", "blue"))(10))
pheatmap(meanMar, scale="row")
dev.off()


cl.x <- "Celltype"
labelCoordinates <- pDat[,.(tSNE_1=median(tSNE_1), tSNE_2=median(tSNE_2)),by=cl.x]
ggplot(pDat, aes_string(x="tSNE_1",y="tSNE_2", color=cl.x)) + geom_point(alpha=0.5) + # ggtitle(sample.x) +
  geom_label(data=labelCoordinates, aes_string(x="tSNE_1", y="tSNE_2", label=cl.x), color="black", alpha=0.5)
ggsave(dirout(out, "Cluster_tSNE_", cl.x, ".pdf"), width=7, height=7)




# Multinomial model -------------------------------------------------------

predMeta <- pDat[sample == "DMSO" & !Alpha_Beta == "DoublePos"]
dat <- t(as.matrix(pbmc@data[-which(rownames(pbmc@data) %in% c("INS", "GCG")),predMeta$rn]))
set.seed(5001)
glmFit <- glmnet(y=factor(predMeta$Celltype), x=dat,family="multinomial",alpha=1,lambda=0.1)
str(coef(glmFit))

treat <- "DMSO"
for(treat in unique(pDat$sample)){
  message(treat)
  predMeta <- pDat[sample == treat]
  dat <- t(as.matrix(pbmc@data[-which(rownames(pbmc@data) %in% c("INS", "GCG")),predMeta$rn]))
  predMeta$class <- predict(glmFit, dat, type="class")
  jac <- matrix(NA, length(unique(predMeta$class)), length(unique(predMeta$Celltype)))
  row.names(jac) <- unique(predMeta$class)
  colnames(jac) <- unique(predMeta$Celltype)
  for(cl1 in unique(predMeta$class)){ for(cl2 in unique(predMeta$Celltype)){
    x1 <- predMeta[class == cl1]$rn 
    x2 <- predMeta[Celltype == cl2]$rn
    jac[cl1, cl2] <- length(intersect(x1, x2))/length(union(x1,x2))
  }}
  pdf(dirout(out, "Prediction_Accuracy_", treat,".pdf"), width=10, height=7, onefile=F)
  pheatmap(jac,main= paste(treat, "\nAccuracy = ", round(sum(predMeta$Celltype == predMeta$class)/nrow(predMeta),3)))
  dev.off()
}

cf <- do.call(cbind, lapply(coef(glmFit), function(m) as.matrix(m)))
colnames(cf) <- names(coef(glmFit))
pdf(dirout(out, "Prediction_Coefficient.pdf"), width=7, height=10, onefile=F)
pheatmap(as.matrix(cf[apply(cf, 1, sum) != 0 & row.names(cf) != "",]))
dev.off()

predMeta <- pDat[sample == "Artemether"]
predMeta[Alpha_Beta == "DoublePos", Celltype := "INS+ GCG+"]
dat <- t(as.matrix(pbmc@data[-which(rownames(pbmc@data) %in% c("INS", "GCG")),predMeta$rn]))
resp <- predict(glmFit, dat, type="response")[,,1]
rowAnnot <- data.frame(predMeta[,"Celltype", with=F], row.names=predMeta$rn)
pdf(dirout(out, "Prediction_Artemether_Probabilities.pdf"), width=7, height=12, onefile=F)
pheatmap(resp[order(predMeta$Celltype),], annotation_row=rowAnnot, cluster_rows=F, show_rownames=F)
dev.off()

