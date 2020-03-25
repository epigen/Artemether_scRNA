require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(gridExtra)
require(glmnet)
require(ROCR)
project.init2("artemether")


out <- "12_INS_GLC/"
dir.create(dirout(out))

load(file=dirout("10_Seurat/", "DropSeq/","DropSeq.RData"))

pDat <- data.table(cell = colnames(pbmc@data), INS=pbmc@data["INS",], GCG=pbmc@data["GCG",], 
                   #NGN3=pbmc@data["NEUROG3",], 
                   Sample=pbmc@meta.data$sample, nGene = pbmc@meta.data$nGene, nUMI=pbmc@meta.data$nUMI, 
                   pbmc@meta.data,
                   pbmc@dr$tsne@cell.embeddings)
table(pDat$Sample)

# DO THE SAME ON THE RAW DATA ---------------------------------------------
# pDat <- data.table(INSraw=pbmc@raw.data["INS",], GCGraw=pbmc@raw.data["GCG",], Sample=gsub("(.+?)\\..+", "\\1", colnames(pbmc@raw.data)), nGene = apply(pbmc@raw.data != 0 , 2, sum), nUMI=apply(pbmc@raw.data, 2, sum))
# pDat[, INS := log(((INSraw)/nUMI) * 1e4)]
# pDat[INS == -Inf, INS := 0]
# pDat[, GCG := log(((GCGraw)/nUMI) * 1e4)]
# pDat[GCG == -Inf, GCG := 0]
# table(pDat$Sample)
# out <- paste0("12_INS_GLC/", "FullData/")
# dir.create(dirout(out))

cutoff <- 5
pDat[,group := "Negative"]
pDat[INS > cutoff,group := "INS+"]
pDat[GCG > cutoff,group := "GCG+"]
pDat[INS > cutoff & GCG > cutoff,group := "INS+ GCG+"]
pDat



# PERCENT OF ALL READS ----------------------------------------------------
# 
# dmsoDat <- pbmc@raw.data[row.names(pbmc@data),pDat[sample == "DMSO"]$cell]
# sum(dmsoDat["INS",])/sum(dmsoDat)
# sum(dmsoDat["GCG",])/sum(dmsoDat)
# 
# arteDat <- pbmc@raw.data[row.names(pbmc@data),pDat[sample == "Artemether"]$cell]
# sum(arteDat["INS",])/sum(arteDat)
# sum(arteDat["GCG",])/sum(arteDat)


# FACS LIKE PLOTS ---------------------------------------------------------

ggplot(pDat, aes(x=INS, y=GCG, color=nUMI, size=nGene)) + geom_point(alpha=0.5) + facet_grid(. ~ Sample) + theme_bw(24)
ggsave(dirout(out, "INS_GCG_nGene_nUMI.pdf"),width=16, height=5)

(p <- ggplot(pDat[INS != 0 & GCG != 0], aes(x=INS, y=GCG)) + geom_hex() + facet_grid(. ~ Sample) + theme_bw(24))
ggsave(dirout(out, "INS_GCG_hex.pdf"),width=16, height=5)

(p <- ggplot(pDat, aes(x=INS, y=GCG)) + geom_point(alpha=0.5) + facet_grid(. ~ Sample) + theme_bw(24))
ggsave(dirout(out, "INS_GCG.pdf"),width=16, height=5)

(p + geom_hline(yintercept=4, color="blue", size=1, alpha=0.5) + 
   geom_vline(xintercept=4, color="blue", size=1, alpha=0.5) +
   geom_hline(yintercept=6, color="blue", size=1, alpha=0.5) + 
   geom_vline(xintercept=6, color="blue", size=1, alpha=0.5)
)
ggsave(dirout(out, "INS_GCG_lines2.pdf"),width=16, height=5)

(p + geom_hline(yintercept=cutoff, color="red", size=2, alpha=0.5) + 
   geom_vline(xintercept=cutoff, color="red", size=2, alpha=0.5))
ggsave(dirout(out, "INS_GCG_lines.pdf"),width=16, height=5)



# DENSITY PLOTS -----------------------------------------------------------

# INSULIn
ggplot(pDat, aes(x=INS, color=sample)) + stat_ecdf(size=1.5) + theme_bw(24)
ggsave(dirout(out, "INS_ECDF.pdf"),width=8, height=7)
ggplot(pDat, aes(x=INS, fill=sample)) + geom_density(alpha=0.3) + theme_bw(24)
ggsave(dirout(out, "INS_Density.pdf"),width=8, height=7)

# GLUCAGON
ggplot(pDat, aes(x=GCG, color=sample)) + stat_ecdf(size=1.5) + theme_bw(24)
ggsave(dirout(out, "GCG_ECDF.pdf"),width=8, height=7)
ggplot(pDat, aes(x=GCG, fill=sample)) + geom_density(alpha=0.3) + theme_bw(24)
ggsave(dirout(out, "GCG_Density.pdf"),width=8, height=7)


# t-SNE PLOTS -----------------------------------------------------------

ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=INS)) + geom_point(alpha=0.5) + 
  facet_grid(. ~ Sample) + scale_color_gradient(low="grey", high="blue")
ggsave(dirout(out, "INS_tSNE.pdf"),width=16, height=5)

ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=GCG)) + geom_point(alpha=0.5) + 
  facet_grid(. ~ Sample) + scale_color_gradient(low="grey", high="blue")
ggsave(dirout(out, "GCG_tSNE.pdf"),width=16, height=5)


# QC PLOTS ----------------------------------------------------------------

ggplot(pDat, aes(x=INS, y=GCG, color=group)) + geom_point()
ggsave(dirout(out, "QC_Assignment_GROUPS.pdf"), height=7, width=7)

ggplot(pDat, aes(x=group, y=nUMI, fill=group)) + geom_violin() + facet_grid(. ~ Sample) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggsave(dirout(out, "QC_nUMI.pdf"),width=16, height=5)

ggplot(pDat, aes(x=group, y=nGene, fill=group)) + geom_violin() + facet_grid(. ~ Sample) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggsave(dirout(out, "QC_nGene.pdf"),width=16, height=5)

ggplot(pDat, aes(x=group, y=percent.mito, fill=group)) + geom_violin() + facet_grid(. ~ Sample) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggsave(dirout(out, "QC_PercentMito.pdf"),width=16, height=5)

ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=group)) + geom_point(alpha=0.5) + 
  facet_grid(. ~ Sample) + scale_color_manual(values=c("turquoise", "gold", "black", "grey"))
ggsave(dirout(out, "group_tSNE.pdf"),width=16, height=5)


# PERCENTAGE --------------------------------------------------------------

pDat2 <- pDat[,.N, by=c("group", "Sample")]
pDat2[, sampleCount := sum(N), by="Sample"]
pDat2[, percentage := N/sampleCount * 100]
ggplot(pDat2, aes(x=group, fill=Sample, y=percentage)) + geom_bar(stat="identity", position="dodge") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggsave(dirout(out, "DoublePositive_Percentage.pdf"), height=7, width=6)


fisher.test(as.matrix(with(pDat[Sample %in% c("Artemether", "DMSO")], table(group != "INS+ GCG+", Sample))))
fisher.test(as.matrix(with(pDat[Sample %in% c("Artemether", "GABA")], table(group != "INS+ GCG+", Sample))))




# SIGNATURES --------------------------------------------------------------

str(bSig <- fread("metadata//Ackermann2016_Bcell.txt", header=FALSE)$V1)
str(aSig <- fread("metadata//Ackermann2016_Acell.txt", header=FALSE)$V1)
sigs <- rbind(data.table(Gene=bSig, cell="beta"), data.table(Gene=aSig, cell="alpha"))
sigs[Gene == "SYNDIG", Gene := "SYNDIG1"]

sigs[!Gene %in% row.names(pbmc@data)]
"ANXA8" %in% row.names(pbmc@raw.data)
"ANX8" %in% row.names(pbmc@raw.data)
"KRTAP" %in% row.names(pbmc@raw.data)
row.names(pbmc@raw.data)[grepl("KRTAP", row.names(pbmc@raw.data))]
"SYNDIG" %in% row.names(pbmc@raw.data)

sigs <- sigs[Gene %in% row.names(pbmc@data)]
annotSample <- data.frame(Group = pDat$group, row.names=pDat$cell)
annotGene <- data.frame(Cell = sigs$cell, row.names=sigs$Gene)


treatment <- "Artemether"
for(treatment in unique(pDat$Sample)){
  cnts <- table(pDat[grepl(treatment, cell)]$group)
  pdf(dirout(out, "Signature_", treatment, ".pdf"), height=10, width=10,onefile=F)
  pheatmap(pbmc@data[sigs$Gene,pDat[grepl(treatment, cell)][order(group)]$cell], 
           cluster_rows=FALSE,cluster_cols=FALSE,
           show_colnames=FALSE,
           annotation_col=annotSample, annotation_row=annotGene,
           color=colorRampPalette(c("black", "yellow"))(30), border_color=NA,
           gaps_col=c(cnts["GCG+"], sum(cnts[c("GCG+", "INS+")]), sum(cnts[c("GCG+", "INS+", "INS+ GCG+")])))
  dev.off()
  
  pdf(dirout(out, "Signature_", treatment, "3.pdf"), height=10, width=10,onefile=F)
  pheatmap(pbmc@data[sigs$Gene,pDat[grepl(treatment, cell) & group == "INS+ GCG+"]$cell], 
           show_colnames=FALSE,
           annotation_col=annotSample, annotation_row=annotGene,
           color=colorRampPalette(c("black", "yellow"))(30), border_color=NA
           )
  dev.off()
  
  meanVals <- lapply(unique(pDat$group), function(grp){
    m <- matrix(apply(pbmc@data[sigs$Gene, pDat[grepl(treatment, cell) & group == grp]$cell], 1, mean))
    colnames(m) <- grp
    rownames(m) <- sigs$Gene
    m
  })
  meanVals <- do.call(cbind, meanVals)
  meanVals <- meanVals - apply(meanVals, 1, min, na.rm=TRUE)
  meanVals <- meanVals / apply(meanVals, 1, max, na.rm=TRUE)
  
  pdf(dirout(out, "Signature_", treatment, "2.pdf"), height=10, width=3,onefile=F)
  pheatmap(meanVals[apply(meanVals, 1, sum,na.rm=TRUE) > 0,colnames(meanVals) != "Negative"],
           annotation_row=annotGene, cluster_rows=FALSE,
           color=colorRampPalette(c("black", "yellow"))(30))
  dev.off()
}




# PREDICTION --------------------------------------------------------------
stopifnot(all(pDat$cell == colnames(pbmc@data)))
predMeta <- pDat[Sample == "DMSO" & group %in% c("GCG+", "INS+")]
dat <- t(as.matrix(pbmc@data[-which(rownames(pbmc@data) %in% c("INS", "GCG")),predMeta$cell]))


# CROSS VALIDATION ON DMSO
trainIndex <- createDataPartition(predMeta$group, p = .75, list = TRUE, times = 3)
trainRes <- data.table()
for(trainX in 1:length(trainIndex)){
  idx <- trainIndex[[trainX]]
  message(trainX)
  for(alpha in seq(0, 1, 0.2)){
    glmFit <- glmnet(y=factor(predMeta[idx]$group), 
                     x=dat[idx,],
                     family="binomial",
                     alpha=alpha
                     )
    pred <- predict(glmFit, dat[-idx,], type="response")
    for(i in 1:ncol(pred)){
      (auc <- performance(prediction(pred[,i], predMeta[-idx]$group),"auc")@y.values[[1]])
      trainRes <- rbind(trainRes, data.table(auc=auc, alpha=alpha, lambda = glmFit$lambda[i], indx=trainX))
    }
  }
}
trainRes[, maxAUC := max(auc), by="indx"]
trainRes[auc == maxAUC]
with(trainRes[auc == maxAUC], table(alpha, indx))
trainRes[indx == 2 & alpha==0]
ggplot(trainRes, aes(x=lambda, y=auc)) + geom_point() + facet_grid(factor(alpha) ~ indx) +
  scale_x_log10() + theme_bw(12)
ggsave(dirout(out, "Prediction_training.pdf"), height=14, width=14)



# TRAIN ONE FINAL MODEL
# glmFit <- glmnet(y=factor(predMeta$group), x=dat,family="binomial",alpha=0,lambda=5)
set.seed(5001)
glmFit <- glmnet(y=factor(predMeta$group), x=dat,family="binomial",alpha=1,lambda=0.1)
save(glmFit, file=dirout(out, "Prediction_LogisticRegression.RData"))

predMeta <- pDat
dat <- t(as.matrix(pbmc@data[-which(rownames(pbmc@data) %in% c("INS", "GCG")),predMeta$cell]))
predMeta$class <- predict(glmFit, dat, type="response")
ggplot(predMeta, aes(x=group, y=class, fill=group)) + geom_violin() +   theme_bw(12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  facet_grid(. ~ sample) + ylab("beta cell probability") + xlab("") + ylim(0,1)
ggsave(dirout(out, "Prediction_ClassProbabilities.pdf"),width=7, height=5)


# ERROR ON DMSO, GABA, and Artemether
treat <- "DMSO"
for(treat in unique(pDat$Sample)){
  predMeta <- pDat[Sample == treat]
  dat <- t(as.matrix(pbmc@data[-which(rownames(pbmc@data) %in% c("INS", "GCG")),predMeta$cell]))
  predMeta$class <- predict(glmFit, dat, type="response")
  (auc <- performance(prediction(
    predMeta[group %in% c("INS+", "GCG+")]$class, 
    predMeta[group %in% c("INS+", "GCG+")]$group),"auc")@y.values[[1]])
  ggplot(predMeta, aes(x=group, y=class)) + geom_violin() + ggtitle(paste("DMSO,",treat,"= ", round(auc,4)))
  ggsave(dirout(out, "Prediction_",treat,".pdf"),height=7,width=7)
}

# # Make sure I actually filtered INS GCG
# c("INS", "GCG") %in% row.names(pbmc@data)
# c("INS", "GCG") %in% colnames(dat)
# c("INS", "GCG") %in% colnames(dat2)
# head(colnames(dat))
# head(colnames(dat2))
#

# Calculate z score for up and down - is this different between double pos and double neg
(feat <- data.table(as.matrix(coef(glmFit)),keep.rownames=TRUE)[order(abs(s0),decreasing=TRUE)][s0 != 0][rn != "(Intercept)"])
stopifnot(all(colnames(pbmc@data) == pDat$cell))

pDat$Beta <- apply(pbmc@data[feat[s0 > 0]$rn,pDat$cell] * feat[s0 > 0]$s0, 2, sum)
pDat$Alpha <- apply(pbmc@data[feat[s0 < 0]$rn,pDat$cell]* feat[s0 < 0]$s0, 2, sum)
ggplot(pDat, aes(x=group, y=Beta, fill=group)) + geom_violin() + theme_bw(12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  facet_grid(. ~ sample) + ylab("beta features") + xlab("")
ggsave(dirout(out, "Features_Beta.pdf"), width=16, height=7)

ggplot(pDat, aes(x=group, y=Alpha, fill=group)) + geom_violin() + theme_bw(12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  facet_grid(. ~ sample) + ylab("alpha features") + xlab("")
ggsave(dirout(out, "Features_Alpha.pdf"), width=16, height=7)

ggplot(pDat, aes(x=Beta, y=Alpha, color=group)) + geom_point(alpha=0.5) +
  facet_grid(. ~ sample)+ theme_bw(24) +
  scale_color_manual(values=c("turquoise", "gold", "black", "grey"))
ggsave(dirout(out, "Features_points.pdf"), width=16, height=5)



# PLOT FEATURES IN HEATMAP ------------------------------------------------

annotGene2 <- data.frame(weight=feat$s0, row.names=feat$rn)
treatment <- "Artemether"
for(treatment in unique(pDat$Sample)){
  pDatX <- pDat[Sample == treatment]
  cnts <- table(pDatX$group)
  pdf(dirout(out, "Features_", treatment, ".pdf"), height=5, width=10,onefile=F)
  pheatmap(pbmc@data[feat$rn,pDatX[order(group)]$cell], 
           cluster_rows=FALSE,cluster_cols=FALSE,
           show_colnames=FALSE,
           annotation_col=annotSample, annotation_row=annotGene2,
           color=colorRampPalette(c("black", "yellow"))(30), border_color=NA,
           gaps_col=c(cnts["GCG+"], sum(cnts[c("GCG+", "INS+")]), sum(cnts[c("GCG+", "INS+", "INS+ GCG+")])))
  dev.off()
}




g <- "PAX4"
g <- "ARX"
g <- "GABARG2"
g <- "GABARB3"
g <- "GPHN"
g <- "INS"
g %in% row.names(pbmc@data)
table(pbmc@data[g,] == 0)
