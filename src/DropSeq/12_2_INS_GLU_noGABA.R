
# WITHOUT GABA SAMPLES

require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(gridExtra)
require(glmnet)
require(ROCR)
project.init2("artemether")


out <- "12_2_INS_GLC_noGABA_DeltaCells/"
dir.create(dirout(out))

load(file=dirout("DropSeq/10_Seurat/", "DropSeq/","DropSeq.RData"))


# QC nUMIs ----------------------------------------------------------------
ggplot(subset(pbmc@meta.data, nUMI > 500), aes(x=nUMI)) + geom_density() + 
  scale_x_log10() + facet_grid(. ~ sample) + theme_bw(16)
ggsave(dirout(out, "Quality_nUMIs_lt500.pdf"), width=15, height=5)


# INS GCG EXPRESSION ------------------------------------------------------
pDat <- data.table(cell = colnames(pbmc@data), INS=pbmc@data["INS",], GCG=pbmc@data["GCG",], SST=pbmc@data["SST",],
                   #NGN3=pbmc@data["NEUROG3",], 
                   Sample=pbmc@meta.data$sample, nGene = pbmc@meta.data$nGene, nUMI=pbmc@meta.data$nUMI, 
                   pbmc@meta.data,
                   pbmc@dr$tsne@cell.embeddings)
pDat <- pDat[Sample != "GABA"]
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
pDat[SST > 6,group := "Negative"]
pDat

save(pDat, file=dirout(out, "Artemether_groups.RData"))

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

pDat$INS_raw <- pmin(pbmc@raw.data["INS", pDat$cell], 100)
ggplot(pDat, aes(x=INS, y=GCG, color=INS_raw)) + geom_point() + facet_grid(. ~ Sample) + theme_bw(24)
ggsave(dirout(out, "INS_GCG_RawINS_Reads_Max100.pdf"),width=11, height=5)

ggplot(pDat, aes(x=INS, y=GCG, color=nUMI, size=nGene)) + geom_point(alpha=0.5) + facet_grid(. ~ Sample) + theme_bw(24)
ggsave(dirout(out, "INS_GCG_nGene_nUMI.pdf"),width=11, height=5)

(p <- ggplot(pDat[INS != 0 & GCG != 0], aes(x=INS, y=GCG)) + geom_hex() + facet_grid(. ~ Sample) + theme_bw(24))
ggsave(dirout(out, "INS_GCG_hex.pdf"),width=11, height=5)

(p <- ggplot(pDat, aes(x=INS, y=GCG)) + geom_point(alpha=0.5) + facet_grid(. ~ Sample) + theme_bw(16))
ggsave(dirout(out, "INS_GCG.pdf"),width=11, height=5)

(p + geom_hline(yintercept=4, color="blue", size=1, alpha=0.5) + 
   geom_vline(xintercept=4, color="blue", size=1, alpha=0.5) +
   geom_hline(yintercept=6, color="blue", size=1, alpha=0.5) + 
   geom_vline(xintercept=6, color="blue", size=1, alpha=0.5)
)
ggsave(dirout(out, "INS_GCG_lines2.pdf"),width=11, height=5)

(p + geom_hline(yintercept=cutoff, color="red", size=2, alpha=0.5) + 
   geom_vline(xintercept=cutoff, color="red", size=2, alpha=0.5))
ggsave(dirout(out, "INS_GCG_lines.pdf"),width=7, height=3.5)



# DENSITY PLOTS -----------------------------------------------------------

# INSULIn
ggplot(pDat, aes(x=INS, color=sample)) + stat_ecdf(size=1.5) + theme_bw(16) + ylab("Fraction of all cells")
ggsave(dirout(out, "INS_ECDF.pdf"),width=4.5, height=3.5)
ggplot(pDat, aes(x=INS, fill=sample)) + geom_density(alpha=0.3) + theme_bw(24) + ylab("Density")
ggsave(dirout(out, "INS_Density.pdf"),width=8, height=7)

# GLUCAGON
ggplot(pDat, aes(x=GCG, color=sample)) + stat_ecdf(size=1.5) + theme_bw(16)+ ylab("Fraction of all cells")
ggsave(dirout(out, "GCG_ECDF.pdf"),width=4.5, height=3.5)
ggplot(pDat, aes(x=GCG, fill=sample)) + geom_density(alpha=0.3) + theme_bw(24)+ ylab("Density")
ggsave(dirout(out, "GCG_Density.pdf"),width=8, height=7)


# t-SNE PLOTS -----------------------------------------------------------

ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=INS)) + geom_point(alpha=0.5) + 
  facet_grid(. ~ Sample) + scale_color_gradient(low="grey", high="blue")
ggsave(dirout(out, "INS_tSNE.pdf"),width=11, height=5)

ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=GCG)) + geom_point(alpha=0.5) + 
  facet_grid(. ~ Sample) + scale_color_gradient(low="grey", high="blue")
ggsave(dirout(out, "GCG_tSNE.pdf"),width=11, height=5)


# QC PLOTS ----------------------------------------------------------------

ggplot(pDat, aes(x=INS, y=GCG, color=group)) + geom_point()
ggsave(dirout(out, "QC_Assignment_GROUPS.pdf"), height=7, width=7)

ggplot(pDat, aes(x=group, y=nUMI, fill=group)) + geom_violin() + facet_grid(. ~ Sample) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + ylab("Number of UMIs")+ xlab("")
ggsave(dirout(out, "QC_nUMI.pdf"),width=11, height=5)

ggplot(pDat, aes(x=group, y=nGene, fill=group)) + geom_violin() + facet_grid(. ~ Sample) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + ylab("Identified genes")+ xlab("")
ggsave(dirout(out, "QC_nGene.pdf"),width=11, height=5)

ggplot(pDat, aes(x=group, y=percent.mito, fill=group)) + geom_violin() + facet_grid(. ~ Sample) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + ylab("Fraction mitochondrial reads") + xlab("")
ggsave(dirout(out, "QC_PercentMito.pdf"),width=11, height=5)

ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=group)) + geom_point(alpha=0.5) + 
  facet_grid(. ~ Sample) + scale_color_manual(values=c("turquoise", "gold", "black", "grey"))
ggsave(dirout(out, "group_tSNE.pdf"),width=11, height=5)


# PERCENTAGE --------------------------------------------------------------

pDat2 <- pDat[,.N, by=c("group", "Sample")]
pDat2[, sampleCount := sum(N), by="Sample"]
pDat2[, percentage := N/sampleCount * 100]
ggplot(pDat2, aes(x=group, fill=Sample, y=percentage)) + geom_bar(stat="identity", position="dodge") + theme_bw(16) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + ylab("Percentage of all cells") + xlab("")
ggsave(dirout(out, "DoublePositive_Percentage.pdf"), height=7, width=4)


grp <- c("INS+ GCG+", "INS+", "GCG+","Negative")
fishRes <- data.table()
for(gg in grp){
  ff <- fisher.test(as.matrix(with(pDat[Sample %in% c("Artemether", "DMSO")], table(group != gg, Sample))))
  fishRes <- rbind(fishRes, data.table(pval = ff$p.value, oddsRatio = ff$estimate, group=gg))
}
write.table(fishRes, file=dirout(out, "FisherTest.tsv"), sep="\t", quote=F, row.names=F)


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
  hmDat <- pbmc@data[sigs$Gene,pDat[grepl(treatment, cell)][order(group)]$cell]

  pdf(dirout(out, "Signature_", treatment, "_Subset.pdf"), height=7, width=7,onefile=F)
  pheatmap(hmDat[apply(hmDat != 0, 1, sum) > 20,], 
           cluster_rows=FALSE,cluster_cols=FALSE,
           show_colnames=FALSE, fontsize=16, 
           annotation_col=annotSample, annotation_row=annotGene, 
           color=colorRampPalette(c("black", "yellow"))(30), border_color=NA,
           gaps_col=c(cnts["GCG+"], sum(cnts[c("GCG+", "INS+")]), sum(cnts[c("GCG+", "INS+", "INS+ GCG+")])))
  dev.off()
  
  pdf(dirout(out, "Signature_", treatment, ".pdf"), height=7, width=5,onefile=F)
  pheatmap(hmDat, 
           cluster_rows=FALSE,cluster_cols=FALSE,
           show_colnames=FALSE,
           annotation_col=annotSample, annotation_row=annotGene,
           color=colorRampPalette(c("black", "yellow"))(30), border_color=NA,
           gaps_col=c(cnts["GCG+"], sum(cnts[c("GCG+", "INS+")]), sum(cnts[c("GCG+", "INS+", "INS+ GCG+")])))
  dev.off()
  
  pdf(dirout(out, "Signature_", treatment, "3.pdf"), height=7, width=5,onefile=F)
  pheatmap(hmDat[,pDat[grepl(treatment, cell) & group == "INS+ GCG+"]$cell], 
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
  
  pdf(dirout(out, "Signature_", treatment, "2.pdf"), height=7, width=3,onefile=F)
  pheatmap(meanVals[apply(meanVals, 1, sum,na.rm=TRUE) > 0,colnames(meanVals) != "Negative"],
           annotation_row=annotGene, cluster_rows=FALSE,fontsize=16, fontsize_row=6,cluster_cols=F,
           color=colorRampPalette(c("black", "yellow"))(30))
  dev.off()
}