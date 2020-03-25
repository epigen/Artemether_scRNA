require("project.init")
require(Seurat)
require(Matrix)
require(methods)
project.init2("artemether")

out <- "EXPORT_01_Averages/"
dir.create(dirout(out))



# Read in Data ------------------------------------------------------------
metaH <- loadHMeta2()
pbmcH <- loadHData()
load(dirout("07_03_CleanHuman_CleanSpikeIns/Uncorrected.Data.RData"))
quantile(Matrix::colSums(data.raw.cells[,1:20]))
dat.uncorr.tpm <- toTPM_fromRaw(data.raw.cells, 1e6)
dat.uncorr.log <- toLogTPM(dat.uncorr.tpm)

rowMeansLog <- function(m){
  log(Matrix::rowMeans(toTPM(m)) + 1)
}


# AVERAGES ACROSS SAMPLES -------------------------------------------------
avDat.uncor <- sapply(split(metaH$rn, factor(metaH$sample)), function(rns) rowMeansLog(dat.uncorr.log[,rns]))
write.csv(avDat.uncor[row.names(pbmcH@data),], dirout(out, "Full_Sample_Averages_uncorrected.csv"))

avDat <- sapply(split(metaH$rn, factor(metaH$sample)), function(rns) rowMeansLog(pbmcH@data[,rns]))
write.csv(avDat, dirout(out, "Full_Sample_Averages.csv"))

# Markers -----------------------------------------------------------------
(load(dirout("20_01_AlphaCells/", "Markers.RData")))
markers <- markers[V4 > 0.5]
markers <- markers[V2.1 == "human"]

genes <- markers[,c("V1", "V1.1")]
colnames(genes) <- c("gene", "group")
genes <- rbind(genes, data.table(fread("metadata/Dediff_DOWN.csv"), group = "Dediff_down"))
genes <- rbind(genes, data.table(fread("metadata/Dediff_UP.csv"), group = "Dediff_up"))

# Get averages for those genes --------------------------------------------
genes <- genes[gene %in% row.names(pbmcH@data)]
for(ct in c("Alpha", "Beta")){
  for(repl in c("I", "II")){
    for(treat in c("DMSO", "A10", "FoxO")){
      genes[[paste(treat, ct, repl, sep="_")]] <- rowMeansLog(pbmcH@data[genes$gene, metaH[celltype2 == ct & replicate == repl & treatment == treat]$rn])      
    }
  }
}
write.tsv(genes, dirout(out, "AverageData.tsv"))



# HEATMAP -----------------------------------------------------------------
genes <- fread("AverageData.tsv")
genes$gene.unique <- make.unique(genes$gene)
mt <- as.matrix(genes[,3:14])
row.names(mt) <- genes$gene.unique
pdf(dirout(out, "Heatmap.pdf"), width=7, height=15, onefile=F)
pheatmap(mt, annotation_row=data.frame(row.names=genes$gene.unique, group=genes$group), cluster_rows=T, fontsize_row=6)
dev.off()



# plot for specific gene groups -------------------------------------------
genes[group %in% c("Alpha", "Beta"), c("gene.unique", "group", "DMSO_Alpha_I", "FoxO_Alpha_I", "DMSO_Alpha_II", "FoxO_Alpha_II")]
genes[group %in% c("Alpha", "Beta"), c("gene.unique", "group", "DMSO_Alpha_I", "A10_Alpha_I", "DMSO_Alpha_II", "A10_Alpha_II")]

for(ct in c("Alpha", "Beta")){
  for(repl in c("I", "II")){
    for(treat in c("A10", "FoxO", "DMSO")){
      val.treat <- genes[group %in% c("Alpha", "Beta")][[paste(treat, ct, repl, sep="_")]]
      val.dmso <- genes[group %in% c("Alpha", "Beta")][[paste("DMSO", ct, repl, sep="_")]]
      qplot(val.treat, val.dmso) +
        xlab(treat) + ylab("DMSO") + ggtitle(paste(treat, ct, repl, round(cor(val.treat, val.dmso, method='spearman'),3), sep=" "))
      ggsave(dirout(out, paste("Cor", treat, ct, repl, sep="_"), ".pdf"))
    }
  }
}