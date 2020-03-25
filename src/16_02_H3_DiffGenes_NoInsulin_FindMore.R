require("project.init")
require(Seurat)
require(Matrix)
require(methods)
project.init2("artemether")

out <- "16_H3_02_DiffGenes/"
dir.create(dirout(out))

list.files(dirout("15_Counts/"))
meta <- loadH3Meta2()
meta[,celltype := celltype2]
meta[treatment == "DMSO", replicate := gsub("III\\d", "III", replicate)]
table(meta$celltype)
pbmc <- loadH3Data2()

marker.genes <- c("INS", "GCG", "SST", "PPY")


# Diff expression of marker genes -----------------------------------------

outTreat <- paste0(out, "VsDMSO_negbinom_onlyIns/")
seurat.genes.use <- marker.genes
seurat.thresh.use <- 0


source("src/12_DiffGenes_SCRIPT_AverageDMSO.R", echo=TRUE)

ff <- list.files(dirout(outTreat), pattern="III\\d_.+_.+_.+\\.tsv$")
res <- data.table(); for(f in ff){print(f);res <- rbind(res, data.table(fread(dirout(outTreat, f)), file=f), fill=T)}
res[, file := gsub("SI_", "SI.", file)]
res[, file := gsub("Acinar_like", "Acinar.like", file)]
res[, file := gsub("\\.tsv$", "", file)]
res <- cbind(res, data.table(do.call(rbind, strsplit(res$file, "_"))))
colnames(res) <- make.unique(colnames(res))
res$qval <- p.adjust(res$p_val, method="BH")
res[, direction := ifelse(avg_diff > 0, "up", "down")]
# res[[7]] <- NULL
# res
# res[V1 == "INS" & V4 == "Alpha"]
# ExpMean(meta[replicate == "III1_72" & treatment == "A10" & celltype2 == "Alpha"]$INS) - ExpMean(meta[replicate == "III1_72" & treatment == "DMSO" & celltype2 == "Alpha"]$INS)
# ggplot(meta[replicate == "III1_72" & celltype2 == "Alpha"], aes(x=INS, color=treatment)) + stat_ecdf()
save(res, file=dirout(outTreat, "AllRes.RData"))

# ggplot(res[!grepl("^\\d", V1)], aes(x=V4, y=V1.1, color=avg_diff, size=pmin(5, -log10(qval)))) + geom_point() + facet_grid(V1 ~ V2) +
#   theme_bw(12) + scale_color_gradient2(low="blue", high="red") + xRot()
# ggsave(dirout(out, "VsDMSO_negbinom_onlyIns.pdf"), width=15, height=8)

rm(list=c("seurat.genes.use", "seurat.thresh.use", "outTreat"))


# Diff expression without marker genes -----------------------------------------
pbmcOrig <- pbmc

gg <- row.names(pbmcOrig@data)
rawDat <- pbmcOrig@raw.data[gg[!gg %in% marker.genes], pbmc@cell.names]
str(rawDat)

pbmc <- CreateSeuratObject(raw.data = rawDat, min.cells = 3, min.genes = -Inf, scale.factor=SCALE.FACTOR, project = "X")
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = SCALE.FACTOR)

# dmso <- log(mean(toTPM(pbmc@data["AC004556.1", meta[treatment == "DMSO" & replicate == "I" & celltype == "Beta"]$rn, drop=F]))+1)
# foxo <- log(mean(toTPM(pbmc@data["AC004556.1", meta[treatment == "FoxO" & replicate == "I" & celltype == "Beta"]$rn, drop=F]))+1)
# a10 <- log(mean(toTPM(pbmc@data["AC004556.1", meta[treatment == "A10" & replicate == "I" & celltype == "Beta"]$rn, drop=F]))+1)
# foxo-dmso
# a10-dmso

rm(list=c("pbmcOrig"))

outTreat <- paste0(out, "VsDMSO_negbinom_noIns/")
seurat.thresh.use <- 0.1
source("src/12_DiffGenes_SCRIPT_AverageDMSO.R", echo=TRUE)


ff <- list.files(dirout(outTreat), pattern="III\\d_.+_.+_.+\\.tsv$")
res <- data.table(); for(f in ff){print(f);res <- rbind(res, data.table(fread(dirout(outTreat, f)), file=f), fill=T)}
res[, file := gsub("SI_", "SI.", file)]
res[, file := gsub("Acinar_like", "Acinar.like", file)]
res[, file := gsub("\\.tsv$", "", file)]
res <- cbind(res, data.table(do.call(rbind, strsplit(res$file, "_"))))
colnames(res) <- make.unique(colnames(res))
res$qval <- p.adjust(res$p_val, method="BH")
res[, direction := ifelse(avg_diff > 0, "up", "down")]
save(res, file=dirout(outTreat, "AllRes.RData"))