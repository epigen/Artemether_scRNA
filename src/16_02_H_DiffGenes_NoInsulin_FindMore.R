require("project.init")
require(Seurat)
require(Matrix)
require(methods)
project.init2("artemether")

out <- "16_H_02_DiffGenes/"
dir.create(dirout(out))

list.files(dirout("15_Counts/"))
meta <- loadHMeta2()
meta[,celltype := celltype2]
table(meta$celltype)
pbmc <- loadHData2()

marker.genes <- c("INS", "GCG", "SST", "PPY")


# Diff expression of marker genes -----------------------------------------

outTreat <- paste0(out, "VsDMSO_negbinom_onlyIns/")
seurat.genes.use <- marker.genes
seurat.thresh.use <- 0.1

source("src/12_DiffGenes_SCRIPT.R", echo=TRUE)

(load(dirout(outTreat, "AllRes.RData")))
# ggplot(res[!grepl("^\\d", V1)], aes(x=V3.1, y=V1.1, color=V3, size=pmin(5, -log10(V2)))) + geom_point() + facet_grid(V1 ~ V2.1) +
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
# 
# dmso <- log(mean(toTPM(pbmc@data["AC004556.1", meta[treatment == "DMSO" & replicate == "I" & celltype == "Beta"]$rn, drop=F]))+1)
# foxo <- log(mean(toTPM(pbmc@data["AC004556.1", meta[treatment == "FoxO" & replicate == "I" & celltype == "Beta"]$rn, drop=F]))+1)
# a10 <- log(mean(toTPM(pbmc@data["AC004556.1", meta[treatment == "A10" & replicate == "I" & celltype == "Beta"]$rn, drop=F]))+1)
# foxo-dmso
# a10-dmso

rm(list=c("pbmcOrig"))

outTreat <- paste0(out, "VsDMSO_negbinom_noIns/")
seurat.thresh.use <- 0.1
source("src/12_DiffGenes_SCRIPT.R", echo=TRUE)
