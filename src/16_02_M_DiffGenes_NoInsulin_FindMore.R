require("project.init")
require(Seurat)
require(Matrix)
require(methods)
project.init2("artemether")

out <- "16_M_02_DiffGenes/"
dir.create(dirout(out))

list.files(dirout("15_Counts/"))
meta <- loadMMeta2()
meta[,celltype := celltype2]
table(meta$celltype)
pbmc <- loadMData2()

marker.genes <- c("Ins1", "Ins2", "Gcg", "Sst", "Ppy")



# Diff expression of marker genes -----------------------------------------

outTreat <- paste0(out, "VsDMSO_negbinom_onlyIns/")
seurat.genes.use <- marker.genes
seurat.thresh.use <- 0.1

source("src/12_DiffGenes_SCRIPT.R", echo=TRUE)

(load(dirout(outTreat, "AllRes.RData")))
# ggplot(res[!grepl("^\\d", V1)], aes(x=V3, y=V1, color=V3, size=pmin(5, -log10(V2)))) + geom_point() + facet_grid(V1 ~ V2.1) +
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

rm(list=c("pbmcOrig"))

outTreat <- paste0(out, "VsDMSO_negbinom_noIns/")
seurat.thresh.use <- 0.1
source("src/12_DiffGenes_SCRIPT.R", echo=TRUE)

# list.files(dirout(outTreat))
#
# x <- loadMMeta2()
# hits <- fread(dirout(outTreat, "3_Agg_Genes.tsv"))
# tar <- hits[V2.1=="A10" & V3.1 == "Alpha" & direction == "up"]$V1
#
# cells = x[celltype2 == "Alpha" & treatment %in% c("DMSO", "A10")]
# cells = cells[sample(1:nrow(cells), 600)][order(replicate, treatment)]$rn
# dat = pbmc@data[tar, cells]
# dat = dat / apply(dat, 1, max)
# pheatmap(dat, cluster_cols=F,
#          annotation_col=data.frame(row.names=x$rn, replicate=x$replicate, treatment=x$treatment),
#          color=colorRampPalette(c("black", "yellow"))(50))
#
# x2 = x[rn %in% cells]
# x2 <- cbind(x2, data.table(as.matrix(t(dat[,x2$rn]))))
# x2 = melt(x2, id.vars=c("replicate", "treatment"),measure.vars=tar)
# ggplot(x2, aes(x=replicate, fill=treatment,y=value)) + geom_violin() +
#   facet_grid(. ~ variable)
# ggsave(dirout(outTreat, "Diff_Expr_example.pdf"), width=15, height=5)
