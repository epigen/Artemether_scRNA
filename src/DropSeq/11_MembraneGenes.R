require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)
require(gridExtra)

project.init2("artemether")


out <- "11_Find_Sort_Markers/"
dir.create(dirout(out))

mem <- fread(dirout(out, "uniprot-goa%3A%28_membrane+part+%5B0044425%5D_%29+AND+reviewed%3Ayes+AND+org--.tab"))
mem <- sapply(strsplit(mem[["Gene names"]], " "), function(v) v[1])

(alpha.markers <- fread(dirout("10_Seurat/DropSeq/Cluster_res.0.5/Markers_Cluster0_version2.tsv"))[V3 %in% mem])
(beta.markers <- fread(dirout("10_Seurat/DropSeq/Cluster_res.0.5/Markers_Cluster1_version2.tsv"))[V3 %in% mem])

(alpha.markers <- fread(dirout("10_Seurat/DropSeq/Cluster_res.0.5/Markers_Cluster0.tsv"))[V1 %in% mem & V3 > 0])
(beta.markers <- fread(dirout("10_Seurat/DropSeq/Cluster_res.0.5/Markers_Cluster1.tsv"))[V1 %in% mem & V3 > 0])


load(file=dirout("10_Seurat/", "DropSeq/","DropSeq.RData"))

pbmc@ident <- factor(pbmc@meta.data$res.0.5)
names(pbmc@ident) <- colnames(pbmc@data)
DoHeatmap(object = pbmc, genes.use = c("GCG", alpha.markers$V1,"INS", beta.markers$V1), slim.col.label = TRUE, remove.key = TRUE)
ggsave(dirout(out, "PotentialMarkers.jpg"), width=7, height=15)

pDat <- data.table(pbmc@dr$tsne@cell.embeddings, sample=pbmc@meta.data$sample)
glist <- list()
for(tr in unique(pDat$sample)){
  glist[[tr]] <- ggplot(pDat, aes(x=tSNE_1, y=tSNE_2)) + 
    geom_point(data=pDat, color="lightgrey") + 
    geom_point(data=pDat[sample == tr], color="blue", alpha=0.5, size=0.5) + ggtitle(tr)
}
(p <- grid.arrange(grobs=glist, ncol=3))
ggsave(dirout(out, "Treatment.pdf"),width=10, height=5,plot=p)


pDat <- data.table(INS=pbmc@data["INS",], GCG=pbmc@data["GCG",], Sample=pbmc@meta.data$sample, nGene = pbmc@meta.data$nGene, nUMI=pbmc@meta.data$nUMI)
ggplot(pDat, aes(x=INS, y=GCG, color=nUMI, size=nGene)) + geom_point(alpha=0.5) + facet_grid(. ~ Sample) + theme_bw(24)
ggsave(dirout(out, "INS_GCG.pdf"),width=16, height=5)


pDat$group <- "Negative"
pDat[INS > 5, group := "INS+"]
pDat[GCG > 5, group := "GCG+"]
pDat[GCG > 5 & INS > 5, group := "INS+ GCG+"]
pbmc@meta.data$Alpha_Beta <- pDat$group
max.nr.of.cores <- 2
extra.genes.to.plot <- c()
outS <- paste0("10_Seurat/", "DropSeq/")
source(paste0(Sys.getenv("CODEBASE"), "10x_datasets/src/FUNC_Seurat2.R"), echo=TRUE)
#save(pbmc, file=dirout("10_Seurat/", "DropSeq/","DropSeq.RData"))


(double.markers <- fread(dirout("10_Seurat/DropSeq/Cluster_Alpha_Beta/Markers_ClusterDoublePos.tsv"))[V1 %in% mem & V3 > 0]$V1)
pbmc@ident <- factor(pbmc@meta.data$Alpha_Beta)
names(pbmc@ident) <- colnames(pbmc@data)
DoHeatmap(object = pbmc, genes.use = c("GCG", "INS", double.markers), slim.col.label = TRUE, remove.key = TRUE)
ggsave(dirout(out, "PotentialMarkers_DoublePos1.jpg"), width=7, height=7)

DoHeatmap(object = pbmc, genes.use = c("GCG", alpha.markers$V1,"INS", beta.markers$V1), slim.col.label = TRUE, remove.key = TRUE)
ggsave(dirout(out, "PotentialMarkers_DoublePos2.jpg"), width=7, height=15)
