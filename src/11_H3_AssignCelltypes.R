require("project.init")
require(Seurat)
require(Matrix)
require(methods)
project.init2("artemether")

(load(dirout("10_H3_01_Seurat/SeuratObject.RData")))
inDirFull <- "10_H3_01_Seurat/"
inDirSub <- paste0(inDirFull, "Subcluster/")


# Cell Type assignment ----------------------------------------------------


# IN ENDOCRINE CELLS
meta <- fread(dirout(inDirSub, "MetaData_tSNE.tsv"))
ggplot(meta, aes(x=tSNE_1, y=tSNE_2, color=sample)) + geom_point()
ggsave(dirout(inDirSub, "tSNE_sample", '.jpg'), width=7, height=7)

data <- pbmc@data
for(gg in c("INS", "GCG", "SST", "PPY", "REG1A")){
  stopifnot(gg %in% row.names(data))
  meta[[gg]] <- data[gg,meta$rn]
  ggplot(meta, aes_string(x=gg, color="sample")) + geom_density()
  ggsave(dirout(inDirSub, "Density_",gg, '.pdf'), width=7, height=5)
}

meta[,nr_celltypes := sum(c(INS, GCG, SST, PPY) > 7.5),by="rn"]
ggplot(meta, aes(x=tSNE_1, y=tSNE_2, color=nr_celltypes)) + geom_point() + scale_color_gradient(low="white", high="red")
ggsave(dirout(inDirSub, "tSNE_nr_celltypes", '.jpg'), width=7, height=7)

x <- 7.5; dpF <- data.table()
for(x in seq(0,15,0.5)){
  dpF <- rbind(dpF, data.table(meta[, sum(GCG > x & INS > x) / .N, by="sample"], cutoff=x))}
dpF <- data.table(dpF, do.call(rbind, strsplit(dpF$sample, "_")))
ggplot(dpF, aes(x=cutoff, y=V1, color=V3)) + geom_line() + facet_grid(V2 ~ .)
ggsave(dirout(inDirSub, "DoublePositives", '.pdf'), width=10, height=7)


meta$celltype <- ""
meta[INS > 7.5 & GCG < 7.5 & SST < 7.5,celltype := "Beta"]
meta[INS < 7.5 & GCG > 7.5 & SST < 7.5,celltype := "Alpha"]
meta[INS < 7.5 & GCG < 7.5 & SST > 7.5,celltype := "Delta"]
meta[PPY > 7.5 & INS < 7.5 & GCG < 7.5 & SST < 7.5,celltype := "Gamma"]
meta[REG1A > 1.25 & tSNE_1 > 25 & tSNE_2 < 4 & tSNE_2 > -10,celltype := "Acinar_like"]
meta[celltype == "", celltype := "Endocrine"]
# tail(sort(meta[celltype == "Beta"]$GCG),20)

ggplot(meta, aes(x=tSNE_1, y=tSNE_2, color=celltype)) + geom_point(alpha=0.3)
ggsave(dirout(inDirSub, "tSNE_celltypes", '.jpg'), width=7, height=7)

ggplot(meta, aes(x=tSNE_1, y=tSNE_2, color=celltype)) + geom_point(alpha=0.3) + facet_grid(.~celltype)
ggsave(dirout(inDirSub, "tSNE_celltypes_wide", '.jpg'), width=20, height=7)

write.tsv(meta, dirout(inDirSub, "MetaData_AnnotatedCells.tsv"))
rm(list="meta")




#  IN ALL CELLTYPES
metaSub <- fread(dirout(inDirSub, "MetaData_AnnotatedCells.tsv"))
metaFull <- fread(dirout(inDirFull, "MetaData_tSNE.tsv"))
metaFull <- merge(metaFull, metaSub[,c("rn", "celltype")], by="rn", all.x=T)
metaFull[is.na(celltype), celltype := res.0.5]

# Spike ins
meta.spikeIns <- fread(list.files(dirout("07_04_CleanHuman3_CleanSpikeIns/"), pattern="Cont.Genes.Data.+tsv", full.names=T))
metaFull <- merge(metaFull, meta.spikeIns[,c("spikeIn", "rn"),with=F], by="rn", all.x=T)
metaFull[is.na(spikeIn), spikeIn := "NA"]
with(metaFull, table(celltype, as.character(spikeIn)))

# assign more celltypes
metaFull[ celltype == "3",celltype := "Ductal"]
metaFull[ celltype == "6",celltype := "Acinar"]
metaFull[ celltype %in% c("8", "12", "10"),celltype := "Endothelial"] # CD3
metaFull[ celltype == "9",celltype := "SI_human"]
metaFull[spikeIn == "mouse", celltype := "SI_mouse"]

write.tsv(metaFull, dirout(inDirFull, "MetaData_AnnotatedCells.tsv"))
