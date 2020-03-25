require("project.init")
require(Seurat)
require(Matrix)
require(methods)
project.init2("artemether")

(load(dirout("10_M_01_Seurat/SeuratObject.RData")))
inDirFull <- "10_M_01_Seurat/"
inDirSub <- paste0(inDirFull, "Subcluster/")


# Cell Type assignment ----------------------------------------------------

# ENDOCRINE
meta <- fread(dirout(inDirSub, "MetaData_tSNE.tsv"))
ggplot(meta, aes(x=tSNE_1, y=tSNE_2, color=sample)) + geom_point()
ggsave(dirout(inDirSub, "tSNE_sample", '.jpg'), width=7, height=7)

data <- pbmc@data
for(gg in c("Ins2", "Gcg", "Sst", "Ppy")){
  stopifnot(gg %in% row.names(data))
  meta[[gg]] <- data[gg,meta$rn]
  ggplot(meta, aes_string(x=gg, color="sample")) + geom_density()
  ggsave(dirout(inDirSub, "Density_",gg, '.pdf'), width=7, height=5)
}

x <- list()
for(gg in c("Ins2", "Gcg", "Sst", "Ppy")){
  x[[gg]] <- meta[get(gg) > 7.5]$rn}
pdf(dirout(inDirSub, "Venn.pdf"))
gplots::venn(x)
dev.off()

meta[,nr_celltypes := sum(c(Ins2, Gcg, Sst, Ppy) > 7.5),by="rn"]
ggplot(meta, aes(x=tSNE_1, y=tSNE_2, color=nr_celltypes)) + geom_point() + scale_color_gradient(low="white", high="red")
ggsave(dirout(inDirSub, "tSNE_nr_celltypes", '.jpg'), width=7, height=7)

meta$celltype <- ""
meta[Ins2 > 7.5 & Gcg < 7.5 & Sst < 7.5,celltype := "Beta"]
meta[Ins2 < 7.5 & Gcg > 7.5 & Sst < 7.5,celltype := "Alpha"]
meta[Ins2 < 7.5 & Gcg < 7.5 & Sst > 7.5,celltype := "Delta"]
meta[Ppy > 7.5 & Ins2 < 7.5 & Gcg < 7.5 & Sst < 7.5,celltype := "Gamma"]
meta[celltype == "", celltype := "Endocrine"]
# tail(sort(meta[celltype == "Beta"]$Gcg),20)

ggplot(meta, aes(x=tSNE_1, y=tSNE_2, color=celltype)) + geom_point(alpha=0.3)
ggsave(dirout(inDirSub, "tSNE_celltypes", '.jpg'), width=7, height=7)

write.tsv(meta, dirout(inDirSub, "MetaData_AnnotatedCells.tsv"))
rm(list="meta")



# ALL CELLTYPES
metaSub <- fread(dirout(inDirSub, "MetaData_AnnotatedCells.tsv"))
metaFull <- fread(dirout(inDirFull, "MetaData_tSNE.tsv"))
metaFull <- merge(metaFull, metaSub[,c("rn", "celltype")], by="rn", all.x=T)
metaFull[is.na(celltype), celltype := res.0.5]

# Spike ins
meta.spikeIns <- fread(list.files(dirout("07_11_CleanMouse_CleanSpikeIns/"), pattern="Cont.Genes.Data.+tsv", full.names=T))
metaFull <- merge(metaFull, meta.spikeIns[,c("spikeIn", "rn"),with=F], by="rn", all.x=T)
metaFull[is.na(spikeIn), spikeIn := "NA"]
with(metaFull, table(celltype, as.character(spikeIn)))

# assign more celltypes
metaFull[ celltype == "7",celltype := "SI_Human"]
metaFull[ celltype == "8",celltype := "Acinar"]
metaFull[ celltype == "9",celltype := "Endothelial1"] # SPARC
metaFull[ celltype == "12",celltype := "Endothelial2"] # SPARC
metaFull[ celltype == "14",celltype := "SI_Mouse"]

write.tsv(metaFull, dirout(inDirFull, "MetaData_AnnotatedCells.tsv"))
