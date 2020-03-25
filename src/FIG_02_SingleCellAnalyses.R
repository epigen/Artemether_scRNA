require("project.init")
project.init2("artemether")

out <- "FIG_02_SingleCellAnalyses/"
dir.create(dirout(out))

cleanTreatment <- function(vec){ gsub("^FoxO$", "FoxOi", gsub("^A10$", "Artemether", vec)) }
cleanTreatment2 <- function(vec){ gsub("FoxO", "FoxOi", gsub("A10", "Artemether", vec)) }

# LOAD DATA ---------------------------------------------------------------
plot.data.files <- dirout(out, "Data.RData")
if(!file.exists(plot.data.files)){
  metaH <- loadHMetaUnfiltered()
  pbmcH <- loadHData()
  metaM <- loadMMetaUnfiltered()()
  pbmcM <- loadMData()
  
  metaH <- assignGenesMeta(c("INS", "GCG"), metaH, pbmcH)
  metaM <- assignGenesMeta(c("Ins2", "Gcg"), metaM, pbmcM)
  metaM[,INS := Ins2]
  metaM[,GCG := Gcg]
  metaH <- metaH[!sample %in% c("hIslets_II_GABA", "hIslets_I_pla2g16")]
  dat.log <- list(human=pbmcH@data[,metaH$rn], mouse=pbmcM@data[,metaM$rn])
  dat.raw <- list(human=pbmcH@raw.data[,metaH$rn], mouse=pbmcM@raw.data[,metaM$rn])
  meta <- list(human=metaH, mouse=metaM)
  save(dat.log, dat.raw, meta, file=plot.data.files)
} else {
  load(plot.data.files)
}
meta <- lapply(meta, function(x){x$treatment <- cleanTreatment(x$treatment);x})
metaH <- meta$human
metaM <- meta$mouse

names(COLORS.TREATMENTS) <- cleanTreatment(names(COLORS.TREATMENTS))

#pancGenes <- intersect(unique(fread("metadata//PancreasMarkers.tsv")$Human_GeneName), row.names(data.raw))
celltypes.of.interest <- c("Alpha", "Beta")
use.treatments <- c("Artemether", "GABA", "FoxOi")



# EXPORT NUMBERS ----------------------------------------------------------
write.tsv(rbind(data.table(metaH[,.N, by=c("replicate", "treatment")], organism = "human"),data.table(metaM[,.N, by=c("replicate", "treatment")], organism = "mouse"))[order(organism, replicate, treatment)],
          dirout(out, "Supp_CellNumbers.tsv"))

# JACCARD OF BARCODES -----------------------------------------------------
x <- gdata::read.xls("metadata/MF153_mouse_human_Islets_10x_sample_annotation_sheet.xlsx",sheet=4)
x <- data.table(x)
x <- x[,c("Library.Name", "Sample.Name"),with=F]
x$Sample.Name <- gsub("_\\d+$", "", x$Sample.Name)
x <- unique(x[Sample.Name != ""])
x$Pool <- as.numeric(factor(x$Library.Name))
name.map <- rbind(fread("metadata/Aggregate_Human1_2.csv"), fread("metadata/Aggregate_Mouse_withLTI.csv"))
name.map[,x := basename(gsub("/outs/molecule_info.h5", "", molecule_h5))]
name.map[,x := gsub("_\\d+?_?S\\d+$", "", x)]
annot <- merge(x, name.map, by.x="Sample.Name", by.y="x")
annot <- data.frame(row.names=annot$library_id, Pool=gsub("MF153_P", "Pool", annot$Library.Name))
ll <- c(split(gsub("-\\d+", "", metaH$rn), factor(metaH$sample)), split(gsub("-\\d+", "", metaM$rn), factor(metaM$sample)))
ll <- ll[names(ll) %in% row.names(annot)]
names(ll) <- cleanTreatment2(names(ll))
row.names(annot) <- cleanTreatment2(row.names(annot))
annot <- annot[names(ll),,drop=F]
annot.col <- list(Pool=RColorBrewer::brewer.pal(n=length(unique(annot$Pool)),name="Accent"))
names(annot.col$Pool) <- unique(annot$Pool)
jac <- jaccard(ll)
jac2 <- jac/max(jac)
cleanDev()
pdf(dirout(out, "Barcode_Overlap.pdf"),width=7, height=6, onefile=F)
pheatmap(jac, clustering_distance_rows= as.dist(1-jac2), clustering_distance_cols=as.dist(1-jac2), annotation_row=annot, annotation_col=annot,
         annotation_colors=annot.col)
dev.off()

# nUMI in samples (including GABAII) --------------------------------------
metaH.old <- loadHMeta()[sample != "hIslets_I_pla2g16"]
metaH.old$treatment <- cleanTreatment(metaH.old$treatment)
ggplot(metaH.old, aes(x=treatment, y=nUMI)) + geom_boxplot(coef=Inf) +
  facet_grid(replicate ~ .) + theme_bw(16) + xRot()
ggsave(dirout(out, "nUMI_Human.pdf"), height=6, width=4)

ggplot(metaM, aes(x=treatment, y=nUMI)) + geom_boxplot(coef=Inf) +
  facet_grid(replicate ~ .) + theme_bw(16) + xRot()
ggsave(dirout(out, "nUMI_Mouse.pdf"), height=6, width=4)



# nGene and class probabilities -------------------------------------------
options(stringsAsFactors=T)
(p <- ggplot(metaH, aes(y=nGene, x=nUMI)) + stat_density_2d(aes(fill=stat(level)), geom="polygon") + 
  theme_bw(16)  + xlab("Number of UMIs") + ylab("Number of genes"))
ggsave(dirout(out, "nGene_nUMI.pdf"), h=4,w=4, plot=p+ guides(fill=F))
ggsave(dirout(out, "nGene_nUMI_Guide.pdf"), h=4,w=4, plot=p)
options(stringsAsFactors=F)

(p <- ggplot(metaH, aes(y=nGene, x=nUMI)) + geom_hex() + 
   theme_bw(16)  + xlab("Number of UMIs") + ylab("Number of genes"))
ggsave(dirout(out, "nGene_nUMI_hex.pdf"), h=4,w=4, plot=p+ guides(fill=F))
ggsave(dirout(out, "nGene_nUMI_hex_Guide.pdf"), h=4,w=4, plot=p)

class.probs <- fread(dirout("14_03_Classifier_moreCelltypes_noEndocrine/human/Prediction.tsv"), check.names=T)
class.probs$maxProb <- rowMax(as.matrix(class.probs[, colnames(class.probs)[15:26], with=F]))
class.probs$maxSecond <- apply(as.matrix(class.probs[, colnames(class.probs)[15:26], with=F]), 1, function(row) sort(row, decreasing=T)[2])
ggplot(class.probs, aes(x=maxProb)) + geom_density() + theme_bw(16) + xlab("Largest class probability") + ylab("Density")
ggsave(dirout(out, "Class_probability_Density.pdf"), h=4,w=4)
(p <- ggplot(class.probs, aes(y=nGene, x=maxProb)) + geom_hex() +
  theme_bw(16) + xlab("Largest class probability") + ylab("Number of genes"))
ggsave(dirout(out, "nGene_maxProb.pdf"), h=4,w=4, plot=p+guides(fill=F))
ggsave(dirout(out, "nGene_maxProb_guide.pdf"), h=4,w=4, plot=p)
metaH.strict <- metaH[celltype != "Endocrine" | rn %in% class.probs[maxProb/2 > maxSecond]$rn]
write.tsv(metaH.strict, dirout(out, "Reassigned_Human_Strict.tsv"))


# T-SNE Plots of celltypes -------------------------------------------------------------

# Celltypes
cl.x <- "celltype"
pDat <- metaH
clPlot <- function(pDat, cl.x, label){
  labelCoordinates <- pDat[,.(tSNE_1=median(tSNE_1), tSNE_2=median(tSNE_2)),by=cl.x]
  plim <- ceiling(max(max(abs(pDat$tSNE_1)), max(abs(pDat$tSNE_2)))/10)*10
  (p <- ggplot(pDat, aes_string(x="tSNE_1",y="tSNE_2", color=cl.x)) + theme_bw(16) + xlim(-plim, plim) + ylim(-plim,plim) +
     xlab("t-SNE dimension 1") + ylab("t-SNE dimension 2"))
  ggsave(dirout(out, label, ".jpg"), height=6, width=6, plot=p + geom_point(alpha=0.3) + guides(color=F))
  ggsave(dirout(out, label, ".pdf"), height=6, width=6, 
         plot=p + geom_point(data=pDat[rn %in% sapply(split(pDat$rn, factor(pDat[[cl.x]])), function(x) x[1])], alpha=0.3))
  ggsave(dirout(out, label, "Labels.pdf"), height=6, width=6, 
         plot=p + geom_label(data=labelCoordinates, aes_string(x="tSNE_1", y="tSNE_2", label=cl.x), color="black", alpha=0.5) + theme(panel.grid=element_blank()))
  ggsave(dirout(out, label, "Done.jpg"), height=6, width=6, 
         plot=p + geom_point(alpha=0.3) + guides(color=F) + geom_label(data=labelCoordinates, aes_string(x="tSNE_1", y="tSNE_2", label=cl.x), color="black", alpha=0.5))
}
clPlot(metaH, "celltype", "tSNE_human")
clPlot(metaM, "celltype", "tSNE_mouse")
metaH.E <- fread(dirout("10_H_01_Seurat/Subcluster/MetaData_AnnotatedCells.tsv"))[!sample %in% c("hIslets_II_GABA", "hIslets_I_pla2g16")]
metaM.E <- fread(dirout("10_M_01_Seurat/Subcluster/MetaData_AnnotatedCells.tsv"))
clPlot(metaH.E, "celltype", "tSNE_human_endocrine")
clPlot(metaM.E, "celltype", "tSNE_mouse_endocrine")

# MARKERS
marker.plot <- function(pDat, genes, dat, label, ncol=3){
  plots <- list()
  gg <- "INS"
  for(gg in genes){
    pDat$E <- dat[gg,pDat$rn]
    pDat[, E := E/max(pDat$E)]
    plots[[gg]]  <- ggplot(pDat, aes_string(x="tSNE_1",y="tSNE_2", color="E")) + theme_bw(16) + geom_point(alpha=0.3) +
      scale_color_gradient(low="lightgrey", high="blue") + ggtitle(gg) +
        guides(color=F) + theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank())
  }
  p <- gridExtra::grid.arrange(grobs=plots, ncol=ncol)
  ggsave(dirout(out, label, "_Markers.jpg"), width=ncol*3, height=ceiling(length(genes)/ncol)*3, plot=p)
}
"Reg1a" %in% row.names(dat.log$mouse)
marker.plot(metaH, c("INS", "GCG", "REG1A", "KRT8", "SPARC", "SPP1"), dat.log$human, "tSNE_human")
marker.plot(metaH.E, c("INS", "GCG", "SST", "PPY"), dat.log$human, "tSNE_human_endocrine_small", ncol=2)
marker.plot(metaH.E, c("INS", "GCG", "SST", "PPY", "TTR", "REG1A"), dat.log$human, "tSNE_human_endocrine")
marker.plot(metaM, c("Ins2", "Gcg", "Sparc", "Cd14", "Ptprc", "Ttr"), dat.log$mouse, "tSNE_mouse")
marker.plot(metaM.E, c("Ins2", "Gcg", "Sst", "Ppy", "Ttr", 'Iapp'), dat.log$mouse, "tSNE_mouse_endocrine")

# MARKER DISTRIBUTIONS USED
pDat <- data.table()
for(gg in c("INS", "GCG", "SST", "PPY", "REG1A")){
  pDat <- rbind(pDat, data.table(meta$human, Expression = dat.log$human[gg, meta$human$rn], gene=gg))}
ggplot(pDat, aes(x=Expression, color=replicate, linetype=treatment)) + geom_density() + facet_wrap(~gene, ncol=1) + theme_bw(16) + xlab("") +
  geom_vline(data=data.table(gene=c("INS","GCG","SST","PPY","REG1A"), value=c(rep(7.5,4), 1.5)), aes(xintercept=value), size=2, color="black")
ggsave(dirout(out, "Marker_cutoffs.human.pdf"), height=15, width=7)
pDat <- data.table()
for(gg in c("Ins2", "Gcg", "Sst", "Ppy")){
  pDat <- rbind(pDat, data.table(meta$mouse, Expression = dat.log$mouse[gg, meta$mouse$rn], gene=gg))}
ggplot(pDat, aes(x=Expression, color=replicate, linetype=treatment)) + geom_density() + facet_wrap(~gene, ncol=1) + theme_bw(16) + xlab("") +
  geom_vline(data=data.table(gene=c("Ins2","Gcg","Sst","Ppy"), value=c(rep(7.5,4))), aes(xintercept=value), size=2, color="black")
ggsave(dirout(out, "Marker_cutoffs.mouse.pdf"), height=15, width=10)
x <- meta$human; x[is.na(spikeIn), spikeIn := "NA"]; ggplot(x, aes(x=factor(res.0.5), fill=spikeIn)) + geom_bar() + theme_bw(16) + xRot() + ylab("Cells") + xlab("Cluster")
ggsave(dirout(out, "Marker_Spikeins_human.pdf"), height=4, width=6)
x <- meta$mouse; x[is.na(spikeIn), spikeIn := "NA"]; ggplot(x, aes(x=factor(res.0.5), fill=spikeIn)) + geom_bar() + theme_bw(16) + xRot() + ylab("Cells") + xlab("Cluster")
ggsave(dirout(out, "Marker_Spikeins_mouse.pdf"), height=4, width=6)


# Celltype TSNEs - individual ----------------------------------------------
xnam <- "30_04_ClusteringExperiments"
for(xnam in c("30_04_ClusteringExperiments", "30_05_ClusteringExperiments")){
  ff <- list.files(dirout(xnam), pattern="MetaData_tSNE.tsv", recursive=T)
  names(ff) <- make.names(dirname(ff))
  pDat <- lapply(names(ff), function(fnam){x <- fread(dirout(xnam, "/", ff[[fnam]])); x$sample <- fnam; return(x)})
  pDat <- do.call(rbind, pDat)
  pDat <- merge(pDat, loadHMetaUnfiltered()()[!sample %in% c("hIslets_I_pla2g16")][,c("rn", "PredictedCell"), with=F], by="rn")
  pDat <- cbind(pDat, do.call(rbind, strsplit(pDat$sample, "_")))
  pDat <- cbind(pDat, do.call(rbind, strsplit(pDat$V1, "\\.")))
  colnames(pDat) <- make.unique(colnames(pDat))
  pDat[,treatment := cleanTreatment(V2.1)]
  p <- ggplot(pDat[PredictedCell %in% c("Alpha", "Beta", "Gamma", "Delta")][V1 != "Aold"], 
         aes(x=tSNE_1, y=tSNE_2, color=PredictedCell)) + 
    scale_color_manual(values=c(Alpha="#6a3d9a", Beta="#33a02c", Gamma="#ff7f00", Delta="#e31a1c")) +
    geom_point(alpha=0.5) + facet_grid(V1.1 + V2 ~ treatment) +
    xlab("t-SNE dimension 1") + ylab("t-SNE dimension 2") + theme_bw(16)
  ggsave(dirout(out, "IndividualTSNEs_guides_",xnam,".pdf"), w=16,h=16, plot=p)
  ggsave(dirout(out, "IndividualTSNEs_",xnam,".jpg"), w=16,h=16, plot=p+guides(color=F))
}

# xx.list <- list(corrected = "30_01_ClusteringExperiments", raw="30_02_ClusteringExperiments_raw", raw2="30_03_ClusteringExperiments_raw_subcluster")
# xx.nam <- "raw2"
# for(xx.nam in c("corrected", "raw", "raw2")){
#   ff <- list.files(dirout(xx.list[[xx.nam]], "/"), pattern="MetaData_tSNE.tsv", recursive=T)
#   names(ff) <- dirname(ff)
#   pDat <- lapply(names(ff), function(fnam){x <- fread(dirout(xx.list[[xx.nam]], "/", ff[[fnam]])); x$sample <- fnam; return(x)})
#   pDat <- do.call(rbind, pDat)
#   pDat <- merge(pDat, loadHMetaUnfiltered()()[!sample %in% c("hIslets_I_pla2g16")][,c("rn", "PredictedCell"), with=F], by="rn")
#   pDat <- cbind(pDat, do.call(rbind, strsplit(pDat$sample, "_")))
#   pDat[,treatment := cleanTreatment(V1)]
#   p <- ggplot(pDat[PredictedCell %in% c("Alpha", "Beta", "Gamma", "Delta")][V1 != "Aold"], 
#          aes(x=tSNE_1, y=tSNE_2, color=PredictedCell)) + 
#     scale_color_manual(values=c(Alpha="#6a3d9a", Beta="#33a02c", Gamma="#ff7f00", Delta="#e31a1c")) +
#     geom_point(alpha=0.5) + facet_grid(V2 ~ treatment) +
#     xlab("t-SNE dimension 1") + ylab("t-SNE dimension 2") + theme_bw(16)
#   ggsave(dirout(out, "IndividualTSNEs_guides_",xx.nam,".pdf"), w=16,h=8, plot=p)
#   ggsave(dirout(out, "IndividualTSNEs_",xx.nam,".jpg"), w=16,h=8, plot=p+guides(color=F))
# }


# CELLTYPE PREDICTION -----------------------------------------------------
# PLOT PREDICTIONS
for(orgx in c("human", "mouse")){
  predH <- fread(dirout("14_03_Classifier_moreCelltypes_noEndocrine/",orgx,"/", "Prediction.tsv"))
  predH <- predH[!sample %in% c("hIslets_II_GABA", "hIslets_I_pla2g16")]
  predH$treatment <- cleanTreatment(predH$treatment)
  classes <- intersect(colnames(predH), predH$celltype)
  meanDat <- melt(predH, id.vars=c("treatment", "replicate", "celltype"), measure.vars=classes)[,.(value = mean(value)), by=c("treatment", "replicate", "celltype","variable")]
  ggplot(meanDat, aes(x=treatment, y=variable, fill=value)) + geom_tile() + 
    theme_bw(12) + xlab("") + ylab("") + 
    facet_grid(replicate ~ celltype) + xRot() +
    scale_fill_gradient(name="Mean Prob", low="white", high="red",limits=c(0,1))
  ggsave(dirout(out, "Predictions_",orgx,".pdf"), width=15, height=5)
  
  # PLOT ACCURACY
  predH$PredictedCell <- apply(as.matrix(predH[,unique(intersect(colnames(predH), predH$celltype)), with=F]),1,function(row) names(row)[which(row == max(row))])
  predH.acc <- predH[celltype %in% classes, sum(celltype == PredictedCell)/.N, by=c("celltype", "replicate", "treatment")]
  predH.acc$Class <- predH.acc$celltype
  ggplot(predH.acc,aes(x=treatment, y=V1,color=Class, shape=Class)) + 
    geom_jitter(width=0.1, height=0) + ylim(0,1) +
    scale_shape_manual(values=rep(c(1,16,2,18,3,4), 20)) + 
    facet_wrap(. ~ replicate, scales="free_x") +theme_bw(16) + xRot() +
    xlab("Treatment") + ylab("Precision")
  ggsave(dirout(out, "Predictions_Precision_",orgx,".pdf"), height=4, width=3*length(unique(predH$replicate)))
}

# Do the initially unassignable cells have more genes (doublets?)
pDat <- meta$human[celltype2 %in% c("Alpha", "Beta", "Gamma", "Delta")]
#pDat$celltype <- factor(pDat$celltype, levels=c("Endocrine", "Alpha", "Beta", "Gamma", "Delta"))
ggplot(pDat, aes(x=(celltype=="Endocrine"), y=nGene)) + geom_boxplot(coef=Inf) + theme_bw(16) + 
  facet_grid(treatment ~ celltype2) + xRot() +
  xlab("Reassigned cells") + ylab("Number of Genes")
ggsave(dirout(out, "Predictions_nGene.pdf"),w=8,h=8)


# COUNTS - Alpha/Beta cell ratio ------------------------------------------
# Human Counts
for(orgx in c("human", "mouse")){
  pDat <- meta[[orgx]][celltype2 != "",.N, by=c("celltype2", "replicate", "treatment", "celltype")]
  pDat[,sum := sum(N), by=c("replicate", "treatment")]
  pDat[,percent := N / sum * 100]
  pDat[celltype == "Endocrine", X := "Predicted"]
  pDat[celltype != "Endocrine", X := "Assigned"]
  ggplot(pDat, aes(x=treatment,y=percent, fill=X)) + geom_bar(stat="identity", position="stack") + facet_grid(replicate ~ celltype2)  +
   scale_fill_manual(name="", values=c("grey", "blue")) + theme_bw(16)+ xRot() + xlab("") + ylab("Percent")
  ggsave(dirout(out, "Predicted_",orgx,".pdf"), height=length(unique(pDat$replicate))* 2,width=18)
}

# Endocrine of all / beta cells
pDat <- rbind(
  data.table(meta$human[, length(rn[celltype=="Endocrine"])/.N*100, by=c("sample", "replicate", "treatment")][order(sample)], type="All", org="human"),
  data.table(meta$mouse[, length(rn[celltype=="Endocrine"])/.N*100, by=c("sample", "replicate", "treatment")][order(sample)], type="All", org="mouse"),
  data.table(meta$human[, length(rn[celltype=="Endocrine" & celltype2=="Beta"])/.N*100, by=c("sample", "replicate", "treatment")][order(sample)], type="Beta", org="human"),
  data.table(meta$mouse[, length(rn[celltype=="Endocrine" & celltype2=="Beta"])/.N*100, by=c("sample", "replicate", "treatment")][order(sample)], type="Beta", org="mouse"),
  data.table(meta$human[, length(rn[celltype=="Endocrine" & celltype2=="Alpha"])/.N*100, by=c("sample", "replicate", "treatment")][order(sample)], type="Alpha", org="human"),
  data.table(meta$mouse[, length(rn[celltype=="Endocrine" & celltype2=="Alpha"])/.N*100, by=c("sample", "replicate", "treatment")][order(sample)], type="Alpha", org="mouse"))
ggplot(pDat, aes(x=treatment, y=V1, fill=treatment)) + geom_bar(stat="identity") + facet_grid(type ~ paste(org, replicate), space="free_x", scales="free") + theme_bw(16) +xRot() +
  xlab("Treatment") + ylab("Percent of cells") + scale_fill_manual(values=COLORS.TREATMENTS)
ggsave(dirout(out, "Endocrine_Percent.pdf"), height=6,width=8)



# Alpha vs Beta
pDat <- rbind(data.table(meta$mouse[,sum(celltype2 == "Alpha")/ sum(celltype2 == "Beta"), by=c("treatment", "replicate")], org = "mouse"),
  data.table(meta$human[,sum(celltype2 == "Alpha")/ sum(celltype2 == "Beta"), by=c("treatment", "replicate")][order(replicate)], org = "human"))
pDat[,grp := paste(org, replicate)]; pDat$grp <- factor(pDat$grp, levels=unique(pDat$grp))
ggplot(pDat, aes(x=treatment, y=V1, fill=treatment)) + geom_bar(stat="identity") + 
  facet_wrap(~grp,scales="free",ncol=3) + 
  theme_bw(16) + xRot() + ylab("Alpha / beta cell ratio") + xlab("") +
  scale_fill_manual(values=COLORS.TREATMENTS) + guides(fill=F)
ggsave(dirout(out, "AlphaBetaRatio.pdf"), height=6,width=5)


# INS / GCG distributions -------------------------------------------------
for(org in c("human", "mouse")){
  for(ct in celltypes.of.interest){
    ttt <- intersect(use.treatments, meta[[org]]$treatment)
    pDat <- do.call(rbind, lapply(ttt, function(treat) data.table(meta[[org]][treatment %in% c("DMSO", treat)], treatx=treat)))
    gg <- if(ct == "Beta") "INS" else "GCG"
    for(gg in c("GCG", "INS")){
      pDat <- pDat[celltype2 == ct]
      ggplot(pDat, aes_string(x=gg, color="treatment")) + stat_ecdf() +
        theme_bw(16) + facet_grid(replicate ~ treatx, scale="free") +
        scale_color_manual(values=COLORS.TREATMENTS)+ guides(color=F) + 
        ylab("Fraction of cells")
      ggsave(dirout(out, "ECDF_", org, "_", ct,"_", gg, ".pdf"),height=2*length(unique(pDat$replicate)),width=length(unique(pDat$treatx)) * 2 + 1)
      }}}


# QUANTILE DIFFERENCES ----------------------------------------------------
ct <- celltypes.of.interest[1]
org <- "human"
for(ct in celltypes.of.interest){
  #for(gg in c("GCG", "INS")){
    res <- data.table()
    for(org in c("human", "mouse")){
      for(replx in unique(meta[[org]]$replicate)){
        for(treatX in use.treatments){
            gg <- if(ct == "Beta") "INS" else "GCG"
            pDat <- meta[[org]][celltype2 == ct & replicate == replx]
            q1 <- quantile(pDat[treatment == treatX & replicate == replx][[gg]], probs=seq(0,1,1e-2))
            q2 <- quantile(pDat[treatment == "DMSO" & replicate == replx][[gg]], probs=seq(0,1,1e-2))
            res <- rbind(res, data.table(diff= q1 - q2,treatment = treatX, replicate=replx, org=org, gene=gg))
          }
        }
      }
      ggplot(res, aes(y=diff, color=replicate, x=org)) + 
        geom_boxplot(coef=Inf) + theme_bw(16) +
        #scale_fill_manual(values=c("white", "grey")) + 
        facet_grid(. ~ treatment, scales="free", space="free") + xRot() +
        xlab("Treatment") + ylab(paste(gg, " difference to DMSO")) + ggtitle(ct)
      ggsave(dirout(out, "ECDF_Quantiles_", ct, "_", gg, ".pdf"),height=6, width=6)
  #}
}


# DIFF GENES SIMILARITIES -------------------------------------------------
file.copy(dirout("17_02_06_Gene_plots_q0.1_ByGroup/HM_Similarity_human.pdf"), dirout(out))
file.copy(dirout("17_02_06_Gene_plots_q0.1_ByGroup/HM_Similarity_mouse.pdf"), dirout(out))

