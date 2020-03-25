require("project.init")
require(SoupX)
project.init2("artemether")

out <- "FIG_01_Contamination/"
dir.create(dirout(out))


# LOAD EVERYTHING ----------------------------------------------------------
cleanTreatment <- function(vec){ gsub("FoxO", "FoxOi", gsub("A10", "Artemether", vec)) }

# LOAD RAW DATA
data.raw <- Read10X(paste0(getOption("PROCESSED.PROJECT"), "results_pipeline/cellranger_count/","Human1_2","/outs/raw_gene_bc_matrices_mex/GRCh38/"))

# LOAD OTHER STUFF
list.files(dirout("07_03_CleanHuman_CleanSpikeIns"), pattern="RData$")
(load(dirout("07_03_CleanHuman_CleanSpikeIns/","Uncorrected.Data.RData")))
#"ml.rho.pred"         "ml.rho.cv"           "ml.rho.pred.with.cv" "data.raw.cells"      "meta.raw.cells"      "spikeIn.soups" 
(load(dirout("07_03_CleanHuman_CleanSpikeIns/","CorrectedData.RData")))
# "corr.export"
(load(dirout("07_03_CleanHuman_CleanSpikeIns/","Soup.results.RData")))
(load(dirout("07_03_CleanHuman_CleanSpikeIns/", "Meta.raw.RData")))

metaH <- loadHMetaUnfiltered()

# CLEAN ANNOTATION --------------------------------------------------------
meta.raw.cells$library_id <- cleanTreatment(meta.raw.cells$library_id)
colnames(spikeIn.soups) <- cleanTreatment(colnames(spikeIn.soups))
names(soupX.contFrac) <- cleanTreatment(names(soupX.contFrac))
names(soupX.rho) <- cleanTreatment(names(soupX.rho))
names(soupX.soup) <- cleanTreatment(names(soupX.soup))
colnames(soupX.soup.tpm) <- cleanTreatment(colnames(soupX.soup.tpm))
names(soupX.strainedCells) <- cleanTreatment(names(soupX.strainedCells))
meta.raw$library_id <- cleanTreatment(meta.raw$library_id)

# remove samples not used
table(meta.raw.cells$library_id)
rm.samples <- c("hIslets_II_Aold", "hIslets_I_pla2g16")
metaH <- metaH[!sample %in% rm.samples]
meta.raw <- meta.raw[!library_id %in% rm.samples]
meta.raw.cells <- meta.raw.cells[!library_id %in% rm.samples]
meta.raw.spikes <- meta.raw[spikeIn != "NA"]
meta.raw.spikes <- meta.raw.spikes[!library_id %in% rm.samples]
for(samplex in rm.samples){
  soupX.soup[[samplex]] <- NULL
  soupX.strainedCells[[samplex]] <- NULL
}

# limit raw data to cells used
data.raw <- data.raw[,meta.raw$rn]



# ASSIGN DATA -------------------------------------------------------------

# Cell data
data.tpm.cells <- t(t(data.raw.cells) / Matrix::colSums(data.raw.cells)) * SCALE.FACTOR
data.log.cells <- toLogTPM(data.tpm.cells)

# Spike in data
data.raw.spike <- data.raw[,meta.raw.spikes$rn]
data.tpm.spike <- t(t(data.raw.spike) / Matrix::colSums(data.raw.spike)) * SCALE.FACTOR
spikeIn.soups.tpm <- t(t(spikeIn.soups)/colSums(spikeIn.soups)) * SCALE.FACTOR

# external spike ins
(load(dirout("05_SpikeIns_Clean/SpikeData_human.RData")))
spike.clean <- add.genes.sparseMatrix(toTPM(dat.log.spike), row.names(data.raw))
spike.clean.meta <- fread(dirout("05_SpikeIns_Clean/SpikeMeta_human.tsv"))
spike.clean.diffOrg <- Matrix::rowMeans(spike.clean[,spike.clean.meta[org == "mouse"]$barcode])
spike.clean.sameOrg <- Matrix::rowMeans(spike.clean[,spike.clean.meta[org == "human"]$barcode])

# Pancreatic genes
pancGenes <- intersect(unique(fread("metadata//PancreasMarkers.tsv")$Human_GeneName), row.names(data.raw))

# Samples
samples <- unique(meta.raw$library_id)

# Contaminating genes
(cont.genes <- intersect(names(tail(sort(rowMeans(spikeIn.soups)), 100)), pancGenes))

# Corrected values tpm (regenerated from exported count values)
nUMIs <- Matrix::colSums(data.raw.cells[,colnames(corr.export)])
ml.corr.values.all.merge <- t(t(corr.export/nUMIs) * SCALE.FACTOR)


# CORRECT DATA FROM SPIKE-INS in mouse and human spike-ins ---------------------------------------------
ml.corr.values <- lapply(samples, function(x){
  rn <- meta.raw.spikes[library_id == x]$rn
  ret <- data.tpm.spike[,rn] - (spikeIn.soups[row.names(data.tpm.spike),paste(x, "spikeIn")] %o% ml.rho.pred.with.cv[rn])
  ret[ret < 0] <- 0
  t(t(ret) / Matrix::colSums(ret)) * SCALE.FACTOR
})
ml.corr.values.merge <- do.call(cbind, ml.corr.values)


# COMPARE ESTIMATED TO SPIKE IN TRUTH -------------------------------------------
(sum.I <- sum(data.raw.cells[,meta.raw.cells[grepl("_I_", library_id)]$rn]))
(ins.I <- Matrix::rowSums(data.raw.cells[c("INS", "GCG"),meta.raw.cells[grepl("_I_", library_id)]$rn]))
ins.I / sum.I

(sum.II <- sum(data.raw.cells[,meta.raw.cells[grepl("_II_", library_id)]$rn]))
(ins.II <- Matrix::rowSums(data.raw.cells[c("INS", "GCG"),meta.raw.cells[grepl("_II_", library_id)]$rn]))
ins.II / sum.II

write.tsv(meta.raw.spikes, file=dirout(out, "Meta.raw.spikes.tsv"))
median(meta.raw.spikes[grepl("_I_", library_id)]$contamination, na.rm=T)
median(meta.raw.spikes[grepl("_II_", library_id)]$contamination, na.rm=T)
max(meta.raw.spikes$contamination, na.rm=T)


# COMPARE ESTIMATED RHO

# Soup
meta.raw.spikes$contSoupX <- soupX.rho.agg[meta.raw.spikes$rn]
(p <- ggplot(meta.raw.spikes[!is.na(contamination)], aes(x=contSoupX, y=contamination, color=library_id)) +  geom_point(alpha=0.5) +
  theme_bw(12) + xlim(0,.7) + ylim(0,.7) + geom_abline(size=2, color="lightgrey", alpha=0.5) + 
  xlab("Contamination soupX") + ylab("Contamination cross-alignment")) +
  ggtitle(with(meta.raw.spikes[!is.na(contamination)], cor(contSoupX, contamination), method="pearson"))
ggsave(dirout(out, "Rho_estimation_SOUP.pdf"), height=4, width=4)
ggsave(dirout(out, "Rho_estimation_SOUP_NoGuid.pdf"), height=4, width=4, plot=p + guides(color=F))

# machine learning
meta.raw.spikes$contMLCV <- ml.rho.cv[meta.raw.spikes$rn]
with(meta.raw.spikes, table(is.na(contMLCV), spikeIn))
(p <- ggplot(meta.raw.spikes[!is.na(contamination)], aes(x=contMLCV, y=contamination, color=library_id)) +  geom_point(alpha=0.5) +
  theme_bw(12) + xlim(0,.4) + ylim(0,.4) + geom_abline(size=2, color="lightgrey", alpha=0.5) +
  xlab("Contamination spike-ins") + ylab("Contamination cross-alignment") +
  ggtitle(with(meta.raw.spikes[!is.na(contamination)], cor(contMLCV, contamination), method="pearson")))
ggsave(dirout(out, "Rho_estimation_ML.pdf"), height=4, width=4)
ggsave(dirout(out, "Rho_estimation_ML_noGuide.pdf"), height=4, width=4, plot=p + guides(color=F))

write.tsv(meta.raw.spikes, dirout(out, "Rho_estimnation_Numbers.tsv"))


# COMPARE ESTIMATED CONTAMINATION (SOUP)
stopifnot(all(round(colSums(spikeIn.soups.tpm)) == 1e6))
stopifnot(all(colSums(soupX.soup.tpm) == 1e6))
for(x in unique(meta.raw.spikes$library_id)){
  qplot(spikeIn.soups.tpm[,paste0(x, " spikeIn")], soupX.soup.tpm[row.names(spikeIn.soups),x]) + 
    geom_abline(size=2, color="grey", alpha=0.5) + scale_x_log10() + scale_y_log10() + theme_bw(12) +
    xlab("Signature spike-ins") + ylab("Signature soupX") +
    annotate("text", x=10, y=10000, label=paste("R =", round(corS(spikeIn.soups.tpm[,paste0(x, " spikeIn")], soupX.soup.tpm[row.names(spikeIn.soups),x]),3)))
  ggsave(dirout(out, "Soup_estimation_",x,".jpg"), height=4, width=4)
}


# COMPARE CORRECTED VALUES
compare.values <- fread(dirout("07_03_CleanHuman_CleanSpikeIns/", "Values_compare.tsv"))
compare.values$V3 <- cleanTreatment(compare.values$V3)
compare.values <- compare.values[!V3 %in% c("pla2g16", "Aold")]
ggplot(compare.values, aes(x=V3, y=value, color=V4)) + geom_boxplot() + facet_grid(V5 ~ V2, space="free_x", scale="free_x") +
  theme_bw(12) + ylab("Correlation to reference") + xlab("Replicate") + xRot()
ggsave(dirout(out, "Values_compare_human.pdf"), height=4, width=8)

# P values
compare.values[,id2 := paste(V2, V3, V5)]
compare.values.p <- sapply(unique(compare.values$id2), function(idx) 
  t.test(compare.values[id2 == idx][V4 == "spikeIns"]$value, compare.values[id2 == idx][V4 == "original"]$value, alternative="greater")$p.value)
compare.values.p.soupX <- sapply(unique(compare.values$id2), function(idx) 
  t.test(compare.values[id2 == idx][V4 == "soupX"]$value, compare.values[id2 == idx][V4 == "original"]$value, alternative="greater")$p.value)
write.tsv(
  rbind(
    data.table(names(compare.values.p), p=compare.values.p, type="SpikeIns"), 
    data.table(names(compare.values.p.soupX), p=compare.values.p.soupX, type="SoupX")
  ),dirout(out, "Values_compare_human_p.values.tsv"))

ggplot(compare.values[V4 != "soupX"], aes(x=V3, y=value, color=V4)) + geom_boxplot() + 
  facet_grid(V5 ~ V2, space="free_x", scale="free_x") +
  theme_bw(12) + guides(color=F) + ylab("Correlation to reference") + xlab("Replicate") +
  scale_color_manual(values=c("black", "red")) + xRot()
ggsave(dirout(out, "Values_compare_NoSoupX.pdf"), height=4, width=5)

# MOUSE
compare.values <- fread(dirout("07_11_CleanMouse_CleanSpikeIns/", "Values_compare.tsv"))
compare.values$V3 <- cleanTreatment(compare.values$V3)
compare.values <- compare.values[!V3 %in% c("A1")]
table(compare.values$V3)
ggplot(compare.values[V4 != "soupX"], aes(x=V3, y=value, color=V4)) + geom_boxplot() + 
  facet_grid(V5 ~ V2, space="free_x", scale="free_x") +
  theme_bw(12) + guides(color=F) + ylab("Correlation to reference") + xlab("Replicate") +
  scale_color_manual(values=c("black", "red")) + xRot()
ggsave(dirout(out, "Values_compare_NoSoupX_mouse.pdf"), height=4, width=5)

# P values
compare.values[,id2 := paste(V2, V3, V5)]
compare.values.p <- sapply(unique(compare.values$id2), function(idx) 
  t.test(compare.values[id2 == idx][V4 == "spikeIns"]$value, compare.values[id2 == idx][V4 == "original"]$value, alternative="greater")$p.value)
compare.values.p.soupX <- sapply(unique(compare.values$id2), function(idx) 
  t.test(compare.values[id2 == idx][V4 == "soupX"]$value, compare.values[id2 == idx][V4 == "original"]$value, alternative="greater")$p.value)
write.tsv(
  rbind(
    data.table(names(compare.values.p), p=compare.values.p, type="SpikeIns"), 
    data.table(names(compare.values.p.soupX), p=compare.values.p.soupX, type="SoupX")
  ),dirout(out, "Values_compare_mouse_p.values.tsv"))

# Plots for individual comparisons
outValCompare <- paste0(out, "Values_compare/")
dir.create(dirout(outValCompare))
for(x in samples){
  for(orgXX in c("human", "mouse")){
    spikeRef <- spike.clean.sameOrg; if(orgXX != "human") spikeRef <- spike.clean.diffOrg
    rn <- meta.raw.spikes[library_id == x & spikeIn == orgXX]$rn
    
    pDat <- data.table(
      Reference = spikeRef[row.names(data.tpm.spike)],
      Corrected = Matrix::rowMeans(ml.corr.values.merge[row.names(data.tpm.spike), rn]),
      Original = Matrix::rowMeans(data.tpm.spike[, rn]))
    ggplot(pDat, aes(x=Reference, y=Original)) + geom_abline(size=2, color="lightgrey", alpha=0.5) + 
      geom_point(alpha=0.3) + geom_point(aes(y=Corrected), color="red", alpha=0.3) + theme_bw(12) + 
      ylab("Expression in spike-ins") +
      annotate("text", x=2000, y=19000, label=paste("R =", round(cor(pDat$Reference, pDat$Original),3)), color="black") +
      annotate("text", x=2000, y=21000, label=paste("R =", round(cor(pDat$Reference, pDat$Corrected),3)), color="red")
    ggsave(dirout(outValCompare, "0_", orgXX, "_", x, ".jpg"), height=4, width=4)
    ggsave(dirout(outValCompare, "0_", orgXX, "_", x, ".pdf"), height=4, width=4)
  }}


gg <- "GCG"
for(gg in cont.genes){
  m <- meta.raw.spikes
  pDat <- rbind (
    data.table(m, Gene=soupX.strainedCells.all[gg,m$rn], type="soupX"),
    data.table(m, Gene=ml.corr.values.merge[gg,m$rn], type="ML"),
    data.table(m, Gene=data.tpm.spike[gg,m$rn],type="raw"))
  ggplot(pDat, aes(x=library_id, y=Gene, color=type)) + geom_boxplot() + theme_bw(16) + 
    facet_grid(spikeIn ~ .) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + ylab(gg)
  ggsave(dirout(out, "GeneAfterCorr_",gg,".pdf"), height=8, width=9)
}



# PLOTS SHOWING CONTAMINATION IN ALL CELLS  ---------------------------------------------------
# Raw data INS vs GCG

# Cleaned data
pDat <-  fread(dirout("07_03_CleanHuman_CleanSpikeIns/", "Cont.Genes.Data_HUMAN.tsv"))
pDat$library_id <- cleanTreatment(pDat$library_id)
pDat <- pDat[!library_id %in% c("hIslets_I_pla2g16", "hIslets_II_Aold")]
table(pDat$library_id)
plot.genes <- gsub("_corrected", "", grep("_corrected", colnames(pDat), value=T))
pDat <- melt(pDat,id.vars="library_id", measure.vars=paste0(plot.genes, rep(c("_corrected", "_raw"), times=length(plot.genes))))
pDat <- cbind(pDat, do.call(rbind, strsplit(as.character(pDat$variable), "_")))
pDat$V2 <- factor(pDat$V2, levels=c("raw", "corrected"))
ggplot(pDat[V1 %in% c("INS", "GCG")], aes(x=gsub(".+Islets\\_(.+)", "\\1", library_id),y=value, fill=V2)) + 
  geom_violin(color=NA) + theme_bw(12) +  xlab("Replicate") + ylab("Log(TPM)") +
  scale_fill_manual(values=c("black", "red")) + 
  facet_grid(. ~ V1) + guides(fill=F, color=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggsave(dirout(out, "Contamination_INSvsGCG.pdf"), height=4, width=8)

ggplot(pDat[V2 == "raw" & V1 %in% c("INS", "GCG")], aes(x=gsub(".+Islets\\_(.+)", "\\1", library_id),y=value, fill=V2)) + 
  geom_violin(color=NA) + theme_bw(12) +  xlab("Replicate") + ylab("Log(TPM)") +
  scale_fill_manual(values=c("black", "red")) + 
  facet_grid(. ~ V1) + guides(fill=F, color=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggsave(dirout(out, "Contamination_INSvsGCG_original.pdf"), height=4, width=4)

# xx <- fread(dirout("EXT_HCA/Hormones.tsv"))
# xx <- data.table(melt(xx[,-"barcode",with=T]), library_id = "Human cell atlas", V2="raw")
# xx$V1 <- xx$variable
# xx <- rbind(xx, pDat, fill=T)
# ggplot(xx[V2 == "raw" & V1 %in% c("INS", "GCG")], 
#        aes(x=gsub(".+Islets\\_(.+)", "\\1", library_id),y=value, fill=V2)) + 
#   geom_violin(color=NA) + theme_bw(12) +  xlab("Replicate") + ylab("Log(TPM)") +
#   scale_fill_manual(values=c("black", "red")) + 
#   facet_grid(. ~ V1) + guides(fill=F, color=F) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))


# CONTAMINATION IN SPIKE INS --------------------------------------------------------

# HEATMAP
pDat <- cbind(
  spike.clean.diffOrg[row.names(data.tpm.cells)],
  spike.clean.sameOrg[row.names(data.tpm.cells)])
colnames(pDat) <- paste(gsub("^(.)", "\\U\\1", c("mouse", "human"), perl=T), "reference")
for(x in samples){
  for(orgX in c("mouse", "human")){
    pDat <- cbind(pDat, Matrix::rowMeans(data.tpm.spike[,meta.raw.spikes[library_id == x & spikeIn == orgX]$rn]))
    colnames(pDat)[ncol(pDat)] <- paste(gsub(".+Islets\\_(.+)", "Rep \\1", x), orgX)
  }
}
pDat <- pDat[,!grepl("pla2g16", colnames(pDat))]
pDat <- log(pDat + 1)
x <- rowMeans(pDat[,grepl("Rep", colnames(pDat))]) - rowMeans(pDat[,grepl("reference", colnames(pDat))])
tail(sort(x))
head(sort(x))
n <- 20
cleanDev()
pdf(dirout(out, "Contamination_vsSpikeIns_HM.pdf"), height=8, width=6, onefile=F)
pheatmap(pDat[c(names(head(sort(x), n)), names(tail(sort(x), n))),], cluster_rows=F)
dev.off()

# BOXPLOT
pDat <- data.table()
genes <- c("INS", "GCG","TTR")
for(orgX in c("human", "mouse")){
  for(gg in genes){
    pDat <- rbind(pDat, data.table(Gene = gg, Expression = spike.clean[gg,spike.clean.meta[org == orgX]$barcode], sample=paste("Reference", gsub("^(.)", "\\U\\1", orgX, perl=T))))}}
for(x in samples){
  for(orgX in c("mouse", "human")){
    for(gg in genes){
      pDat <- rbind(pDat, data.table(Gene = gg, Expression = data.tpm.spike["INS",meta.raw.spikes[library_id == x & spikeIn == orgX]$rn], sample=paste(gsub(".+Islets\\_(.+)", "\\1", x), orgX)))}}}
pDat[,Expression := log(Expression + 1)]
ggplot(pDat, aes(x=sample, y=Expression))+ 
  geom_boxplot(alpha=0.4) + 
  #geom_violin(color="white", fill="lightblue") + geom_jitter(height=0, width=0.1, alpha=0.1) +
  theme_bw(12) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + facet_grid(Gene ~ .) + xlab("")
ggsave(dirout(out, "Contamination_vsSpikeIns_INS_Boxplot.pdf"), height=8, width=5)

# FRACTION OF CONTAMINATION

# Contamination cross-alignment
ggplot(meta.raw.spikes[spikeIn == "mouse"], aes(x=gsub(".+Islets\\_(.+)", "\\1", library_id), y=contamination * 100)) + geom_boxplot() +
  theme_bw(12) + xlab("Replicate") + ylab("Contamination (%)")+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggsave(dirout(out, "Contamination_CrossAlignment.pdf"), height=4, width=2)

# Correlate this to reads in cells
meta.raw.spikes[spikeIn == "mouse"]
seqstats <- fread(dirout("04_SequencingStats/SequencingStats_..tsv"),check.names=T)
stopifnot(all(unique(meta.raw.spikes$old_id) %in% seqstats$dataset))
meta.raw.spikes <- cbind(meta.raw.spikes, do.call(rbind, strsplit(meta.raw.spikes$library_id, "_")))
pDat <- merge(meta.raw.spikes[spikeIn == "mouse"][,median(contamination), by=c("V3", "V2", "old_id")], 
  seqstats[dataset %in% meta.raw.spikes$old_id][,c("dataset", "Fraction.Reads.in.Cells")],
  by.x="old_id", by.y="dataset")
pDat[,rc := as.numeric(gsub("%", "", Fraction.Reads.in.Cells))]
ggplot(pDat, aes(x=V1*100, y=rc)) + geom_point(shape=1) + theme_bw(12) +
  xlab("Contamination (%)") + ylab("Reads in cells (%)") + xlim(5,20)
ggsave(dirout(out, "Contamination_CrossAlignment_vs_fractionReadsInCells.pdf"), w=2,h=2)


# correlation of contamination
write("sample,spearman",dirout(out, "Contamination_CrossAlignment_Correlation.tsv"), append=F)
for(x in unique(meta.raw.spikes$old_id)){
  cont <- fread(dirout("06_SpikeIns_Islets/",x, "_hmG/Contamination.tsv"))
  labx <- gsub("^.Islets_(I+)_(.+)", "\\2 Replicate \\1", meta.raw.spikes[old_id == x]$library_id[1])
  corX <- round(corS(cont$inHuman, cont$inMouse), 3)
  ggplot(cont, aes(x=inHuman, y=inMouse)) + geom_point(alpha=0.2) + theme_bw(12) + 
    scale_x_log10() + scale_y_log10() +
    annotate("text", x=0.01, y=100, label=paste("R =",corX )) + 
    xlab("Expression in human") + ylab("Contamination in mouse") + 
    ggtitle(labx)
  ggsave(dirout(out, "Contamination_CrossAlignment_Correlation_",x,".pdf"), height=4, width=4)
  ggsave(dirout(out, "Contamination_CrossAlignment_Correlation_",x,".jpg"), height=4, width=4)
  write(paste(x,corX, sep=","),dirout(out, "Contamination_CrossAlignment_Correlation.tsv"), append=T)
  
  pDat.org <- fread(dirout("06_SpikeIns_Islets/",x, "_hmG/IdentifyOrganism_full.tsv"))
  max.reads <- max(c(pDat.org$rh, pDat.org$rm))
  min.reads <- max(min(c(pDat.org$rh, pDat.org$rm)), 1)
  (p <- ggplot(pDat.org, aes(x=rh, y=rm, color=org, size=nGene)) + geom_point(alpha=0.3) + geom_abline()+ ggtitle(labx) +
     scale_y_log10(limits=c(min.reads,max.reads)) + scale_x_log10(limits=c(min.reads,max.reads)) + xlab("Human reads") + ylab("Mouse reads") + theme_bw(12))
  ggsave(dirout(out, "IdentifyOrganism_",x,".jpg"), height=4, width=5)
  ggsave(dirout(out, "IdentifyOrganism_",x,".pdf"), height=4, width=5)
  ggsave(dirout(out, "IdentifyOrganism2_",x,".jpg"), height=4, width=4, plot=p+guides(color=F, fill=F, size=F))
}

# Numbers for Identify Organism plots
ff <- list.files(dirout("06_SpikeIns_Islets"))
fx <- ff[1]
res <- data.table()
for(fx in ff){
  pDat.org <- fread(dirout("06_SpikeIns_Islets/",fx,"/IdentifyOrganism_full.tsv"))
  res <- rbind(res, data.table(pDat.org[,.N, by="org"], sample=fx))
}
res[, Percent := N/sum(N)*100, by="sample"]
res[,sample := gsub("_S\\d+$", "", gsub("_hmG$", "", sample))]
res <- cbind(res, do.call(rbind, strsplit(res$sample, "_")))
res[is.na(org), org := "Ambiguous"]
write.tsv(res, dirout(out, "IdentifyOrganism_AmbiguousCellNumbers.tsv"))


# NUMBER OF GENES BEFORE / AFTER CORRECTION -------------------------------
x <- as(object=corr.export, "dgCMatrix")
metaH$nGene.corr <- diff(x[,metaH$rn]@p)
x <- as(object=data.raw.cells, "dgCMatrix")
metaH$nGene.raw <- diff(x[,metaH$rn]@p)
metaH[,treatment := cleanTreatment(treatment)]
pDat <- melt(metaH[!grepl("^\\d+$", celltype2)][!grepl("_", celltype2)],measure.vars=c("nGene.corr", "nGene.raw"), id.vars=c("celltype2", "replicate", "treatment"))
ggplot(pDat, aes(x=celltype2, color=variable, y=value)) + geom_boxplot(coef=Inf) + facet_grid(treatment ~ replicate) +
  theme_bw(12) + xlab("Cell type") + ylab("Number of Genes") + scale_color_manual(values=c("red", "black")) +
  xRot() + guides(color=F)
ggsave(dirout(out, "nGene_cor_vs_raw.pdf"), height=6, width=6)



# EXAMPLE PLOT ------------------------------------------------------------
set.seed(1219)
expr.cont <- runif(5) * 1:10
for(i in c(1219, 3433)){
  set.seed(i*2)
  nGenes <- 10
  expr <- runif(5) * 1:10
  
  prop <- 0.8
  pDat <- rbind(
    data.table(value=prop * expr, gene=1:nGenes, type="true", label = "raw"),
    data.table(value=(1-prop) * expr.cont, gene=1:nGenes, type="cont", label = "raw"))
  
  #pDat$label <- factor(pDat$label, levels=c("raw", "corrected","clean"))
  ggplot(pDat, aes(x=as.factor(gene), y=value,fill=type)) + geom_bar(stat="identity") + theme_bw(12) + 
    scale_fill_manual(values=c("blue", "grey")) + ylab("TPM") + xlab("Gene") + guides(fill=F)
  ggsave(dirout(out, "Example_",i,".pdf"), width=3, height=2)
}



# LINEAR REGRESSION EXAMPLE -----------------------------------------------

x <- meta.raw.spikes$library_id[1]
pp <- list()
for(gg in c("INS", "GCG", "TTR")){
  m <- meta.raw.spikes[library_id == x & spikeIn == "mouse"]
  m$Expression <- data.tpm.spike[gg,m$rn]
  fit <- lm(contamination ~ Expression + 0,data=m)
  a = coef(fit)[1]
  pp[[gg]] <- ggplot(m, aes(x=Expression, y=contamination)) + #geom_point(alpha=.1, color="black") +
    geom_abline(intercept=0, slope=a) + ylab("") + theme_classic() +
    theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) + 
    theme(axis.text.x=element_blank(), axis.text.y=element_blank(),panel.grid=element_blank()) + 
    xlim(0, max(m$Expression, na.rm=T)) + ylim(0, max(m$contamination,na.rm=T)) + 
    xlab(gg)
}
(p <- gridExtra::grid.arrange(grobs=pp,ncol=3))
ggsave(dirout(out, "LinearRegression.pdf"), plot=p, height=1, width=3)


