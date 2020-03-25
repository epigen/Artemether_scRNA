require("project.init")
require(Seurat)
require(Matrix)
require(methods)
require(SoupX)
project.init2("artemether")

out <- "07_04_CleanHuman3_CleanSpikeIns/"
dir.create(dirout(out))


# LOAD DATA ---------------------------------------------------------------
data.raw <- Read10X(paste0(getOption("PROCESSED.PROJECT"), "results_pipeline/cellranger_count/","Human3","/outs/raw_gene_bc_matrices_mex/GRCh38/"))

# external spike ins
(load(dirout("05_SpikeIns_Clean/SpikeData_human.RData")))
spike.clean <- add.genes.sparseMatrix(toTPM(dat.log.spike), row.names(data.raw))
spike.clean.meta <- fread(dirout("05_SpikeIns_Clean/SpikeMeta_human.tsv"))
spike.clean.diffOrg <- Matrix::rowMeans(spike.clean[,spike.clean.meta[org == "mouse"]$barcode])
spike.clean.sameOrg <- Matrix::rowMeans(spike.clean[,spike.clean.meta[org == "human"]$barcode])

pancGenes <- intersect(unique(fread("metadata//PancreasMarkers.tsv")$Human_GeneName), row.names(data.raw))

# PREPARE CELL ANNOTATION -------------------------------------------------
sample.annot <- fread("metadata/Aggregate_Human3.csv")
sample.annot$nr <- as.character(1:nrow(sample.annot))
sample.annot[,old_id := gsub(".+cellranger_count\\/(.+)\\/outs.+", "\\1", molecule_h5)]

meta.raw.file <- dirout(out, "Meta.raw.RData")
if(!file.exists(meta.raw.file)){
  meta.raw <- data.table(rn = colnames(data.raw))
  meta.raw[,barcode := gsub("(.+)\\-(\\d+)", "\\1", rn)]
  meta.raw[,nr := gsub("(.+)\\-(\\d+)", "\\2", rn)]
  meta.raw <- merge(meta.raw, sample.annot, by="nr", all.x=T)
  meta.raw$molecule_h5 <- NULL
  
  islet.spike.annot <- do.call(rbind, lapply(list.files(dirout("06_SpikeIns_Islets/"), pattern="^hIslets"), function(x){
    data.table(fread(dirout("06_SpikeIns_Islets/", x, "/IdentifyOrganism.tsv")), old_id = gsub("_hmG$", "", x))}))
  
  meta.raw[,old_id := gsub("_human$", "", old_id)]
  table(unique(meta.raw$old_id) %in% islet.spike.annot$old_id)
  meta.raw <- merge(meta.raw, islet.spike.annot, by=c("barcode", "old_id"), all.x=T)
  meta.raw$nUMI <- Matrix::colSums(data.raw[,meta.raw$rn])
  meta.raw$nGene <- Matrix::colSums(data.raw[,meta.raw$rn] != 0)
  save(meta.raw, file=meta.raw.file)
} else {
  load(meta.raw.file)  
}
# meta.raw <- meta.raw[grepl("DMSO", library_id)]
data.raw <- data.raw[,meta.raw$rn]
table(meta.raw$library_id)
str(data.raw)

# Define cells
meta.raw.cells <- meta.raw[nUMI >= MIN.UMIS & nGene >= MIN.GENES]
data.raw.cells <- data.raw[,meta.raw.cells$rn]
data.tpm.cells <- t(t(data.raw.cells) / Matrix::colSums(data.raw.cells)) * SCALE.FACTOR
data.log.cells <- toLogTPM(data.tpm.cells)

# Define spike Ins
meta.raw.spikes <- meta.raw[spikeIn != "NA"]
with(meta.raw.spikes, table(library_id, spikeIn))
data.raw.spike <- data.raw[,meta.raw.spikes$rn]
data.tpm.spike <- t(t(data.raw.spike) / Matrix::colSums(data.raw.spike)) * SCALE.FACTOR

# Samples
samples <- unique(meta.raw$library_id)


# SOUP ANALYSIS -----------------------------------------------------------
soupX.file <- dirout(out, "Soup.results.RData")
if(!file.exists(soupX.file)){
  #estimate soup
  soupX.soup <- lapply(samples, function(x) estimateSoup(tod=data.raw[,meta.raw[library_id == x]$rn]))
  names(soupX.soup) <- samples
  soupX.soup.tpm <- do.call(cbind, lapply(soupX.soup, function(dt){ret <- dt$cnts/dt$total * SCALE.FACTOR; names(ret) <- row.names(dt); ret}))
  
  # estimate contaminated fraction
  soupX.contFrac <- lapply(samples, function(x){
    calculateContaminationFraction(
      toc=data.raw[,meta.raw.spikes[library_id == x]$rn], 
      soupProfile=soupX.soup[[x]], 
      nonExpressedGeneList=pancGenes,
      excludeMethod="qCut")})
  names(soupX.contFrac) <- samples
  
  # estimate rho
  soupX.rho <- lapply(samples, function(x){
    interpolateCellContamination(
      rhos=soupX.contFrac[[x]], 
      nUMIs=Matrix::colSums(data.raw[,meta.raw.spikes[library_id == x]$rn]))})
  names(soupX.rho) <- samples
  str(soupX.rho)
  soupX.rho.agg <- do.call(c, lapply(samples, function(x){names(soupX.rho[[x]]) <- meta.raw.spikes[library_id == x]$rn; soupX.rho[[x]]}))
  
  # clean cell
  soupX.strainedCells <- lapply(samples, function(x){
    strainCells(
      toc=data.raw[,meta.raw.spikes[library_id == x]$rn], 
      rhos=soupX.rho[[x]],
      soupProfiles=soupX.soup[[x]], 
      retainSparsity=TRUE)})
  names(soupX.strainedCells) <- samples
  soupX.strainedCells.all <- do.call(cbind, soupX.strainedCells)
  soupX.strainedCells.all <- t(t(soupX.strainedCells.all)/Matrix::colSums(soupX.strainedCells.all))*SCALE.FACTOR
  stopifnot(all(round(Matrix::colSums(soupX.strainedCells.all)) == SCALE.FACTOR))
  
  # save
  save(list=ls()[grepl("^soupX\\.", ls())], file=soupX.file)
} else {
  (load(soupX.file))
}


# CALCULATE ESTIMATED SOUP FROM SPIKE INS ---------------------------------
spikeIn.soups <- do.call(cbind, lapply(samples, function(x){
  m <- meta.raw.spikes[library_id == x & spikeIn == "mouse"]
  cont <- t(t(data.tpm.spike[,m$rn]) * 1/(1-m$contamination)) - spike.clean.diffOrg[row.names(data.tpm.spike)]
  cont[cont < 0] <- 0
  soup <- Matrix::rowMeans(t(t(cont) / m$contamination))
  #soup <- Matrix::rowMeans(t(t(cont) / Matrix::colSums(cont) * SCALE.FACTOR))
  return(soup)
}))
colnames(spikeIn.soups) <- paste0(samples, " spikeIn")
spikeIn.soups.tpm <- t(t(spikeIn.soups)/colSums(spikeIn.soups)) * SCALE.FACTOR


# ESTIMATED RHOS FROM MACHINE LEARNING ------------------------------------------------------- 
(cont.genes <- intersect(names(tail(sort(rowMeans(spikeIn.soups)), 100)), pancGenes))

# Cross validation
gg <- "INS"
pred.cv <- list()
for(gg in cont.genes){
  message(gg)
  pred.cv[[gg]] <- c()
  sampleX <- meta.raw.spikes$library_id[1]
  for(sampleX in unique(meta.raw.spikes$library_id)){
    m <- meta.raw.spikes[library_id == sampleX & spikeIn == "mouse"]
    folds <- caret::createFolds(m$contamination, k=3)
    f <- folds[[1]]
    for(f in folds){
      lm.dat <- data.table(cont=m$contamination, gene=data.tpm.spike[gg,m$rn], rn=m$rn)
      fit <- lm(cont ~ gene + 0, data=lm.dat[-f])
      res <- predict(fit, lm.dat[f])
      names(res) <- lm.dat[f]$rn
      pred.cv[[gg]] <- c(pred.cv[[gg]], res)
    }}}
ml.rho <- do.call(cbind, lapply(pred.cv, function(vec) vec[meta.raw.spikes$rn]))
ml.rho[ml.rho == 0] <- NA
ml.rho <- ml.rho[!is.na(row.names(ml.rho)),]
str(ml.rho)
ml.rho.cv <- rowMedians(ml.rho, na.rm=T)
names(ml.rho.cv) <- row.names(ml.rho)

# Prediction
pred.final <- list()
gg <- "INS"
for(gg in cont.genes){
  message(gg)
  pred.final[[gg]] <- c()
  sampleX <- meta.raw.cells$library_id[1]
  for(sampleX in unique(meta.raw.cells$library_id)){
    # Prediction
    m <- meta.raw.cells[library_id == sampleX]
    m$gene <- data.tpm.cells[gg,m$rn]
    fit <- lm(contamination ~ gene + 0, data=m[spikeIn == "mouse"])
    res <- predict(fit, m);names(res) <- m$rn
    pred.final[[gg]] <- c(pred.final[[gg]], res)
  }}

ml.rho.final <- do.call(cbind, lapply(pred.final, function(vec) vec[meta.raw.cells$rn]))
ml.rho.final[ml.rho.final == 0] <- NA
ml.rho.final <- ml.rho.final[!is.na(row.names(ml.rho.final)),]
str(ml.rho.final)
ml.rho.pred <- rowMedians(ml.rho.final, na.rm=T)
names(ml.rho.pred) <- row.names(ml.rho.final)
ml.rho.pred.with.cv <- ml.rho.pred
ml.rho.pred.with.cv[names(ml.rho.cv)] <- ml.rho.cv

ggplot(data.table(contamination=ml.rho.pred), aes(x=contamination)) + geom_density() + xlim(0,1) + geom_vline(xintercept=0.2, color='red', size=1) +theme_bw(12)
ggsave(dirout(out, "Capped_Contamination.pdf"), width=4, heigh=4)
#ml.rho.pred.with.cv[ml.rho.pred.with.cv > 0.2] <- 0.2
ml.rho.pred[ml.rho.pred > 0.2] <- 0.2

# CORRECT DATA FROM SPIKE-INS in mouse and human spike-ins ---------------------------------------------
ml.corr.values <- lapply(samples, function(x){
  rn <- meta.raw.spikes[library_id == x]$rn
  ret <- data.tpm.spike[,rn] - (spikeIn.soups[row.names(data.tpm.spike),paste(x, "spikeIn")] %o% ml.rho.pred.with.cv[rn])
  ret[ret < 0] <- 0
  t(t(ret) / Matrix::colSums(ret)) * SCALE.FACTOR
})
ml.corr.values.merge <- do.call(cbind, ml.corr.values)


# CORRECT ALL DATA FROM SPIKE-INS ---------------------------------------------
ml.corr.values.all <- lapply(samples, function(x){
  rn <- meta.raw.cells[library_id == x]$rn
  ret <- data.tpm.cells[,rn] - (spikeIn.soups[row.names(data.tpm.cells),paste(x, "spikeIn")] %o% ml.rho.pred[rn])
  ret[ret < 0] <- 0
  t(t(ret) / Matrix::colSums(ret)) * SCALE.FACTOR
})
ml.corr.values.all.merge <- do.call(cbind, ml.corr.values.all)


corr.export <- ml.corr.values.all.merge
nUMIs <- Matrix::colSums(data.raw.cells[,colnames(corr.export)])
corr.export <- t(t(corr.export/SCALE.FACTOR) * nUMIs)
stopifnot(all(round(Matrix::colSums(corr.export)) == round(nUMIs)))
save(corr.export, ml.corr.values.all.merge, file=dirout(out, "CorrectedData.RData"))
save(ml.rho.pred, ml.rho.cv, ml.rho.pred.with.cv, data.raw.cells, meta.raw.cells, spikeIn.soups, file=dirout(out, "Uncorrected.Data.RData"))


# COMPARE ESTIMATED TO SPIKE IN TRUTH -------------------------------------------

# COMPARE ESTIMATED RHO

# Soup
meta.raw.spikes$contSoupX <- soupX.rho.agg[meta.raw.spikes$rn]
ggplot(meta.raw.spikes[!is.na(contamination)], aes(x=contSoupX, y=contamination, color=library_id)) +  geom_point(alpha=0.5) +
  theme_bw(12) + xlim(0,.7) + ylim(0,.7) + geom_abline(size=2, color="lightgrey", alpha=0.5) + 
  xlab("Contamination soupX") + ylab("Contamination cross-alignment") + guides(color=F)
ggsave(dirout(out, "Rho_estimation_SOUP.pdf"), height=4, width=4)

# machine learning
meta.raw.spikes$contMLCV <- ml.rho.cv[meta.raw.spikes$rn]
with(meta.raw.spikes, table(is.na(contMLCV), spikeIn))
ggplot(meta.raw.spikes[!is.na(contamination)], aes(x=contMLCV, y=contamination, color=library_id)) +  geom_point(alpha=0.5) +
  theme_bw(12) + xlim(0,.4) + ylim(0,.4) + geom_abline(size=2, color="lightgrey", alpha=0.5) +
  xlab("Contamination spike-ins") + ylab("Contamination cross-alignment") + guides(color=F)
ggsave(dirout(out, "Rho_estimation_ML.pdf"), height=4, width=4)



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
pdf(dirout(out, "Soup_estimation_Correlation.pdf"), width=11, height=10, onefile=F)
pheatmap(corS(cbind(spikeIn.soups.tpm, soupX.soup.tpm[row.names(spikeIn.soups),])),cluster_rows=F, cluster_cols=F); 
dev.off()


# COMPARE CORRECTED VALUES
compare.values <- data.table()
for(dataType in c("original", "soupX", "spikeIns")){
  for(orgXX in c("human", "mouse")){
    datX <- data.tpm.spike; 
    if(dataType == "soupX") datX <- soupX.strainedCells.all; 
    if(dataType == "spikeIns") datX <- ml.corr.values.merge
    spikeRef <- spike.clean.sameOrg; if(orgXX != "human") spikeRef <- spike.clean.diffOrg
    corRes <- lapply(samples, function(x){
      apply(as.matrix(datX[, meta.raw.spikes[library_id == x & spikeIn == orgXX]$rn]), 2, function(col) cor(col, spikeRef))
    }); names(corRes) <- paste(samples, dataType, orgXX, sep="_")
    compare.values <- rbind(compare.values, do.call(rbind, lapply(names(corRes), function(nam) data.table(Group=nam, value=corRes[[nam]], rn=names(corRes[[nam]])))))
  }
}
compare.values <- cbind(compare.values, data.table(do.call(rbind, strsplit(compare.values$Group, "_"))))
ggplot(compare.values, aes(x=paste(V3,V4), y=value, color=V6)) + geom_boxplot() + facet_grid(V5 ~ ., space="free_x", scale="free_x") +
  theme_bw(12) + ylab("Correlation to reference") + xlab("Replicate") + xRot()
ggsave(dirout(out, "Values_compare.pdf"), height=4, width=8)

ggplot(compare.values[V6 != "soupX"], aes(x=paste(V3,V4), y=value, color=V6)) + geom_boxplot() + facet_grid(V5 ~ ., space="free_x", scale="free_x") +
  theme_bw(12) + guides(color=F) + ylab("Correlation to reference") + xlab("Replicate") +
  scale_color_manual(values=c("black", "red")) + xRot()
ggsave(dirout(out, "Values_compare_NoSoupX.pdf"), height=4, width=8)
write.tsv(compare.values, file=dirout(out, "Values_compare.tsv"))

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



# CONTAMINATION FIGURES ---------------------------------------------------
# Raw data INS vs GCG

# Cleaned data
pDat <- meta.raw.cells
plot.genes <- cont.genes
for(gg in plot.genes){
  pDat[[paste0(gg, "_corrected")]] <- log(ml.corr.values.all.merge[gg,][pDat$rn] + 1)
  pDat[[paste0(gg, "_raw")]] <- data.log.cells[gg,][pDat$rn]
}
write.tsv(pDat,dirout(out, "Cont.Genes.Data.tsv"))
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


# SPIKE IN FIGURES --------------------------------------------------------

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
pDat <- log(pDat + 1)
x <- rowMeans(pDat[,grepl("Rep", colnames(pDat))]) - rowMeans(pDat[,grepl("reference", colnames(pDat))])
tail(sort(x))
head(sort(x))
n <- 20
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


# Contamination cross-alignment
ggplot(meta.raw.spikes[spikeIn == "mouse"], aes(x=gsub(".+Islets\\_(.+)", "\\1", library_id), y=contamination * 100)) + geom_boxplot() +
  theme_bw(12) + xlab("Replicate") + ylab("Contamination (%)")+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggsave(dirout(out, "Contamination_CrossAlignment.pdf"), height=4, width=2)
for(x in unique(meta.raw.spikes$old_id)){
  cont <- fread(dirout("06_SpikeIns_Islets/",x, "_hmG/Contamination.tsv"))
  ggplot(cont, aes(x=inHuman, y=inMouse)) + geom_point(alpha=0.2) + theme_bw(16) + 
    scale_x_log10() + scale_y_log10() +
    annotate("text", x=0.001, y=100, label=paste("R =", round(corS(cont$inHuman, cont$inMouse), 3))) + 
    xlab("Expression in human") + ylab("Contamination in mouse")
  ggsave(dirout(out, "Contamination_CrossAlignment_Correlation_",x,".pdf"), height=4, width=4)
  ggsave(dirout(out, "Contamination_CrossAlignment_Correlation_",x,".jpg"), height=4, width=4)
  
  pDat.org <- fread(dirout("06_SpikeIns_Islets/",x, "_hmG/IdentifyOrganism_full.tsv"))
  max.reads <- max(c(pDat.org$rh, pDat.org$rm))
  min.reads <- max(min(c(pDat.org$rh, pDat.org$rm)), 1)
  (p <- ggplot(pDat.org, aes(x=rh, y=rm, color=org, size=nGene)) + geom_point(alpha=0.3) + geom_abline()+
     scale_y_log10(limits=c(min.reads,max.reads)) + scale_x_log10(limits=c(min.reads,max.reads)) + xlab("Human reads") + ylab("Mouse reads") + theme_bw(12))
  ggsave(dirout(out, "IdentifyOrganism_",x,".jpg"), height=4, width=5)
  ggsave(dirout(out, "IdentifyOrganism2_",x,".jpg"), height=4, width=4, plot=p+guides(color=F, fill=F, size=F))
}




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
