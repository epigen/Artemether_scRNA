require("project.init")
project.init2("artemether")

out <- "FIG_04_SpikeInDetails/"
dir.create(dirout(out))

cleanTreatment <- function(vec){ gsub("^FoxO$", "FoxOi", gsub("^A10$", "Artemether", vec)) }
cleanTreatment2 <- function(vec){ gsub("FoxO", "FoxOi", gsub("A10", "Artemether", vec)) }

# Base contamination ------------------------------------------------------
csa <- fread(dirout("05_SpikeIns_Clean/CrossSpecies_Alignment.tsv"))
(load(dirout("05_SpikeIns_Clean/scData.RData")))
rn <- split(row.names(pbmc@data), factor(gsub("(.+?)_.+", "\\1", row.names(pbmc@data))))
cnts <- sapply(rn, function(gg) Matrix::colSums(pbmc@raw.data[gg,]))
cnts <- data.table(cnts,keep.rownames=T)
cnts[,sum := hg19 + mm10]
cnts <- merge(cnts, csa, by.x="rn", by.y="barcode")
cnts[org == "human", cont := mm10/sum]
cnts[org == "mouse", cont := hg19/sum]

ggplot(cnts, aes(x=org, y=cont * 100)) + 
  geom_violin(aes(fill=org)) +
  geom_boxplot(color="black", fill=NA,coef=Inf) + 
  theme_bw(12) + guides(color=F, fill=F) + xRot() +
  xlab("") + ylab("Contamination (%)")
ggsave(dirout(out, "Mouse_Human_ReferenceCont.pdf"), w=2,h=3)



# Base contamination external ---------------------------------------------
ff <- list.files(dirout("EXT_04_SpeciesMixtures/", "/"), full.names=T, pattern="x$")
fx <- ff[1]
res <- data.table()
for(fx in ff){
  x <- Read10X(fx)
  outfile <- dirout(out, "CrossSpecies_", gsub("_filter.+$", "", basename(fx)), ".tsv")
  if(file.exists(outfile)) next
  rn <- split(row.names(x), factor(gsub("(.+?)_.+", "\\1", row.names(x))))
  cnts <- sapply(rn, function(gg) Matrix::colSums(x[gg,]))
  cnts <- data.table(cnts,keep.rownames=T)
  cnts[,sum := hg19 + mm10]
  cnts[,ratio := log2(hg19/mm10)]
  cnts[ratio > 2, org := "human"]
  cnts[ratio < -2, org := "mouse"]
  cnts[org == "human", cont := mm10/sum]
  cnts[org == "mouse", cont := hg19/sum]
  write.tsv(cnts, outfile)
}

ff <- list.files(dirout(out), pattern="CrossSpecies", full.names=T)
xx <- lapply(ff, function(fx){
  x <- fread(fx)
  x[,sample := gsub("CrossSpecies_(.+).tsv", "\\1", basename(fx))]
  })
xx <- do.call(rbind, xx)

ggplot(xx[!is.na(org)], aes(x=org, y=cont * 100)) + 
  geom_violin(aes(fill=org)) + facet_grid(. ~ sample) + 
  geom_boxplot(color="black", fill=NA,coef=Inf) + 
  theme_bw(12) + guides(color=F, fill=F) + xRot() +
  xlab("") + ylab("Contamination (%)")
ggsave(dirout(out, "Mouse_Human_ExternalCont.pdf"), w=5,h=3)



# Correlation of contamination in single cells -----------------------------------------------
nam.map <- fread("metadata/Aggregate_Human1_2.csv")
res <- data.table()
fx <- "hIslets_I_DMSO_1_S33929_hmG"
for(fx in list.files(dirout("06_SpikeIns_Islets/"))){
  nam.clean <- nam.map[grepl(gsub("_hmG", "", fx), molecule_h5)]$library_id
  message(fx)
  print(nam.clean)
  if(length(nam.clean) == 0) next
  if(nam.clean %in% res$ds) next
  (load(dirout("06_SpikeIns_Islets/", fx, "/scData.RData")))
  org.assign <- fread(dirout("06_SpikeIns_Islets/",fx, "/IdentifyOrganism.tsv"))
  print(nrow(org.assign[org != orgX]))
  
  org.idx <- substr(row.names(pbmc@data), 0,4) == "hg19"
  orgX <- if(grepl("hIslets", fx)) "human" else "mouse"
  genes.use <- if(orgX == "human") row.names(pbmc@data)[org.idx] else row.names(pbmc@data)[!org.idx]
  m <- pbmc@data[genes.use,org.assign[org != orgX]$barcode]
  m <- toTPM_fromRaw(m,scale.factor=1e6)
  corMT <- cor(as.matrix(m)) #[,sample(colnames(m), 100)]))
  res <- rbind(res, data.table(cor=corMT[upper.tri(corMT)], ds=nam.clean))
}

res[,treatment := gsub(".Islets_I+_", "", cleanTreatment2(ds))]
res[,Donor := paste("Donor", gsub(".Islets_(I+)_.+", "\\1", ds))]
ggplot(res[!treatment %in% c("Aold", "pla2g16")], aes(x=treatment, y=cor)) + geom_boxplot(coef=Inf) + facet_grid(. ~ Donor) + theme_bw(12) + xRot() +
  xlab("") + ylab("Correlation")
ggsave(dirout(out, "Contamination_Correlation_singleCells.pdf"), width=4,h=4)



# INS / GCG plots ---------------------------------------------------------
metaH <- loadHMeta2()
pDat <-  fread(dirout("07_03_CleanHuman_CleanSpikeIns/", "Cont.Genes.Data_HUMAN.tsv"))
metaH <- merge(metaH, pDat[,c("rn", "INS_raw", "GCG_raw")], by="rn")
metaH <- metaH[treatment == "DMSO"]
pDat <- melt(metaH, measure.vars=c("INS_raw", "GCG_raw"))
pDat[,variable := gsub("_raw", "", variable)]
ggplot(pDat, aes(x=value, color=replicate, fill=replicate)) +
  xlim(-1, max(pDat$value)) + 
  geom_density() + facet_wrap(~variable, ncol=1, strip.position="right") +
  theme_bw(12) +
  theme(panel.grid=element_blank()) +
  xlab("Expression log(TPM)") + ylab("Density") + guides(linetype=F, fill=F, color=F) +
  scale_fill_manual(values=c(I="#b2df8a99", II="#cab2d699")) +
  scale_color_manual(values=c(I="#b2df8a", II="#cab2d6")) +
  geom_vline(xintercept=0, color="#e31a1c", linetype=2)
ggsave(dirout(out, "INS_GCG.pdf"), width=4,h=4)
