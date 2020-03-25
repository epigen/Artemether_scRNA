require("project.init")
project.init2("artemether")

out <- "FIG_03_DiffEXPR/"
dir.create(dirout(out))



# CLEAN FUNCTION ----------------------------------------------------------
cleanTreatment <- function(vec){ gsub("^FoxO$", "FoxOi", gsub("^A10$", "Artemether", vec)) }

# LOAD DATA ---------------------------------------------------------------
plot.data.files <- dirout(out, "Data.RData")
if(!file.exists(plot.data.files)){
  metaH <- loadHMeta2()
  pbmcH <- loadHData2()
  metaM <- loadMMeta2()
  pbmcM <- loadMData2()
  metaH3 <- loadH3Meta2()
  pbmcH3 <- loadH3Data2()
  
  metaM[,INS := Ins2]
  metaM[,GCG := Gcg]
  metaH <- metaH[!sample %in% c("hIslets_II_GABA", "hIslets_I_pla2g16")]
  dat.log <- list(human=pbmcH@data[,metaH$rn], mouse=pbmcM@data[,metaM$rn], human3=pbmcH3@data[,metaH3$rn])
  dat.raw <- list(human=pbmcH@raw.data[,metaH$rn], mouse=pbmcM@raw.data[,metaM$rn], human3=pbmcH3@raw.data[,metaH3$rn])
  meta <- list(human=metaH, mouse=metaM, human3=metaH3)
  #   dat.log$human3 <- pbmcH3@data[,metaH3$rn]
  #   dat.raw$human3 <- pbmcH3@raw.data[,metaH3$rn]
  #   meta$human3 <- metaH3
  save(dat.log, dat.raw, meta, file=plot.data.files)
} else {
  load(plot.data.files)
}

meta <- lapply(meta, function(x){x$treatment <- cleanTreatment(x$treatment);x})
save(meta, file=dirout(out, "Meta.RData"))
celltypes.of.interest <- c("Alpha", "Beta")
use.treatments <- c("Artemether", "GABA", "A1", "FoxOi")



# MetaData export ---------------------------------------------------------
meta.export <- meta
gg <- c("Ins1", "Ins2", "Gcg", "Sst", "Ppy","INS", "GCG", "SST", "PPY", "PDX1", "MAFA")
for(mnam in names(meta.export)){
  for(g in gg){
    if(g %in% row.names(dat.log[[mnam]])){
      meta.export[[mnam]][[g]] <- dat.log[[mnam]][g, meta[[mnam]]$rn]
    }
  }
  write.tsv(meta.export[[mnam]], dirout(out, "MetaData_Export_", mnam, ".tsv"))
}


# Clean log data (normalized after removing insulin, glucagon,...)
genes.exclude <- c("Ins1", "Ins2", "Gcg", "Sst", "Ppy","INS", "GCG", "SST", "PPY")
dat.log.clean <- list()
for(org in names(dat.raw)){
  dat.log.clean[[org]] <- toLogTPM(toTPM_fromRaw(dat.raw[[org]][!row.names(dat.raw[[org]]) %in% genes.exclude,], scale.factor=1e6))
}


# Load diff genes
diffGenes.file <- dirout(out, "DiffGenes.RData")
if(!file.exists(diffGenes.file)){
  
  (load(dirout("16_H_02_DiffGenes/VsDMSO_negbinom_noIns/AllRes.RData")))
  resH <- res[V2 != "pla2g16" & !(V1.1=="II" & V2=="GABA")]
  resH[,org:="human"]
  
  (load(dirout("16_M_02_DiffGenes/VsDMSO_negbinom_noIns/AllRes.RData")))
  resM <- res
  resM[,org:="mouse"]
  
  (load(dirout("16_H3_02_DiffGenes/VsDMSO_negbinom_noIns/AllRes.RData")))
  resH3 <- res
  resH3Ag <- copy(resH3)
  resH3Ag[,direction:=ifelse(avg_diff > 0, "up", "down")]
  resH3Ag <- resH3Ag[,.(p_val = fishers.method(p_val), dir = length(unique(direction)), avg_diff = mean(avg_diff), pct.1 = mean(pct.1), pct.2 = mean(pct.2)), by=c("V1", "V2", "V3", "V4")]
  resH3Ag[,file := paste(V2, V3, V4, sep="_")] 
  resH3Ag <- resH3Ag[dir == 1][,c("V1", "p_val", "avg_diff", "pct.1", "pct.2", "file", "V2", "V3", "V4")]
  resH3Ag[,direction:=ifelse(avg_diff > 0, "up", "down")]
  resH3Ag[,org:="human3aggregated"]
  resH3[,org:="human3"]
  
  resH3[,timepoint := V2]
  resH3$V2 <- NULL
  resH3Ag[,timepoint := V2]
  resH[,timepoint := NA]
  resM[,timepoint := NA]
  
  resH[,qval := p.adjust(p_val, method="BH")]
  resM[,qval := p.adjust(p_val, method="BH")]
  resH3[,qval := p.adjust(p_val, method="BH")]
  resH3Ag[,qval := p.adjust(p_val, method="BH")]
  
  resX <- list(resH, resM, resH3, resH3Ag)
  resX <- lapply(resX, function(x){
    colnames(x)[1:9] <- c("gene", "pval", "logFC", "perc1", "perc2", "group", "replicate", "treatment", "celltype")
    return(x)
    })
  resX <- do.call(rbind, resX)
  source("src/FUNCTIONS_Synonyms_etc.R")
  resX <- merge(resX, homol, by.x="gene", by.y="human", all.x=T)
  resX <- merge(resX, homol, by.x="gene", by.y="mouse", all.x=T)
  resX[org == "mouse", mouse := gene]
  resX[grepl("human", org), human := gene]
  diffGenes <- resX
  diffGenes[,direction:=ifelse(logFC > 0, "up", "down")]
  save(diffGenes, file=diffGenes.file)
  write.tsv(homol, dirout(out, "Homology_unique.tsv"))
  write.tsv(diffGenes, dirout(out, "DiffGenes_byReplicate.tsv"))
  write.tsv(diffGenes[celltype %in% c("Beta", "Alpha")], dirout(out, "DiffGenes_byReplicate_AlphaBetaOnly.tsv"))
} else {
  (load(diffGenes.file))
}
diffGenes$treatment <- cleanTreatment(diffGenes$treatment)
table(diffGenes$treatment)



# HORMONE DIFF RES --------------------------------------------------------
(load(dirout("16_H_02_DiffGenes/VsDMSO_negbinom_noIns/AllRes.RData")))
resH <- res[V2 != "pla2g16" & !(V1.1=="II" & V2=="GABA")]
resH[,org:="human"]

(load(dirout("16_M_02_DiffGenes/VsDMSO_negbinom_noIns/AllRes.RData")))
resM <- res
resM[,org:="mouse"]

(load(dirout("16_H3_02_DiffGenes/VsDMSO_negbinom_noIns/AllRes.RData")))
resH3 <- res

# Supplementary table
diffGenes[,direction := ifelse(logFC > 0, "up", "down")]
write.tsv(diffGenes[!treatment %in% c("A1", "Aold")][qval < 0.1][,-"group"], dirout(out, "DiffGenes_SuppTable_q0.1.tsv"))



# DIFF GENES CHERRY PICKED EXAMPLES -----------------------------------------------------
# x <- data.table(gdata::read.xls("metadata/TopHits_v2_Brenda_2018_09_17.xlsx",sheet=1))[,1:4,with=F]
# names(x) <- paste(rep(c("human", "mouse"),each=2), x[1,])
# x <- as.list(x[-1,])
# x$"human beta 2" <- data.table(gdata::read.xls("metadata/TopHits_v2_Brenda_2018_09_17.xlsx",sheet=2))$Gene
# x$"mouse beta 2" <- data.table(gdata::read.xls("metadata/TopHits_v2_Brenda_2018_09_17.xlsx",sheet=2))$Gene
# xnam <- names(x)[1]
# for(xnam in names(x)){
#   xnam.split <- strsplit(xnam, " ")[[1]]
#   pDat <- diffGenes[toupper(gene) %in% c(x[[xnam]])][org == xnam.split[1]][toupper(celltype) == toupper(xnam.split[2])]
#   pDat <- hierarch.ordering(pDat, "gene", "group","logFC")
#   ggplot(pDat, aes(x=replicate, y=gene, color=logFC, size=pmin(5, -log10(qval)))) + 
#     geom_point() + facet_grid(celltype ~ treatment, scales="free", space="free") + 
#     scale_color_gradient2(name="logFC", high="indianred", low="navyblue") + theme_bw(16) + 
#     scale_size_continuous(name=expression(-log[10](q))) +
#     xlab("Replicate") + ylab("")
#   ggsave(dirout(out, "DiffGenes_", xnam, ".pdf"), width=length(unique(pDat$group)) * 0.35 + 3.5, height=length(unique(pDat$gene))*0.25+1)
# }

# mouse vs human
# pDat <- diffGenes[!is.na(mouse) & !is.na(human)][human %in% do.call(c, x[c("human beta", "mouse beta")]) | toupper(mouse) %in% do.call(c, x[c("human beta", "mouse beta")])][treatment == "Artemether" & celltype == "Beta"]
# pDat[, id := paste(group, org)]
# pDat <- hierarch.ordering(pDat, "human", "id","logFC")
# ggplot(pDat, aes(x=replicate, y=human, color=logFC, size=pmin(5, -log10(qval)))) + 
#   geom_point() + facet_grid(celltype ~ org, scales="free", space="free") + 
#   scale_color_gradient2(name="logFC", high="indianred", low="navyblue") + theme_bw(16) + 
#   scale_size_continuous(name=expression(-log[10](q))) +
#   xlab("Replicate") + ylab("")
# ggsave(dirout(out, "DiffGenes_MouseVsHuman_Artemether_Beta", ".pdf"), width=length(unique(pDat$id)) * 0.35 + 2.5, height=length(unique(pDat$human))*0.25+1)



# LOG FOLD CHANGES -----------------------------------------
avg.logFC.file <- dirout(out, "Avg.logfc.RData")
if(!file.exists(avg.logFC.file)){
  avg <- list()
  org <- names(meta)[1]
  for(org in c("human", "mouse")){
    metaX <- meta[[org]][celltype2 %in% c("Beta", "Alpha")]
    avg[[org]] <- sapply(split(metaX$rn, factor(with(metaX, paste(org, celltype2, treatment, replicate)))), function(rns){
      log(Matrix::rowMeans(toTPM(dat.log.clean[[org]][,rns])) + 1)
      })
  }
  logfcs <- list()
  for(org in names(meta)){
    x <- avg[[org]]
    for(i in colnames(x)){
      x[,i] <- avg[[org]][,i] - avg[[org]][,gsub("(.+) (.+) (.+) (.+)", "\\1 \\2 DMSO \\4",i)]
    }
    logfcs[[org]] <- x[,!grepl("DMSO", colnames(x))]
  }
  source("src/FUNCTIONS_Synonyms_etc.R")
  
  # Averages mapped to homologs
  avg.homol <- avg
  row.names(avg.homol$mouse) <- homol[match(row.names(avg.homol$mouse), homol$mouse)]$human
  avg.homol$mouse <- avg.homol$mouse[!is.na(row.names(avg.homol$mouse)) & row.names(avg.homol$mouse) %in% row.names(avg.homol$human),]
  avg.homol$human <- avg.homol$human[row.names(avg.homol$mouse),]
  stopifnot(all(row.names(avg.homol$human) == row.names(avg.homol$mouse)))
  sapply(avg.homol, dim)
  
  # logfc mapped to homologs
  logfcs.homol <- logfcs
  row.names(logfcs.homol$mouse) <- homol[match(row.names(logfcs.homol$mouse), homol$mouse)]$human
  logfcs.homol$mouse <- logfcs.homol$mouse[!is.na(row.names(logfcs.homol$mouse)) & row.names(logfcs.homol$mouse) %in% row.names(logfcs.homol$human),]
  logfcs.homol$human <- logfcs.homol$human[row.names(logfcs.homol$mouse),]
  stopifnot(all(row.names(logfcs.homol$human) == row.names(logfcs.homol$mouse)))
  sapply(logfcs.homol, dim)
  
  # write files
  for(org in names(avg)){
    write.csv(avg[[org]], dirout(out, "Averages_", org, ".csv"))
    write.csv(avg.homol[[org]], dirout(out, "Averages_Homol_", org, ".csv"))
    write.csv(logfcs[[org]], dirout(out, "LogFC_", org, ".csv"))
    write.csv(logfcs.homol[[org]], dirout(out, "LogFC_Homol_", org, ".csv"))
  }
  save(avg, avg.homol, logfcs, logfcs.homol, file=avg.logFC.file)
} else {
  (load(avg.logFC.file))
}


# Artemether vs FoxOi replicates -------------------------------------------------------------------

# plot Artemether vs FoxOi for each replicate
org <- "human"
repl <- "I"
for(org in c("human", "mouse")){
  for(repl in unique(meta[[org]]$replicate)){
    pDat <- data.table(
      FoxOi = logfcs[[org]][,paste(org, "Beta", "FoxOi", repl)], 
      Artemether = logfcs[[org]][,paste(org, "Beta", "Artemether", repl)],
      Gene=row.names(logfcs[[org]]))
    
    pDat[Gene %in% diffGenes[celltype == "Beta" & treatment == "FoxOi" & replicate == repl & qval < 0.05]$gene, FoxOisig := "FoxOi"]
    pDat[Gene %in% diffGenes[celltype == "Beta" & treatment == "Artemether" & replicate == repl & qval < 0.05]$gene, AmSig := "Am"]
    pDat[,sig := gsub("_NA$", "", gsub("^NA_", "", paste0(FoxOisig, "_", AmSig)))]
    
    ggplot(pDat, aes(x=FoxOi, y=Artemether)) + 
      geom_point(data=pDat[sig == "NA"], alpha=0.2) +
      geom_point(data=pDat[sig != "NA"], aes(color=sig), alpha=0.2) +
      theme_bw(16) + ggtitle(round(cor(pDat$FoxOi, pDat$Artemether, method="spearman"),3))
    ggsave(dirout(out, "Diff_FoxOi_Artemether_", org, "_", repl, ".pdf"), height=4, width=5)
    ggsave(dirout(out, "Diff_FoxOi_Artemether_", org, "_", repl, ".jpg"), height=4, width=5) 
  }
}


# Artemether vs FoxOi averages ----------------------------------------------------
logfcs.comb <- cbind(logfcs.homol$human, logfcs.homol$mouse)
logfcs.comb <- sapply(split(colnames(logfcs.comb), factor(gsub(" I+", "", colnames(logfcs.comb)))), function(cc) rowMeans(logfcs.comb[,cc,drop=F]))
pDT <- data.table(logfcs.comb)
colnames(pDT) <- make.names(colnames(pDT))
ggplot(pDT, aes(y=human.Beta.Artemether,x =human.Beta.FoxOi)) + geom_point(alpha=0.1) +
  theme_bw(16) + ggtitle(round(cor(pDT$human.Beta.Artemether, pDT$human.Beta.FoxOi, method="spearman"),3))
ggsave(dirout(out, "Diff_FoxOi_Artemether_human_comb.pdf"), height=4, width=4)
ggsave(dirout(out, "Diff_FoxOi_Artemether_human_comb.jpg"), height=4, width=4)

ggplot(pDT, aes(y=mouse.Beta.Artemether,x =mouse.Beta.FoxOi)) + geom_point(alpha=0.1) +
  theme_bw(16) + ggtitle(round(cor(pDT$mouse.Beta.Artemether, pDT$mouse.Beta.FoxOi, method="spearman"),3))
ggsave(dirout(out, "Diff_FoxOi_Artemether_mouse_comb.pdf"), height=4, width=4)
ggsave(dirout(out, "Diff_FoxOi_Artemether_mouse_comb.jpg"), height=4, width=4)


# mouse vs human -------------------------------------------------------------------

# plot mouse vs human Artemether
pDat <- data.table(
  human= rowMeans(logfcs.homol$human[,grepl("Beta Artemether", colnames(logfcs.homol$human))]),
  mouse= rowMeans(logfcs.homol$mouse[,grepl("Beta Artemether", colnames(logfcs.homol$mouse))]),
  gene= row.names(logfcs.homol$human))
ggplot(pDat, aes(x=human, y=mouse)) + geom_point(alpha=0.2) +
  ggtitle(round(cor(pDat$mouse, pDat$human, method="spearman"),3))
ggsave(dirout(out, "Diff_human_mouse_Artemether.pdf"), height=4, width=5)
ggsave(dirout(out, "Diff_human_mouse_Artemether.jpg"), height=4, width=5) 

# plot mouse vs human FoxOi
pDat <- data.table(
  human= rowMeans(logfcs.homol$human[,grepl("Beta FoxOi", colnames(logfcs.homol$human))]),
  mouse= rowMeans(logfcs.homol$mouse[,grepl("Beta FoxOi", colnames(logfcs.homol$mouse))]),
  gene= row.names(logfcs.homol$human))
ggplot(pDat, aes(x=human, y=mouse)) + geom_point(alpha=0.2) + 
  ggtitle(round(cor(pDat$mouse, pDat$human, method="spearman"),3))
ggsave(dirout(out, "Diff_human_mouse_FoxOi.pdf"), height=4, width=5)
ggsave(dirout(out, "Diff_human_mouse_FoxOi.jpg"), height=4, width=5)


# HM OF CORRELATIONS ------------------------------------------------------
# Beta
x <- cbind(logfcs.homol$human, logfcs.homol$mouse)
x <- x[,grepl("Beta Artemether", colnames(x)) | grepl("Beta FoxOi", colnames(x)) | grepl("Beta GABA", colnames(x))]
cleanDev()
# pdf(dirout(out, "Diff_human_mouse_Beta_HM_all.pdf"),width=5, height=5)
# pheatmap(cor(x, method="spearman"))
# dev.off()
x <- x[rownames(x) %in% diffGenes[celltype == "Beta" & treatment %in% c("Artemether", "FoxOi", "GABA") & qval < 0.05]$gene,]
cMT <- cor(x, method="spearman")
diag(cMT) <- NA
pdf(dirout(out, "Diff_human_mouse_Beta_HM_sig.pdf"),width=6, height=6)
pheatmap(cMT,clustering_distance_rows=dist(1-cMT), clustering_distance_cols=dist(1-cMT), 
         color=colorRampPalette(c("#edf8fb", "#8c96c6", "#810f7c"))(50))
dev.off()
pdf(dirout(out, "Diff_human_mouse_Beta_Clustering_sig.pdf"),width=6, height=6)
plot(hclust(dist(1-cMT)))
dev.off()



# Alpha
x <- cbind(logfcs.homol$human, logfcs.homol$mouse)
x <- x[,grepl("Alpha Artemether", colnames(x)) | grepl("Alpha FoxOi", colnames(x)) | grepl("Alpha GABA", colnames(x))]
cleanDev()
# pdf(dirout(out, "Diff_human_mouse_Alpha_HM_all.pdf"),width=5, height=5)
# pheatmap(cor(x, method="spearman"))
# dev.off()
x <- x[rownames(x) %in% diffGenes[celltype == "Alpha" & treatment %in% c("Artemether", "FoxOi", "GABA") & qval < 0.05]$gene,]
cMT <- cor(x, method="spearman")
diag(cMT) <- NA
pdf(dirout(out, "Diff_human_mouse_Alpha_HM_sig.pdf"),width=6, height=6)
pheatmap(cMT,clustering_distance_rows=dist(1-cMT), clustering_distance_cols=dist(1-cMT), 
         color=colorRampPalette(c("#edf8fb", "#8c96c6", "#810f7c"))(50))
dev.off()
pdf(dirout(out, "Diff_human_mouse_Alpha_Clustering_sig.pdf"),width=6, height=6)
plot(hclust(dist(1-cMT)))
dev.off()



# INS pos ALPHA cells -------------------------------------------------------------
ab.cor <- list()
(load(dirout("21_01_CellCorrelation/Markers/human.RData")))
ab.cor[["human"]] <- meta2[!sample %in% c("hIslets_I_pla2g16", "hIslets_II_GABA")]
(load(dirout("21_01_CellCorrelation/Markers/mouse.RData")))
ab.cor[["mouse"]] <- meta2
(load(dirout("21_01_CellCorrelation/Markers/human3.RData")))
ab.cor[["human3"]] <- meta2

orgx <- "human3"
for(orgx in c("human", "mouse", "human3")){
  pDat <- ab.cor[[orgx]]
  pDat <- pDat[celltype2 == "Alpha"]
  ins.gene <-  if(grepl("human", orgx)) "INS" else "Ins2"
  pDat$ins <-  c("INS-", "INS+")[(pDat[[ins.gene]] > 0) + 1]
  ggplot(pDat, aes(x=alphaCor, y=betaCor, color=ins)) + geom_point() +
    scale_color_manual(name="", values=c("#00000020", "#0000ff80")) + theme_bw(14) +
    xlab("Correlation to alpha cells") + ylab("Correlation to beta cells")
  ggsave(dirout(out, "INS_POS_ALPHA_", orgx, ".pdf"), width=5, height=4)
  ggsave(dirout(out, "INS_POS_ALPHA_", orgx, ".jpg"), width=5, height=4)
  
  ggplot(pDat, aes(x=alphaCor, color=ins)) + geom_density() +
    scale_color_manual(name="", values=c("#00000080", "#0000ff80")) + theme_bw(14) +
    xlab("Correlation to alpha cells") + ylab("Density")
  ggsave(dirout(out, "INS_POS_ALPHAcor_", orgx, ".pdf"), width=5, height=4)
  
  ggplot(pDat, aes(x=betaCor, color=ins)) + geom_density() +
    scale_color_manual(name="", values=c("#00000080", "#0000ff80")) + theme_bw(14) +
    xlab("Correlation to beta cells") + ylab("Density")
  ggsave(dirout(out, "INS_POS_BETAcor_", orgx, ".pdf"), width=5, height=4)
  
  if(orgx == "human3"){
    ggplot(pDat, aes(x=alphaCor, y=betaCor, color=ins)) + geom_point() + facet_wrap(~replicate) +
      scale_color_manual(name="", values=c("#00000020", "#0000ff80")) + theme_bw(14) +
      xlab("Correlation to alpha cells") + ylab("Correlation to beta cells")
    ggsave(dirout(out, "INS_POS_ALPHA_", orgx, "_byGroup.pdf"), width=15, height=8)
    ggsave(dirout(out, "INS_POS_ALPHA_", orgx, "_byGroup.jpg"), width=15, height=8)
  }
}

# Enrichment tests for INS+ -----------------------------
(load(dirout(out, "Meta.RData")))
meta$mouse$Ins1 <- dat.log$mouse["Ins1", meta$mouse$rn]
meta$mouse$INS <- NULL
res <- data.table()
for(orgx in names(meta)){
  for(insx in c("Ins1", "Ins2", "INS")){
    metaO <- meta[[orgx]]
    if(!insx %in% colnames(metaO)) next
    metaO[get(insx) > 0, insulin := "INS+"]
    metaO[!get(insx) > 0, insulin := "INS-"]
    for(treatx in c("Artemether", "FoxOi")){
      for(repx in c("I", "II", "III")){
        if(!repx %in% unique(metaO$replicate)) next
        m <- with(metaO[treatment %in% c("DMSO", treatx)][replicate == repx][celltype == "Alpha"], table(insulin, treatment))
        m <- m[c("INS+", "INS-"),c(treatx, "DMSO")]
        f <- fisher.test(m)
        res <- rbind(res, data.table(organism=orgx, treatment = treatx, replicate=repx, pvalue=f$p.value, oddsRatio=f$estimate, insulin=insx))
      }
      m <- with(metaO[treatment %in% c("DMSO", treatx)][celltype == "Alpha"], table(insulin, treatment))
      m <- m[c("INS+", "INS-"),c(treatx, "DMSO")]
      f <- fisher.test(m)
      res <- rbind(res, data.table(organism=orgx, treatment = treatx, replicate="All", pvalue=f$p.value, oddsRatio=f$estimate, insulin=insx))
    }
  }
}
res
write.tsv(res, dirout(out, "INS_Pos_FisherTest.tsv"))
(load(dirout(out, "Meta.RData")))


# INS+ Alpha cells - Alpha cell probability cutoff -----------------------------------------------------------------
#list.files(dirout("14_03_Classifier_moreCelltypes_noEndocrine/human/"))
orgx <- "human"
for(orgx in c("human", "mouse")){
  probs <- fread(dirout("14_03_Classifier_moreCelltypes_noEndocrine/",orgx,"/Prediction.tsv"))
  probs <- merge(meta[[orgx]], probs[,c("rn", "Alpha", "Beta"),with=F], by="rn")[celltype == "Alpha" | (celltype == "Endocrine" & PredictedCell == "Alpha")]
  probs[,Alpha2 := Alpha]
  probs[celltype == "Alpha", Alpha2 := 1]
  alpha2.co <- 0.5
  res <- data.table()
  for(alpha2.co in unique(sort(floor(probs$Alpha2*100)/100))){
    res <- rbind(res, data.table(probs[Alpha2 >= alpha2.co, .(INSpos = sum(INS > 0), INSneg = sum(INS <= 0)), by=c("treatment", "replicate")], cutoff=alpha2.co))
  }
  res[,percentPos := INSpos/(INSpos + INSneg) * 100]
  ggplot(res, aes(x=cutoff, y=percentPos, color=treatment)) + facet_grid(. ~ replicate) + 
    geom_line() + theme_bw(16) + xlab("Alpha cell probability") + ylab("Percent of INS+ cells")
  ggsave(dirout(out, "INSpos_",orgx,"_AlphaProbCutoff.pdf"), width=length(unique(res$replicate)) * 4 + 1, height=4)
  write.tsv(res, dirout(out, "INSpos_",orgx,"_AlphaProbCutoff.tsv"))
  
  ggplot(probs, aes(x=Alpha2, color=treatment)) + stat_ecdf() + facet_grid(. ~ replicate) + theme_bw(16) + xlab("Alpha cell probability") + ylab("ECDF")
  ggsave(dirout(out, "INSpos_",orgx,"_AlphaProbDistr_ECDF.pdf"), width=length(unique(res$replicate)) * 4 + 1, height=4)
  
  ggplot(probs, aes(x=Alpha2, color=treatment)) + geom_density() + facet_grid(. ~ replicate) + theme_bw(16) + xlab("Alpha cell probability")
  ggsave(dirout(out, "INSpos_",orgx,"_AlphaProbDistr_Density.pdf"), width=length(unique(res$replicate)) * 4 + 1, height=4)
}



# INS Differential expression ---------------------------------------------
# without replicates
res <- data.table()
meta.de <- meta
meta.de$human3_36 <- meta.de$human3[grepl("_36", replicate)]
meta.de$human3_72 <- meta.de$human3[grepl("_72", replicate)]
meta.de$human3 <- NULL
for(orgx in names(meta.de)){
  orgx.dat <- if(grepl("human3_", orgx)) "human3" else orgx
  x <- CreateSeuratObject(raw.data=dat.raw[[orgx.dat]], min.cells=0, min.genes=0)
  x@data <- dat.log[[orgx.dat]][row.names(x@raw.data), colnames(x@raw.data)]
  for(ctx in c("Alpha", "Beta")){
    mx <- meta.de[[orgx]][celltype == ctx]
    for(treatx in c("FoxOi", "Artemether")){
      if(treatx == "FoxOi" & grepl("human3_", orgx)) next
      bar.dmso <- mx[treatment == "DMSO"]$rn  
      bar.treat <- mx[treatment == treatx]$rn
      genes.to.use <- if(grepl("human", orgx)) c("INS", "GCG", "SST", "PPY") else c("Ins1", "Ins2", "Gcg", "Sst", "Ppy")
      stopifnot(all(colnames(x@raw.data) == colnames(x@data)))
      stopifnot(all(row.names(x@raw.data) == row.names(x@data)))
      stopifnot(all(bar.dmso %in% colnames(x@raw.data)))
      stopifnot(all(bar.treat %in% colnames(x@raw.data)))
      stopifnot(all(genes.to.use %in% row.names(x@raw.data)))
      dtest <- data.table(FindMarkers(object = x, min.pct=0, ident.1=bar.treat, ident.2 = bar.dmso, genes.use = genes.to.use, thresh.use=0, test.use="negbinom"), keep.rownames=T)
      #p <- data.table(Seurat:::NegBinomDETest(object = x, cells.1 = bar.treat, cells.2 = bar.dmso, genes.use = genes.to.use, latent.vars = "nUMI", print.bar = TRUE, min.cells = 10), keep.rownames=T)
      res <- rbind(res, data.table(dtest, organism=orgx, celltype=ctx, treatment=treatx))
    }
  }
}
write.tsv(res, dirout(out, "DE_NoReplicates_Insulin.tsv"))

# with replicates
resx <- data.table()
for(x in c("H", "H3", "M")){
  (load(dirout("16_",x,"_02_DiffGenes/VsDMSO_negbinom_onlyIns/AllRes.RData")))
  res <- res[,1:6]
  colnames(res) <- paste0("V",1:6)
  resx <- rbind(resx, data.table(res, dataset=x), fill=T)
}
write.tsv(resx[V1 %in% c("INS", "Ins1", "Ins2")][(grepl("Fox", V6) | grepl("A10", V6))][grepl("Alpha", V6) | grepl("Beta", V6)], dirout(out, "DE_Replicates_Insulin.tsv"))
write.tsv(resx[(grepl("Fox", V6) | grepl("A10", V6))][grepl("Alpha", V6) | grepl("Beta", V6)], dirout(out, "DE_Replicates_Hormones.tsv"))

# FoxOi Literature ----------------------------------------------------------
foxoList <- fread("metadata/FoxO_Genes_Literature.csv")
foxoList <- merge(foxoList, homol, by.x="gene", by.y="mouse")
colnames(foxoList) <- gsub(" \\(.+$", "", colnames(foxoList))
foxoList <- foxoList[human %in% row.names(avg.homol$human)]
foxoList.cor <- cor(cbind(logfcs.homol$human[foxoList$human,], logfcs.homol$mouse[foxoList$human,]), foxoList$LogFC, method="spearman")
pDat <- data.table(Correlation=foxoList.cor[,1], do.call(rbind, strsplit(row.names(foxoList.cor), " ")))
ggplot(pDat, aes(x=V3, y=Correlation, color=V4)) + 
  facet_grid(V2 ~ V1, scales="free_x") + 
  geom_jitter(height=0, width=0.1,size=2, alpha=0.7)+ 
  theme_bw(16) + xRot()
ggsave(dirout(out, "FoxoLiterature_LogFC_Correlation.pdf"), width=5, height=6)



# Enrichments -------------------------------------------------------------
diffGenes.cnt <- diffGenes[celltype %in% c("Beta", "Alpha") & treatment %in% c("Artemether", "FoxOi")][,.N, by=c("org", "treatment", "celltype", "direction", "gene")]
diffGenes.cnt <- diffGenes.cnt[(org == "mouse" & N == 3) | (org == "human" & N == 2)]
diffGenes.cnt <- dcast.data.table(diffGenes.cnt, org + celltype + gene ~ treatment,value.var="direction")
diffGenes.cnt[,Artemether_FoxOi := gsub("NA", "x", paste0(Artemether, "_", FoxOi))]
enr.lists <- split(diffGenes.cnt$gene, factor(with(diffGenes.cnt, paste(celltype, Artemether_FoxOi, org))))
enr.file <- dirout(out, "Enr.RData")
source("src/FUNC_Enrichr.R")
if(!file.exists(enr.file)){
  enr <- enrichGeneList.oddsRatio.list(geneLists=enr.lists, enrichrDBs=c("WikiPathways_2016", "NCI-Nature_2016","Reactome_2016","KEGG_2016"))
  save(enr, diffGenes.cnt, file=enr.file)
} else {
  load(enr.file)
}
enr2 <- cbind(enr, do.call(rbind, strsplit(enr$grp, " ")))
enrichr.plot.many(enrichRes=enr2, out=dirout(out),label="Enrichr")



# Contamination percent in spike-ins --------------------------------------
(load(dirout("07_03_CleanHuman_CleanSpikeIns/","Uncorrected.Data.RData")))
#"ml.rho.pred"         "ml.rho.cv"           "ml.rho.pred.with.cv" "data.raw.cells"      "meta.raw.cells"      "spikeIn.soups" 
(load(dirout("07_03_CleanHuman_CleanSpikeIns/","CorrectedData.RData")))
# "corr.export"
metaH <- meta$human
data.raw.log <- data.raw.cells[row.names(dat.log$human),colnames(dat.log$human)]
data.raw.log <- toTPM_fromRaw(data.raw.log, scale.factor=SCALE.FACTOR)
data.raw.log <- toLogTPM(data.raw.log)
metaH$INS_Raw <- data.raw.log["INS", metaH$rn]


# RAW INSULIN IN SPIKE-INS ------------------------------------------------
ggplot(metaH[celltype2 == "SI_mouse" & treatment %in% c("Artemether", "DMSO")], aes(x=treatment, y=INS_Raw, color=treatment)) + 
  geom_boxplot(coef=Inf) + geom_jitter(height=0, shape=1) + 
  facet_grid(. ~ replicate) + theme_bw(14) +
  xlab("") + ylab("INS Log(TPM)") + guides(color=F)
ggsave(dirout(out, "INS_RAW_SPIKE_INS.pdf"), width=4, height=4)

ggplot(metaH[celltype2 == "SI_mouse" & treatment %in% c("Artemether", "DMSO")], aes(x=treatment, y=exp(INS_Raw)-1, color=treatment)) + 
  geom_boxplot(coef=Inf) + geom_jitter(height=0, shape=1) + 
  facet_grid(. ~ replicate) + theme_bw(14) +
  xlab("") + ylab("INS TPM") + guides(color=F)
ggsave(dirout(out, "INS_RAW_SPIKE_INS_TPM.pdf"), width=4, height=4)
