require("project.init")
project.init2("artemether")

source("src/FUNC_Enrichr.R")

# min 2 replicates

out <- "17_02_06_Gene_plots_q0.1_ByGroup/"
dir.create(dirout(out))

res.file <- dirout(out, "Diff_Genes.tsv")
if(!file.exists(res.file)){
  # Load data
  (load(dirout("16_H_02_DiffGenes/VsDMSO_negbinom_noIns/AllRes.RData")))
  resH <- res
  
  (load(dirout("16_M_02_DiffGenes/VsDMSO_negbinom_noIns/AllRes.RData")))
  resM <- res
  
  resH <- resH[V1.1 != "II" | V2 != "GABA"]
  resH <- resH[V2 != "pla2g16"]
  
  resH[,qval := p.adjust(p_val, method="BH")]
  resM[,qval := p.adjust(p_val, method="BH")]
  
  # Combine results
  resH[,org:="human"]
  resM[,org:="mouse"]
  resX <- rbind(resH, resM)
  write.tsv(resX, dirout(out, "AllResults.tsv"))
  res.agg <- resX[qval < 0.1]
  res.agg <- res.agg[, count := .N, by=c("V1", "V2", "V3", "direction", "org")]
  res.agg <- res.agg[, logFC := mean(avg_diff), by=c("V1", "V2", "V3", "direction", "org")]
  res.agg[,totalCnt := length(unique(V1.1)), by=c("V2", "org")]
  res.agg[,.(totalCnt = length(unique(V1.1))), by=c("V2", "org")]
  res.agg <- res.agg[count >= totalCnt]
  res <- res.agg[, .(count = .N, logFC = mean(avg_diff)), by=c("V1", "V2", "V3", "direction", "org")]
  table(res$count, paste(res$V2, res$org))
  write.tsv(res, res.file)
} else {
  res <- fread(res.file)
}

# Plot numbers
ggplot(res, aes(x=V3, fill=direction)) + geom_bar(position="dodge") + theme_bw(12) + 
  facet_grid(org ~ V2) + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
ggsave(dirout(out, "Numbers.pdf"), width=20, height=10)



# MAP TO HOMOLOGS ---------------------------------------------------------
source("src/FUNCTIONS_Synonyms_etc.R")
# homologs.file <- dirout(out, "Homologs.tsv")
# if(!file.exists(homologs.file)){
#   source("src/FUNCTIONS_Synonyms_etc.R")
#   resHomologs <- res
#   resHomologs <- merge(resHomologs, homol, by.x="V1", by.y="mouse", all.x=T)
#   resHomologs <- merge(resHomologs, homol, by.x="V1", by.y="human", all.x=T)
#   resHomologs[org == "human", human := V1]
#   resHomologs[org == "mouse", mouse := V1]
#   write.tsv(resHomologs, file=homologs.file)
# } else {
#   resHomologs <- fread(homologs.file)
# }
# resHomologs <- resHomologs[!is.na(human) & !is.na(mouse)]
# resHomologs[,id := paste(V2, V3, org, direction)]
# x <- resHomologs[((V2 == "FoxO" & direction == "down") | 
#                     (V2 == "A10" & org == "human") | 
#                     (V2 == "A10" & org == "mouse" & direction == "down")) & V3 == "Beta"]
# gplots::venn(split(x$human, factor(x$id)))


# # Enrichments
# ll2 <- with(res, split(V1, factor(paste(V2, V3, direction, org, sep="_"))))
# enr <- enrichGeneList.oddsRatio.list(ll2, enrichrDBs=ENRICHR.DBS)
# enrichr.plot.many(enrichRes=enr, out=dirout(out), label="Enrichr_")
# enr <- fread(dirout(out, "Enrichr_.tsv"))
# 
# # ENRICHR OVERLAPS --------------------------------------------------------
# enr <- cbind(enr, do.call(rbind, strsplit(enr$grp, "_")))
# enr2 <- enr[V1 %in% c("A10", "FoxO") & V2 %in% c("Alpha", "Beta") & hitLength > 2]
# enrichr.plot.many(enrichRes=enr2, out=dirout(out), label="EnrFocus_")
# enr2 <- enr[V1 %in% c("A10") & V2 %in% c("Alpha", "Beta") & V3 == "up" & hitLength > 2]
# enrichr.plot.many(enrichRes=enr2, out=dirout(out), label="EnrFocus2_")
# 
# 
# 
# # MOUSE VS HUMAN ----------------------------------------------------------
# enrX <- enr[database %in% c("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "TRANSFAC_and_JASPAR_PWMs")]
# enrB <- enrX[V3 == "up" & V1 == "A10" & V2 == "Beta"]
# enrB <- enrB[category %in% enrB[,.N,by="category"][N == 1]$category]
# enrA <- enrX[V3 == "up" & V1 == "A10" & V2 == "Alpha"]
# enrA <- enrA[category %in% enrA[,.N,by="category"][N == 1]$category]
# enrC <- rbind(enrA, enrB)
# enrCC <- enr[V3 == "up" & V1 == "A10" & V2 %in% c("Beta", "Alpha") & category %in% enrC$category]
# enrCC <- hierarch.ordering(dt=enrCC,toOrder="category",orderBy="grp",value.var="oddsRatio")
# ggplot(enrCC, aes(x=category, y=V4, color=oddsRatio, size=pmin(5, -log10(qval)))) + 
#   geom_point()  + 
#   geom_point(data=enrCC[qval < 0.05], color="black", size=1) +
#   theme_bw(12) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   facet_grid(V2 ~ database, space="free_x", scales="free_x") +
#   scale_color_gradient(low= "grey", high="red")
# ggsave(dirout(out,"MvsH.pdf"), width=14, height=4)
# 
# write.tsv(enrCC, dirout(out, "MvsH.tsv"))


# Gene HITS overlaps - JACCARD AND MDS ------------------------------------------------------
source("src/FUNCTIONS_Synonyms_etc.R")
hits <- copy(res)
hits[,grp := paste(V2, V3, direction, org, sep="_")]
hits[org == "human",gene := V1]
hits[org == "mouse"]$gene <- homol[match(hits[org == "mouse"]$V1, homol$mouse)]$human
write.tsv(hits, dirout(out, "Hits_Mapped.tsv"))
hits <- hits[!is.na(gene)]

x <- toMT(dt=hits,row="gene", col="grp", val="logFC")
x[is.na(x)] <- 0
write.table(x, dirout(out, "Hits_Overlaps.tsv"), sep="\t", row.names=TRUE, quote=FALSE)

ll <- split(hits$gene, factor(hits$grp))
cleanDev()
pdf(dirout(out, "Hits_Jaccard.pdf"), height=15, width=15,onefile=F)
pheatmap(jaccard(ll),cluster_rows=T, cluster_cols=T,color=colorRampPalette(c("grey", "blue"))(50))
dev.off()
# 
# xx <- list(human=hits[org=="human"], mouse=hits[org=="mouse"], comb=hits)
# for(nam in names(xx)){
#   hitsX <- xx[[nam]]
#   MDS.scale=cmdscale(as.dist(1-jaccard(split(hitsX$gene, factor(hitsX$grp)))), eig = TRUE, k = 2)
#   MDS = data.table(sample_name =rownames(MDS.scale$points), MDS1 = MDS.scale$points[,1], MDS2 = MDS.scale$points[, 2])
#   MDS <- cbind(MDS, do.call(rbind, strsplit(MDS$sample_name, "_")))
#   hitsX <- xx[[nam]]
#   MDS.scale=cmdscale(as.dist(1-jaccard(split(hitsX$gene, factor(hitsX$grp)))), eig = TRUE, k = 2)
#   MDS = data.table(sample_name =rownames(MDS.scale$points), MDS1 = MDS.scale$points[,1], MDS2 = MDS.scale$points[, 2])
#   MDS <- cbind(MDS, do.call(rbind, strsplit(MDS$sample_name, "_")))
#   ggplot(MDS, aes(x=MDS1, y=MDS2)) + 
#     geom_point(aes(color=V1), alpha=0.5) +
#     geom_text(aes(label=V2,color=V1), data=MDS[V2 %in% c("Alpha", "Beta")]) + 
#     theme_bw(12) + facet_grid(V3 ~ V4) +
#     scale_color_manual(values=COLORS.TREATMENTS)
#   ggsave(dirout(out, "Hits_MDS_",nam,".pdf"), width=length(unique(MDS$V4)) * 4 + 1, height=9)   
# }



# Correlation in logfold changes ------------------------------------------
hm.col <- colorRampPalette(c("lightgrey", "khaki", "darkorchid4"))(51)
orgx <- "human"
for(orgx in c("human", "mouse")){
  resX <- fread(dirout(out, "AllResults.tsv"))
  resX[,id := gsub("I+_", "", file)]
  resX.H <- resX[V3 %in% c("Alpha", "Beta")][org == orgx][,.(logFC = mean(avg_diff)), by=c("V1", "id")]
  resX.H.MT <- toMT(resX.H, row="V1", col="id", val="logFC")
  #ggplot(data.table(resX.H.MT), aes(x=A10_Beta, y=FoxO_Beta)) + geom_point()
  corMT <- cor(resX.H.MT, use="pairwise.complete.obs", method="spearman")
  diag(corMT) <- NA
  cleanDev()
  pdf(dirout(out, "HM_Similarity_",orgx,".pdf"), width=6, height=6, onefile=F)
  pheatmap(corMT, color=hm.col, na_col="lightgrey")#, breaks=seq(0,1,0.02))
  dev.off()
}