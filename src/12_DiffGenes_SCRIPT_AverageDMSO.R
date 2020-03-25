# Diff vs DMSO ------------------------------------------------------------
if(!exists("outTreat")) stop("MISSING outTreat")
if(!exists("pbmc")) stop("MISSING pbmc")
if(!exists("meta")) stop("MISSING meta")
if(!exists("nrCores")) nrCores <- 5
if(!exists("seurat.thresh.use")) seurat.thresh.use <- 0.25
if(!exists("seurat.genes.use")) seurat.genes.use <- row.names(pbmc@data)

dir.create(dirout(outTreat))

treat <- "FoxO"
cl <- "Beta"
repl <- "I"
meta <- meta[match(row.names(pbmc@meta.data), meta$rn)]
stopifnot(all(meta$rn == row.names(pbmc@meta.data)))
require(doMC); registerDoMC(cores=nrCores)

for(repl in unique(meta$replicate)){
  meta.rpl <- meta[replicate == repl]
  for(treat in unique(meta.rpl[treatment != "DMSO"]$treatment)){
    ctxx <- unique(meta.rpl[treatment == treat][!is.na(celltype)]$celltype)
    foreach(cl = ctxx) %dopar% {
      message(paste(repl, treat, cl, sep="_"))
      (ff <- dirout(outTreat, repl, "_", treat, "_", cl, ".tsv"))
      print(ff)
      if(!file.exists(ff)){
        print('starting analysis')
        bar.treat <- meta[replicate == repl & treatment == treat & celltype == cl]$rn
        bar.dmso <- meta[replicate == gsub("III\\d", "III", repl) & treatment == "DMSO" & celltype == cl]$rn
        #print(paste(length(bar.treat), " ", treat, " ", length(bar.dmso), " DMSO"))
        if(length(bar.treat) > 10 & length(bar.dmso) > 10){
          cluster.markers <- FindMarkers(pbmc,
                                         thresh.use=seurat.thresh.use,
                                         genes.use=seurat.genes.use,
                                         test.use="negbinom",
                                         ident.1 = bar.treat,
                                         ident.2 = bar.dmso)
          write.table(cluster.markers, ff, sep="\t", quote=F, row.names=TRUE)
        }
      }
    }
  }
}
source("src/FUNC_Enrichr.R")

# AGGREGATE
ff <- list.files(dirout(outTreat), pattern="I_.+_.+\\.tsv$")
res <- data.table(); for(f in ff){print(f);res <- rbind(res, data.table(fread(dirout(outTreat, f)), file=f), fill=T)}
res[, file := gsub("SI_", "SI.", file)]
res[, file := gsub("Acinar_like", "Acinar.like", file)]
res[, file := gsub("\\.tsv$", "", file)]
res <- cbind(res, data.table(do.call(rbind, strsplit(res$file, "_"))))
colnames(res) <- make.unique(colnames(res))
res$qval <- p.adjust(res$V2, method="BH")
res[, direction := ifelse(V3 > 0, "up", "down")]
save(res, file=dirout(outTreat, "AllRes.RData"))

# INDIVIDUAL LISTS
# res0 <- res[qval < 0.05 & abs(V3) > log(3)]
# ggplot(res0, aes(x=V3.1, fill=direction)) + geom_bar(position="dodge") + facet_grid(V1.1 ~ V2.1) + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)); ggsave(dirout(outTreat, "0_Numbers.pdf"), width=20, height=10)
# write.tsv(res0, dirout(outTreat, "0_IndividualHits.tsv"))
# # analyze individual lists
# ll <- with(res0, split(V1, factor(paste0(file, "_", direction))))
# # enr <- enrichGeneList.oddsRatio.list(ll, enrichrDBs=ENRICHR.DBS)
# # enrichr.plot.many(enrichRes=enr, out=dirout(outTreat), label="0_Enrichr_vsDMSO")
# pdf(dirout(outTreat, "0_Jaccard.pdf"), height=29, width=29,onefile=F)
# pheatmap(jaccard(ll),cluster_rows=T, cluster_cols=T,color=colorRampPalette(c("grey", "blue"))(50))
# dev.off()
#
# # AGGREGATED LISTS ACROSS REPLICATES
# res.agg <- res0[, .(count = .N, logFC = mean(V3)), by=c("V1", "V2.1", "V3.1", "direction")][count == length(unique(res0$V1.1))]
# ggplot(res.agg, aes(x=V3.1, fill=direction)) + geom_bar(position="dodge") + facet_grid(. ~ V2.1) + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)); ggsave(dirout(outTreat, "2_Agg_Numbers.pdf"), width=20, height=5)
# write.tsv(res.agg, dirout(outTreat, "2_Agg_Genes.tsv"))
# res.agg <- res.agg[!grepl("SI\\.", V3.1)]
# ll2 <- with(res.agg, split(V1, factor(paste(V2.1, V3.1, direction, sep="_"))))
# pdf(dirout(outTreat, "2_Agg_Jaccard.pdf"), height=6, width=6,onefile=F)
# pheatmap(jaccard(ll2),cluster_rows=T, cluster_cols=T,color=colorRampPalette(c("grey", "blue"))(50))
# dev.off()
# enr <- enrichGeneList.oddsRatio.list(ll2, enrichrDBs=ENRICHR.DBS)
# enrichr.plot.many(enrichRes=enr, out=dirout(outTreat), label="2_Enrichr_vsDMSO")
# res.agg[,x := paste0(V2.1, V3.1)]
# res.agg <- hierarch.ordering(res.agg, toOrder="V1", orderBy="x", value.var="logFC")
# ggplot(res.agg[!grepl("SI\\.", V3.1)], aes(x=V3.1, y=V1, color=logFC)) + geom_point() + facet_grid(. ~ V2.1) + theme_bw(12) +
#   scale_color_gradient2(low="blue", high="red") + xRot()
# ggsave(dirout(outTreat, "2_Agg_Genes.pdf"), height=20, width=7)
#
#
#
#
# # AGGREGATE WITH LESS STRICTNESS ------------------------------------------
# res.agg <- res[qval < 0.05, .(count = .N, logFC = mean(V3)), by=c("V1", "V2.1", "V3.1", "direction")][count == length(unique(res0$V1.1))]
# ggplot(res.agg, aes(x=V3.1, fill=direction)) + geom_bar(position="dodge") + facet_grid(. ~ V2.1) + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)); ggsave(dirout(outTreat, "3_Agg_Numbers.pdf"), width=20, height=5)
# write.tsv(res.agg, dirout(outTreat, "3_Agg_Genes.tsv"))
# res.agg <- res.agg[!grepl("SI\\.", V3.1)]
# ll2 <- with(res.agg, split(V1, factor(paste(V2.1, V3.1, direction, sep="_"))))
# pdf(dirout(outTreat, "3_Agg_Jaccard.pdf"), height=15, width=15,onefile=F)
# pheatmap(jaccard(ll2),cluster_rows=T, cluster_cols=T,color=colorRampPalette(c("grey", "blue"))(50))
# dev.off()
# enr <- enrichGeneList.oddsRatio.list(ll2, enrichrDBs=ENRICHR.DBS)
# enrichr.plot.many(enrichRes=enr, out=dirout(outTreat), label="3_Agg_Enrichr_vsDMSO")
#
# res.agg[,grp := paste(V2.1, V3.1, direction, sep="_")]
# x <- toMT(dt=res.agg,row="V1", col="grp",val="logFC")
# x[is.na(x)] <- 0
# pdf(dirout(outTreat, "3_Agg_Overlaps.pdf"), height=25, width=5, onefile=F)
# pheatmap(x[rowSums(x != 0) > 1,], fontsize_row=3, fontsize_col=5, color=colorRampPalette(c("blue", "white","red"))(50))
# dev.off()
