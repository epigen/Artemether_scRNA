require("project.init")
require(Seurat)
require(Matrix)
require(methods)
project.init2("artemether")

inDir <- "20_01_AlphaCells/"
out <- "20_02_Signatures/"
dir.create(dirout(out))

metaH <- loadHMeta2()
metaM <- loadMMeta2()

pbmcM <- loadMData2()
pbmcH <- loadHData2()


load(dirout(inDir, "Markers.RData"))

# Calculate marker based scores ----------------------------------------------------------------
scores.file <- dirout(out, "Scores.RData")
if(!file.exists(scores.file)){
  org <- "mouse"
  for(org in c("mouse", "human")){
    pbmc <- pbmcM
    meta <- metaM
    if(org == "human"){
      pbmc <- pbmcH
      meta <- metaH
    }
    mark <- markers[V2 == org & myAUC > 0.75]
    ct <- "Alpha"
    for(ct in unique(mark$V1.1)){
      x <- as.matrix(pbmc@data[mark[V1.1 == ct]$V1,meta$rn])
      dr <- sapply(split(meta$rn, factor(meta$celltype2)),function(rns) rowMeans(x[,rns, drop=F]))
      dr <- dr[,-which(colnames(dr) == ct)] - dr[,ct]
      x <- x[names(which(
        apply(dr < -log(2), 1, sum) >= ncol(dr) * 1 &
          apply(dr < 0, 1, sum) == ncol(dr)
        )),,drop=F]
      x <- x - rowMin(x)
      x <- x / rowMax(x)
      meta[[paste0(ct, "_score")]] <- colMeans(x)
      h <- sapply(split(meta$rn, factor(with(meta, paste(celltype2, treatment, replicate, sep="_")))),function(rns) rowMeans(x[,rns, drop=F]))
      pdf(dirout(out, "GeneAverages_", org, "_", ct, ".pdf"), width=25, height=8, onefile=F)
      pheatmap(h, cluster_cols=F)
      dev.off()
    }
    if(org == "mouse") metaM <- meta
    if(org == "human") metaH <- meta
  }
  save(metaH, metaM, file=scores.file)
} else { load(scores.file)}

# Plot scores ----------------------------------------------------------------
org <- "human"
score <- "Alpha_score"
for(org in c("mouse", "human")){
  meta <- if(org == "human") metaH  else metaM
  meta$treatment <- factor(meta$treatment, levels = c("DMSO", sort(unique(meta[treatment != "DMSO"]$treatment))))
  for(score in c("Beta_score", "Alpha_score")){
#     for(treat in unique(meta$treatment)){
#       ggplot(meta[treatment %in% c(treat, "DMSO") & celltype2 %in% c("Alpha", "Beta")], aes_string(x=score, color="treatment")) + theme_bw(12) + 
#         stat_ecdf() + facet_grid(replicate ~ celltype2)+ scale_color_manual(values=c("black", "red"))
#       ggsave(dirout(out, "Density_", org, "_", treat, "_", score, ".pdf"), width=8, height=12)
#     }
#     meta[celltype2 %in% c("Alpha", "Beta")][,score := mean(get(score)), by=c("treatment", "replicate", "celltype2")]
#     ggplot(meta[celltype2 %in% c("Alpha", "Beta")], aes_string(y=score, x="treatment")) + theme_bw(12) + 
#       geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) + facet_grid(replicate ~ celltype2)
#     ggsave(dirout(out, "Violin_", org, "_", "_", score, ".pdf"), width=8, height=12)
#     
    ggplot(meta[celltype == "Endocrine" & celltype2 %in% c("Alpha", "Beta")], aes_string(y=score, x="treatment")) + theme_bw(12) + 
      geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) + facet_grid(replicate ~ celltype2)
    ggsave(dirout(out, "ViolinEndocrine_", org, "_", "_", score, ".pdf"), width=8, height=12)
  }
#   ggplot(meta[celltype2 %in% c("Alpha", "Beta", "Gamma", "Delta") & treatment %in% c("DMSO", "A10")],
#          aes(x=Alpha_score, y=Beta_score, color=treatment)) + stat_ellipse() + facet_grid(replicate ~ .)
#   ggsave(dirout(out, "Ellipses_", org, ".pdf"), width=4, height=8)
#   
#   ggplot(meta[celltype2 %in% c("Alpha", "Beta", "Gamma", "Delta") & treatment %in% c("DMSO", "A10")],
#          aes(x=Alpha_score, y=Beta_score,color=treatment)) + geom_point(alpha=0.1) + facet_grid(replicate ~ .) +
#     scale_color_manual(values=c("blue", "red"))
#   ggsave(dirout(out, "Points_", org, ".pdf"), width=4, height=8)
#   
#   ggplot(meta[celltype2 %in% c("Alpha", "Beta", "Gamma", "Delta") & treatment %in% c("DMSO", "A10")],
#          aes(x=Alpha_score, y=Beta_score)) + geom_hex() + facet_grid(replicate ~ treatment) +
#     scale_color_manual(values=c("blue", "red"))
#   ggsave(dirout(out, "Hex_", org, ".pdf"), width=8, height=8)
}


