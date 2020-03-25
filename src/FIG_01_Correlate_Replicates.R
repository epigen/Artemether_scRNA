require("project.init")
project.init2("artemether")

out <- "FIG_01_Correlate_Replicates/"
dir.create(dirout(out))

orgx <- "human"
cleanTreatment <- function(vec){ gsub("^FoxO$", "FoxOi", gsub("^A10$", "Artemether", vec)) }


# Get average expression --------------------------------------------------
metaH <- if(orgx == "mouse") loadMMetaUnfiltered() else loadHMetaUnfiltered()

gene.means.raw.file <- dirout(out, "gene.means.raw.RData")
if(!file.exists(gene.means.raw.file)){
  if(orgx == "mouse") (load(dirout("07_11_CleanMouse_CleanSpikeIns/","Uncorrected.Data.RData"))) else (load(dirout("07_03_CleanHuman_CleanSpikeIns/","Uncorrected.Data.RData")))
  str(data.raw.cells)
  data.raw.cells <- toTPM_fromRaw(data.raw.cells, scale.factor=SCALE.FACTOR)
  
  gene.means.raw <- sapply(with(metaH, split(rn, factor(paste(replicate, treatment, celltype2)))), function(rns) Matrix::rowMeans(data.raw.cells[,rns,drop=F]))
  str(gene.means.raw)
  rm(list=c("data.raw.cells", "meta.raw.cells", "spikeIn.soups", "ml.rho.pred", "ml.rho.pred.with.cv"))
  save(gene.means.raw, file=gene.means.raw.file)
} else {
  load(gene.means.raw.file)
}

gene.means.cor.file <- dirout(out, "gene.means.cor.RData")
if(!file.exists(gene.means.cor.file)){
  pbmcH <- if(orgx == "mouse") loadMData() else loadHData()
  
  data.cor.cells <- toTPM(pbmcH@data)
  rm(list="pbmcH")
  gene.means.cor <- sapply(with(metaH, split(rn, factor(paste(replicate, treatment, celltype2)))), function(rns) Matrix::rowMeans(data.cor.cells[,rns,drop=F]))
  str(gene.means.cor)
  save(gene.means.cor, file=gene.means.cor.file)
} else {
  load(gene.means.cor.file)
}



# Correlate cell types before after ---------------------------------------
rep <- "I"
res <- data.table()
for(rep in c("I", "II")){
  treatx <- "DMSO"
  for(treatx in c("DMSO", "A10", "GABA", "FoxO")){
    g1 <- paste(rep, treatx, "Beta", sep=" ")
    g2 <- paste(rep, treatx, "Alpha", sep=" ")
    res <- rbind(res, data.table(replicate=rep, treatment=treatx, 
                                 corrected=cor(gene.means.cor[,g1], gene.means.cor[,g2]),
                                 raw=cor(gene.means.raw[,g1], gene.means.raw[,g2])))
  }
}
res$treatment <- cleanTreatment(res$treatment)
pDat <- melt(res)
pDat$data.type <- factor(pDat$variable, levels=c("raw", "corrected"))
(p <- ggplot(pDat[!(treatment == "GABA" & replicate == "II")], aes(x=data.type,y=value, color=replicate, group=replicate)) + geom_point() + geom_line() + facet_wrap(~treatment, ncol=2) +
  theme_bw(16) + xRot() +
  ylab("Correlation of\nalpha and beta cells") + xlab("Data"))
ggsave(dirout(out, "Alpha_Beta_Correlation_guide.pdf"), w=4,h=4)
ggsave(dirout(out, "Alpha_Beta_Correlation.pdf"), w=4,h=5, plot=p+guides(color=F))


# Prepare analysis --------------------------------------------------------
# gg <- names(tail(sort(apply(gene.means.cor, 1, min)), 5000))
# #gg <- rownames(gene.means.cor)
# xx <- list(corrected=gene.means.cor[gg,], raw=gene.means.raw[gg,])
# #xx <- list(corrected=gene.means.cor, raw=gene.means.raw)


# # Correlate replicates ----------------------------------------------------
# res <- data.table()
# for(xnam in names(xx)){
#   x <- xx[[xnam]]
#   grps <- unique(gsub("^I+ ", "X ", colnames(x)))
#   i <- grps[1]
#   for(i in grps){
#     repl1 <- gsub("^X ", "I ", i)
#     repl2 <- gsub("^X ", "II ", i)
#     if(repl1 %in% colnames(x) & repl2 %in% colnames(x))
#     res <- rbind(res, data.table(grp=gsub("^X ", "", i), type=xnam, value=cor(x[,repl1], x[,repl2], method="spearman")))
#   }
# }
# 
# res <- cbind(res, do.call(rbind, strsplit(res$grp, " ")))
# res[,treatment := cleanTreatment(V1)]
# ct.use <- c("Beta", "Alpha") #, "Delta", "Gamma")
# ggplot(res[V2 %in% ct.use & V1 != "GABA"], aes(x=type, y=value, color=V2, group=paste(V2, V1))) + 
#   #geom_jitter(height=0, width=0.1) + 
#   geom_line() + geom_point() +
#   facet_grid(. ~ treatment) + theme_bw(12) + xRot()
# ggsave(dirout(out, "Expr_Alpha_Beta_",orgx,".pdf"),w=6,h=4)
# 
# res.wide <- dcast.data.table(res, treatment + V2 ~ type, value.var="value")
# res.wide[,diff := corrected - raw]
# ggplot(res.wide, aes(x=treatment, y=V2, fill=diff)) + geom_tile() + scale_fill_gradient2() + xRot()
# ggsave(dirout(out, "Expr_AllCelltypes_",orgx,".pdf"),w=4,h=6)
# 
# 
# 
# # Correlate log fold changes ----------------------------------------------
# res <- data.table()
# xnam <- names(xx)[1]
# for(xnam in names(xx)){
#   x <- xx[[xnam]]
#   treatments <- unique(gsub(".+ (.+) .+", "\\1", colnames(x)))
#   treatments <- treatments[treatments != "DMSO"]
#   treatx <- treatments[1]
#   for(treatx in treatments[treatments != "DMSO"]){
#     x.treat <- x[,grepl(treatx, colnames(x))]
#     grps <- unique(gsub("^I+ ", "X ", colnames(x.treat)))
#     i <- grps[1]
#     for(i in grps){
#       repl1 <- gsub("^X ", "I ", i)
#       dmso1 <- gsub("^X ", "I ", gsub(treatx, "DMSO", i))
#       repl2 <- gsub("^X ", "II ", i)
#       dmso2 <- gsub("^X ", "II ", gsub(treatx, "DMSO", i))
#       if(sum(c(repl1, repl2, dmso1, dmso2) %in% colnames(x)) == 4){
#         fc1 <- x[,repl1] - x[,dmso1]
#         fc2 <- x[,repl2] - x[,dmso2]
#         res <- rbind(res, data.table(grp=gsub("^X ", "", i), type=xnam, value=cor(fc1, fc2, method="spearman")))
#       }
#     }
#   }
# }
# 
# res <- cbind(res, do.call(rbind, strsplit(res$grp, " ")))
# res[,treatment := cleanTreatment(V1)]
# ct.use <- c("Beta", "Alpha") #, "Delta", "Gamma")
# ggplot(res[V2 %in% ct.use & V1 != "GABA"], aes(x=type, y=value, color=V2, group=paste(V2, V1))) + 
#   #geom_jitter(height=0, width=0.1) + 
#   geom_line() + geom_point() +
#   facet_grid(. ~ treatment) + theme_bw(12) + xRot()
# ggsave(dirout(out, "logFC_Alpha_Beta_", orgx,".pdf"),w=6,h=4)
# 
# res.wide <- dcast.data.table(res, treatment + V2 ~ type, value.var="value")
# res.wide[,diff := corrected - raw]
# ggplot(res.wide, aes(x=treatment, y=V2, fill=diff)) + geom_tile() + scale_fill_gradient2() + xRot()
# ggsave(dirout(out, "logFC_AllCelltypes_", orgx,".pdf"),w=4,h=6)