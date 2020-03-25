require("project.init")
require(Matrix)
require(methods)
project.init2("artemether")

out <- "EXT_02_GEO_Analysis/"
inDir <- "EXT_01_GEO/"
dir.create(dirout(out))

list.files(dirout(inDir))


# GSE73727
file.GSE73727 <- dirout(out, "GSE73727.RData")
if(!file.exists(file.GSE73727)){
  pathGSE73727 <- dirout(inDir, "/GSE73727/")
  data.GSE73727 <- do.call(cbind, lapply(list.files(pathGSE73727, pattern="RPKM"), function(f){
    x <- read.table(paste0(pathGSE73727, f), sep="\t", header=TRUE)
    x <- data.table(x)[,max(RPKM), by="ensG"]
    data.frame(data=x$V1, row.names=x$ensG)
  }))
  mt.GSE73727 <- as.matrix(data.GSE73727)
  mt.GSE73727 <- toSparseMT(mt.GSE73727)
  colnames(mt.GSE73727) <- gsub("(GSM\\d+).+", "\\1", list.files(pathGSE73727, pattern="RPKM"))
  str(mt.GSE73727)
  #design
  des <- readLines(dirout(inDir, "/GSE73727_series_matrix.txt.gz"))
  substr(des, 0,40)
  des <- data.table(
    sample=gsub('\\"', "", strsplit(des[36], "\t")[[1]]),
    cell_type=gsub("assigned cell type: ", "", gsub('\\"', "", strsplit(des[44], "\t")[[1]])))
  colnames(mt.GSE73727) <- make.unique(des[match(colnames(mt.GSE73727), des$sample)]$cell_type)
  save(mt.GSE73727, file=file.GSE73727)
}else { load(file.GSE73727)}


# GSE81547
file.GSE81547 <- dirout(out, "GSE81547.RData")
if(!file.exists(file.GSE81547)){
  pathGSE81547 <- dirout(inDir, "/GSE81547/")
  data.GSE81547 <- do.call(cbind, lapply(list.files(pathGSE81547), function(f){
    read.table(paste0(pathGSE81547, f),row.names=1)
  }))
  mt.GSE81547 <- as.matrix(data.GSE81547)
  mt.GSE81547 <- toSparseMT(mt.GSE81547)
  colnames(mt.GSE81547) <- gsub("(GSM\\d+).+", "\\1", list.files(pathGSE81547))
  
  #design
  des <- readLines(dirout(inDir, "/GSE81547_series_matrix.txt.gz"))
  substr(des, 0,40)
  des <- data.table(
    sample=gsub('\\"', "", strsplit(des[32], "\t")[[1]]),
    cell_type=gsub("inferred_cell_type: ", "", gsub('\\"', "", strsplit(des[42], "\t")[[1]])))
  colnames(mt.GSE81547) <- make.unique(des[match(colnames(mt.GSE81547), des$sample)]$cell_type)
  save(mt.GSE81547, file=file.GSE81547)
} else { load(file.GSE81547)}


# GSE81608
file.GSE81608 <- dirout(out, "GSE81608.RData")
if(!file.exists(file.GSE81608)){
  data <- read.table(dirout(inDir, "/GSE81608_human_islets_rpkm.txt.gz"), sep="\t", header=T)
  str(data)
  
  # Design
  des <- readLines(dirout(inDir, "/GSE81608_series_matrix.txt.gz"))
  #substr(des, 0,40)
  des <- data.table(
    sample=strsplit(des[28], "\t")[[1]],
    cell_type=strsplit(des[43], "\t")[[1]])
  des <- des[-1]
  des[,sample := gsub("Pancreatic islet cell s", "S", sample)]
  des[,cell_type := gsub("cell subtype: ", "", cell_type)]
  row.names(data) <- paste0("Gene_", data$gene.id)
  data$gene.id <- NULL
  des[,sample := gsub(" ", "_", sample)]
  des[,sample := gsub('\\"', "", sample)]
  des[,cell_type := gsub('\\"', "", cell_type)]
  head(des)
  colnames(data) <- make.unique(des[match(colnames(data), des$sample)]$cell_type)
  mt.GSE81608 <- as.matrix(data)
  mt.GSE81608 <- toSparseMT(mt.GSE81608)
  save(mt.GSE81608, file=file.GSE81608)
}else { load(file.GSE81608)}


# GSE83139
file.GSE83139 <- dirout(out, "GSE83139.RData")
if(!file.exists(file.GSE83139)){
  data.GSE83139 <- read.table(dirout(inDir, "GSE83139_tbx-v-f-norm-ntv-cpms.csv.gz"), sep="\t", header=T)
  mt.GSE83139 <- data.GSE83139
  row.names(mt.GSE83139) <- make.names(mt.GSE83139$gene)
  mt.GSE83139 <- as.matrix(mt.GSE83139[,8:ncol(mt.GSE83139)])
  mt.GSE83139 <- toSparseMT(mt.GSE83139)
  str(mt.GSE83139)
  save(mt.GSE83139, file=file.GSE83139)
}else { load(file.GSE83139)}


# GSE84133
file.GSE84133 <- dirout(out, "GSE84133.RData")
if(!file.exists(file.GSE84133)){
  pathGSE84133 <- dirout(inDir, "/GSE84133/")
  data.GSE84133 <- lapply(list.files(pathGSE84133), function(f){
    message(f)
    x <- read.csv(paste0(pathGSE84133, f))
    row.names(x) <- make.unique(paste0(gsub(".+(lib\\d+).+", "\\1", x$X), "_", x$assigned_cluster))
    x1 <- t(as.matrix(x[,4:ncol(x)]))
    print(str(x1))
    x2 <- toSparseMT(x1)
    print(str(x2))
    return(x2)
  })
  lapply(data.GSE84133, str)
  names(data.GSE84133) <- gsub("GSM\\d+_(.+)_umifm_counts.csv.gz", "\\1", list.files(pathGSE84133))
  data.GSE84133 <- lapply(data.GSE84133, function(x){
    x <- x[,Matrix::colSums(x) > MIN.UMIS]
    x <- x[,Matrix::colSums(x != 0) > MIN.GENES]
    x <- toTPM_fromRaw(x, SCALE.FACTOR)
  })
  save(data.GSE84133, file=file.GSE84133)
} else { { load(file.GSE84133)}}


# GSE85241
file.GSE85241 <- dirout(out, "GSE85241.RData")
if(!file.exists(file.GSE85241)){
  data.GSE85241 <- read.table(dirout(inDir, "GSE85241_cellsystems_dataset_4donors_updated.csv.gz"), sep="\t", header=T)
  mt.GSE85241 <- data.GSE85241
  row.names(mt.GSE85241) <- make.unique(gsub("__chr.+$", "", row.names(mt.GSE85241)))
  mt.GSE85241 <- as.matrix(mt.GSE85241)
  mt.GSE85241 <- toSparseMT(mt.GSE85241)
  mt.GSE85241 <- t(t(mt.GSE85241)/Matrix::colSums(mt.GSE85241)) * SCALE.FACTOR
  save(mt.GSE85241, file=file.GSE85241)
}else { load(file.GSE85241)}


# E-MTAB-5061
file.MTAB_5061 <- dirout(out, "MTAB_5061.RData")
if(!file.exists(file.MTAB_5061)){
  data.MTAB_5061 <- read.table(dirout(inDir, "E-MTAB-5061.pancreas_refseq_rpkms_counts_3514sc.txt"), sep="\t", header=F)
  row.names(data.MTAB_5061) <- make.unique(data.MTAB_5061$V1)
  he <- readLines(dirout(inDir, "E-MTAB-5061.pancreas_refseq_rpkms_counts_3514sc.txt"))
  he1 <- strsplit(he[1], "\t")[[1]]
  he1 <- he1[-1]
  mt.MTAB_5061 <- data.MTAB_5061[,3:(length(he1) + 2)]
  colnames(mt.MTAB_5061) <- he1
  mt.MTAB_5061 <- as.matrix(mt.MTAB_5061)
  mt.MTAB_5061 <- toSparseMT(mt.MTAB_5061)

  #design
  des <- fread(dirout(inDir, "/E-MTAB-5061.sdrf.txt"))
  des <- des[match(colnames(mt.MTAB_5061), des[["Source Name"]])]
  colnames(mt.MTAB_5061) <- make.unique(des[["Factor Value[cell type]"]])
  stopifnot(all(des[["Source Name"]] == colnames(mt.MTAB_5061)))
  colnames(mt.MTAB_5061) <- make.unique(des[match(colnames(mt.MTAB_5061), des[["Source Name"]])][["Factor Value[cell type]"]])
  save(mt.MTAB_5061, file=file.MTAB_5061)
} else { (load(file.MTAB_5061))}



getExprGene <- function(mt, gene){ 
  if(gene == "INS") ids=c("INS", "Ins2", "ENSG00000254647")
  if(gene == "GCG") ids=c("GCG", "Gcg", "ENSG00000115263")
  mt[ids[which(ids %in% row.names(mt))],]
}

pDat <- data.table()
for(gg in c("INS", "GCG")){
  pDat <- rbind(pDat, data.table(
    Gene=gg, Expr=getExprGene(mt.GSE73727, gg), Sample=colnames(mt.GSE73727),
    Name="Li et al.",Type="Sorted"),fill=T)
  
  pDat <- rbind(pDat, data.table(
    Gene=gg, Expr=getExprGene(mt.GSE81547, gg), Sample=colnames(mt.GSE81547),
    Name="Enge et al.",Type="Sorted"),fill=T)  
  
  # pDat <- rbind(pDat, data.table(Gene=gg, Expr=getExprGene(mt.GSE81608, gg), Name="Xin et al.",Sample=colnames(mt.GSE81608)),fill=T)
  
  pDat <- rbind(pDat, data.table(
    Gene=gg, Expr=getExprGene(mt.GSE83139, gg), Sample=colnames(mt.GSE83139),
    Name="Wang et al.",Type="Sorted"),fill=T)
  
  pDat <- rbind(pDat, do.call(rbind, lapply(names(data.GSE84133), function(x) data.table(
    Gene=gg, 
    Expr=getExprGene(data.GSE84133[[x]], gg), 
    Name=paste("Baron et al.", x), 
    Type="Droplet",
    Sample=paste0(x, ".", colnames(data.GSE84133[[x]]))))),fill=T)
    
  pDat <- rbind(pDat, data.table(
    Gene=gg, Expr=getExprGene(mt.GSE85241, gg), Sample=colnames(mt.GSE85241),
    Name="Muraro et al.",Type="Sorted"),fill=T)

  pDat <- rbind(pDat, data.table(
    Gene=gg, Expr=getExprGene(mt.MTAB_5061, gg), Sample=colnames(mt.MTAB_5061),
    Name="Segerstolpe et al.",Type="Sorted"),fill=T)
}
pDat[grepl("alpha", Sample), cell_type := "alpha"]
pDat[grepl("beta", Sample), cell_type := "beta"]
pDat[is.na(cell_type), cell_type := "other"]

write.tsv(pDat, dirout(out, "INS_GCG_Data.tsv"))

pDat$Type <- factor(pDat$Type, levels=c("Sorted", "Droplet"))
pDat$cell_type <- factor(pDat$cell_type, levels=c("alpha", "beta", "other"))

ggplot(pDat, aes(x=Name, y=log(Expr+1))) + geom_violin(fill="lightblue") + 
  facet_grid(Gene ~ Type, space="free_x", scale="free") +
  theme_bw(12) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + 
  xlab('') + ylab("Log expression")
ggsave(dirout(out, "AllExpression.pdf"), height=4, width=4)

ggplot(pDat, aes(x=Name, y=log(Expr+1), fill=cell_type)) + geom_violin() + 
  facet_grid(Gene ~ Type, space="free_x", scale="free") +
  theme_bw(12) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + 
  xlab('') + ylab("Log expression") +
  scale_fill_manual(values=c("red", "blue", "black"))
ggsave(dirout(out, "CellTypes.pdf"), height=4, width=8)

lapply(data.GSE84133, function(x){
  ret <- colnames(x)
  tail(sort(table(gsub("\\.\\d+", "", lapply(strsplit(ret, "_"), function(xx) xx[2])))),4)
})