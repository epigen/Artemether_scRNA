hierarch.ordering <- function(dt, toOrder, orderBy, value.var, aggregate = FALSE){
  if(!aggregate){orMT <- t(as.matrix(dcast.data.table(dt, get(orderBy) ~ get(toOrder), value.var=value.var)[,-"orderBy",with=F]))}
  else{orMT <- t(as.matrix(dcast.data.table(dt, get(orderBy) ~ get(toOrder), value.var=value.var, fun.aggregate=mean)[,-"orderBy",with=F]))}
  orMT[is.na(orMT)] <- 1
  orMT[orMT == Inf] <- max(orMT[orMT != Inf])
  hclustObj <- hclust(dist(orMT))
  dt[[toOrder]] <- factor(dt[[toOrder]], levels=hclustObj$labels[hclustObj$order])
  return(dt)
}

jaccard <- function(ll){
  jacc <- matrix(NA, ncol=length(ll), nrow=length(ll))
  row.names(jacc) <- names(ll)
  colnames(jacc) <- names(ll)
  for(g1 in names(ll)){
    for(g2 in names(ll)){
      g1x <- ll[[g1]]
      g2x <- ll[[g2]]
      jacc[g1, g2] <- length(intersect(g1x, g2x))/length(union(g1x,g2x))}}
  diag(jacc) <- 0
  return(jacc)
}

toMT <- function(dt, row, col, val){
  retDT <- dcast.data.table(dt, get(row) ~ get(col), value.var=val)
  retMT <- as.matrix(retDT[,-"row"])
  row.names(retMT) <- retDT$row
  return(retMT)
}

toGR <- function(charVec){
  if(length(charVec) == 0) return(NULL)
  gr <- strsplit(charVec, "_"); 
  gr <- data.table(do.call(rbind, gr[sapply(gr, length) == 3]));
  colnames(gr) <- c("CHR","START","END")
  gr[, START := as.numeric(START)]
  gr[, END := as.numeric(END)]
  as(gr, "GRanges")
}

grToDT <- function(gr){
  data.table(chr=as.character(seqnames(gr)),
             start=start(gr),
             end=end(gr))}

grToVec <- function(gr){
  x <- grToDT(gr)
  paste(x$chr, x$start, x$end, sep="_")
}

goodGrStrings <- function(charVec){
  x <- strsplit(charVec, "_")
  sapply(x, length) == 3
}

write.tsv <- function(...){
  write.table(..., sep="\t", row.names=FALSE, quote=FALSE);
}

mapNumericToColors <- function(x, cols=c("royalblue1", "gray80", "indianred"), alpha=1){
  x <- x/max(abs(x))
  res <- c()
  res[which(x >= 0)] <- apply(colorRamp(c(cols[2], cols[3]))(x[which(x >= 0)]), 1, function(row) rgb(row[1]/255, row[2]/255, row[3]/255, alpha))
  res[which(x < 0)] <- apply(colorRamp(c(cols[2], cols[1]))(abs(x[which(x < 0)])), 1, function(row) rgb(row[1]/255, row[2]/255, row[3]/255, alpha))
  res
}

xRot <- function(){theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))}

gg.removeGrids <- function(){theme(panel.grid=element_blank())}

cleanDev <- function(n=2){sapply(1:n, function(i){try({dev.off()}, silent=TRUE)}); return(TRUE)}

corS <- function(...){cor(..., method="spearman")}

minMax <- function(x){ (x-min(x))/(max(x) - min(x))}

fishers.method <- function(p){pchisq((sum(log(p))*-2), df=length(p)*2, lower.tail=F)}
