require(enrichR)
require(data.table)



# Get Enrichr gene sets ---------------------------------------------------
enrichrGetGenesets <- function(databases){
  setNames(lapply(databases, function(dbx){
    fpath <- paste0("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=",dbx)
    fhandle <- file(fpath)
    dblines <- tryCatch({
      readLines(con=fhandle)
    }, error=function(e){
      message(e, "\nFailed reading database: ", dbx)
      NULL
    })
    close(fhandle)
    if(is.null(dblines)){
      return(list())
    }else {
      res <- strsplit(dblines, "\t")
      names(res) <- sapply(res, function(x) x[1])
      res <- lapply(res, function(x) x[3:length(x)])
      return(res)
    }
  }), databases)
}


# Call Enrichr --------------------------------------------------
enrichGeneList.oddsRatio <- function(gene.list, databases = "KEGG_2016", fdr.cutoff = 0.1, genome.size=20000) {
  res <- enrichr(gene.list, databases)
  res <- lapply(res, data.table)
  res <- res[sapply(res, nrow) > 0]
  res <- do.call(rbind, lapply(names(res), function(rnam) data.table(res[[rnam]], database=rnam)))
  if(is.null(res)) return(data.table())
  res[,Genes := gsub(";", ",", Genes)]
  name.map <- list(
    c("Term", "category"), 
    c("P.value", "pval"),
    c("Adjusted.P.value", "qval"),
    c("Combined.Score", "combinedScore"),
    c("Genes", "genes"),
    c("Odds.Ratio", "oddsRatio"))
  for(i in 1:length(name.map)){
    colnames(res)[colnames(res) == name.map[[i]][1]] <- name.map[[i]][2]
  }
  
  res[,listLength := length(gene.list)]
  res[,dbLength := gsub("(\\d+)\\/(\\d+)", "\\2", Overlap)]
  res[,hitLength := gsub("(\\d+)\\/(\\d+)", "\\1", Overlap)]
  res <- res[,colnames(res) %in% c("database", "category", "pval", "zScore", "combinedScore", "genes", "qval", "dbLength", "listLength", "hitLength", "oddsRatio"),with=F]
  return(res)
}

enrichGeneList.oddsRatio.list <- function(geneLists, enrichrDBs="KEGG_2016"){
  enrichRes <- data.table()
  for(grp.x in names(geneLists)){
    print(paste0("Enrichment of list: ", grp.x))
    ret=enrichGeneList.oddsRatio(geneLists[[grp.x]],databases = enrichrDBs)
    if(nrow(ret) > 0){
      enrichRes <- rbind(enrichRes, data.table(ret, grp = grp.x))
    }
  }
  return(enrichRes)
}



# PLOT ENRICHR ------------------------------------------------------------
enrichr.plot <- function(enrichRes, qval.cap = 4, qval.cutoff = 0.05){
  enrichRes$category <- gsub("\\_(\\w|\\d){8}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){4}-(\\w|\\d){12}", "", enrichRes$category)
  enrichRes$category <- abbreviate(enrichRes$category, minlength=50) # substr(enrichRes$category2, 0, 50)
  enrichRes[,term := paste0(category, "_", dbLength)]
  enrichRes <- enrichRes[term %in% enrichRes[,.(min(qval)), by="term"][V1 < qval.cutoff]$term]
  if(!is.na(qval.cap)) enrichRes[, mLog10Q := pmin(-log10(qval),qval.cap)]
  
  # order terms by similarity (of OR)
  if(length(unique(enrichRes$term)) >= 2){
    try({
      orMT <- t(as.matrix(dcast.data.table(enrichRes, grp ~ term, value.var="oddsRatio")[,-"grp",with=F]))
      orMT[is.na(orMT)] <- 1
      orMT[orMT == Inf] <- max(orMT[orMT != Inf])
      hclustObj <- hclust(dist(orMT))
      enrichRes$term <- factor(enrichRes$term, levels=hclustObj$labels[hclustObj$order])
    },silent=T)
  }
  
  # order groups by similarity (of OR)
  if(length(unique(enrichRes$grp)) >= 2){
    try({
      orMT <- t(as.matrix(dcast.data.table(enrichRes, term ~ grp, value.var="oddsRatio")[,-"term",with=F]))
      orMT[is.na(orMT)] <- 1
      orMT[orMT == Inf] <- max(orMT[orMT != Inf])
      hclustObj <- hclust(dist(orMT))
      enrichRes$grp <- factor(enrichRes$grp, levels=hclustObj$labels[hclustObj$order])
    }, silent=T)
  }
  
  # plot
  p <- ggplot(enrichRes, aes(x=grp, y=term, size=log10(oddsRatio), color=mLog10Q)) + 
    geom_point() + scale_color_gradient(low="grey", high="red") + theme_bw(12) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  if(!is.na(qval.cap)) p <- p + ggtitle(paste0("-log10(q) capped at ", qval.cap))
}

enrichr.plot.many <- function(enrichRes, out, label="Enrichr", hitLengthCutoff = 3,...){
  enrichRes <- enrichRes[hitLength >= hitLengthCutoff]
  if(nrow(enrichRes) > 0){
    write.table(enrichRes, file=paste0(out, label, ".tsv"), sep="\t", quote=F, row.names=F)
    enrClass <- enrichRes$database[1]
    for(enrClass in unique(enrichRes$database)){
      enrichResX <- enrichRes[database == enrClass]
      if(nrow(enrichResX) > 0){
        enrichr.plot(enrichResX, ...)
        ggsave(paste0(out, label, "_", enrClass, ".pdf"), width=min(29, 6+ length(unique(enrichResX$grp))*0.3), height=min(29, length(unique(enrichResX$category))*0.3 + 4))
      }
    }
  }
}



# Fisher exact test with defined background -------------------------------------------------------
fisher.test.enrichment <- function(geneSets, gene.list, bg=NULL){
  #   geneSets <- list(KEGG=list(x=LETTERS[1:8], y=LETTERS[1:4]), GO=list(x=LETTERS[5:14]))
  #   gene.list <- list(One=LETTERS[5:10], Two=LETTERS[10:20])
  #   bg <- LETTERS[1:24]
  if(is.null(bg)) bg <- do.call(c, gene.list)
  bg <- unique(bg)
  res <- data.table()
  for(gl in names(gene.list)){
    for(db in names(geneSets)){
      for(gs in names(geneSets[[db]])){
        geneset.x <- unique(geneSets[[db]][[gs]])
        genelist.x <- unique(gene.list[[gl]])
        cont.table <- table((bg %in% geneset.x) + 1, (bg %in% genelist.x) + 10)
        if(!all(dim(cont.table) == c(2,2))){
          ft <- list(p.value = 1, estimate=1)
          cont.table <- matrix(0, 2,2)
        } else {
          ft <- fisher.test(cont.table)
        }
        res <- rbind(res, data.table(
          database=db, 
          geneset=gs, 
          list=gl, 
          pval=ft$p.value, 
          oddsRatio=ft$estimate, 
          genes=paste(intersect(geneset.x,genelist.x), collapse=","), 
          geneset.length=sum(cont.table[2,]),
          genelist.length=sum(cont.table[,2]),
          overlap.length=cont.table[2,2]
        ))
      }
    }
  }
  return(res)
}

