require("project.init")
require(Seurat)
require(Matrix)
require(methods)
require(doMC)
require(glmnet)
require(ROCR)

project.init2("artemether")

out <- "14_03_Classifier_moreCelltypes_noEndocrine/"
dir.create(dirout(out))

metaH <- loadHMeta()
pbmcH <- loadHData()
metaH3 <- loadH3Meta()
pbmcH3 <- loadH3Data()
metaM <- loadMMeta()
pbmcM <- loadMData()

with(metaH[treatment == "DMSO"], table(celltype, sample))
with(metaM[treatment == "DMSO"], table(celltype, sample))
with(metaH3[treatment == "DMSO"], table(celltype, sample))
  
org <- "human"
# c("mouse", "human", "human3")
for(org in c("mouse", "human", "human3")){

  meta <- metaH3
  data <- pbmcH3@data
  predGeneExclude <- c("INS", "GCG", "SST", "PPY")
  if(org == "mouse"){
    meta <- metaM
    data <- pbmcM@data
    predGeneExclude <- c("Ins2", "Gcg", "Sst", "Ppy")
  }
  if(org == "human"){
    meta <- metaH
    data <- pbmcH@data
  }
  celltypes.to.use <- names(which(table(meta[treatment == "DMSO"]$celltype) > 50))
  celltypes.to.use <- celltypes.to.use[celltypes.to.use != "Endocrine"]
  
  
  replicates <- length(unique(meta$replicate))
  
  outOrg <- paste0(out, org, "/")
  dir.create(dirout(outOrg))
  
  fit.file <- dirout(outOrg, "ModelFit.RData")
  if(!file.exists(fit.file)){
    trainData <- meta[treatment == "DMSO" & celltype %in% celltypes.to.use]
    set.seed(2121)
    cells.to.use <- unique(do.call(c, lapply(split(trainData$rn, factor(paste0(trainData$celltype, trainData$replicate))), function(x) sample(x, 100, replace=T))))
    trainData <- trainData[rn %in% cells.to.use]
    ggplot(trainData, aes(x=celltype, fill=replicate)) + geom_bar() + xRot()
    ggsave(dirout(outOrg, "TrainingCounts.pdf"), width=5, height=5)
    
    
    dat <- t(as.matrix(data[-which(rownames(data) %in% predGeneExclude),trainData$rn]))
    registerDoMC(cores=3)
    glmFit <- cv.glmnet(y=factor(trainData$celltype), x=dat[trainData$rn,],family="multinomial",alpha=1, parallel=TRUE, keep=TRUE, nfolds=5)
    str(coef(glmFit))
    trainData$fold <- glmFit$foldid
    classes <- glmFit$glmnet.fit$classnames
    lambda <- glmFit$lambda.1se
    lambda.idx <- which(glmFit$lambda == lambda)
    save(glmFit, trainData, lambda, lambda.idx, classes, file=fit.file)
  } else { load(fit.file)}
  
  with(trainData, table(fold, celltype))
  str(glmFit$fit.preval)
  
  
  # GET CV PERFORMANCE ------------------------------------------------------
  clx <- classes[1]
  cvRoc <- data.table()
  for(clx in classes){
    clx.i <- which(classes == clx)
    labels <- (trainData$celltype == clx) + 0
    predX <- prediction(predictions=glmFit$fit.preval[,clx.i,lambda.idx], labels=labels)
    auc <- performance(prediction.obj=predX, measure="auc")@y.values[[1]]
    perf <- performance(prediction.obj=predX, measure="tpr", x.measure="fpr")
    cvRoc <- rbind(cvRoc, data.table(celltype=clx, tpr=perf@y.values[[1]], fpr=perf@x.values[[1]], auc=auc, random=NA))
    for(rand.i in 1:50){
      predX <- prediction(predictions=glmFit$fit.preval[,clx.i,lambda.idx], labels=sample(labels))
      auc <- performance(prediction.obj=predX, measure="auc")@y.values[[1]]
      perf <- performance(prediction.obj=predX, measure="tpr", x.measure="fpr")
      cvRoc <- rbind(cvRoc, data.table(celltype=clx, tpr=perf@y.values[[1]], fpr=perf@x.values[[1]], auc=auc, random=rand.i))
    }
  }
  save(cvRoc, file=dirout(outOrg, "cvRoc.RData"))
  
  
  # CV
  cvRoc[,realData := is.na(random)]
  p <- ggplot(cvRoc, aes(x=fpr, y=tpr)) +  
    theme_bw(16) +
    facet_grid(celltype ~ .) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
  ggsave(dirout(outOrg, "ROC_CV.pdf"), width=4, height=length(classes) * 2 + 1,plot=p + geom_path(data=cvRoc[realData == TRUE], alpha=1))
  ggsave(dirout(outOrg, "ROC_CV_rand.jpg"), width=4, height=length(classes) * 2 + 1, plot=p+geom_path(data=cvRoc[realData == FALSE], alpha=0.3, color="grey", aes(group=random)))
  write.tsv(unique(cvRoc[realData==TRUE, c("celltype", "auc")]), file=dirout(outOrg, "AUC.tsv"))
  
  
  # COEFFICIENTS
  cfs <- sapply(classes, function(clx) coef(glmFit$glmnet.fit)[[clx]][,lambda.idx])
  cfs <- cfs[apply(cfs, 1, sum) > 0.1,]
  pdf(dirout(outOrg, "Coefficients.pdf"), width=length(classes) + 1, height=nrow(cfs) * 0.2 + 1, onefile=F); pheatmap(cfs); dev.off()
  
  # PREDICT
  pred.file <- dirout(outOrg, "Prediction.tsv")
  if(!file.exists(pred.file)){
    predData <- meta[!rn %in% trainData$rn]
    dat <- t(as.matrix(data[-which(rownames(data) %in% predGeneExclude),predData$rn]))
    pred <- predict(glmFit, dat, type = "response")
    predData <- cbind(predData,data.table(pred[predData$rn,,1]))
    write.tsv(predData, pred.file)
  } else {
    predData <- fread(pred.file)
  }
  
  # PREDICT
  pred.all.file <- dirout(outOrg, "Prediction_allCells.tsv")
  if(!file.exists(pred.all.file)){
    predDataAll <- copy(meta)
    dat <- t(as.matrix(data[-which(rownames(data) %in% predGeneExclude),predDataAll$rn]))
    pred <- predict(glmFit, dat, type = "response")
    predDataAll <- cbind(predDataAll,data.table(pred[predDataAll$rn,,1]))
    write.tsv(predDataAll, pred.all.file)
  } else {
   predDataAll <- fread(pred.all.file)
  }
  
  
  accData <- predData[treatment == "DMSO" & celltype %in% classes]
  accData$predCell <- apply(as.matrix(accData[,classes, with=F]),1,function(row) classes[which(row == max(row))])
  write.tsv(sum(accData$celltype == accData$predCell)/nrow(accData), dirout(outOrg, "Accuracy_test.tsv"))
  
  meanDat <- melt(predData, id.vars=c("treatment", "replicate", "celltype"), measure.vars=classes)[,.(value = mean(value)), by=c("treatment", "replicate", "celltype","variable")]
  ggplot(meanDat, aes(x=treatment, y=variable, fill=value)) + geom_tile() + facet_grid(replicate ~ celltype) + xRot() +
    scale_fill_gradient(low="white", high="red")
  ggsave(dirout(outOrg, "Predictions.pdf"), width=20, height=3*replicates)
  
  
  clx <- classes[1]
  for(clx in c("Alpha", "Beta")){
    ggplot(predData[celltype %in% c("Beta", "Alpha")], aes_string(x=clx, color="treatment", linetype="treatment")) + 
      facet_grid(celltype ~ replicate, scales="free") + stat_ecdf() + theme_bw(12)
    ggsave(dirout(outOrg, "Predict_", clx, ".pdf"), width=8, height=4*replicates)
  }
}

#ggplot(predData[treatment == "A10"], aes(x=Alpha, y=Beta)) + geom_point(alpha=0.3) + facet_grid(replicate ~ celltype) + theme_bw(12)
