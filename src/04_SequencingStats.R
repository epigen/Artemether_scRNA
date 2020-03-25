require(data.table)
require(project.init)

project.init2("artemether")

out <- "04_SequencingStats/"
dir.create(dirout(out))


datasets <- c(".")

for(dataset.x in datasets){
  ds <- list.dirs(path=paste0(Sys.getenv("PROCESSED"), "artemether/results_pipeline/cellranger_count/"), recursive=F)
  ds <- ds[grepl(dataset.x, ds)]
  ff <- list.files(paste0(ds, "/outs/"), pattern="^metrics_summary.csv$", recursive=F, full.names=T)
  st <- lapply(ff, fread)
  names(st) <- basename(gsub("\\/outs.+", "", ff))
  
  cc <- data.table()
  for(i in names(st)){
    cc <- rbind(cc, data.table(st[[i]], dataset=i), fill=T)
  }
  write.tsv(cc, dirout(out, "SequencingStats_", dataset.x, ".tsv"))
}
