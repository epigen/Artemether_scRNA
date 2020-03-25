require(project.init)

project.init2("artemether") 
out <- "GEO/"
dir.create(dirout(out))
outFTP <- paste0(out, "FTP/")
dir.create(dirout(outFTP))



# FIND ALL FASTQs ---------------------------------------------------------
bsf.dirs <- unique(c(
  "/data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX/",
  "/data/groups/lab_bsf/sequences/BSF_0410_HNF7TBBXX_l1_10x/fastq_path/HNF7TBBXX/",
  "/data/groups/lab_bsf/sequences/BSF_0412_HNJH5BBXX_l1_10x/fastq_path/HNJH5BBXX/",
  "/data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX/",
  "/data/groups/lab_bsf/sequences/BSF_0409_HNCYFBBXX_l1_l2_l3_l4_l5_l6_l7_10x/fastq_path/HNCYFBBXX/",
  "/data/groups/lab_bsf/sequences/BSF_0410_HNF7TBBXX_l1_10x/fastq_path/HNF7TBBXX/",
  "/data/groups/lab_bsf/sequences/BSF_0412_HNJH5BBXX_l1_10x/fastq_path/HNJH5BBXX/",
  "/data/groups/lab_bsf/sequences/BSF_0429_HNTVGBBXX_l5_l6_l7_l8_10x/fastq_path/HNTVGBBXX/",
  "/data/groups/lab_bsf/sequences/BSF_0460_HV235BBXX_l1_l7_l8_10x/fastq_path/HV235BBXX/",
  "/data/groups/lab_bsf/sequences/BSF_0525_H2K53BBXY_l2_l3_l4_l5_l6_l7_l8_10x/fastq_path/H2K53BBXY/"
))
ff <- do.call(c, lapply(bsf.dirs, function(x){list.files(x, recursive=T, full.names=T)}))
names(ff) <- basename(ff)
ff <- ff[grepl("hIslet", names(ff)) | grepl("mIslet", names(ff)) | grepl("MF179", names(ff))]


# LIST OF RAW FASTQ FILES -----------------------------------------------------
if(!file.exists(dirout(out, "RawAnnotation.tsv"))){
  raw.ann <- data.table(path = ff)
  raw.ann[,file := basename(path)]
  raw.ann[,single := ifelse(grepl("I1", file), "single", "paired")]
  raw.ann[grepl("I1", file), length := 8]
  raw.ann[grepl("R1", file), length := 26]
  raw.ann[grepl("R2", file), length := 58]
  raw.ann[, sample := gsub("_S\\d{1,3}_.+$", "", file)]
  with(raw.ann, table(sample, length))
  for(i in 1:nrow(raw.ann)){
    raw.ann[i, checksum := md5sum(path)]
  }
  raw.ann <- raw.ann[!grepl("pla2g16", sample)]
  raw.ann <- raw.ann[!grepl("LT1", sample)]
  stopifnot(nrow(raw.ann) == length(unique(raw.ann$file)))
  write.tsv(raw.ann, dirout(out, "RawAnnotation.tsv"))
}

# LIST OF PAIRED-END RAW FILES
raw.ann <- fread(dirout(out, "RawAnnotation.tsv"))
paired.ann <- raw.ann[grepl("R1", file), "file", with=F]
paired.ann[,file2 := gsub("R1", "R2", file)]
write.tsv(paired.ann, dirout(out, "PairedAnnotation.tsv"))


# COLLECT PROCESSED DATA ----------------------------------------------------------
outFTP_processed <- paste0(out, "FTP_processed/")
dir.create(dirout(outFTP_processed))
outZIP <- paste0(out, "ZIP/")
dir.create(dirout(outZIP))

all.samples <- c("Human1_2", "Human3", "MouseLTI", "MF179")
samplex <- all.samples[1]
sample.to.run <- data.table()
for(samplex in all.samples){
  
  # Clean annotation
  meta <- NA
  if(samplex == "Human1_2") meta <- loadHMetaUnfiltered()
  if(samplex == "MouseLTI") meta <- loadMMetaUnfiltered()
  if(samplex == "Human3") meta <- loadH3MetaUnfiltered()
  if(samplex == "MF179") meta <- fread(dirout("05_SpikeIns_Clean/", "CrossSpecies_Alignment.tsv"))
  if("treatment" %in% colnames(meta)) meta <- meta[treatment != "pla2g16"]
  write.tsv(meta, dirout(outFTP_processed, samplex, "_CellAnnotation_final.tsv"))
  
  # Raw annotation
  if(samplex != "MF179"){
    merge.csv <- fread(paste0(Sys.getenv("PROCESSED"), '/artemether/results_pipeline/cellranger_count/',samplex,'/outs/aggregation_csv.csv'), header=T)
    sample.to.run <- rbind(sample.to.run, data.table(sample=samplex, run=sapply(strsplit(merge.csv$molecule_h5, "/"), function(x) rev(x)[3])))
    #     cell.ann <- fread(paste0(Sys.getenv("PROCESSED"), '/artemether/results_pipeline/cellranger_count/',samplex,'/outs/cell_barcodes.csv'), header=F)
    #     cell.ann[,i := as.numeric(gsub(".+-(\\d+)", "\\1", V2))]
    #     cell.ann$sample <- gsub("^\\s*(.+)\\s*$", "\\1", merge.csv[cell.ann$i]$library_id)
    #     cell.ann <- cell.ann[,c("V2", "sample"), with=F]
    #     colnames(cell.ann) <- c("Cell barcode", "Sample")
    #     write.tsv(cell.ann, dirout(outFTP_processed, samplex, "_CellAnnotation_raw.tsv"))
  } else {
    sample.to.run <- rbind(sample.to.run, data.table(sample=samplex, run=samplex))
  }
  
  # Corrected reads
  data <- NA
  if(samplex == "Human1_2") data <- loadHData()
  if(samplex == "MouseLTI") data <- loadMData()
  if(samplex == "Human3") data <- loadH3Data()
  if(!is.na(data)){
    data.export <- data@raw.data[,data.table(data@meta.data, keep.rownames=T)[!grepl("pla2g16", sample)]$rn]
    writeMM(data@raw.data, dirout(outFTP_processed, samplex, "_correctedTPM_Matrix.mtx"))
    write.table(data@raw.data@Dimnames[[1]], dirout(outFTP_processed, samplex, "_correctedTPM_Matrix_genes.txt"), col.names=F, quote=F, row.names=F)
    write.table(data@raw.data@Dimnames[[2]], dirout(outFTP_processed, samplex, "_correctedTPM_Matrix_barcodes.txt"), col.names=F, quote=F, row.names=F)
    #system(paste("zip", dirout(outFTP_processed, samplex, "_correctedTPM_Matrix.zip"), dirout(outZIP, samplex, "_correctedTPM_Matrix*")))
  }
}
sample.to.run <- sample.to.run[!grepl("pla2g16", run)]
write.tsv(sample.to.run, dirout(out, "Samples2run.mapping.tsv"))



# Generate Table of processed data ----------------------------------------------------------
processed.ann <- data.table(file.path=list.files(dirout(outFTP_processed), full.names=T))
# for(samplex in all.samples){
#   fx <- paste0(Sys.getenv("PROCESSED"), '/artemether/results_pipeline/cellranger_count/',samplex, ifelse(samplex == "MF179", "_hmG", ""),'/outs/raw_gene_bc_matrices_h5.h5')
#   stopifnot(file.exists(fx))
#   processed.ann <- rbind(processed.ann, data.table(file.path = fx))
# }
for(samplex in all.samples){
  processed.ann[grepl(samplex, file.path), sample := samplex]
}
processed.ann[, file := basename(file.path)]
processed.ann[grepl("raw_gene_bc_matrices_h5.h5$", file.path), file := paste0(sample, "_", basename(file.path))]
for(i in 1:nrow(processed.ann)){
  processed.ann[i, checksum := md5sum(file.path)]
}
stopifnot(length(unique(processed.ann$file))== nrow(processed.ann))
write.tsv(processed.ann, dirout(out, "ProcessedAnnotation.tsv"))


# List of samples and raw / processed data --------------------------------
processed.ann <- fread(dirout(out, "ProcessedAnnotation.tsv"))
processed.by.sample <- split(processed.ann$file, factor(processed.ann$sample))
processed.by.sample.max <- max(sapply(processed.by.sample,length))
processed.by.sample <- lapply(processed.by.sample, function(x){return(c(x, rep("", times=processed.by.sample.max-length(x))))})

raw.ann <- fread(dirout(out, "RawAnnotation.tsv"))
raw.by.run <- split(raw.ann$file, factor(raw.ann$sample))
raw.by.run.max <- max(sapply(raw.by.run,length))
raw.by.run <- lapply(raw.by.run, function(x){return(c(x, rep("", times=raw.by.run.max-length(x))))})

run.ann <- fread(dirout(out, "Samples2run.mapping.tsv"))
run.ann[sample == "Human3", run := gsub("_human$", "", run)]
run.ann[,molecule := "mRNA"]
run.ann[grepl("^hI", run),organism := "Human"]
run.ann[grepl("^mI", run),organism := "Mouse"]
run.ann[grepl("^MF179", run),organism := "Human/Mouse"]
run.ann[,source_name := "Pancreatic Islet cells"]
run.ann[grepl("DMSO", run), treament:= "DMSO (Control)"]
run.ann[grepl("GABA", run), treament:= "GABA"]
run.ann[grepl("FoxO", run), treament:= "FOXO inhibitor"]
run.ann[is.na(treament), treament:= paste("Artemether (10muM)")]
run.ann[grepl("36hr", run), treament_time := "36hr"]
run.ann[grepl("72hr", run), treament_time := "72hr"]
run.ann[is.na(treament_time), treament_time := "72hr"]
run.ann[is.na(treament_time), treament_time := "72hr"]
run.ann <- cbind(run.ann, do.call(rbind, processed.by.sample[run.ann$sample]))
run.ann <- cbind(run.ann, data.table(do.call(rbind, raw.by.run[run.ann$run])))
write.tsv(run.ann, dirout(out, "Run_Annotation.tsv"))


# Commands to upload files ------------------------------------------------
write("", dirout(out, "Upload_commands.sh"))
processed.ann <- fread(dirout(out, "ProcessedAnnotation.tsv"))
for(i in 1:nrow(processed.ann)){
  stopifnot(file.exists(processed.ann[i]$file.path))
  write(paste("put -z", processed.ann[i]$file.path, processed.ann[i]$file), dirout(out, "Upload_commands.sh"), append=T)
}
raw.ann <- fread(dirout(out, "RawAnnotation.tsv"))
for(i in 1:nrow(raw.ann)){
  stopifnot(file.exists(raw.ann[i]$path))
  write(paste("put -z", raw.ann[i]$path, raw.ann[i]$file), dirout(out, "Upload_commands.sh"), append=T)
}




