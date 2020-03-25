# TO DOWNLOAD THE RELEVANT DATA:
# ###
# ### HomoloGene
# ###
# cd ~/resources_nfortelny
# mkdir HomoloGene
# cd HomoloGene
# wget ftp://ftp.ncbi.nih.gov/pub/HomoloGene/build68/homologene.data -O build68.tsv
# cd ../
#
# ###
# ### Human and Mouse gene data and synonyms
# ###
# cd ~/resources_nfortelny
# wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz -O NCI_mouse_gene_Info.tsv
# wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz -O NCI_human_gene_Info.tsv



# SYNONYMS ----------------------------------------------------------------
syn <- readLines("~/resources_nfortelny/NCI_human_gene_Info.tsv")
syn[1]
syn <- syn[-1]
syn <- t(sapply(strsplit(syn, "\t"), function(x) c(x[3], x[5])))
syn.split <- strsplit(syn[,2], "\\|")
synDT.human <- data.table(gene=do.call(c,sapply(1:nrow(syn), function(x) rep(syn[x,1], length(syn.split[[x]])))), syn=do.call(c,syn.split))
synDT.human <- synDT.human[!syn %in% synDT.human$gene] #syn %in% synDT[,.N, by="syn"][N == 1]$syn & 
rm(list=c("syn"))

# MOUSE
syn <- readLines("~/resources_nfortelny/NCI_mouse_gene_Info.tsv")
syn[1]
syn <- syn[-1]
syn <- t(sapply(strsplit(syn, "\t"), function(x) c(x[3], x[5])))
syn.split <- strsplit(syn[,2], "\\|")
synDT.mouse <- data.table(gene=do.call(c,sapply(1:nrow(syn), function(x) rep(syn[x,1], length(syn.split[[x]])))), syn=do.call(c,syn.split))
synDT.mouse <- synDT.mouse[!syn %in% synDT.mouse$gene] #syn %in% synDT[,.N, by="syn"][N == 1]$syn & 
rm(list=c("syn"))

# FUNCTION
cleanGenes3 <- function(v, synDT){
  ret <- v
  for(i in unique(ret[!ret %in% synDT$gene])){
    gxx <- synDT[syn == i]$gene
    if(length(gxx) == 1){
      ret[ret == i] <- gxx
    } else{
      message("Can't map ", i, ", it maps to ", paste(gxx, collapse=" "))
    }
  }
  return(ret)
}


# HOMOLOGENE --------------------------------------------------------------
# mm.genes <- readLines("~/resources_nfortelny/NCI_mouse_gene_Info.tsv")
HOMOLOGENE <- fread("~/resources_nfortelny/HomoloGene/build68.tsv")
homologene <- HOMOLOGENE[V2 %in% c(9606, 10090)]
table(table(homologene$V1))
# Filter to only one to one mapping
homologene1to1 <- homologene[V1 %in% homologene[,.(count=.N, orgs=length(unique(V2))), by="V1"][count==2 & orgs==2]$V1]
homol <- data.table(mouse=homologene1to1[V2==10090][order(V1)]$V4, human=homologene1to1[V2==9606][order(V1)]$V4)
homol[toupper(mouse) != human]
