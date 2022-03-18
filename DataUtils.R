
##################################################
## R scripts for OmicsNet 
## Description: Data IO functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

.on.public.web <<- TRUE;

# init resources for analysis
Init.Data<-function(){
  print("init");
  .on.public.web <<- TRUE;
  if(exists("dataSet")){
    return(1)
  }

  initNet();
  net.info <-list()
  net.info$gene.ids <- vector();
  net.info$protein.ids <- vector();
  net.info$int.ids <- vector();
  net.info <<- net.info;
  partialToBeSaved <<- c("Rload.RData", "Rhistory.R")
  dataSet <- list();    
  dataSet$toremove <- c("Metabolic pathways",
                 "Biosynthesis of secondary metabolites",
                 "Microbial metabolism in diverse environments",
                 "Biosynthesis of antibiotics",
                 "Carbon metabolism",
                 "2-Oxocarboxylic acid metabolism",
                 "Fatty acid metabolism",
                 "Biosynthesis of amino acids",
                 "Degradation of aromatic compounds");
  uploadedGraph <<- FALSE; 
  netbuild.opt <<- "first"
  isKo <<- FALSE

  options(stringsAsFactors = FALSE)
  graph_name <<- NULL;
  lib.path <<- "../../data/";
  data.org <<- NULL; 
  module.count <<- 0; 
  dataSet$seeds.proteins <- vector();
  dataSet$mat <- dataSet$exp.mat <- dataSet$seed <- dataSet$inv <- list();
  dataSet$mic.thresh <- 0.8;
  current.msg <<- "";
  if(file.exists("/home/glassfish/sqlite/")){
    sqlite.path <<- "/home/glassfish/sqlite/";  #public server + qiang local
  }else if(file.exists("/Users/xia/Dropbox/sqlite/")){
    sqlite.path <<- "/Users/xia/Dropbox/sqlite/"; #xia local
  }else if(file.exists("/Users/jeffxia/Dropbox/sqlite/")){
    sqlite.path <<- "/Users/jeffxia/Dropbox/sqlite/"; #xia local2
  }else if(file.exists("/home/zgy/sqlite/")){
    sqlite.path <<- "/home/zgy/sqlite/"; #zgy local
  }else if(file.exists("/Users/jessicaewald/sqlite")){
    sqlite.path <<- "/Users/jessicaewald/sqlite/"; #ewald local
  }else if(file.exists("/Users/lzy/sqlite")){
    sqlite.path <<- "/Users/lzy/sqlite/"; #luyao local
  }else{
    sqlite.path <<-"/media/zzggyy/disk/sqlite/"; #zgy local
  }
  return(.set.nSet(dataSet));
}

SetUploadedGraph<- function(bool){
  uploadedGraph <<- bool;
}

SetOrganism <- function(org){
  data.org <<- org;
}

SetNetType <- function(netType){
  net.type <<- netType;
}

SetCurrentDataMulti <- function(){
  dataSet$type = nms.vec;
  .set.nSet(dataSet);
}

SetCurrentInteractionMulti <- function(){
  dataSet$interaction.type = interaction.vec;
  .set.nSet(dataSet);
}


SetFileType <- function(fileType){
  file.type <<- fileType;
}

CheckQueryTypeMatch <- function(qvec, type){
  if(type == "snp"){
    queryType <- "rsid";
  }else if(type %in% c("gene","tf")){
    queryType <- "entrez";
  }else if(type %in% c("met", "peak", "m2m")){
    queryType <- "kegg";
  }else if(type == "mic"){
    queryType <- mic.taxa;
  }else if(type == "mir"){
    queryType <- "mir_acc";
  }
  
  GeneAnotDB <-doProteinIDMapping(qvec, queryType);
  if(is.null(nrow(GeneAnotDB)) && GeneAnotDB == 0){
    return(FALSE);
  }else{
    na.inx <- is.na(GeneAnotDB[,1]) | is.na(GeneAnotDB[,2]);
    return(sum(!na.inx)>0);
  }
}

# geneIDs is text one string, need to make to vector
PrepareInputList <- function(dataSetObj="NA", inputList, org, type, queryType, netInvType){
  if(!exists("dataSet")){
    Init.Data();
  }
  dataSet <- .get.nSet(dataSetObj); 
  dataSet$name <- type;
  dataSet$orig <- inputList;
  if(type == "ko"){
    isKo <<- TRUE;
  }
  
  current.msg <<- NULL;
  data.org <<- org;
  if(type == "gene" & exists("signature.gene")){
    if(nrow(signature.gene)>200){
      signature.gene <- as.matrix(signature.gene[1:200,]);
    }
    prot.mat <- gene.mat <- signature.gene;
  }else if(type == "met" & exists("signature.cmpds")){
    if(nrow(signature.cmpds)>50){
      signature.cmpds <- as.matrix(signature.cmpds[1:50,]);
    }
    prot.mat <- gene.mat <- signature.cmpds;
    queryType <- signature.type;
  } else if(type == "peak"){
    prot.mat <- gene.mat <- inputList;
  } else{
    prot.mat <- gene.mat <- .parseListInput(inputList, queryType);
  }
  
  gene.mat <<- gene.mat
  GeneAnotDB <-doProteinIDMapping(rownames(gene.mat), queryType);
  na.inx <- is.na(GeneAnotDB[,1]) | is.na(GeneAnotDB[,2]);

  if(type == "ko"){
    dataSet$ko <- GeneAnotDB;
  }
  if(sum(!na.inx) < 2){
    current.msg <<- "Less than two hits found in uniprot database. ";
    print(current.msg);
  }
  if(type == "met" || type == "mir" ){
    col1 <- 1;
    if(queryType == "mir_acc"){
    col1 <- 2;
    }
    hit.inx <- match(rownames(prot.mat), GeneAnotDB[,col1]);
    if(is.null(dim(prot.mat))){
      prot.mat <- matrix(prot.mat);
    }
    
    rownames(prot.mat) <- GeneAnotDB[,"gene_id"][hit.inx];
    na.inx1 <- is.na(rownames(prot.mat));
    prot.mat = as.matrix(prot.mat[!na.inx1,])
  }else if(type == "peak"){
    hit.inx <- match(rownames(prot.mat), GeneAnotDB[,1]);
    if(is.null(dim(prot.mat))){
      prot.mat <- matrix(prot.mat);
    }
    
    rownames(prot.mat)[!is.na(hit.inx)] <- as.character(GeneAnotDB[,2][!is.na(hit.inx)]);
  } else{
    rownames(prot.mat) <- GeneAnotDB[,"gene_id"];
    prot.mat <- prot.mat[!na.inx, , drop=F];
  } 
  # now merge duplicates
  prot.mat <- RemoveDuplicates(prot.mat, "mean", quiet=T); 
  
  if(is.null(dim(prot.mat))){
    prot.mat <- matrix(prot.mat);
  }
  
  seed.proteins <- rownames(prot.mat);
  dataSet$GeneAnotDB <- GeneAnotDB;

  if(type == "mic" && sum(gene.mat) == 0){
    gene.mat[gene.mat == 0,] <- 1;
    prot.mat <- gene.mat
  }

  if(type == "tf"){
    dataSet$mat[[type]] <- gene.mat;
    dataSet$exp.mat[[type]] <- prot.mat;
    dataSet$seed[[type]] <- prot.mat;
    dataSet$inv[[type]] <- netInvType;
  }else if(type %in% c("gene", "geneonly", "protein")){
    if(length(dataSet$mat[["gene"]]) == 0){
      dataSet$mat[["gene"]] <- gene.mat;
      dataSet$exp.mat[["gene"]] <- prot.mat;
      dataSet$seed[["gene"]] <- prot.mat;
      dataSet$inv[["gene"]] <- netInvType;
      if(type == "geneonly"){
        dataSet$gene_type_vec = rep(2, nrow(prot.mat))
      }else if (type == "protein"){
        dataSet$gene_type_vec = rep(3, nrow(prot.mat))
      }else if (type == "snp"){
        dataSet$gene_type_vec = rep(4, nrow(prot.mat))
      }else{
        dataSet$gene_type_vec = rep(1, nrow(prot.mat))
      }
    }else{
      dataSet$mat[["gene"]] <- rbind(dataSet$mat[["gene"]], gene.mat);
      dataSet$exp.mat[["gene"]] <- rbind(dataSet$exp.mat[["gene"]], prot.mat);
      dataSet$seed[["gene"]] <- rbind(dataSet$seed[["gene"]], prot.mat);
      dataSet$inv[["gene"]] <- netInvType;
      if(type == "geneonly"){
        dataSet$gene_type_vec = c(dataSet$gene_type_vec, rep(2, nrow(prot.mat)));
      }else if (type == "protein"){
        dataSet$gene_type_vec = c(dataSet$gene_type_vec, rep(3, nrow(prot.mat)));
      }else if (type == "snp"){
        dataSet$gene_type_vec = c(dataSet$gene_type_vec, rep(4, nrow(prot.mat)));
      }else{
        dataSet$gene_type_vec = c(dataSet$gene_type_vec, rep(1, nrow(prot.mat)));
      }
    }
  }else{
    dataSet$mat[[type]] <- gene.mat;
    dataSet$exp.mat[[type]] <- prot.mat;
    dataSet$seed[[type]] <- prot.mat;
    dataSet$inv[[type]] <- netInvType;
  }
  
  if(length(dataSet$type) == 1){
    dataSet$prot.mat = rbind(dataSet$prot.mat, prot.mat);
    dataSet$seeds.proteins <- c(dataSet$seeds.proteins, seed.proteins); 
  } else{
    for(i in 1:length(dataSet$exp.mat)){
      dataSet$prot.mat = c(dataSet$prot.mat, dataSet$exp.mat[[i]]);
      dataSet$seeds.proteins <- c(dataSet$seeds.proteins, rownames(dataSet$exp.mat[[i]]));
    }
  }
  dataSet$seeds.proteins <- unique(dataSet$seeds.proteins)
  dataSet$seeds.expr = as.matrix(dataSet$prot.mat);

  if(.on.public.web){
  .set.nSet(dataSet);
  return (seed.proteins);
  }else{
  return (.set.nSet(dataSet));
  }
}


# given a text from input area: one or two-col entries (geneID and logFC)
# parse out return the a matrix containing the logFc, with rows named by gene ID
# if no 2nd col (logFC), the value will be padded by 0s
.parseListInput <- function(geneIDs, IDtype){
  lines <- unlist(strsplit(geneIDs, "\r|\n|\r\n")[1]);
  if(substring(lines[1],1,1)=="#"){
    lines <- lines[-1];
  }
  if(IDtype == "meta" || IDtype == "name"){
    gene.lists <- lines;
  }else{
    gene.lists <- strsplit(lines, "\\s+");
  }
  gene.mat <- do.call(rbind, gene.lists);
  
  if(dim(gene.mat)[2] == 1){ # add 0
    gene.mat <- cbind(gene.mat, rep(0, nrow(gene.mat)));
    current.msg <- "Only one column found in the list. Abundance will be all zeros. ";
  }else if(dim(gene.mat)[2] > 2){
    gene.mat <- gene.mat[,1:2];
    current.msg <- "More than two columns found in the list. Only first two columns will be used. ";
  }
  
  rownames(gene.mat) <- gene.mat[,1];
  gene.mat <- gene.mat[,-1, drop=F];
  gene.mat <- RemoveDuplicates(gene.mat, "mean", quiet=T); 
  good.inx <- !is.na(gene.mat[,1]);
  gene.mat <- gene.mat[good.inx, , drop=F];
  
  if(IDtype %in% c("class", "family","order" ,"genus","species","strain", "phylum")){
    rownames(gene.mat) <- gsub("_", " ", rownames(gene.mat));
  }

  return(gene.mat);
}


SetNetInv <- function(netInv){
  net.inv <<- netInv;
}

SetDbType <- function(dbType){
  dbu.type <<- dbType;
}

SetNetwType <- function(netwType){
  netw.type <<- netwType;
}

SetMirBo <- function(mirBo){
  requireMirExp <<- mirBo
}

GetInputListStat <- function(listNm){
  geneList <- rownames(dataSet$exp.mat[["gene"]]);
  tfList <- rownames(dataSet$exp.mat[["tf"]]);
  mirList <- rownames(dataSet$exp.mat[["mir"]]);
  metList <- rownames(dataSet$exp.mat[["met"]]);
  peakList <- rownames(dataSet$seed[["peakRaw"]]);
  micList <- rownames(dataSet$exp.mat[["mic"]]);
  snpList <- rownames(dataSet$seed[["snp"]]);
  
  koList <- rownames(dataSet$exp.mat[["ko"]]);
  return(c(length(rownames(dataSet$exp.mat[["gene"]])[dataSet$gene_type_vec == 1]),length(rownames(dataSet$exp.mat[["gene"]])[dataSet$gene_type_vec == 3]),length(unique(tfList)),length(unique(mirList)),length(unique(metList)), length(unique(koList)), length(unique(peakList)), length(unique(micList)), length(unique(snpList))));
}

GetFilesToBeSaved <-function(naviString){
  return(unique(partialToBeSaved));
}

GetCurrentJson <-function(type){
  return(dataSet$jsonNms[[type]]);
}

SetSnpTissueType <- function(type){
  dataSet$snpTissueType <<- type;
}

#'Export R Command History
#'@param mSetObj Input the name of the created dataSetObj (see Init.Data)
#'@export
GetRCommandHistory <- function(dataSetObj=NA){
  dataSet <- .get.nSet(dataSetObj); 
  if(length(dataSet$cmdSet) == 0){
    return("No commands found");
  }
  return(dataSet$cmdSet);
}

#'Record R Commands
#'@param mSetObj Input the name of the created dataSetObj (see Init.Data)
#'@param cmd Commands 
#'@export
RecordRCommand <- function(dataSetObj=NA, cmd){
  dataSet <- .get.nSet(dataSetObj); 
  dataSet$cmdSet <- c(dataSet$cmdSet, cmd);
  return(.set.nSet(dataSet));
}


SaveRCommands <- function(dataSetObj=NA){
  dataSet <- .get.nSet(dataSetObj);  
  cmds <- paste(dataSet$cmdSet, collapse="\n");
  pid.info <- paste0("# PID of current job: ", Sys.getpid());
  cmds <- c(pid.info, cmds);
  write(cmds, file = "Rhistory.R", append = FALSE);
}
