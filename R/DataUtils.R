##################################################
## R scripts for OmicsNet
## Description: Data I/O functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

.onAttach <- function (libname, pkgname){
  .on.public.web <<- FALSE;
  k1 <- paste("OmicsNetR",
              utils::packageVersion( "OmicsNetR"),
              "initialized Successfully !")
  k0 <- "\n";
  packageStartupMessage(c(k1,k0));
}

#' Initialize dataSet object for downstream functions
#'
#' @export
Init.Data<-function(){

  #if(exists("dataSet")){
  #  return(1)
  #}
  if(!exists(".on.public.web")){
    .on.public.web <<- FALSE;
  }

  dataSet <- list();
  dataSet$imgSet <- list();
  dataSet$snpRegion <- F;
  dataSet$allOrgs <- "";
  dataSet$viewTable <- list();
  dataSet$require.exp <- T;
  dataSet$ppiZero <- F;
  dataSet$min.score <- 900;
  dataSet$seeds.proteins <- vector();
  dataSet$mat <- dataSet$exp.mat <- dataSet$seed <- list();
  dataSet$mic.thresh <- 0.8;
  dataSet$toremove <- c("Metabolic pathways",
                 "Biosynthesis of secondary metabolites",
                 "Microbial metabolism in diverse environments",
                 "Biosynthesis of antibiotics",
                 "Carbon metabolism",
                 "2-Oxocarboxylic acid metabolism",
                 "Fatty acid metabolism",
                 "Biosynthesis of amino acids",
                 "Degradation of aromatic compounds");

  dataSet$mic.exclude.opt <- "currency";
  dataSet$currExclude <- "true";
  dataSet$orphExclude <- dataSet$uniExclude <- "TRUE";
  dataSet$snp2gene$phesc.opt <- "eqtl";
  dataSet$snp2gene$vep.opt <- "dis";
  dataSet$snp2gene$vep.num <- 1;
  dataSet$snp2gene$vep.dis <- 5;


  dataSet <<- dataSet;

  # set up other global vars
  net.info <<- list(gene.ids=vector(),protein.ids=vector());
  initNetwork();
  uploadedGraph <<- FALSE;
  data.org <<- NULL;
  module.count <<- 0;
  max.row <<- 2000; # maximal rows to display interaction table in web page
  current.msg <<- "";
  partialToBeSaved <<- c("Rload.RData", "Rhistory.R")
  lib.path <<- "../../data/";
  if(file.exists("/home/glassfish/sqlite/")){
    sqlite.path <<- "/home/glassfish/sqlite/";  #public server + qiang local
  }else if(file.exists("/Users/xialab/Dropbox/sqlite/")){
    sqlite.path <<- "/Users/xialab/Dropbox/sqlite/"; #xia local
  }else if(file.exists("/Users/jeffxia/Dropbox/sqlite/")){
    sqlite.path <<- "/Users/jeffxia/Dropbox/sqlite/"; #xia local2
  }else if(file.exists("/home/zgy/sqlite/")){
    sqlite.path <<- "/home/zgy/sqlite/"; #zgy local
  }else if(file.exists("/home/ly/sqlite/")){
    sqlite.path <<- "/home/ly/sqlite/"; #ly local
  }else if(file.exists("/home/fiona/sqlite")){
    sqlite.path <<- "/home/fiona/sqlite/"; #fiona local
  }else if(file.exists("/Users/lzy/sqlite")){
    sqlite.path <<- "/Users/lzy/sqlite/"; #luyao local
  }else{
    sqlite.path <<-"/media/zzggyy/disk/sqlite/"; #zgy local
  }

  if(!.on.public.web) {
    sqlite.path <<- paste0(getwd(), "/");
    lib.path <<- "https://www.omicsnet.ca/OmicsNet/resources/data/";
    message("A dataset has been initiated.");
    return(dataSet);
  } else {
    .on.public.web <<- TRUE;
  }
  return(1)
}

SetUploadedGraph<- function(bool){
  uploadedGraph <<- bool;
}

SetOrganism <- function(org){
  data.org <<- org;
}

SetCurrentDataMulti <- function(){
  dataSet$type <- nms.vec;
  rm('nms.vec', envir = .GlobalEnv);
  return(.set.nSet(dataSet));
}


#' Process input list
#'
#' @param dataSetObj Input the name of the created dataSetObj (see Init.Data)
#' @param inputList Tab-delimited String of input omics list
#' @param org organism code: hsa, mmu, microbiome, rno, cel, dme, dre, sce, bta, gga
#' @param type character, type
#' @param queryType character, queryType
#'
#' @export
PrepareInputList <- function(dataSetObj="NA", inputList, org, type, queryType){

  if(!exists("dataSet")){
    Init.Data();
  }
  dataSet <- .get.nSet(dataSetObj);
  dataSet$name <- type;
  dataSet$orig <- inputList;

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

  if(type == "snp"){
    if(queryType == "region"){
          dataSet$snpRegion <- TRUE
     }else{
          dataSet$snpRegion <- FALSE
     }
  }

  if(sum(!na.inx) < 2){
    current.msg <<- "Less than two hits found in uniprot database. ";
    print(current.msg);
  }
  if(type == "met" || type == "mir" ){
    col1 <- 1;
    if(queryType == "mir_acc" || (type == "met" && queryType != "kegg")){
        col1 <- 2;
    }
    hit.inx <- match(tolower(rownames(prot.mat)), tolower(GeneAnotDB[,col1]));
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

  if(type %in% c("gene", "geneonly", "protein1")){
    if(length(dataSet$mat[["gene"]]) == 0){
      dataSet$mat[["gene"]] <- gene.mat;
      dataSet$exp.mat[["gene"]] <- prot.mat;
      dataSet$seed[["gene"]] <- prot.mat;
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
  }

  if(length(dataSet$type) == 1){
    dataSet$prot.mat <- rbind(dataSet$prot.mat, prot.mat);
    dataSet$seeds.proteins <- c(dataSet$seeds.proteins, seed.proteins);
  } else{
    for(i in 1:length(dataSet$exp.mat)){
      dataSet$prot.mat <- c(dataSet$prot.mat, dataSet$exp.mat[[i]]);
      dataSet$seeds.proteins <- c(dataSet$seeds.proteins, rownames(dataSet$exp.mat[[i]]));
    }
  }
  dataSet$seeds.proteins <- unique(dataSet$seeds.proteins)
  dataSet$seeds.expr <- as.matrix(dataSet$prot.mat);
  message(paste0("Preparing input list with the type of ", type, " completed."))

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

SetDbType <- function(dbType){
  dbu.type <<- dbType;
}


SetMirBo <- function(mirBo){
  requireMirExp <<- mirBo
}

SetAllOrgs <- function(orgs){
  dataSet$allOrgs <<- orgs
}

GetAllOrgs <- function(){
  return(dataSet$allOrgs);
}

GetInputTypes <- function(dataSetObj=NA){
  dataSet <- .get.nSet(dataSetObj);
  types.vec <- c(names(dataSet$exp.mat));
  return(types.vec);
}

GetInputNames <- function(dataSetObj=NA){
  dataSet <- .get.nSet(dataSetObj);
  types.vec <- c(names(dataSet$exp.mat));
  length.vec <- sapply(dataSet$exp.mat,function(x) {
            length(unique(rownames(x)))
                    })
  names.vec <- convertInputTypes2Names(types.vec);
  nms.and.size.vec <- vector();

  containsProtein <- F;
  if(length(rownames(dataSet$exp.mat[["gene"]])[dataSet$gene_type_vec == 3]) > 0){
    containsProtein <- T;
  }

  containsGene <- F;
  if(length(rownames(dataSet$exp.mat[["gene"]])[dataSet$gene_type_vec == 1]) > 0){
    containsGene <- T;
  }
  for(i in 1:length(names.vec)){
    nm <- names.vec[i]
    if(!containsGene && nm == "Gene"){
        nm <- "Protein";
    }else if(containsProtein && nm == "Gene"){
        nm <- "mRNA/Protein";
    }
    nms.and.size.vec[i] <- paste0(nm, " (", length.vec[i], ")");
  }
  dataSet$inputInfo <- nms.and.size.vec;
  dataSet <<- dataSet
  return(nms.and.size.vec);
}

convertInputTypes2Names <- function(types){
  typesCol <- c("gene","protein", "tf", "mir", "met", "peak", "mic", "snp", "ko");
  namesCol <- c("mRNA","Protein", "Transcription factor", "miRNA", "Metabolite", "LC-MS peak", "Microbial taxa", "SNP", "KO");
  db.map <- data.frame(type=typesCol, name=namesCol);
  hit.inx <- match(types, db.map[, "type"]);
  names <- db.map[hit.inx, "name"];
  mode(names) <- "character";
  na.inx <- is.na(names);
  names[na.inx] <- types[na.inx];
  return(names);
}

convertInteraction2Names <- function(types){
  typesCol <- c("gene", "tf", "mir", "met", "peak", "mic", "snp", "ko", "m2m");
  namesCol <- c("PPI", "TF-gene", "miRNA-gene", "Metabolite-protein", "Peak-metabolite", "Taxon-metabolite", "SNP Annotation", "ko", "Metabolite-metabolite");
  db.map <- data.frame(type=typesCol, name=namesCol);
  hit.inx <- match(types, db.map[, "type"]);
  #print(hit.inx);
  names <- db.map[hit.inx, "name"];
  mode(names) <- "character";
  return(names);
}

doEmblGene2EntrezMapping <- function(emblgene.vec){
  db.path <- paste(lib.path, data.org, "/entrez_embl_gene.rds", sep="");
  db.map <-  readRDS(db.path);
  hit.inx <- match(emblgene.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}

GetInputListStat <- function(listNm){
  geneList <- rownames(dataSet$exp.mat[["gene"]]);
  tfList <- rownames(dataSet$exp.mat[["tf"]]);
  mirList <- rownames(dataSet$exp.mat[["mir"]]);
  metList <- rownames(dataSet$exp.mat[["met"]]);
  peakList <- rownames(dataSet$exp.mat[["peak"]]);
  micList <- rownames(dataSet$exp.mat[["mic"]]);
  snpList <- rownames(dataSet$seed[["snp"]]);
  koList <- rownames(dataSet$exp.mat[["ko"]]);
  return(c(length(rownames(dataSet$exp.mat[["gene"]])),length(rownames(dataSet$exp.mat[["protein"]])),length(unique(tfList)),length(unique(mirList)),length(unique(metList)), length(unique(koList)), length(unique(peakList)), length(unique(micList)), length(unique(snpList))));
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
#'@param dataSetObj Input the name of the created dataSetObj (see Init.Data)
#'@export
GetRCommandHistory <- function(dataSetObj=NA){
  dataSet <- .get.nSet(dataSetObj);
  if(length(dataSet$cmdSet) == 0){
    return("No commands found");
  }
  return(dataSet$cmdSet);
}

#'Record R Commands
#'@param dataSetObj Input the name of the created dataSetObj (see Init.Data)
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

RemoveEdgeEntry <- function(tblnm, row.id) {
    inx <- which(rownames(dataSet$viewTable[tblnm][[1]]) == row.id);
    if(length(inx) > 0){
        dataSet$viewTable[tblnm][[1]] <- dataSet$viewTable[tblnm][[1]][-inx,];
        edgeu.res.list[[tblnm]] <- edgeu.res.list[[tblnm]][-inx,];
    }
    edgeu.res.list <<- edgeu.res.list;
    dataSet<<-dataSet
    return(1)
}

GetEdgeResRowNames <- function(netType){
    resTable <- dataSet$viewTable[[netType]]
    if(nrow(resTable) > max.row){
        resTable <- resTable[1:max.row, ];
        current.msg <<- "Due to computational constraints, only top 2000 rows will be displayed.";
    }
    rownames(resTable);
}

GetEdgeResCol <- function(netType, colInx){

    resTable <- dataSet$viewTable[[netType]]

    if(nrow(resTable) > max.row){
        resTable <- resTable[1:max.row, ];
        #current.msg <<- "Due to computational constraints, only top 1000 rows will be displayed.";
    }
    res <- resTable[, colInx];
    hit.inx <- is.na(res) | res == ""; # note, must use | for element-wise operation
    res[hit.inx] <- "N/A";
    return(res);
}


# batch remove based on
UpdateEdgeTableEntries <- function(table.nm,table.type, col.id, method, value, action) {

    col <- dataSet$viewTable[[table.nm]][,col.id];

    if(method == "contain"){
        hits <- grepl(value, col, ignore.case = TRUE);
    }else if(method == "match"){
        hits <- tolower(col) %in% tolower(value);
    }

    if(action == "keep"){
        hits = !hits;
    }

    if(sum(hits) > 0){
        row.ids <- rownames(dataSet$viewTable[[table.nm]])[hits];
        dataSet$viewTable[[table.nm]] <- dataSet$viewTable[[table.nm]][!hits,];
        edgeu.res.list[[table.nm]]$table <-  dataSet$viewTable[[table.nm]]
        UpdateModifiedTable(edgeu.res.list);
        fast.write.csv(dataSet$viewTable[[table.nm]], file="omicsnet_edgelist.csv", row.names=FALSE);
        .set.nSet(dataSet);
        return(row.ids);
    }else{
        return("NA");
    }
}

UpdateModifiedTable <- function(edgeu.res.list) {
  edge.res <- data.frame();

  if(length(edgeu.res.list) > 0){
    for(i in 1:length(edgeu.res.list)){
      edge.res <- rbind(edge.res, edgeu.res.list[[i]]$table);
    }
  }

  omics.net$edge.data <- unique(edge.res);
  omics.net <<- omics.net;
  #print(omics.net$edge.data);

  return(1)
}

GetOrg <- function(){
  return(data.org);
}

PrepareSqliteDB <- function(sqlite_Path, onweb = TRUE) {
  if(onweb) {return(TRUE)};
  if(file.exists(sqlite_Path)) {return(TRUE)};

  dbNM <- basename(sqlite_Path);
  DonwloadLink <- paste0("https://www.xialab.ca/resources/sqlite/", dbNM);
  download.file(DonwloadLink, sqlite_Path);
  return(TRUE)
}

#importFrom("grDevices", "col2rgb", "colorRampPalette", "dev.off","hsv", "rainbow")
#importFrom("stats", "aggregate", "complete.cases", "dnorm", "filter","formula", "lm", "median", "na.omit", "p.adjust", "phyper","quantile", "wilcox.test")
#importFrom("utils", "capture.output", "install.packages","installed.packages", "object.size", "read.csv","read.table", "sessionInfo", "write.csv")


SetFunctionByDbVersion <- function(db_version) {
  # Defined according to DB update
  func_names <- list(
    'mir2gene.sqlite' = c('doProteinIDMapping', 'QueryMirSQLite'),
    'omicsnet_met.sqlite' = c('QueryMetSQLiteNet', 'QueryM2mSQLiteNet', 'QueryMicM2mSQLiteNet'),
    'ppi.sqlite' = c('QueryPpiSQLite'),
    'tf2gene.sqlite' = c('QueryTFSQLite')
  )
  version_suffix <- '_2022'

  for (db_name in names(func_names)) {
    funcs <- func_names[[db_name]]
    for (func in funcs) {
      if (db_version == "previousDB") {
        SetToCurrent(func, db_name, version_suffix)
      } else {
        SetToPrevious(func, db_name, version_suffix)
      }
    }
  }
}

SetToCurrent <- function(func_name, db_name, version_suffix) {
  db_name_versioned <- sub("(\\.sqlite)$", paste0(version_suffix, "\\1"), db_name)

  original_func <- get(func_name, envir = .GlobalEnv)
  original_body <- body(original_func)

  modified_body <- lapply(original_body, function(line) {
    if (is.call(line) && any(grepl(db_name_versioned, deparse(line)))) {
      modified_line <- gsub(db_name_versioned, db_name, deparse(line))
      return(parse(text = modified_line)[[1]])
    } else {
      return(line)
    }
  })

  modified_func <- original_func
  body(modified_func) <- as.call(modified_body)
  assign(func_name, modified_func, envir = .GlobalEnv)
}

SetToPrevious <- function(func_name, db_name, version_suffix) {
  db_name_versioned <- sub("(\\.sqlite)$", paste0(version_suffix, "\\1"), db_name)

  original_func <- get(func_name, envir = .GlobalEnv)
  original_body <- body(original_func)

  modified_body <- lapply(original_body, function(line) {
    if (is.call(line) && any(grepl(db_name, deparse(line)))) {
      modified_line <- gsub(db_name, db_name_versioned, deparse(line))
      return(parse(text = modified_line)[[1]])
    } else {
      return(line)
    }
  })

  modified_func <- original_func
  body(modified_func) <- as.call(modified_body)
  assign(func_name, modified_func, envir = .GlobalEnv)
}