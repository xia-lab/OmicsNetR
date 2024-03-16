##################################################
## R script for OmicsNet
## Description: ID mapping, GO/Pathway enrichment analysis
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

LoadLib <- function(fun.type){
  fun.type <- tolower(fun.type);
  if(fun.type %in% c("kegg", "reactome", "motif", "bp", "cc", "mf","panth_mf", "go_bp", "go_cc", "go_panthbp", "go_panthcc", "go_panthmf","mirfamily", "mircluster", "mirdisease","mirtissue", "mirfunction")){
    LoadEnrLib("standard", fun.type);
  }else if(fun.type == 'integ'){
    LoadKEGGLibOther("integ");
  }else if(fun.type == 'keggm'){
    LoadKEGGLibOther("keggm");
  }else if(fun.type %in% c("pathway", "blood", "urine", "csf", "predicted", "location", "drug")){
    LoadEnrLib("mset", fun.type);
  }else if(fun.type %in% c("mir", "jasper", "drug", "encode")){
    LoadEnrLib("other", fun.type);
  }else{
    print(paste("Unknown lib option:", fun.type));
    return(0);
  }
}

LoadEnrLib <- function(type, subtype){
  if (type == "msets"){
    path <- paste(lib.path, "msets", "/", tolower(subtype), ".rds", sep="");
  }else if (type == "other"){
    path <- paste(lib.path, data.org, "/", tolower(paste0(subtype, "_enr")), ".rds", sep="");
  }else{
    path <- paste(lib.path, data.org, "/", tolower(subtype), ".rds", sep="");
  }
  res <- readRDS(path);
  if(!.on.public.web){
    nmdb <- basename(go.path);
    download.file(go.path, destfile = nmdb, method="libcurl", mode = "wb");
    path <- nmdb;
  }
  if(is.null(names(res))){ # new go lib does not give names
    names(res) <- c("link", "term", "sets");
  }
  set <- readRDS(path)

  #some sets list do not have names
  if("sets" %in% names(set)){
    setsInx <- which(names(set) == "sets");
    termInx <- which(names(set) == "term");
    linkInx <- which(names(set) == "link");
  }else{
    setsInx <- 3;
    termInx <- 2;
    linkInx <- 1;
  }
  current.geneset <- set[[setsInx]];

  set.ids<- names(current.geneset);
  names(current.geneset) <- set[[termInx]];
  current.geneset <<- current.geneset

  current.setlink <<- set[[linkInx]];
  current.setids <<- set.ids;
  current.universe <<- unique(unlist(current.geneset));
}


SearchReg <- function(file.nm, fun.type, IDs, count){
  regCount <<- count;
  regBool <<- "true";
  res <- PerformNetEnrichment(file.nm, fun.type, IDs);
}

regEnrichment <- function(file.nm, fun.type, IDs, netInv){
  regBool <<- "false";
  res <- PerformNetEnrichment(file.nm, fun.type, IDs);
}

#' Perform gene enrichment analysis or identify gene regulatory targets
#'
#' @param file.nm File name of result table to be exported in csv format, do not include file extension
#' @param fun.type Enrichment database type
#' @param IDs String of ids to be tested separated by "; "
#'
#' @export
#'
PerformNetEnrichment <- function(file.nm, fun.type, IDs){
  # note: hit.query, resTable must synchronize
  # prepare query
  ora.vec <- NULL;
  idtype <- "entrez"
  ora.vec <- unlist(strsplit(IDs, "; "));
  names(ora.vec) <- ora.vec;
  if(fun.type %in% c("trrust", "encode", "jaspar", "mirnet", "met", "drugbank")){
    netInv <- "inverse";
    res <- PerformRegEnrichAnalysis(file.nm, fun.type, ora.vec, netInv, idtype);
  } else{
    res <- PerformEnrichAnalysis(file.nm, fun.type, ora.vec);
  }
  return(res);
}

PerformRegEnrichAnalysis <- function(file.nm, fun.type, ora.vec, netInv, idType){
  require("plyr")
  ora.nms <- names(ora.vec);
  # prepare for the result table
  set.size<-100;

  if (fun.type %in% c("chea", "encode", "jaspar", "trrust")){
    table.nm <- paste(data.org, fun.type, sep="_");
    res <- QueryTFSQLite(table.nm, ora.vec, netInv);
    # no hits
    if(nrow(res)==0){ return(c(0,0)); }
    edge.res <- data.frame(gene=res[,"entrez"], symbol=res[,"symbol"],id=res[,"tfid"], name=res[,"tfname"]);
    node.ids <- c(res[,"entrez"], res[,"tfid"]);
    node.nms <- c(res[,"symbol"], res[,"tfname"]);

  }else if(fun.type == "mirnet"){ # in miRNA, table name is org code, colname is id type
    table.nm <- data.org
    res <- QueryMirSQLite(table.nm, "entrez", ora.vec, "inverse", data.org);
    if(nrow(res)==0){ return(c(0,0)); }
    edge.res <- data.frame(gene=res[,"entrez"], symbol=res[,"symbol"], id=res[,"mir_acc"], name=res[,"mir_id"] );
    node.ids <- c(res[,"entrez"], res[,"mir_acc"])
    node.nms <- c(res[,"symbol"], res[,"mir_id"]);

  }else if(fun.type == "met"){
    table.nm <- paste(data.org, "kegg", sep="_");
    res <- QueryMetSQLiteNet(table.nm, ora.vec, "inverse");
    if(nrow(res)==0){ return(c(0,0)); }
    edge.res <- data.frame(gene=res[,"entrez"], symbol=res[,"symbol"], id=res[,"kegg"], name=res[,"met"] );
    node.ids <- c(res[,"entrez"], res[,"kegg"])
    node.nms <- c(res[,"symbol"], res[,"met"]);

  }else{
    table.nm <- paste("drug", data.org, sep="_");
    ora.vec <- doEntrez2UniprotMapping(ora.vec);
    res <- QueryDrugSQLite(ora.vec);
    if(nrow(res)==0){ return(c(0,0)); }
    edge.res <- data.frame(gene=doUniprot2EntrezMapping(res[,"upid"]), symbol=res[,"symbol"], id=res[,"dbid"], name=res[,"dbname"] );
    node.ids <- c(doUniprot2EntrezMapping(res[,"upid"]), res[,"dbid"])
    node.nms <- c(res[,"symbol"], res[,"dbname"]);
  }

  edge.res$mix <- paste0(edge.res[,1], edge.res[,3]);
  edge.res <- edge.res[!duplicated(edge.res$mix),];

  row.names(edge.res) <- 1:nrow(edge.res);
  fast.write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);

  hit.freq <- count(edge.res, 'id')
  hit.freq <- hit.freq[order(hit.freq$freq, decreasing = TRUE),]
  hit.freqnm <- count(edge.res, 'name')
  hit.freqnm <- hit.freqnm[order(hit.freqnm$freq, decreasing = TRUE),]
  hit.freq$name <- hit.freqnm$name
  hits.gene <- list()

  if(idType == "uniprot"){
    for(i in 1:nrow(hit.freq)){
      df <- edge.res[which(edge.res$id == hit.freq$id[i]),];
      gene.vec <- as.vector(df$gene);
      gene.vec <- doEntrez2UniprotMapping(gene.vec)
      hits.gene[[i]] <- gene.vec;
    }
  } else if(idType == "entrez") {
    for(i in 1:nrow(hit.freq)){
      df <- edge.res[which(edge.res$id == hit.freq$id[i]),];
      gene.vec <- as.vector(df$gene);
      hits.gene[[i]] <- gene.vec;
    }
  }else {
    for(i in 1:nrow(hit.freq)){
      df <- edge.res[which(edge.res$id == hit.freq$id[i]),];
      gene.vec <- as.vector(df$gene);
      gene.vec <- doEntrez2EmblProteinMapping(gene.vec)
      hits.gene[[i]] <- gene.vec;
    }
  }
  names(hits.gene) = hit.freq$name;

  resTable1 = data.frame(hit.freq$id, hit.freq$name, hit.freq$freq)
  colnames(resTable1) = c("Id", "Name", "Hits")
  resTable1 = resTable1[order(-resTable1$Hits),]
  current.msg <<- "Regulatory element enrichment analysis was completed";

  if(regBool == "true"){
    resTable1 <- resTable1[which(resTable1$Hits == regCount),]
  }
  # write json
  require("RJSONIO");
  fun.ids <- resTable1[,1];
  fun.nms <- resTable1[,2];
  fun.hits <- resTable1[,3];
  json1.res <- list(
    fun.ids = fun.ids,
    fun.nms = fun.nms,
    fun.hits = fun.hits,
    fun.genes = hits.gene
  );
  json.mat <- toJSON(json1.res);
  json.nm <- paste(file.nm, ".json", sep="");
  sink(json.nm)
  cat(json.mat);
  sink();

  resTable1$Features = hits.gene
  # write csv
  csv.nm <- paste(file.nm, ".csv", sep="");
  fast.write.csv(resTable1, file=csv.nm, row.names=F);

  #record table for report
  type = "regNetwork";
  dataSet$imgSet$enrTables[[type]] <- list()
  dataSet$imgSet$enrTables[[type]]$table <- resTable1;
  dataSet$imgSet$enrTables[[type]]$library <- fun.type
  dataSet$imgSet$enrTables[[type]]$algo <- "Overrepresentation Analysis"
  dataSet <<- dataSet

  return(1);
}

# note: hit.query, resTable must synchronize
# ora.vec should contains entrez ids, named by their gene symbols
PerformEnrichAnalysis <- function(file.nm, fun.type, ora.vec){
  # prepare lib
  LoadLib(fun.type);

  # prepare query
  ora.nms <- names(ora.vec);

  # prepare for the result table
  set.size<-length(current.geneset);
  res.mat<-matrix(0, nrow=set.size, ncol=5);
  rownames(res.mat)<-names(current.geneset);
  colnames(res.mat)<-c("Total", "Expected", "Hits", "P.Value", "FDR");

  # need to cut to the universe covered by the pathways, not all genes
  if(tolower(fun.type) %in% c("chea", "jaspar", "encode", "mir")){
    current.universe <- toupper(current.universe)
  }
  hits.inx <- ora.vec %in% current.universe;
  ora.vec <- ora.vec[hits.inx];
  ora.nms <- ora.nms[hits.inx];

  q.size<-length(ora.vec);

  # get the matched query for each pathway
  hits.query <- lapply(current.geneset,
                       function(x) {
                         ora.nms[ora.vec%in%unlist(x)];
                       }
  );
  hits.query <- lapply(hits.query, function(x){unique(x)});
  names(hits.query) <- names(current.geneset);
  hit.num<-unlist(lapply(hits.query, function(x){length(x)}), use.names=FALSE);

  # total unique gene number
  uniq.count <- length(current.universe);

  # unique gene count in each pathway
  set.size <- unlist(lapply(current.geneset, length));

  res.mat[,1]<-set.size;
  res.mat[,2]<-q.size*(set.size/uniq.count);
  res.mat[,3]<-hit.num;

  # use lower.tail = F for P(X>x)
  raw.pvals <- phyper(hit.num-1, set.size, uniq.count-set.size, q.size, lower.tail=F);
  res.mat[,4]<- raw.pvals;
  res.mat[,5] <- p.adjust(raw.pvals, "fdr");

  # now, clean up result, synchronize with hit.query
  res.mat <- res.mat[hit.num>0,,drop = F];
  hits.query <- hits.query[hit.num>0];

  if(nrow(res.mat)> 1){
    # order by p value
    ord.inx<-order(res.mat[,4]);
    res.mat <- signif(res.mat[ord.inx,],3);
    hits.query <- hits.query[ord.inx];

    imp.inx <- res.mat[,4] <= 0.05;
    if(sum(imp.inx) < 10){ # too little left, give the top ones
      topn <- ifelse(nrow(res.mat) > 10, 10, nrow(res.mat));
      res.mat <- res.mat[1:topn,];
      hits.query <- hits.query[1:topn];
    }else{
      res.mat <- res.mat[imp.inx,];
      hits.query <- hits.query[imp.inx];
      if(sum(imp.inx) > 120){
        # now, clean up result, synchronize with hit.query
        res.mat <- res.mat[1:120,];
        hits.query <- hits.query[1:120];
      }
    }
  }

  #get gene symbols
  resTable <- data.frame(Pathway=rownames(res.mat), res.mat);
  current.msg <<- "Functional enrichment analysis was completed";

  # write json
  require("RJSONIO");
  fun.anot <- hits.query;
  fun.padj <- resTable[,6]; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
  fun.pval <- resTable[,5]; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
  hit.num <- resTable[,4]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
  fun.ids <- as.vector(current.setids[names(fun.anot)]);
  if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };
  json.res <- list(
    fun.link = current.setlink[1],
    fun.anot = fun.anot,
    fun.ids = fun.ids,
    fun.pval = fun.pval,
    fun.padj = fun.padj,
    hit.num = hit.num
  );
  json.mat <- toJSON(json.res);
  json.nm <- paste(file.nm, ".json", sep="");
  sink(json.nm)
  cat(json.mat);
  sink();

  #record table for report
  type = "network";
  dataSet$imgSet$enrTables[[type]] <- list()
  dataSet$imgSet$enrTables[[type]]$table <- resTable;
  dataSet$imgSet$enrTables[[type]]$library <- fun.type
  dataSet$imgSet$enrTables[[type]]$algo <- "Overrepresentation Analysis"
  dataSet <<- dataSet

  # write csv
  csv.nm <- paste(file.nm, ".csv", sep="");
  a <- lapply(fun.anot,function(x){paste(x,collapse = "; ")})
  resTable$ids <- unlist(a);
  fast.write.csv(resTable, file=csv.nm, row.names=F);

  return(1);
}

doProteinIDMapping <- function(q.vec, type, dbType = "NA"){
  if(type %in% c("rsid")){
    hit.inx <- startsWith(q.vec, "rs");
    q.vec <- q.vec[hit.inx]
    res <- data.frame(gene_id=q.vec, accession=q.vec);
    entrezs <- res;
  }else if(type %in% c("region")){
    hit.inx <- grepl("\\:",q.vec);
    q.vec <- q.vec[hit.inx]
    res <- data.frame(gene_id=q.vec, accession=q.vec);
    entrezs <- res;
  }else if(type %in% c("class", "family","order" ,"genus","species","strain", "phylum")){
    require('RSQLite');
    mic.taxa <<- type;
    db.path <- paste(lib.path, "microbiome", "/taxaInfo.rds", sep="");

    if(!.on.public.web){
      nm <- basename(db.path);
      download.file(db.path, destfile = nm, method="libcurl", mode = "wb");
      db.path <- nm;
    }

    taxlist <- readRDS(db.path);
    db.map <- taxlist[[type]]
    db.map[, type] <- gsub("\\[|\\]","", db.map[, type])
    db.map[, type] <- gsub("_"," ", db.map[, type])
    uncls.idx =  c(grep("unclassified",q.vec), grep("Unclassified",q.vec))
    q.vec[uncls.idx] = gsub("unclassified","", q.vec[uncls.idx])
    q.vec[uncls.idx] = gsub("Unclassified","", q.vec[uncls.idx])
    q.vec[uncls.idx] = gsub("^\\s+|\\s+$", "", q.vec[uncls.idx])
    q.vec[uncls.idx] = paste0("unclassified ",q.vec[uncls.idx])
    db.map <- db.map[db.map[, type] %in% q.vec,];
    db.map <- data.frame(a=db.map[, type], b=db.map[,type]);

    hit.inx <- match(q.vec, db.map[, "a"]);
    entrezs <- db.map[hit.inx, ]
    entrezs <- entrezs[c(1,2)]
    colnames(entrezs) = c("gene_id", "accession");
  }else if(type %in% c("entrez", "ctd", "drugbank","ko")){
    # need to get only our data
    if(type == "entrez"){
      db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");
    }else if(type == "ctd"){
      db.path <- paste(lib.path, "chem", "/ctd.rds", sep="")
    }else if(type == "drugbank"){
      db.path <- paste(lib.path, "drug", "/drugbank.rds", sep="")
    }else if(type == "ko"){
      db.path <- paste(lib.path, "microbiome", "/ko.rds", sep="")
    }

    if(!.on.public.web){
      nm <- basename(db.path);
      download.file(db.path, destfile = nm, method="libcurl", mode = "wb");
      db.path <- nm;
    }

    db.map <-  readRDS(db.path);
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
    entrezs <- db.map[hit.inx, ]
    entrezs <- entrezs[c(1,2)]
    colnames(entrezs) = c("gene_id", "accession");
  }else if(type %in% c("mir_acc", "mir_id", "mirnet")){
    require('RSQLite');
    path <- paste0(sqlite.path, "mir2gene.sqlite");
    if(!PrepareSqliteDB(path, .on.public.web)){
      stop("Sqlite database is missing, please check your internet connection!");
    }

    mir.db <- dbConnect(SQLite(), path);
    query <- paste (shQuote(q.vec),collapse=",");
    table.nm <- data.org
    statement <- paste("SELECT * FROM ", table.nm, " WHERE ", type," COLLATE NOCASE IN (",query,")", sep="");

    mirtable <- dbSendQuery(mir.db, statement);
    mir.dic <- fetch(mirtable, n=-1);
    if(nrow(mir.dic) == 0){
      return(0);
    }
    entrezs <- mir.dic[c("mir_id", "mir_acc")];
    entrezs <- entrezs[!duplicated(entrezs[,"mir_id"]),]
    rownames(entrezs) = seq.int(nrow(entrezs));
    entrezs <- data.frame(lapply(entrezs, as.character), stringsAsFactors=FALSE)
    colnames(entrezs) = c("gene_id", "accession");
  }else if(type == "symbol"){
    db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");

    if(!.on.public.web){
      nm <- basename(db.path);
      download.file(db.path, destfile = nm, method="libcurl", mode = "wb");
      db.path <- nm;
    }

    gene.map <- readRDS(db.path);
    if(data.org == "hsa"){
      q.vec = toupper(q.vec);
    }
    hit.inx <- match(q.vec, gene.map[, "symbol"]);
    entrezs <- gene.map[hit.inx, ];
    entrezs = entrezs[c(1,2)];
    colnames(entrezs) <- c("gene_id", "accession");
  }else if(type %in% c("meta", "kegg", "chebi", "name", "bigg", "pubchem", "hmdb")){

    db.path <- paste(lib.path, "lib/compound_db.rds", sep="");

    if(!.on.public.web){
      nm <- basename(db.path);
      download.file(db.path, destfile = nm, method="libcurl", mode = "wb");
      db.path <- nm;
    }

    cmpd.map <- readRDS(db.path)
    q.type <- type;
    cmpd.vec <- q.vec

    if(q.type == "hmdb"){
      hit.inx <- match(tolower(cmpd.vec), tolower(cmpd.map$hmdb_id));
    }else if(q.type == "kegg"){
      hit.inx <- match(tolower(cmpd.vec), tolower(cmpd.map$kegg_id));
    }else if(q.type == "pubchem"){
      hit.inx <- match(tolower(cmpd.vec), tolower(cmpd.map$pubchem_id));
    }else if(q.type == "chebi"){
      hit.inx <- match(tolower(cmpd.vec), tolower(cmpd.map$chebi_id));
    }else if(q.type == "reactome"){
      hit.inx <- match(tolower(cmpd.vec), tolower(cmpd.map$reactome));
    }else{
      print("No support for this compound database");
      return(0);
    }

    if(q.type != "kegg"){
      typ = paste0(q.type, "_id")
      res_entrez <-  cmpd.map[hit.inx, c("kegg_id", typ)];
    }else{
      res_entrez <-  cmpd.map[hit.inx, c("kegg_id", "kegg_id")];
    }
    colnames(res_entrez) <- c("gene_id", "accession")
    entrezs <- res_entrez
  }else if(type == "tf"){
    table.nm <- paste(data.org, dbType, sep="_");
    mir.dic <- QueryTfSQLite(table.nm, q.vec, "inverse");
    res <- mir.dic[ , c("tfid", "tfid")];
    res = res[!duplicated(res),]

    colnames(res) <- c("accession", "gene_id")
    hit.inx <- match(q.vec, res[, "accession"]);
    entrezs <- res[hit.inx, ];
    entrezs = res[c(1,2)];
  }else {
    if(type == "gb"){
      # note, some ID can have version number which is not in the database
      # need to strip it off NM_001402.5 => NM_001402
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      db.path <- paste(lib.path, data.org, "/entrez_gb.rds", sep="");
    }else if(type == "refseq"){
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      db.path <- paste(lib.path, data.org, "/entrez_refseq.rds", sep="");
    }else if(type == "embl_gene"){
      db.path <- paste(lib.path, data.org, "/entrez_embl_gene.rds", sep="");
    }else if(type == "embl_transcript"){
      db.path <- paste(lib.path, data.org, "/entrez_embl_transcript.rds", sep="");
    }else if(type == "embl_protein"){
      db.path <- paste(lib.path, data.org, "/entrez_embl_protein.rds", sep="");
    }else if(type == "orf"){ # only for yeast
      db.path <- paste(lib.path, data.org, "/entrez_orf.rds", sep="");
    }else if(type == "string"){
      db.path <- paste(lib.path, data.org, "/entrez_string.rds", sep="")
    }else if(type == "ecogene"){ # only for ecoli
      db.path <- paste(lib.path, data.org, "/entrez_ecogene.rds", sep="")
    }else if(type == "uniprot"){
      db.path <- paste(lib.path, data.org, "/entrez_uniprot.rds", sep="")
    }else if(type == "flybase"){
      db.path <- paste(lib.path, data.org, "/entrez_flybase.rds", sep="")
    }else{
      print(type);
      print("Unknown data type1");
      return(0);
    }

    if(!.on.public.web){
      nm <- basename(db.path);
      download.file(db.path, destfile = nm, method="libcurl", mode = "wb");
      db.path <- nm;
    }

    db.map <-  readRDS(db.path);
    hit.inx <- match(q.vec, db.map[, "accession"]);
    entrezs <- db.map[hit.inx, ];
  }

  return(entrezs);
}


# mapping between genebank, refseq and entrez
doGeneIDMapping <- function(q.vec, type){
  if(is.null(q.vec)){
    db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");
    db.map <-  readRDS(db.path);
    q.vec <- db.map[, "gene_id"];
    type = "entrez";
  }

  if(type == "symbol"){
    db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");

    if(!.on.public.web){
      nm <- basename(db.path);
      download.file(db.path, destfile = nm, method="libcurl", mode = "wb");
      db.path <- nm;
    }

    db.map <-  readRDS(db.path);
    hit.inx <- match(q.vec, db.map[, "symbol"]);
  }else if(type == "entrez"){
    db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");

    if(!.on.public.web){
      nm <- basename(db.path);
      download.file(db.path, destfile = nm, method="libcurl", mode = "wb");
      db.path <- nm;
    }

    db.map <-  readRDS(db.path);
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
  }else{
    # note, some ID can have version number which is not in the database
    # need to strip it off NM_001402.5 => NM_001402
    q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
    q.vec <- q.mat[,1];
    if(type == "gb"){
      db.path <- paste(lib.path, data.org, "/entrez_gb.rds", sep="");
    }else if(type == "refseq"){
      db.path <- paste(lib.path, data.org, "/entrez_refseq.rds", sep="");
    }else if(type == "embl_gene"){
      db.path <- paste(lib.path, data.org, "/entrez_embl_gene.rds", sep="");
    }else if(type == "embl_transcript"){
      db.path <- paste(lib.path, data.org, "/entrez_embl_transcript.rds", sep="");
    }else if(type == "embl_protein"){
      db.path <- paste(lib.path, data.org, "/entrez_embl_protein.rds", sep="");
    }else if(type == "orf"){ # only for yeast
      db.path <- paste(lib.path, data.org, "/entrez_orf.rds", sep="");
    }else{
      print("Unknown data type2");
      return(0);
    }

    if(!.on.public.web){
      nm <- basename(db.path);
      download.file(db.path, destfile = nm, method="libcurl", mode = "wb");
      db.path <- nm;
    }

    db.map <-  readRDS(db.path);
    hit.inx <- match(q.vec, db.map[, "accession"]);
  }
  entrezs=db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  rm(db.map, q.vec); gc();
  return(entrezs);
}

doEntrez2SymbolMapping <- function(entrez.vec){
  db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");

  if(!.on.public.web){
    nm <- basename(db.path);
    download.file(db.path, destfile = nm, method="libcurl", mode = "wb");
    db.path <- nm;
  }

  gene.map <- readRDS(db.path);

  hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
  symbols <- gene.map[hit.inx, "symbol"];

  # if not gene symbol, use id by itself
  na.inx <- is.na(symbols);
  symbols[na.inx] <- entrez.vec[na.inx];
  return(symbols);
}


doKegg2NameMapping <- function(entrez.vec){
  db.path <- paste(lib.path, "lib/compound_db.rds", sep="");

  if(!.on.public.web){
    nmdb <- basename(db.path);
    download.file(db.path, destfile = nmdb, method="libcurl", mode = "wb");
    db.path <- nmdb;
  }

  gene.map <- readRDS(db.path);

  hit.inx <- match(entrez.vec, gene.map[, "kegg_id"]);
  symbols <- gene.map[hit.inx, "name"];
  symbols <- as.character(symbols)
  # if not gene symbol, use id by itself
  na.inx <- is.na(symbols);
  symbols[na.inx] <- entrez.vec[na.inx];
  return(symbols);
}

doPubchem2NameMapping <- function(entrez.vec){
  db.path <- paste(lib.path, "/lib/pubchem_lib.qs", sep="");

  if(!.on.public.web){
    nmdb <- basename(db.path);
    download.file(db.path, destfile = nmdb, method="libcurl", mode = "wb");
    db.path <- nmdb;
  }

  full.map <- qs::qread(db.path);

  hit.inx <- match(entrez.vec, full.map[, "accession"]);
  symbols <- full.map[hit.inx, "name"];
  symbols <- as.character(symbols)
  # if not gene symbol, use id by itself
  na.inx <- is.na(symbols);
  symbols[na.inx] <- entrez.vec[na.inx];
  return(symbols);
}

doHMDB2NameMapping <- function(entrez.vec){
  db.path <- paste(lib.path, "/lib/hmdb_lib.qs", sep="");

  if(!.on.public.web){
    nmdb <- basename(db.path);
    download.file(db.path, destfile = nmdb, method="libcurl", mode = "wb");
    db.path <- nmdb;
  }

  full.map <- qs::qread(db.path);

  hit.inx <- match(entrez.vec, full.map[, "accession"]);
  symbols <- full.map[hit.inx, "name"];
  symbols <- as.character(symbols)
  # if not gene symbol, use id by itself
  na.inx <- is.na(symbols);
  symbols[na.inx] <- entrez.vec[na.inx];
  return(symbols);
}

doHMDB2KEGGMapping <- function(entrez.vec){
  db.path <- paste(lib.path, "lib/compound_db.rds", sep="");

  if(!.on.public.web){
    nmdb <- basename(db.path);
    download.file(db.path, destfile = nmdb, method="libcurl", mode = "wb");
    db.path <- nmdb;
  }

  gene.map <- readRDS(db.path);

  hit.inx <- match(entrez.vec, gene.map[, "hmdb_id"]);
  symbols <- gene.map[hit.inx, "kegg_id"];
  symbols <- as.character(symbols);
  # if not gene symbol, use id by itself
  symbols[is.na(symbols)] <- entrez.vec[is.na(symbols)];
  symbols[symbols == ""] <- entrez.vec[symbols == ""];
  return(symbols)
}

doPubchem2KEGGMapping <- function(entrez.vec){
  db.path <- paste(lib.path, "/lib/pubchem_lib.qs", sep="");

  if(!.on.public.web){
    nmdb <- basename(db.path);
    download.file(db.path, destfile = nmdb, method="libcurl", mode = "wb");
    db.path <- nmdb;
  }

  full.map <- qs::qread(db.path);

  db.path <- paste(lib.path, "lib/compound_db.rds", sep="");

  if(!.on.public.web){
    nmdb <- basename(db.path);
    download.file(db.path, destfile = nmdb, method="libcurl", mode = "wb");
    db.path <- nmdb;
  }

  gene.map <- readRDS(db.path);

  hit.inx <- match(entrez.vec, full.map[, "accession"]);
  symbols1 <- full.map[hit.inx, "KEGG"];
  symbols1 <- as.character(symbols1);
  if(any(symbols1 == "")){
    pcms <- entrez.vec[symbols1 == ""];
    entrez.vec2 <- gsub("PBCM0+","",pcms);
    hit.inx <- match(entrez.vec2, gene.map[, "pubchem_id"]);
    symbols2 <- gene.map[hit.inx, "kegg_id"];
    symbols2 <- as.character(symbols2);
  }
  symbols <- c(symbols2, symbols1)
  # if not gene symbol, use id by itself
  symbols[is.na(symbols)] <- entrez.vec[is.na(symbols)];
  symbols[symbols == ""] <- entrez.vec[symbols == ""];
  return(symbols);
}

# note, entrez.vec could contain NA/null, cannot use rownames
doEntrezIDAnot <- function(entrez.vec){
  db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");

  if(!.on.public.web){
    nmdb <- basename(db.path);
    download.file(db.path, destfile = nmdb, method="libcurl", mode = "wb");
    db.path <- nmdb;
  }

  gene.map <- readRDS(db.path);

  hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
  anot.mat <- gene.map[hit.inx, c("gene_id", "symbol", "name")];

  na.inx <- is.na(hit.inx);
  anot.mat[na.inx, "symbol"] <- entrez.vec[na.inx];
  anot.mat[na.inx, "name"] <- 'NA';
  return(anot.mat);
}

doUniprot2EntrezMapping <- function(uniprot.vec){
  db.path <- paste(lib.path, data.org, "/entrez_uniprot.rds", sep="");
  if(!.on.public.web){
    nmdb <- basename(db.path);
    download.file(db.path, destfile = nmdb, method="libcurl", mode = "wb");
    db.path <- nmdb;
  }
  db.map <-  readRDS(db.path);
  hit.inx <- match(uniprot.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}

doEntrez2UniprotMapping <- function(entrez.vec){
  db.path <- paste(lib.path, data.org, "/entrez_uniprot.rds", sep="");
  if(!.on.public.web){
    nmdb <- basename(db.path);
    download.file(db.path, destfile = nmdb, method="libcurl", mode = "wb");
    db.path <- nmdb;
  }
  db.map <-  readRDS(db.path);
  hit.inx <- match(entrez.vec, db.map[, "gene_id"]);
  unips <- db.map[hit.inx, "accession"];
  mode(unips) <- "character";
  return(unips);
}

doString2EntrezMapping <- function(string.vec){
  db.path <- paste(lib.path, data.org, "/entrez_string.rds", sep="");
  if(!.on.public.web){
    nmdb <- basename(db.path);
    download.file(db.path, destfile = nmdb, method="libcurl", mode = "wb");
    db.path <- nmdb;
  }
  db.map <-  readRDS(db.path);
  hit.inx <- match(string.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}

doEmblGene2EntrezMapping <- function(emblgene.vec){
  db.path <- paste(lib.path, data.org, "/entrez_embl_gene.rds", sep="");
  if(!.on.public.web){
    nmdb <- basename(db.path);
    download.file(db.path, destfile = nmdb, method="libcurl", mode = "wb");
    db.path <- nmdb;
  }
  db.map <-  readRDS(db.path);
  hit.inx <- match(emblgene.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}

doSymbol2EntrezMapping <- function(symbol.vec){
  db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");
  if(!.on.public.web){
    nmdb <- basename(db.path);
    download.file(db.path, destfile = nmdb, method="libcurl", mode = "wb");
    db.path <- nmdb;
  }
  db.map <-  readRDS(db.path);
  hit.inx <- match(symbol.vec, db.map[, "symbol"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}

doEmblProtein2EntrezMapping <- function(emblprotein.vec){
  db.path <- paste(lib.path, data.org, "/entrez_embl_protein.rds", sep="");
  if(!.on.public.web){
    nmdb <- basename(db.path);
    download.file(db.path, destfile = nmdb, method="libcurl", mode = "wb");
    db.path <- nmdb;
  }
  db.map <-  readRDS(db.path);
  hit.inx <- match(emblprotein.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}

doEntrez2EmblProteinMapping <- function(entrez.vec){
  db.path <- paste(lib.path, data.org, "/entrez_embl_protein.rds", sep="");
  if(!.on.public.web){
    nmdb <- basename(db.path);
    download.file(db.path, destfile = nmdb, method="libcurl", mode = "wb");
    db.path <- nmdb;
  }
  db.map <-  readRDS(db.path);
  hit.inx <- match(entrez.vec, db.map[, "gene_id"]);
  entrezs <- db.map[hit.inx, "accession"];
  entrezs <- as.vector(entrezs)
  return(entrezs);
}

.readDataTable <- function(fileName){
  dat <- try(data.table::fread(fileName, header=TRUE, check.names=FALSE, blank.lines.skip=TRUE, data.table=FALSE));
  if(class(dat) == "try-error" || any(dim(dat) == 0)){
    print("Using slower file reader ...");
    formatStr <- substr(fileName, nchar(fileName)-2, nchar(fileName))
    if(formatStr == "txt"){
      dat <- try(read.table(fileName, header=TRUE, comment.char = "", check.names=F, as.is=T));
    }else{ # note, read.csv is more than read.table with sep=","
      dat <- try(read.csv(fileName, header=TRUE, comment.char = "", check.names=F, as.is=T));
    }
  }
  return(dat);
}

LoadKEGGLibOther<-function(type){
  if(type == "integ"){
    kegg.path <- paste(lib.path, data.org, "/integ.rds", sep="");
  }else if(type == "ginteg"){
    kegg.path <- paste(lib.path,"microbiome", "/integ.rds", sep="");
  }else if(type == "gkeggm"){
    kegg.path <- paste(lib.path,"microbiome",  "/keggm.rds", sep="");
  }else{
    kegg.path <- paste(lib.path, data.org,"/keggm.rds", sep="");
  }
  if(!.on.public.web){
    nmdb <- basename(kegg.path);
    download.file(kegg.path, destfile = nmdb, method="libcurl", mode = "wb");
    kegg.path <- nmdb;
  }
  kegg.anot <- readRDS(kegg.path)
  current.setlink <- " "
  current.mset <- kegg.anot;
  current.setlink <<- current.setlink;
  current.setids <<- names(kegg.anot);
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}

#' Perform metabolite or gene/metabolite enrichment analysis
#'
#' @param dataSetObj Input the name of the created dataSetObj (see Init.Data)
#' @param file.nm File name of result table to be exported in csv format, do not include file extension
#' @param fun.type Enrichment type: keggm (kegg metabolite), integ (joint gene/metabolite)
#' @param ids String of ids to be tested separated by "; "
#'
#' @export
#'
PerformMetEnrichment <- function(dataSetObj=NA, file.nm, fun.type, ids){
  dataSet <- .get.nSet(dataSetObj);
  if(ids=="Not_applicable"){
    ora.vec <- keggp.allfeatures;
    names(ora.vec) <- ora.vec;
  } else {
    ora.vec <- unlist(strsplit(ids, "; "));
    if(length(which(grepl("HMDB", ora.vec)))/length(ora.vec) > 0.6){
      ora.vec <- doHMDB2KEGGMapping(ora.vec);
    }
    if(length(which(grepl("PBCM", ora.vec)))/length(ora.vec) > 0.6){
      ora.vec <- doPubchem2KEGGMapping(ora.vec);
    }
    names(ora.vec) <- ora.vec;
    ora.vecu <<- ora.vec;
  }

  if(fun.type == "integ"){
    res1 = PerformEnrichAnalysisKegg(dataSet, file.nm, "kegg", ora.vec)
    res2 = PerformEnrichAnalysisKegg(dataSet, file.nm, "keggm", ora.vec)
    res3 = PerformEnrichAnalysisKegg(dataSet, file.nm, "integ", ora.vec)
  }else{
    res1 = PerformEnrichAnalysisKegg(dataSet, file.nm, fun.type, ora.vec)
    res.mat <- as.data.frame(res1)
  }
  if(fun.type == "integ" || fun.type == "ginteg"){
    inx = which(rownames(res1) %in% rownames(res2))
    subres1 = as.data.frame(res1[inx,])
    inx = which(rownames(res2) %in% rownames(res1))
    subres2 = as.data.frame(res2[inx,])
    inx = which(rownames(res3) %in% rownames(subres2))
    subres3 = as.data.frame(res3[inx,])

    ord = order(rownames(subres1));
    subres1 = subres1[ord,]
    ord = order(rownames(subres2));
    subres2 = subres2[ord,]
    ord = order(rownames(subres3));
    subres3 = subres3[ord,]
    integ=data.frame(hitsG = subres1$Hits,hitsM = subres2$Hits,hitsTotal = subres3$Hits, P.ValueG=subres1$P.Value, P.ValueM=subres2$P.Value, P.ValueMerge=subres3$P.Value, P.ValueJoint=subres2$P.Value)

    rownames(integ) = rownames(subres1)
    for(i in 1:nrow(integ)){
      if(integ$P.ValueG[i] != 1 && integ$P.ValueM[i] != 1){
        #integ$P.ValueJoint[i] = metap::sumlog(c(integ$P.ValueG[i], integ$P.ValueM[i]))$p
        integ$P.ValueJoint[i] = metap::sumz(p=c(integ$P.ValueG[i], integ$P.ValueM[i]), weight=c(stouffer_gene_percent, stouffer_compound_percent))$p
      }else{
        integ$P.ValueJoint[i]=1
      }
    }

    inxM = integ[,"hitsG"] == 0
    inxG = integ[,"hitsM"] == 0
    inxJ = integ[,"hitsM"] != 0 & integ[,"hitsM"] != 0
    integP = integ[,"P.ValueJoint"]
    integP[inxM] = integ[,"P.ValueM"][inxM]
    integP[inxG] = integ[,"P.ValueG"][inxG]
    integ$integP = integP
    res.mat <- integ
    hits.query<<-hits.query
    SaveIntegEnr(file.nm,res.mat);
  }else{
    SaveSingleOmicsEnr(file.nm,res.mat);
    if(fun.type == "kegg"){
      kg.query <<-hits.query
    }else if(fun.type == "keggm"){
      km.query <<-hits.query
    }
  }

  return(.set.nSet(dataSet));
}

PerformEnrichAnalysisKegg <- function(dataSetObj=NA, file.nm, fun.type, ora.vec){
  dataSet <- .get.nSet(dataSetObj);
  # prepare lib
  LoadLib(fun.type);

  # prepare query
  ora.nms <- names(ora.vec);
  lens <- lapply(current.geneset,
                 function(x) {
                   length(unique(x))
                 }
  );
  keepInx <- lens > 2

  current.geneset = current.geneset[keepInx]
  current.geneset = current.geneset[which(!names(current.geneset) %in% dataSet$toremove)]

  # prepare for the result table
  set.size<-length(current.geneset);
  res.mat<-matrix(0, nrow=set.size, ncol=5);
  rownames(res.mat)<-names(current.geneset);
  colnames(res.mat)<-c("Total", "Expected", "Hits", "P.Value", "FDR");

  # need to cut to the universe covered by the pathways, not all genes

  hits.inx <- ora.vec %in% current.universe;

  #calculate weight for stouffer
  if(fun.type == "kegg"){
    stouffer_gene_percent <<- length(hits.inx)/length(current.universe)
  }else if(fun.type == "keggm"){
    stouffer_compound_percent <<- length(hits.inx)/length(current.universe)
  }

  ora.vec <- ora.vec[hits.inx];
  ora.nms <- ora.nms[hits.inx];

  q.size<-length(ora.vec);

  # get the matched query for each pathway
  hits.query <- lapply(current.geneset,
                       function(x) {
                         ora.nms[ora.vec%in%unlist(x)];
                       }
  );
  hits.query <- lapply(hits.query, function(x){unique(x)})

  names(hits.query) <- names(current.geneset);
  hit.num<-unlist(lapply(hits.query, function(x){length(x)}), use.names=FALSE);

  # total unique gene number
  uniq.count <- length(current.universe);

  # unique gene count in each pathway
  set.size <- unlist(lapply(current.geneset, length));

  res.mat[,1]<-set.size;
  res.mat[,2]<-q.size*(set.size/uniq.count);
  res.mat[,3]<-hit.num;

  # use lower.tail = F for P(X>x)
  raw.pvals <- phyper(hit.num-1, set.size, uniq.count-set.size, q.size, lower.tail=F);
  res.mat[,4]<- raw.pvals;
  res.mat[,5] <- p.adjust(raw.pvals, "fdr");
  all.res.mat <<- res.mat
  hits.query <<- hits.query
  if(fun.type == "integ"){
    integ.query<- list()
    integ.query<- hits.query
    integ.query<<- integ.query
  }
  # now, clean up result, synchronize with hit.query
  return(all.res.mat);
}

SaveSingleOmicsEnr <- function(file.nm,res.mat){
  inx = res.mat[,3]>0
  res.mat <- res.mat[inx,];
  hits.query <- hits.query[inx];
  if(nrow(res.mat)> 1){
    # order by p value
    ord.inx<-order(res.mat[,4]);
    res.mat <- signif(res.mat[ord.inx,],3);
    hits.query <- hits.query[ord.inx];

    imp.inx <- res.mat[,4] < 0.05
    if(sum(imp.inx) < 10){ # too little left, give the top ones
      topn <- ifelse(nrow(res.mat) > 10, 10, nrow(res.mat));
      res.mat <- res.mat[1:topn,];
      hits.query <- hits.query[1:topn];
    }else{
      res.mat <- res.mat[imp.inx,];
      hits.query <- hits.query[imp.inx];
      if(sum(imp.inx) > 120){
        # now, clean up result, synchronize with hit.query
        res.mat <- res.mat[1:120,];
        hits.query <- hits.query[1:120];
      }
    }
  }

  #get gene symbols
  resTable <- data.frame(Pathway=rownames(res.mat), res.mat);
  current.msg <<- "Functional enrichment analysis was completed";

  # write json
  require("RJSONIO");
  fun.anot <- hits.query;
  fun.padj <- resTable[,6]; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj)};
  fun.pval <- resTable[,5]; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval)};
  hit.num <- resTable[,4]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num)};
  fun.ids <- as.vector(current.setids[names(fun.anot)]);
  if(length(fun.ids) ==1) {fun.ids <- matrix(fun.ids)};

  json.res <- list(
    fun.link = current.setlink[1],
    fun.anot = fun.anot,
    fun.ids = fun.ids,
    fun.pval = fun.pval,
    fun.padj = fun.padj,
    hit.num = hit.num
  );
  json.mat <- toJSON(json.res);
  json.nm <- paste(file.nm, ".json", sep="");

  sink(json.nm)
  cat(json.mat);
  sink();
  hitss <- lapply(hits.query, function(x){paste(x, collapse=" ")})
  hitss <- unlist(hitss)
  resTable$Features <- hitss
  csv.nm <- paste(file.nm, ".csv", sep="");
  fast.write.csv(resTable, file=csv.nm, row.names=F);
}

SaveIntegEnr <- function(file.nm,res.mat){
  inx = res.mat$hitsTotal > 0
  res.mat <- res.mat[inx,];
  if(nrow(res.mat)> 1){
    # order by p value
    ord.inx<-order(res.mat[,"integP"]);
    res.mat <- signif(res.mat[ord.inx,],3);
    hits.query <- hits.query[ord.inx];

    imp.inx <- res.mat[,"integP"] > 0
    if(sum(imp.inx) < 10){ # too little left, give the top ones
      topn <- ifelse(nrow(res.mat) > 10, 10, nrow(res.mat));
      res.mat <- res.mat[1:topn,];
      hits.query <- hits.query[1:topn];
    }else{
      res.mat <- res.mat[imp.inx,];
      hits.query <- hits.query[imp.inx];
      if(sum(imp.inx) > 120){
        # now, clean up result, synchronize with hit.query
        res.mat <- res.mat[1:120,];
        hits.query <- hits.query[1:120];
      }
    }
  }


  #get gene symbols
  resTable <- data.frame(Pathway=rownames(res.mat), res.mat);
  current.msg <<- "Functional enrichment analysis was completed";

  # write json
  require("RJSONIO");
  integ.query <- integ.query[resTable$Pathway]
  fun.anot <- integ.query;
  fun.nms <- resTable[,"Pathway"];  if(length(fun.nms) ==1) { fun.nms <- matrix(fun.nmsl) };
  integ.pval <- resTable[,"integP"];  if(length(integ.pval) ==1) { integ.pval <- matrix(integ.pval) };
  fun.pval1 <- resTable[,"P.ValueMerge"]; if(length(fun.pval1) ==1) { fun.pval1 <- matrix(fun.pval1) };
  fun.pval2 <- resTable[,"P.ValueJoint"]; if(length(fun.pval2) ==1) { fun.pval2 <- matrix(fun.pval2) };
  hit.num1 <- resTable[,"hitsM"]; if(length(hit.num1) ==1) { hit.num1 <- matrix(hit.num1) };
  hit.num2 <- resTable[,"hitsG"]; if(length(hit.num2) ==1) { hit.num2 <- matrix(hit.num2) };
  fun.ids <- as.vector(current.setids[names(fun.anot)]);
  if(length(fun.ids) == 1) { fun.ids <- matrix(fun.ids) };

  json.res <- list(
    fun.nms = fun.nms,
    fun.link = current.setlink[1],
    fun.anot = fun.anot,
    fun.ids = fun.ids,
    fun.pval1 = fun.pval1,
    fun.pval2 = fun.pval2,
    integ.pval = integ.pval,
    hit.num1 = hit.num1,
    hit.num2 = hit.num2
  );
  json.mat <- toJSON(json.res);
  json.nm <- paste(file.nm, ".json", sep="");
  integ.query <<-integ.query
  sink(json.nm)
  cat(json.mat);
  sink();

  hitss = lapply(integ.query, function(x){paste(x, collapse=" ")})
  hitss = unlist(hitss)
  resTable$Features = hitss
  # write csv
  csv.nm <- paste(file.nm, ".csv", sep="");
  fast.write.csv(resTable, file=csv.nm, row.names=F);
}


.pathwayMetGenePair <- function(met.vec, gene.vec){

  if(data.org == "microbiome"){
    table.nm <- "microbiome_gem"
  }else{
    table.nm <- paste(data.org, "kegg", sep="_");
  }
  resMet <- QueryMetSQLiteNet(table.nm, rownames(met.vec), "direct");
  resGene <- QueryMetSQLiteNet(table.nm, rownames(gene.vec), "inverse");
  res <- rbind(resMet, resGene)

  geneDf <- data.frame(entrez=rownames(gene.vec), exprG=unname(gene.vec))
  metDf <- data.frame(kegg=rownames(met.vec), exprM=unname(met.vec))
  res <- merge(res, geneDf, by="entrez", all=T)
  res <- merge(res, metDf, by="kegg", all=T)
  res <- unique(res)
  na.inx <- is.na(res)
  res[na.inx] <- 0
  res <- res[(res$symbol!=0),]
  res <- res[(res$met!=0),]
  res <- res[order(-abs(res[,"exprG"]), -abs(res[,"exprM"])),]
  inx <- res$exprG !=0 | res$exprM !=0
  topres <- res[inx,]
  res <- res[-inx,]
  res <- rbind(topres,res)
  restoremove <- c("H2O")
  res <- res[which(res$met != restoremove),]
  res <- res[,c("met","kegg","exprM","symbol","entrez","exprG")]
  colnames(res) <- colnames(res) <- c("Cmpd_name", "ID", "expression", "Gene_name", "ID", "expression")
  return(res)
}


# note, last two par only for STRING database
QueryPpiSQLite <- function(table.nm, q.vec, requireExp, min.score){
  require('RSQLite')
  path <- paste0(sqlite.path, "ppi.sqlite");
  if(!PrepareSqliteDB(path, .on.public.web)){
    stop("Sqlite database is missing, please check your internet connection!");
  }
  ppi.db <- dbConnect(SQLite(), path);
  query <- paste(shQuote(q.vec),collapse=",");
  if(grepl("string$", table.nm)){
    if(requireExp){
      statement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))  AND combined_score >=", min.score, " AND experimental > 0", sep="");
    }else{
      statement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))  AND combined_score >=", min.score, sep="");
    }
  }else{
    statement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))", sep="");
  }

  ppi.res <- .query.sqlite(ppi.db, statement);

  # remove dupliated edges
  # ppi.res <- ppi.res[!duplicated(ppi.res$row_id),]
  return(ppi.res);
}

# table name is org code, id.type is column name
QueryMirSQLite <- function(table.nm, id.type, q.vec, inv, db.nm){
  require('RSQLite');
  path <- paste0(sqlite.path, "mir2gene.sqlite")
  if(!PrepareSqliteDB(path, .on.public.web)){
    stop("Sqlite database is missing, please check your internet connection!");
  }
  mir.db <- dbConnect(SQLite(), path);
  query <- paste (shQuote(q.vec),collapse=",");

  if(inv == "inverse"){
    statement <- paste("SELECT * FROM ", table.nm, " WHERE entrez IN (",query,")", " AND ", db.nm," == 1 ", sep="");
  }else{
    statement <- paste("SELECT * FROM ", table.nm, " WHERE ", id.type," IN (",query,")", " AND ", db.nm," == 1 ", sep="");
  }

  mir.dic <- .query.sqlite(mir.db, statement);
  mir.dic <- mir.dic[complete.cases(mir.dic),]; # get all records
  return(mir.dic);
}

QueryDrugSQLite <- function(q.vec, regsearch){
  require('RSQLite');
  path <- paste0(sqlite.path, "drug.sqlite")
  if(!PrepareSqliteDB(path, .on.public.web)){
    stop("Sqlite database is missing, please check your internet connection!");
  }
  mir.db <- dbConnect(SQLite(), path);
  table.nm = "human"
  query <- paste (shQuote(q.vec),collapse=",");

  statement <- paste("SELECT * FROM ", table.nm, " WHERE upid IN (",query,")", sep="");

  mir.dic <- .query.sqlite(mir.db, statement);
  mir.dic <- mir.dic[complete.cases(mir.dic),]; # get all records
  return(mir.dic);
}

QueryMicSQLite <- function(q.vec, table.nm, sql.nm, min.score, currExclude=T, uniExclude=T, orphExclude = T){

  require('RSQLite');
  path <- paste0(sqlite.path, sql.nm);
  if(!PrepareSqliteDB(path, .on.public.web)){
    stop("Sqlite database is missing, please check your internet connection!");
  }
  mir.db <- dbConnect(SQLite(), path);

  query <- paste (shQuote(q.vec),collapse=",");

  if(uniExclude == "TRUE"){
    statement <- paste("SELECT * FROM ", table.nm, " WHERE ", table.nm, " IN (",query,")"," AND potential >=", min.score, " AND ExcluUni == 0", sep="");
  }else{
    statement <- paste("SELECT * FROM ", table.nm, " WHERE ", table.nm, " IN (",query,")"," AND potential >=", min.score, sep="");
  }

  mir.dic <- .query.sqlite(mir.db, statement);
  mir.dic <- mir.dic[complete.cases(mir.dic),]; # get all records

  met.ids <- mir.dic$metID

  if(orphExclude ==  "TRUE"){
    met.path <- paste(lib.path, "microbiome", "/met4path.rds", sep="");
  }else{
    met.path <- paste(lib.path, "microbiome", "/metInfo.rds", sep="");
  }

  if(!.on.public.web){
    nmdb <- basename(met.path);
    download.file(met.path, destfile = nmdb, method="libcurl", mode = "wb");
    met.path <- nmdb;
  }

  met.info <- readRDS(met.path);
  conv <- data.frame(metID=met.info$metID, KEGG=met.info$KEGG)
  mir.dic <- merge(mir.dic, conv, by="metID");
  mir.dic$KEGG[is.na(mir.dic$KEGG)] <- mir.dic$metID[is.na(mir.dic$KEGG)]

  return(mir.dic);
}

QueryMetSQLiteNet <- function(table.nm, q.vec, inv){
  require('RSQLite');
  path <- paste0(sqlite.path, "omicsnet_met.sqlite")
  if(!PrepareSqliteDB(path, .on.public.web)){
    stop("Sqlite database is missing, please check your internet connection!");
  }
  lnc.db <- dbConnect(SQLite(), path);
  query <- paste (shQuote(q.vec),collapse=",");

  if(inv == "inverse"){
    statement <- paste("SELECT * FROM ", table.nm, " WHERE entrez IN (",query,")", sep="");
  }else{
    statement <- paste("SELECT * FROM ", table.nm, " WHERE kegg IN (",query,")", sep="");
  }

  my.dic <- .query.sqlite(lnc.db, statement);
  return(my.dic);
}

QueryM2mSQLiteNet <- function(table.nm, q.vec, inv){

  require('RSQLite');
  path <- paste0(sqlite.path, "omicsnet_met.sqlite")
  if(!PrepareSqliteDB(path, .on.public.web)){
    stop("Sqlite database is missing, please check your internet connection!");
  }
  lnc.db <- dbConnect(SQLite(), path);
  query <- paste (shQuote(q.vec),collapse=",");

  if(inv == "inverse"){
    statement <- paste("SELECT * FROM ", table.nm, " WHERE productID IN (",query,")", sep="");
  }else{
    statement <- paste("SELECT * FROM ", table.nm, " WHERE sourceID IN (",query,")", sep="");
  }
  my.dic <- .query.sqlite(lnc.db, statement);
  return(my.dic);
}

QueryTFSQLite <- function(table.nm, q.vec, inv){
  require('RSQLite');
  path <- paste0(sqlite.path, "tf2gene.sqlite")
  if(!PrepareSqliteDB(path, .on.public.web)){
    stop("Sqlite database is missing, please check your internet connection!");
  }
  chem.db <- dbConnect(SQLite(), path);

  query <- paste (shQuote(q.vec),collapse=",");

  if(inv == "inverse"){
    statement <- paste("SELECT * FROM ", table.nm, " WHERE entrez IN (",query,")", sep="");
  }else{
    statement <- paste("SELECT * FROM ", table.nm, " WHERE tfid IN (",query,")", sep="");
  }
  chem.dic <- .query.sqlite(chem.db, statement);
  return(chem.dic);
}


Query.snpDB <- function(db.path, q.vec, table.nm, col.nm){

  require('RSQLite');
  db.path <<- db.path;

  db.path <- paste0(db.path, ".sqlite");
  if(!PrepareSqliteDB(db.path, .on.public.web)){
    stop("Sqlite database is missing, please check your internet connection!");
  }
  mir.db <- dbConnect(SQLite(), db.path);

  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", table.nm, " WHERE ", col.nm," IN (", query, ")", sep="");
  mir.dic <- .query.sqlite(mir.db, statement);

  # remove duplicates
  if("mgwas" %in% colnames(mir.dic)){
    dup.inx <- duplicated(mir.dic$mgwas);
    mir.dic <- mir.dic[!dup.inx, ];
  }
  dataSet$snpTable <- mir.dic;

  return(mir.dic);
}


Query.PhenoScanner <- function(snpquery=NULL, genequery=NULL, regionquery=NULL, catalogue="GWAS", pvalue=1E-5, proxies="None", r2=0.8, build=37){
  if(!exists("my.query.phenoscanner")){ # public web on same user dir
    compiler::loadcmp("../../rscripts/OmicsNetR/R/utils_phenoscanner.Rc");
  }
  return(my.query.phenoscanner(snpquery, genequery, regionquery, catalogue, pvalue, proxies, r2, build));
}


QueryVEP <- function(q.vec,vepDis,queryType,snpRegion,content_type="application/json" ){

  require("httr")
  #library(jsonlite)
  #library(xml2)
  server <- "http://rest.ensembl.org"
  if(snpRegion==T){
    ext <- "/vep/human/region/"
  }else{
    ext <- "/vep/human/id/"
  }
  r=list()
  resvep = list()
  vepDis = as.numeric(vepDis)*1000

  for(i in 1:length(q.vec)){
    qr <- gsub("^chr","",q.vec[i])
    r[[i]] <- GET(paste(server, ext, qr,"?distance=",vepDis,sep = ""),  accept(content_type))
    stop_for_status(r[[i]])
    resvep[[i]] = content(r[[i]])[[1]][["transcript_consequences"]]
  }
  names(resvep)=q.vec

  resvep2 = lapply(resvep,function(s) {lapply(s, function(x) {
    list(
      gene_symbol=ifelse(length(x[["gene_symbol"]])!=0,x[["gene_symbol"]],"NA"),
      gene_id=ifelse(length(x[["gene_id"]]) !=0, x[["gene_id"]],'NA'),
      hgnc_id=ifelse(length(x[["hgnc_id"]]) !=0, x[["hgnc_id"]],'NA'),
      transcript_id=ifelse(length(x[["transcript_id"]]) !=0, x[["transcript_id"]],'NA'),
      consequence_terms=ifelse(length(x[["consequence_terms"]]) !=0, paste(unlist(x[["consequence_terms"]]),collapse = ";"),'NA'),
      distance=ifelse(length(x[["distance"]])!=0, x[["distance"]],'NA'))
  })
  }
  )
  resvep3 = do.call(rbind,lapply(resvep2,function(s) {
    do.call(rbind.data.frame,s)
  } )
  )
  resvep3$rsid = gsub("\\.[0-9]*","",rownames(resvep3))

  row.names(resvep3)=NULL
  return(resvep3)
}

QueryMicM2mSQLiteNet <- function(table.nm, q.vec){

  require('RSQLite');
  path <- paste0(sqlite.path, "omicsnet_met.sqlite")
  if(!PrepareSqliteDB(path, .on.public.web)){
    stop("Sqlite database is missing, please check your internet connection!");
  }

  lnc.db <- dbConnect(SQLite(), path);
  query <- paste (shQuote(q.vec),collapse=",");

  statement1 <- paste("SELECT * FROM ", table.nm, " WHERE productID IN (",query,")", sep="");

  statement2 <- paste("SELECT * FROM ", table.nm, " WHERE sourceID IN (",query,")", sep="");


  my.dic1 <- .query.sqlite(lnc.db, statement1);
  lnc.db <- dbConnect(SQLite(), path);

  my.dic2 <- .query.sqlite(lnc.db, statement2);

  my.dic <- unique(rbind(my.dic1,my.dic2))
  return(my.dic);
}

getApiResult <- function(url="NA", init=TRUE){
  json_data = tryCatch({
    rjson::fromJSON(file=url);
  }, warning = function(w) {
    print(w);
    if(init){
      print("API call failed, calling again!");
      Sys.sleep(5);
      getApiResult(url, F);
    }else{
      return(0);
    }
  }, error = function(e) {
    if(init){
      print(e);
      print("API call failed, calling again2!");
      Sys.sleep(5);
      getApiResult(url, F);
    }else{
      return(0);
    }
  }, finally = {

  })
  return(json_data)
}
