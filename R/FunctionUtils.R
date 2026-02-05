##################################################
## R script for OmicsNet
## Description: ID mapping, GO/Pathway enrichment analysis
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

##################################################
## Compound Name Matching Utilities
## Fuzzy matching with synonym database support
##################################################

#' Perform approximate compound name matching using synonym database
PerformCompoundApproxMatch <- function(q, cmpd.db, syn.db) {
  # Check if databases match in size
  if(length(syn.db$syns.vec) != nrow(cmpd.db)) {
    min.size <- min(length(syn.db$syns.vec), nrow(cmpd.db));
    com.nms <- cmpd.db$name[1:min.size];
    syns.vec <- syn.db$syns.vec[1:min.size];
    syns.list <- syn.db$syns.list[1:min.size];
  } else {
    if("lipid" %in% colnames(cmpd.db)) {
      nonLipidInx <- cmpd.db$lipid == 0;
      com.nms <- cmpd.db$name[nonLipidInx];
      syns.vec <- syn.db$syns.vec[nonLipidInx];
      syns.list <- syn.db$syns.list[nonLipidInx];
    } else {
      com.nms <- cmpd.db$name;
      syns.vec <- syn.db$syns.vec;
      syns.list <- syn.db$syns.list;
    }
  }

  q.length <- nchar(q);
  s <- c(0, 0.1, 0.2);
  candidates <- NULL;

  for (j in s) {
    new.q <- q;
    if(q.length > 32){
      new.q <- substr(q, 1, 32);
    }

    matched.inx <- agrep(new.q, syns.vec, ignore.case=TRUE, max.distance=j, useBytes=TRUE);

    if(length(matched.inx) > 0) {
      candidates <- data.frame(
        index = vector(mode = "numeric", length=length(matched.inx)),
        value = vector(mode = "character", length=length(matched.inx)),
        score = vector(mode = "numeric", length=length(matched.inx)),
        stringsAsFactors = FALSE
      );

      for(n in 1:length(matched.inx)){
        nm.vec <- syns.list[[matched.inx[n]]];
        hit3.inx <- agrep(q, nm.vec, ignore.case=TRUE, max.distance=j, useBytes=TRUE);

        if(length(hit3.inx) > 0){
          hit3.nm <- vector(mode = "character", length=length(hit3.inx));
          hit3.score <- vector(mode = "numeric", length=length(hit3.inx));

          for(k in 1:length(hit3.inx)){
            idx <- hit3.inx[k];
            hit3.nm[k] <- nm.vec[idx];
            hit3.score[k] <- j + abs(nchar(nm.vec[idx]) - nchar(q)) / (10 * nchar(q));
          }

          matches2 <- c();
          if(length(grep("^[1-9a-z]{2}", q, ignore.case=TRUE)) > 0){
            matches2 <- grep(paste("^", substr(q, 1, 2), sep=""), hit3.nm, ignore.case=TRUE);
          } else if (length(grep("^[1-9a-z]", q, ignore.case=TRUE)) > 0){
            matches2 <- grep(paste("^", substr(q, 1, 1), sep=""), hit3.nm, ignore.case=TRUE);
          }

          if(length(matches2) > 0){
            hit3.score[matches2] <- hit3.score[matches2] - 0.05;
          }

          best.inx <- which(hit3.score == min(hit3.score))[1];
          candidates[n, 1] <- matched.inx[n];
          candidates[n, 2] <- com.nms[matched.inx[n]];
          candidates[n, 3] <- hit3.score[best.inx];
        }
      }

      rm.inx <- is.na(candidates[,2]) | candidates[,2] == "NA" | candidates[,2] == "";
      candidates <- candidates[!rm.inx, , drop=FALSE];
      candidates <- candidates[order(candidates[,3], decreasing=FALSE), , drop=FALSE];

      if(nrow(candidates) > 10){
        candidates <- candidates[1:10, , drop=FALSE];
      }

      if(nrow(candidates) > 0) {
        return(candidates);
      }
    }
  }

  return(NULL);
}

#' Match compound names with fuzzy matching
MatchCompoundNames <- function(query.vec, cmpd.db, syn.db) {
  n <- length(query.vec);
  hit.inx <- rep(0, n);
  hit.values <- rep("", n);
  match.state <- rep(0, n);

  for(i in 1:n) {
    q <- query.vec[i];
    candidates <- PerformCompoundApproxMatch(q, cmpd.db, syn.db);

    if(!is.null(candidates) && nrow(candidates) > 0) {
      best.idx <- candidates[1, 1];

      if(best.idx > 0 && best.idx <= nrow(cmpd.db)) {
        hit.inx[i] <- best.idx;
        hit.values[i] <- cmpd.db$name[best.idx];
        match.state[i] <- 1;
      }
    }
  }

  return(list(
    hit.inx = hit.inx,
    hit.values = hit.values,
    match.state = match.state
  ));
}

LoadLib <- function(fun.type){
  fun.type <- tolower(fun.type);
  if(fun.type %in% c("kegg", "reactome", "motif", "bp", "cc", "mf","panth_mf", "go_bp", "go_cc","go_mf",  "go_panthbp", "go_panthcc", "go_panthmf","mirfamily", "mircluster", "mirdisease","mirtissue", "mirfunction")){
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
  names(set.ids) <-  set$term;

  names(current.geneset) <- set[[termInx]];
  current.geneset <<- current.geneset

  current.setlink <<- set[[linkInx]];
  current.setids <<- set.ids;
  current.universe <<- unique(unlist(current.geneset));
}


SearchReg <- function(file.nm, fun.type, IDs, count, sourceView="2d"){
  regCount <<- count;
  regBool <<- "true";
  res <- PerformNetEnrichment(file.nm, fun.type, IDs, sourceView);
}

regEnrichment <- function(file.nm, fun.type, IDs, netInv, sourceView="2d"){
  regBool <<- "false";
  res <- PerformNetEnrichment(file.nm, fun.type, IDs, sourceView);
}

#' Perform gene enrichment analysis or identify gene regulatory targets
#'
#' @param file.nm File name of result table to be exported in csv format, do not include file extension
#' @param fun.type Enrichment database type
#' @param IDs String of ids to be tested separated by "; "
#'
#' @export
#'
PerformNetEnrichment <- function(file.nm, fun.type, IDs, sourceView="2d"){
  # note: hit.query, resTable must synchronize
  # prepare query
  ora.vec <- NULL;
  idtype <- "entrez"
  save.type <- ifelse(is.null(sourceView) || sourceView == "", "network", sourceView)
  ora.vec <- unlist(strsplit(IDs, "; "));
  names(ora.vec) <- ora.vec;
  if(fun.type %in% c("trrust", "encode", "jaspar", "mirnet", "met", "drugbank")){
    netInv <- "inverse";
    res <- PerformRegEnrichAnalysis(file.nm, fun.type, ora.vec, netInv, idtype, sourceView);
  } else{
    res <- PerformEnrichAnalysis(file.nm, fun.type, ora.vec, sourceView=sourceView);
  }
  return(res);
}

PerformRegEnrichAnalysis <- function(file.nm, fun.type, ora.vec, netInv, idType, sourceView="2d"){
    if(!exists("my.reg.enrich")){ # public web on same user dir
        compiler::loadcmp("../../rscripts/OmicsNetR/R/utils_reg_enrich.Rc");
  }
    res <- my.reg.enrich(file.nm, fun.type, ora.vec, netInv, idType, sourceView);
    return(res);
    }

# note: hit.query, resTable must synchronize
# ora.vec should contains entrez ids, named by entrez ids ASWELL
PerformEnrichAnalysis <- function(file.nm, fun.type, ora.vec, save.type="network", sourceView="2d"){

  # prepare lib
  LoadLib(fun.type);

  # prepare query
  ora.nms <- names(ora.vec);

  # prepare for the result table
  set.size<-length(current.geneset);
  res.mat<-matrix(0, nrow=set.size, ncol=5);
  rownames(res.mat)<-names(current.geneset);
  colnames(res.mat)<-c("Total", "Expected", "Hits", "Pval", "FDR");

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

  # add a new column containing pathway names
  ###########
  #### note the column position shifts by one as new columns introduced
  ###########
  resTable <- data.frame(Pathway=rownames(res.mat),IDs=as.vector(current.setids[names(hits.query)]), res.mat);
  current.msg <<- "Functional enrichment analysis was completed";

  # write json
  require("RJSONIO");
  fun.anot <- hits.query;
  fun.padj <- resTable$FDR; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
  fun.pval <- resTable$Pval; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
  hit.num <- paste0(resTable$Hits,"/",resTable$Total); if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
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

  gene.vec <- current.universe;
  sym.vec <- doEntrez2SymbolMapping(gene.vec);
  gene.nms <- sym.vec;

  current.geneset.symb <- lapply(current.geneset, 
                       function(x) {
                         gene.nms[gene.vec%in%unlist(x)];
                       }
  );


  #record table for report
  type <- save.type;
  dataSet$imgSet$enrTables[[type]] <- list()
  dataSet$imgSet$enrTables[[type]]$table <- resTable;
  dataSet$imgSet$enrTables[[type]]$library <- fun.type
  dataSet$imgSet$enrTables[[type]]$algo <- "Overrepresentation Analysis"
  dataSet$imgSet$enrTables[[type]]$sourceView <- sourceView


  dataSet$imgSet$enrTables[[type]]$current.geneset <- current.geneset;
  dataSet$imgSet$enrTables[[type]]$hits.query <- hits.query;
  dataSet$imgSet$enrTables[[type]]$current.setids <- current.setids;
  dataSet$imgSet$enrTables[[type]]$res.mat<- res.mat;
  dataSet$imgSet$enrTables[[type]]$current.geneset.symb <- current.geneset.symb;
  #print(paste("saveType==" ,save.type))
  dataSet <<- dataSet

  # write csv
  # print(file.nm);
  csv.nm <- paste(file.nm, ".csv", sep="");
  a <- lapply(fun.anot,function(x){paste(x,collapse = "; ")})
  resTable$ids <- unlist(a);
  fast.write.csv(resTable, file=csv.nm, row.names=F);
  fast.write.csv(resTable, file=paste0(type, "_enr_table.csv"), row.names=F);
  Sys.sleep(0.15);  # CRITICAL: Prevent race condition - allow file system to sync before Java reads

  return(1);
}

# Helper function to query gene database from SQLite
queryGeneDB <- function(table.nm, org){
  if(org == "custom"){
    db.map <- qs::qread("anot_table.qs");
    return(db.map);
  }

  require('RSQLite');
  db.path <- paste0(sqlite.path, org, "_genes.sqlite");

  if(!PrepareSqliteDB(db.path, .on.public.web)){
    stop("Sqlite database is missing, please check your internet connection!");
  }

  conv.db <- dbConnect(SQLite(), db.path);
  tbls <- dbListTables(conv.db);

  if(!table.nm %in% tbls){
    dbDisconnect(conv.db);
    return(NULL);
  }

  db.map <- dbReadTable(conv.db, table.nm);
  dbDisconnect(conv.db);

  return(db.map);
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
      db.map <- queryGeneDB("entrez", data.org);
      if(is.null(db.map)){
        return(data.frame(gene_id=character(0), accession=character(0)));
      }
      hit.inx <- match(q.vec, db.map[, "gene_id"]);
      entrezs <- db.map[hit.inx, ]
      entrezs <- entrezs[c(1,2)]
      colnames(entrezs) = c("gene_id", "accession");
    }else if(type == "ctd"){
      db.path <- paste(lib.path, "chem", "/ctd.rds", sep="")
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
    }else if(type == "drugbank"){
      db.path <- paste(lib.path, "drug", "/drugbank.rds", sep="")
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
    }else if(type == "ko"){
      db.path <- paste(lib.path, "microbiome", "/ko.rds", sep="")
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
    }
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
    gene.map <- queryGeneDB("entrez", data.org);
    if(is.null(gene.map)){
      return(data.frame(gene_id=character(0), accession=character(0)));
    }
    if(data.org == "hsa"){
      q.vec = toupper(q.vec);
    }
    hit.inx <- match(q.vec, gene.map[, "symbol"]);
    entrezs <- gene.map[hit.inx, ];
    entrezs = entrezs[c(1,2)];
    colnames(entrezs) <- c("gene_id", "accession");
  }else if(type %in% c("meta", "kegg", "chebi", "name", "bigg", "pubchem", "hmdb")){

    # For name type with synonym matching, use the .qs database
    if(type == "name"){
      db.path.qs <- paste(lib.path, "lib/compound_db.qs", sep="");

      if(!.on.public.web){
        nm <- basename(db.path.qs);
        download.file(db.path.qs, destfile = nm, method="libcurl", mode = "wb");
        db.path.qs <- nm;
      }

      if(file.exists(db.path.qs)){
        cmpd.map <- qs::qread(db.path.qs);
      }else{
        # Fallback to .rds if .qs not found
        db.path <- paste(lib.path, "lib/compound_db.rds", sep="");
        if(!.on.public.web){
          nm <- basename(db.path);
          download.file(db.path, destfile = nm, method="libcurl", mode = "wb");
          db.path <- nm;
        }
        cmpd.map <- readRDS(db.path);
      }
    }else{
      # For other types, use the .rds database
      db.path <- paste(lib.path, "lib/compound_db.rds", sep="");

      if(!.on.public.web){
        nm <- basename(db.path);
        download.file(db.path, destfile = nm, method="libcurl", mode = "wb");
        db.path <- nm;
      }

      cmpd.map <- readRDS(db.path);
    }

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
    }else if(q.type == "name"){
       # Use synonym-based fuzzy matching for compound names
       syn.db.path <- paste(lib.path, "lib/syn_nms.qs", sep="");

       if(!.on.public.web){
         syn.nm <- basename(syn.db.path);
         download.file(syn.db.path, destfile = syn.nm, method="libcurl", mode = "wb");
         syn.db.path <- syn.nm;
       }

       if(file.exists(syn.db.path)){
         # Load synonym database
         syn.db <- qs::qread(syn.db.path);

         # Perform fuzzy matching using synonyms
         match.results <- MatchCompoundNames(cmpd.vec, cmpd.map, syn.db);
         hit.inx <- match.results$hit.inx;

         # Print matching statistics
         matched.count <- sum(match.results$match.state == 1);
         total.count <- length(cmpd.vec);
         cat(sprintf("Compound name mapping (fuzzy): %d/%d matched (%.1f%%)\n",
                     matched.count, total.count, 100*matched.count/total.count));
       }else{
         # Fallback to exact matching if synonym database not found
         cat("Warning: Synonym database not found. Using exact matching.\n");
         hit.inx <- match(tolower(cmpd.vec), tolower(cmpd.map$name));
       }

    }else{
      print("No support for this compound database");
      return(0);
    }
 
    if(q.type == "name"){
       # For name type, return KEGG ID and original query name
       # (not the matched database name, to ensure proper matching in DataUtils.R)
       res_entrez <- data.frame(
         kegg_id = cmpd.map$kegg[hit.inx],
         query_name = cmpd.vec,
         stringsAsFactors = FALSE
       );
       rownames(res_entrez) <- cmpd.vec;
    }else if(q.type != "kegg"){
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
  }else if(type == "NA"){

    entrezs <- data.frame(gene_id=q.vec,accession=q.vec)
    rownames(entrezs) <- q.vec;
  }else {
    # Determine table name based on type
    table.nm <- NULL;
    if(type == "gb"){
      # note, some ID can have version number which is not in the database
      # need to strip it off NM_001402.5 => NM_001402
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      table.nm <- "entrez_gb";
    }else if(type == "refseq"){
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      table.nm <- "entrez_refseq";
    }else if(type == "embl_gene"){
      table.nm <- "entrez_embl_gene";
    }else if(type == "embl_transcript"){
      table.nm <- "entrez_embl_transcript";
    }else if(type == "embl_protein"){
      table.nm <- "entrez_embl_protein";
    }else if(type == "orf"){ # only for yeast
      table.nm <- "entrez_orf";
    }else if(type == "string"){
      table.nm <- "entrez_string";
    }else if(type == "ecogene"){ # only for ecoli
      table.nm <- "entrez_ecogene";
    }else if(type == "uniprot"){
      table.nm <- "entrez_uniprot";
    }else if(type == "flybase"){
      table.nm <- "entrez_flybase";
    }else{
      table.nm <- paste0("entrez_", type);
    }

    db.map <- queryGeneDB(table.nm, data.org);
    if(is.null(db.map)){
      return(data.frame(gene_id=character(0), accession=character(0)));
    }

    hit.inx <- match(q.vec, db.map[, "accession"]);
    entrezs <- db.map[hit.inx, ];
  }
 
  return(entrezs);
}


# mapping between genebank, refseq and entrez
doGeneIDMapping <- function(q.vec, type){
  if(is.null(q.vec)){
    db.map <- queryGeneDB("entrez", data.org);
    if(is.null(db.map)){
      return(character(0));
    }
    q.vec <- db.map[, "gene_id"];
    type = "entrez";
  }

  if(type == "symbol"){
    db.map <- queryGeneDB("entrez", data.org);
    if(is.null(db.map)){
      return(character(0));
    }
    hit.inx <- match(q.vec, db.map[, "symbol"]);
  }else if(type == "entrez"){
    db.map <- queryGeneDB("entrez", data.org);
    if(is.null(db.map)){
      return(character(0));
    }
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
  }else{
    # note, some ID can have version number which is not in the database
    # need to strip it off NM_001402.5 => NM_001402
    q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
    q.vec <- q.mat[,1];

    table.nm <- NULL;
    if(type == "gb"){
      table.nm <- "entrez_gb";
    }else if(type == "refseq"){
      table.nm <- "entrez_refseq";
    }else if(type == "embl_gene"){
      table.nm <- "entrez_embl_gene";
    }else if(type == "embl_transcript"){
      table.nm <- "entrez_embl_transcript";
    }else if(type == "embl_protein"){
      table.nm <- "entrez_embl_protein";
    }else if(type == "orf"){ # only for yeast
      table.nm <- "entrez_orf";
    }else{
      print("Unknown data type2");
      return(0);
    }

    db.map <- queryGeneDB(table.nm, data.org);
    if(is.null(db.map)){
      return(character(0));
    }
    hit.inx <- match(q.vec, db.map[, "accession"]);
  }
  entrezs=db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  rm(db.map, q.vec); gc();
  return(entrezs);
}

doEntrez2SymbolMapping <- function(entrez.vec){
  gene.map <- queryGeneDB("entrez", data.org);

  if(is.null(gene.map)){
    return(entrez.vec);
  }

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
  if(any(symbols1 == "" | is.na(symbols1))){
    pcms <- entrez.vec[(symbols1 == "" | is.na(symbols1))];
    entrez.vec2 <- gsub("PBCM0+","",pcms);
    hit.inx <- match(entrez.vec2, gene.map[, "pubchem_id"]);
    symbols2 <- gene.map[hit.inx, "kegg_id"];
    symbols2 <- as.character(symbols2);
    idxx <- which(symbols1 == "" | is.na(symbols1))
    symbols1[idxx] <- symbols2
  }
  # if not gene symbol, use id by itself

  symbols1[(is.na(symbols1) | symbols1=="")] <- entrez.vec[(is.na(symbols1) | symbols1=="")];
  return(symbols1);
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
PerformMetEnrichment <- function(dataSetObj=NA, file.nm, fun.type, ids, sourceView="2d", save.type="network"){
  dataSet <- .get.nSet(dataSetObj);
  dataSet$currentSourceView <- sourceView;
  normalizeEnrRes <- function(res){
    df <- as.data.frame(res);
    if(nrow(df) == 0){
      return(df);
    }
    if(!"Hits" %in% colnames(df)){
      if(ncol(df) >= 3){
        df$Hits <- df[,3];
      }
    }
    if(!"P.Value" %in% colnames(df)){
      if("Pval" %in% colnames(df)){
        df$P.Value <- df$Pval;
      }else if(ncol(df) >= 4){
        df$P.Value <- df[,4];
      }
    }
    return(df);
  }
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
    print(paste0("[PerformMetEnrichment] integ sizes: kegg=", nrow(as.data.frame(res1)),
                 " keggm=", nrow(as.data.frame(res2)),
                 " integ=", nrow(as.data.frame(res3))))
  }else{
    res1 = PerformEnrichAnalysisKegg(dataSet, file.nm, fun.type, ora.vec)
    res.mat <- as.data.frame(res1)
  }
  if(fun.type == "integ" || fun.type == "ginteg"){
    print(paste0("[PerformMetEnrichment] integ rownames: kegg=", length(rownames(res1)),
                 " keggm=", length(rownames(res2)),
                 " integ=", length(rownames(res3))))
    # Guard against missing metabolite results (or empty integ results)
    if(nrow(res1) == 0 || nrow(res2) == 0 || nrow(res3) == 0){
      if(nrow(res1) == 0){
        res.mat <- data.frame(
          hitsG = integer(0),
          hitsM = integer(0),
          hitsTotal = integer(0),
          P.ValueG = numeric(0),
          P.ValueM = numeric(0),
          P.ValueMerge = numeric(0),
          P.ValueJoint = numeric(0),
          integP = numeric(0)
        )
      }else{
        subres1 <- as.data.frame(res1)
        ord = order(rownames(subres1))
        subres1 = subres1[ord,]
        integ <- data.frame(
          hitsG = subres1$Hits,
          hitsM = rep(0, nrow(subres1)),
          hitsTotal = subres1$Hits,
          P.ValueG = subres1$P.Value,
          P.ValueM = rep(1, nrow(subres1)),
          P.ValueMerge = subres1$P.Value,
          P.ValueJoint = rep(1, nrow(subres1))
        )
        integ$integP <- integ$P.ValueG
        rownames(integ) <- rownames(subres1)
        res.mat <- integ
      }
      SaveIntegEnr(file.nm, res.mat, save.type);
      return(1);
    }
    # Normalize and filter NA rownames from each result
    res1n <- normalizeEnrRes(res1);
    if (nrow(res1n) > 0) {
      rn1 <- trimws(rownames(res1n))
      valid1 <- !is.na(rn1) & rn1 != "" & tolower(rn1) != "na"
      res1n <- res1n[valid1, , drop = FALSE]
    }

    res2n <- normalizeEnrRes(res2);
    if (nrow(res2n) > 0) {
      rn2 <- trimws(rownames(res2n))
      valid2 <- !is.na(rn2) & rn2 != "" & tolower(rn2) != "na"
      res2n <- res2n[valid2, , drop = FALSE]
    }

    res3n <- normalizeEnrRes(res3);
    if (nrow(res3n) > 0) {
      rn3 <- trimws(rownames(res3n))
      valid3 <- !is.na(rn3) & rn3 != "" & tolower(rn3) != "na"
      res3n <- res3n[valid3, , drop = FALSE]
    }

    all_paths <- unique(na.omit(c(rownames(res1n), rownames(res2n), rownames(res3n))))
    if (length(all_paths) > 0) {
      all_paths <- trimws(all_paths)
      all_paths <- all_paths[!is.na(all_paths) & all_paths != "" & tolower(all_paths) != "na"]
    }
    print(paste0("[PerformMetEnrichment] integ total pathways: ", length(all_paths)))
    if(length(all_paths) == 0){
      res.mat <- data.frame(
        hitsG = integer(0),
        hitsM = integer(0),
        hitsTotal = integer(0),
        P.ValueG = numeric(0),
        P.ValueM = numeric(0),
        P.ValueMerge = numeric(0),
        P.ValueJoint = numeric(0),
        integP = numeric(0)
      )
      SaveIntegEnr(file.nm, res.mat, save.type);
      return(1);
    }

    get_col <- function(df, col, default) {
      if (length(default) > 1) {
        if (length(default) != length(all_paths)) {
          default <- rep(default[1], length(all_paths))
        }
        out <- as.numeric(default)
        names(out) <- all_paths
      } else {
        out <- rep(default, length(all_paths))
        names(out) <- all_paths
      }
      if (is.null(df) || nrow(df) == 0 || is.null(rownames(df))) {
        return(out)
      }
      vec_raw <- df[[col]]
      if (is.null(vec_raw)) {
        return(out)
      }
      if (is.data.frame(vec_raw) || is.matrix(vec_raw)) {
        vec_raw <- vec_raw[, 1]
      }
      vec <- as.vector(vec_raw)
      vec <- as.numeric(vec)
      names(vec) <- rownames(df)
      idx <- match(names(vec), all_paths)
      ok <- !is.na(idx)
      if (any(ok)) {
        out[idx[ok]] <- vec[ok]
      }
      return(out)
    }

    hitsG <- get_col(res1n, "Hits", 0)
    hitsM <- get_col(res2n, "Hits", 0)
    hitsTotal <- get_col(res3n, "Hits", hitsG + hitsM)
    pG <- get_col(res1n, "P.Value", 1)
    pM <- get_col(res2n, "P.Value", 1)
    pMerge <- get_col(res3n, "P.Value", pmin(pG, pM, na.rm = TRUE))

    integ <- data.frame(
      hitsG = as.numeric(hitsG),
      hitsM = as.numeric(hitsM),
      hitsTotal = as.numeric(hitsTotal),
      P.ValueG = as.numeric(pG),
      P.ValueM = as.numeric(pM),
      P.ValueMerge = as.numeric(pMerge),
      P.ValueJoint = as.numeric(pM)
    )
    integ <- as.data.frame(integ, stringsAsFactors = FALSE)
    if (length(all_paths) == nrow(integ)) {
      rownames(integ) <- all_paths
    } else {
      message(sprintf("[PerformMetEnrichment] WARNING: pathway count mismatch (paths=%d, rows=%d) - using sequential rownames",
                      length(all_paths), nrow(integ)))
      rownames(integ) <- seq_len(nrow(integ))
    }
    # Drop any NA/empty pathway names from output
    rn <- trimws(rownames(integ))
    bad_nm <- is.na(rn) | rn == "" | tolower(rn) == "na"
    if (any(bad_nm)) {
      integ <- integ[!bad_nm, , drop = FALSE]
    }

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
    SaveIntegEnr(file.nm,res.mat, save.type);
  }else{
    SaveSingleOmicsEnr(file.nm,res.mat, save.type);
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

  # Filter out NA/empty pathway names BEFORE creating result matrix
  geneset_names <- names(current.geneset)
  if (length(geneset_names) > 0) {
    geneset_names <- trimws(geneset_names)
    valid_idx <- !is.na(geneset_names) & geneset_names != "" & tolower(geneset_names) != "na"
    current.geneset <- current.geneset[valid_idx]
  }

  # prepare for the result table
  set.size<-length(current.geneset);
  res.mat<-matrix(0, nrow=set.size, ncol=5);
  rownames(res.mat)<-names(current.geneset);
  colnames(res.mat)<-c("Total", "Expected", "Hits", "Pval", "FDR");

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

SaveSingleOmicsEnr <- function(file.nm,res.mat, save.type="network"){
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
  if (nrow(resTable) > 0) {
    rn <- trimws(as.character(resTable$Pathway))
    keep <- !is.na(rn) & rn != "" & tolower(rn) != "na"
    if (!all(keep)) {
      resTable <- resTable[keep, , drop = FALSE]
      res.mat <- res.mat[keep, , drop = FALSE]
    }
  }
  current.msg <<- "Functional enrichment analysis was completed";

  # write json
  require("RJSONIO");
  fun.anot <- hits.query;
  fun.padj <- resTable$FDR; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj)};
  fun.pval <- resTable$Pval; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval)};
  hit.num <- paste0(resTable$Hits,"/",resTable$Total); if(length(hit.num) ==1) { hit.num <- matrix(hit.num)};
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

SaveIntegEnr <- function(file.nm,res.mat, save.type="network"){
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
  fun.nms <- resTable[,"Pathway"];  if(length(fun.nms) ==1) { fun.nms <- matrix(fun.nms) };
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

  # Determine which SQLite file to use based on database type
  if(grepl("string$", table.nm)){
    db.path <- paste(sqlite.path, "ppi_string.sqlite", sep="");
  }else if(grepl("intact$", table.nm)){
    db.path <- paste(sqlite.path, "ppi_intact.sqlite", sep="");
  }else if(grepl("huri$", table.nm)){
    db.path <- paste(sqlite.path, "ppi_huri.sqlite", sep="");
  }else if(grepl("innate$", table.nm)){
    # Check if separate InnateDB file exists, otherwise fall back to ppi.sqlite
    innate.path <- paste(sqlite.path, "ppi_innate.sqlite", sep="");
    if(file.exists(innate.path)){
      db.path <- innate.path;
    }else{
      db.path <- paste(sqlite.path, "ppi.sqlite", sep="");
    }
  }else{
    # For other databases (rolland, irefinx, interactome, tissue, etc.), use original ppi.sqlite
    db.path <- paste(sqlite.path, "ppi.sqlite", sep="");
  }

  # DEBUG: Print database path and query info
  #print(paste("DEBUG QueryPpiSQLite: Database file =", db.path));
  #print(paste("DEBUG QueryPpiSQLite: File exists =", file.exists(db.path)));
  #print(paste("DEBUG QueryPpiSQLite: Table name =", table.nm));
  #print(paste("DEBUG QueryPpiSQLite: Number of query genes =", length(q.vec)));
  #print("DEBUG QueryPpiSQLite: First 10 query gene IDs:");
  #print(head(q.vec, 10));
  #print(paste("DEBUG QueryPpiSQLite: requireExp =", requireExp));
  #print(paste("DEBUG QueryPpiSQLite: min.score =", min.score));

  ppi.db <- .connect.sqlite(db.path);

  # DEBUG: Check what IDs are actually in the database
  test_query <- paste("SELECT DISTINCT id1 FROM", table.nm, "LIMIT 10");
  sample_ids <- tryCatch({
    .query.sqlite(ppi.db, test_query, FALSE);
  }, error = function(e) {
    data.frame(id1=c("error"));
  });
  #print("DEBUG QueryPpiSQLite: Sample id1 values from database:");
  #print(sample_ids$id1);

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

  # DEBUG: Print SQL query
  #print("DEBUG QueryPpiSQLite: SQL query:");
  #print(statement);

  # Try to execute query with error handling for missing tables
  ppi.res <- tryCatch({
    .query.sqlite(ppi.db, statement);
  }, error = function(e) {
    dbDisconnect(ppi.db);
    if(grepl("no such table", e$message, ignore.case = TRUE)){
      stop(paste0("Database table '", table.nm, "' not found in ", basename(db.path),
                  ". Please verify that the selected database (", table.nm,
                  ") is available for your organism."));
    }else{
      stop(paste0("Database query failed: ", e$message));
    }
  });

  # DEBUG: Print query results
  #print(paste("DEBUG QueryPpiSQLite: Rows returned (before dedup) =", nrow(ppi.res)));
  #print("DEBUG QueryPpiSQLite: Column names in result:");
  #print(colnames(ppi.res));
  if(nrow(ppi.res) > 0){
    #print("DEBUG QueryPpiSQLite: First row sample:");
    #print(ppi.res[1,]);
  }

  # remove duplicated edges
  # Check if row_id has valid values (not all NA/empty)
  if("row_id" %in% colnames(ppi.res) && any(!is.na(ppi.res$row_id) & ppi.res$row_id != "")){
    # row_id has valid values, use it for deduplication
    ppi.res <- ppi.res[!duplicated(ppi.res$row_id),]
  }else{
    # row_id is empty/NULL, deduplicate by id1+id2 pair
    ppi.res <- ppi.res[!duplicated(paste(ppi.res$id1, ppi.res$id2, sep="_")),]
  }

  #print(paste("DEBUG QueryPpiSQLite: Rows returned (after dedup) =", nrow(ppi.res)));

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



.query.sqlite <- function(db.con, statement, offline=TRUE){
  rs <- dbSendQuery(db.con, statement);
  res <- fetch(rs, n=-1); # get all records
  dbClearResult(rs);
  if(offline){
    dbDisconnect(db.con);
  }
  cleanMem();
  return(res);
}

# OmniPath PPI query function
QueryOmnipathPpiSQLite <- function(sqlite.path, data.org, q.vec){
  require('RSQLite');
  db.path <- paste(sqlite.path, "omnipath_ppi.sqlite", sep="");
  con <- .connect.sqlite(db.path);
  query <- paste(shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", data.org, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))", sep="");
  return(.query.sqlite(con, statement));
}

# OmniPath TF-gene query function
QueryOmnipathTfSQLite <- function(sqlite.path, data.org, q.vec){
  require('RSQLite');
  db.path <- paste(sqlite.path, "omnipath_tf2gene.sqlite", sep="");
  con <- .connect.sqlite(db.path);
  query <- paste(shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", data.org, " WHERE entrez IN (",query,")", sep="");
  return(.query.sqlite(con, statement));
}

# OmniPath miRNA-gene query function
QueryOmnipathMirSQLite <- function(sqlite.path, data.org, q.vec){
  require('RSQLite');
  db.path <- paste(sqlite.path, "omnipath_mir2gene.sqlite", sep="");
  con <- .connect.sqlite(db.path);
  query <- paste(shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", data.org, " WHERE entrez IN (", query, ")", sep="");
  return(.query.sqlite(con, statement));
}

.connect.sqlite <- function(db.path){
  return(dbConnect(SQLite(), db.path));
}
