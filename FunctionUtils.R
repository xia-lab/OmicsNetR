##################################################
## R script for OmicsNet
## Description: ID mapping, GO/Pathway enrichment analysis
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

LoadMotifLib<-function(){
  
  motif.path <- paste(lib.path, data.org, "/motif_set.rds", sep="");
  motif_set<-readRDS(motif.path);
  current.mset <- motif_set$set;
  set.ids<- names(current.mset); 
  names(set.ids) <- names(current.mset) <- motif_set$term;
  current.setlink <<- motif_set$link;
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}


LoadKEGGKO_lib<-function(category){
  
  kegg.rda <- paste0(lib.path, "microbiome/ko_pathways.rda");
  load(kegg.rda);
  current.setlink <- kegg.anot$link;
  current.mset <- kegg.anot$sets$Metabolism;
  
  if(!exists("ko.edge.map")){
    ko.edge.path <- paste0(lib.path, "microbiome/ko_edge.csv");
    ko.edge.map <<- .readDataTable(ko.edge.path); 
    
  } 
  
  kos.01100 <- ko.edge.map$gene[ko.edge.map$net == "ko01100"];
  current.mset <- lapply(current.mset, 
                         function(x) {
                           as.character(unique(x[x %in% kos.01100]));
                         }
  );
  # remove those empty ones
  mset.ln <- lapply(current.mset, length);
  current.mset <- current.mset[mset.ln > 0];
  set.ids<- names(current.mset); 
  names(set.ids) <- names(current.mset) <- kegg.anot$term[set.ids];
  
  
  current.setlink <<- current.setlink;
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}


LoadKEGGLib<-function(){
  kegg.path = ""
  if(isKo){
    kegg.path <- paste(lib.path, "microbiome/koset.rds", sep="");
    kegg.anot <- readRDS(kegg.path)
    current.setlink <- " "
    current.mset <- kegg.anot;
    current.setlink <<- current.setlink;
    current.setids <<- names(kegg.anot);
    current.geneset <<- current.mset;
    current.universe <<- unique(unlist(current.geneset));
  }else{
    kegg.path <- paste(lib.path, data.org, "/kegg1.rds", sep="");
    kegg.anot <- readRDS(kegg.path)
    current.setlink <- kegg.anot$link;
    current.mset <- kegg.anot$sets;
    set.ids<- names(current.mset); 
    names(set.ids) <- names(current.mset) <- kegg.anot$term;
    current.setlink <<- current.setlink;
    current.setids <<- set.ids;
    current.geneset <<- current.mset;
    current.universe <<- unique(unlist(current.geneset));
  }
  
}

LoadREACTOMELib<-function(){
  reactome.path <- paste(lib.path, data.org, "/reactome.rds", sep="");
  reactome.anot <- readRDS(reactome.path)
  current.mset <- reactome.anot$sets;
  set.ids<- names(current.mset); 
  names(set.ids) <- names(current.mset) <- reactome.anot$term;
  current.setlink <<- reactome.anot$link;
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}

LoadMirLib<-function(name){
  reactome.path <- paste(lib.path, data.org, "/", name, ".rds", sep="");
  reactome.anot <- readRDS(reactome.path)
  current.mset <- reactome.anot$sets;
  set.ids<- names(current.mset); 
  names(set.ids) <- names(current.mset) <- reactome.anot$term;
  current.setlink <<- reactome.anot$link;
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}

LoadMsetLib<-function(name){
  reactome.path <- paste(lib.path, "msets", "/", name, ".rds", sep="");
  reactome.anot <- readRDS(reactome.path)
  current.mset <- reactome.anot$sets;
  set.ids<- names(current.mset); 
  names(set.ids) <- names(current.mset) <- reactome.anot$term;
  current.setlink <<- reactome.anot$link;
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}


LoadOtherLib<-function(onto){
  go.path <- paste(lib.path, data.org, "/", tolower(onto),"_enr.rds", sep="");
  if(tolower(onto) == "mir"){
    mir_enr <- readRDS(go.path);
    if(is.null(names(mir_enr))){ # new go lib does not give names
      names(mir_enr) <- c("link", "term", "sets");
    }
    current.link <- mir_enr$link;           
    current.term <- mir_enr$term;
    current.mset = mir_enr$sets
    set.ids<- names(current.mset); 
    names(set.ids) <- names(current.mset) <- current.term;
  }else if(tolower(onto) == "drug"){
    drug_enr <- readRDS(go.path);
    if(is.null(names(drug_enr))){
      names(drug_enr) <- c("link", "term", "sets");
    }
    current.link <- drug_enr$link;
    current.term <- drug_enr$term
    current.mset = drug_enr$sets
    set.ids<- names(current.mset); 
    names(set.ids) <- names(current.mset) <- current.term;
  }else if(tolower(onto) == "jaspar"){
    jaspar_enr <- readRDS(go.path);
    if(is.null(names(jaspar_enr))){
      names(jaspar_enr) <- c("link", "term", "sets");
    }
    current.link <- jaspar_enr$link;
    current.mset <- jaspar_enr$sets;
    set.ids<- names(current.mset); 
    names(set.ids) <- names(current.mset) <- jaspar_enr$term;
  }else{
    encode_enr <- readRDS(go.path);
    if(is.null(names(encode_enr))){
      names(encode_enr) <- c("link", "term", "sets");
    }
    current.link <- encode_enr$link;
    current.mset <- encode_enr$sets;
    set.ids<- names(current.mset); 
    names(set.ids) <- names(current.mset) <- encode_enr$term;
  }
  current.setlink <<- current.link;
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}

LoadGOLib<-function(onto){
  go.path <- paste(lib.path, data.org, "/", tolower(onto), ".rds", sep="");
  if(tolower(onto) %in% c("go_panth","go_bp")){
    go_bp <- readRDS(go.path);
    
    if(is.null(names(go_bp))){ # new go lib does not give names
      names(go_bp) <- c("link", "term", "sets");
    }
    current.link <- go_bp$link;
    current.mset <- go_bp$sets;
    set.ids<- names(current.mset); 
    names(set.ids) <- names(current.mset) <- go_bp$term;
  }else if(tolower(onto) == "go_mf"){
    go_mf <- readRDS(go.path);
    if(is.null(names(go_mf))){
      names(go_mf) <- c("link", "term", "sets");
    }
    current.link <- go_mf$link;
    current.mset <- go_mf$sets;
    set.ids<- names(current.mset); 
    names(set.ids) <- names(current.mset) <- go_mf$term;
  }else{
    go_cc <- readRDS(go.path);
    if(is.null(names(go_cc))){
      names(go_cc) <- c("link", "term", "sets");
    }
    current.link <- go_cc$link;
    current.mset <- go_cc$sets;
    set.ids<- names(current.mset); 
    names(set.ids) <- names(current.mset) <- go_cc$term;
  }
  current.setlink <<- current.link;
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}

SearchReg <- function(file.nm, fun.type, IDs, count){
  regCount <<- count;
  regBool <<- "true";
  res <- PerformNetEnrichment(file.nm, fun.type, IDs);
}

regEnrichment <- function(file.nm, fun.type, IDs, netInv){
  regBool <<- "false";
  res <- PerformNetEnrichment(file.nm, fun.type, IDs, netInv);
}

# note: hit.query, resTable must synchronize
PerformNetEnrichment <- function(file.nm, fun.type, IDs, netInv){
  # prepare query
  ora.vec <- NULL;
  if(omics.net$db.type == 'ppi'){
    if(data.org == "ath"){
      idtype <- "embltranscript"
    }else if(data.org %in% c("bsu", "tbr", "cel", "dme", "eco", "pfa") & net.type == "string"){
      idtype <- "string"
    }else if(data.org %in% c("bta","dre","rno","gga","hsa","mmu") & net.type == "string"){
      idtype <- "emblprotein"
    }else if(data.org %in% c("hsa","mmu", "cel", "dme","sce") & net.type %in% c("innate", "irefinx", "rolland")){
      idtype <- "entrez"
    }else if(data.org == "sce" & net.type == "string"){ # only for yeast
      idtype <- "emblgene";
    }
    if(idtype=="uniprot"){
      uniprot.vec <- unlist(strsplit(IDs, "; "));
      ora.vec <- doUniprot2EntrezMapping(uniprot.vec);
      names(ora.vec) <- uniprot.vec;
    }else if(idtype=="emblprotein"){
      emblprotein.vec <- unlist(strsplit(IDs, "; "))
      ora.vec <- doEmblProtein2EntrezMapping(emblprotein.vec);
      names(ora.vec) <- emblprotein.vec;
    }else if(idtype=="string"){
      string.vec <- unlist(strsplit(IDs, "; "))
      ora.vec <- doString2EntrezMapping(string.vec);
      names(ora.vec) <- string.vec;
    }else if(idtype=="emblgene"){
      emblgene.vec <- unlist(strsplit(IDs, "; "))
      ora.vec <- doEmblGene2EntrezMapping(emblgene.vec);
      names(ora.vec) <- emblgene.vec;
    }else{
      ora.vec <- unlist(strsplit(IDs, "; "));
      names(ora.vec) <- ora.vec;
    }  
  }else{
    ora.vec <- unlist(strsplit(IDs, "; "));
    names(ora.vec) <- ora.vec;
  }
  if(fun.type %in% c("trrust", "encode", "jaspar", "mirnet", "met", "drugbank")){
    res <- PerformRegEnrichAnalysis(file.nm, fun.type, ora.vec, netInv, idType);
  } else{
    res <- PerformEnrichAnalysis(file.nm, fun.type, ora.vec);
  }
  return(res);   
}

PerformRegEnrichAnalysis <- function(file.nm, fun.type, ora.vec, netInv, idType){
  require(plyr)
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
    res <- QueryMirSQLite(table.nm, "entrez", ora.vec, "inverse", mir.type);
    if(nrow(res)==0){ return(c(0,0)); }
    edge.res <- data.frame(gene=res[,"entrez"], symbol=res[,"symbol"], id=res[,"mir_acc"], name=res[,"mir_id"] );
    node.ids <- c(res[,"entrez"], res[,"mir_acc"])
    node.nms <- c(res[,"symbol"], res[,"mir_id"]);
    
  }else if(fun.type == "met"){ # in miRNA, table name is org code, colname is id type
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
    resTable1 = resTable1[which(resTable1$Hits == regCount),]
  }
  # write json
  require(RJSONIO);
  fun.ids <- resTable1[,1]; 
  fun.nms <- resTable1[,2];
  fun.hits <- resTable1[,3]; 
  json1.res <- list(
    fun.ids = fun.ids, 
    fun.nms = fun.nms,
    fun.hits = fun.hits,
    fun.genes = hits.gene 
  );
  json.mat <- toJSON(json1.res, .na='null');
  json.nm <- paste(file.nm, ".json", sep="");
  sink(json.nm)
  cat(json.mat);
  sink();
  
  # write csv
  csv.nm <- paste(file.nm, ".csv", sep="");
  fast.write.csv(resTable1, file=csv.nm, row.names=F);
  
  return(1);
}

# note: hit.query, resTable must synchronize
# ora.vec should contains entrez ids, named by their gene symbols
PerformEnrichAnalysis <- function(file.nm, fun.type, ora.vec){
  # prepare lib
  if(tolower(fun.type) == 'kegg'){ 
    LoadKEGGLib();
  }else if(tolower(fun.type) == 'ko'){ 
    LoadKEGGKO_lib();
  }else if(tolower(fun.type) == 'integ'){ 
    LoadKEGGLibOther("integ");
  }else if(tolower(fun.type) == 'keggm'){ 
    LoadKEGGLibOther("keggm");
  }else if(tolower(fun.type) == 'reactome'){ 
    LoadREACTOMELib();
  }else if(tolower(fun.type) == 'mirfamily'){ # when user choose to perform miRNA family enrichment analysis.
    LoadmiRFamLib();
  }else if(tolower(fun.type) == 'motif'){ 
    LoadMotifLib();
  }else if(tolower(fun.type) %in% c("bp", "cc", "mf","panth_mf", "go_bp", "go_cc", "go_panthbp", "go_panthcc", "go_panthmf")){ 
    LoadGOLib(fun.type);
  }else if(tolower(fun.type) %in% c("mirfamily", "mircluster", "mirdisease","mirtissue", "mirfunction")){
    LoadMirLib(fun.type);
  }else if(tolower(fun.type) %in% c("pathway", "blood", "urine", "csf", "predicted", "location", "drug")){
    LoadMsetLib(fun.type);
  }else if(tolower(fun.type) %in% c("mir", "jasper", "drug", "encode")){
    LoadOtherLib(fun.type);
  }else{
    print(paste("Unknown lib option:", fun.type));
    return(0);
  }

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
  require(RJSONIO);
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
  json.mat <- toJSON(json.res, .na='null');
  json.nm <- paste(file.nm, ".json", sep="");
  sink(json.nm)
  cat(json.mat);
  sink();
  
  # write csv
  csv.nm <- paste(file.nm, ".csv", sep="");
  a <- lapply(fun.anot,function(x){paste(x,collapse = "; ")})
  resTable$ids <- unlist(a);
  fast.write.csv(resTable, file=csv.nm, row.names=F);
  
  return(1);
}

doProteinIDMapping <- function(q.vec, type){
  if(type %in% c("rsid")){
    hit.inx <- startsWith(q.vec, "rs");
    q.vec <- q.vec[hit.inx]
    res <- data.frame(gene_id=q.vec, accession=q.vec);
    entrezs <- res;
  }else if(type %in% c("class", "family","order" ,"genus","species","strain", "phylum")){
    require('RSQLite');
    mic.taxa <<- type;
    db.path <- paste(lib.path, "microbiome", "/taxaInfo.rds", sep="");
    taxlist <- readRDS(db.path);
    db.map <- taxlist[[type]]
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
    db.map <-  readRDS(db.path);
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
    entrezs <- db.map[hit.inx, ]
    entrezs <- entrezs[c(1,2)]
    colnames(entrezs) = c("gene_id", "accession");
  }else if(type %in% c("mir_acc", "mir_id", "mirnet")){
    require('RSQLite');
    path <- paste0(sqlite.path, "mir2gene.sqlite")
    mir.db <- dbConnect(SQLite(), path);
    query <- paste (shQuote(q.vec),collapse=",");
    table.nm <- data.org
    statement <- paste("SELECT * FROM ", table.nm, " WHERE ", type," IN (",query,")", sep="");
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
    cmpd.map = readRDS(db.path)
    q.type = type;
    cmpd.vec = q.vec
    
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
    table.nm <- paste(data.org, tf.type, sep="_");
    mir.dic <- QueryTfSQLite(table.nm, q.vec, "inverse");
    res <- mir.dic[ , c("tfid", "tfid")];
    res = res[!duplicated(res),]
    
    colnames(res) <- c("accession", "gene_id")
    hit.inx <- match(q.vec, res[, "accession"]);
    entrezs <- res[hit.inx, ];
    entrezs = res[c(1,2)];  
  }else {
    if(type == "genbank"){
      # note, some ID can have version number which is not in the database
      # need to strip it off NM_001402.5 => NM_001402
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      db.path <- paste(lib.path, data.org, "/entrez_gb.rds", sep="");
    }else if(type == "refseq"){
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      db.path <- paste(lib.path, data.org, "/entrez_refseq.rds", sep="");
    }else if(type == "emblgene"){
      db.path <- paste(lib.path, data.org, "/entrez_embl_gene.rds", sep="");
    }else if(type == "embltranscript"){
      db.path <- paste(lib.path, data.org, "/entrez_embl_transcript.rds", sep="");
    }else if(type == "emblprotein"){
      db.path <- paste(lib.path, data.org, "/entrez_embl_protein.rds", sep="");
    }else if(type == "orfid"){ # only for yeast
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
    db.map <-  readRDS(db.path);
    hit.inx <- match(q.vec, db.map[, "symbol"]);
  }else if(type == "entrez"){
    db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");
    db.map <-  readRDS(db.path);
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
  }else{
    # note, some ID can have version number which is not in the database
    # need to strip it off NM_001402.5 => NM_001402
    q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
    q.vec <- q.mat[,1];
    if(type == "genbank"){
      db.path <- paste(lib.path, data.org, "/entrez_gb.rds", sep="");
    }else if(type == "refseq"){
      db.path <- paste(lib.path, data.org, "/entrez_refseq.rds", sep="");
    }else if(type == "emblgene"){
      db.path <- paste(lib.path, data.org, "/entrez_embl_gene.rds", sep="");
    }else if(type == "embltranscript"){
      db.path <- paste(lib.path, data.org, "/entrez_embl_transcript.rds", sep="");
    }else if(type == "emblprotein"){
      db.path <- paste(lib.path, data.org, "/entrez_embl_protein.rds", sep="");
    }else if(type == "orfid"){ # only for yeast
      db.path <- paste(lib.path, data.org, "/entrez_orf.rds", sep="");
    }else{
      print("Unknown data type2");
      return(0);
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
  gene.map <- readRDS(db.path);
  
  hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
  symbols <- gene.map[hit.inx, "symbol"];
  
  # if not gene symbol, use id by itself
  na.inx <- is.na(symbols);
  symbols[na.inx] <- entrez.vec[na.inx];
  return(symbols);
}


doKegg2NameMapping <- function(entrez.vec){
  db.path <- paste(lib.path, "microbiome/keggConv.rds", sep="");
  gene.map <- readRDS(db.path);
  
  hit.inx <- match(entrez.vec, gene.map[, "kegg"]);
  symbols <- gene.map[hit.inx, "name"];
  
  # if not gene symbol, use id by itself
  na.inx <- is.na(symbols);
  symbols[na.inx] <- entrez.vec[na.inx];
  return(symbols);
}

# note, entrez.vec could contain NA/null, cannot use rownames
doEntrezIDAnot <- function(entrez.vec){
  db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");
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
  db.map <-  readRDS(db.path);
  hit.inx <- match(uniprot.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];    
  mode(entrezs) <- "character";
  return(entrezs);
}

doEntrez2UniprotMapping <- function(entrez.vec){
  db.path <- paste(lib.path, data.org, "/entrez_uniprot.rds", sep="");
  db.map <-  readRDS(db.path);
  hit.inx <- match(entrez.vec, db.map[, "gene_id"]);
  unips <- db.map[hit.inx, "accession"];    
  mode(unips) <- "character";
  return(unips);
}

doString2EntrezMapping <- function(string.vec){
  db.path <- paste(lib.path, data.org, "/entrez_string.rds", sep="");
  db.map <-  readRDS(db.path);
  hit.inx <- match(string.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}

doEmblGene2EntrezMapping <- function(emblgene.vec){
  db.path <- paste(lib.path, data.org, "/entrez_embl_gene.rds", sep="");
  db.map <-  readRDS(db.path);
  hit.inx <- match(emblgene.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}

doSymbol2EntrezMapping <- function(symbol.vec){
  db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");
  db.map <-  readRDS(db.path);
  hit.inx <- match(symbol.vec, db.map[, "symbol"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}

doEmblProtein2EntrezMapping <- function(emblprotein.vec){
  db.path <- paste(lib.path, data.org, "/entrez_embl_protein.rds", sep="");
  db.map <-  readRDS(db.path);
  hit.inx <- match(emblprotein.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}

doEntrez2EmblProteinMapping <- function(entrez.vec){
  db.path <- paste(lib.path, data.org, "/entrez_embl_protein.rds", sep="");
  db.map <-  readRDS(db.path);
  hit.inx <- match(entrez.vec, db.map[, "gene_id"]);
  entrezs <- db.map[hit.inx, "accession"];
  entrezs <- as.vector(entrezs)
  return(entrezs);
}

doPpiIDMapping <- function(q.vec, direction = "id") {
  if (data.org == "ath") {
    db.path <-paste(lib.path, data.org, "/entrez_embl_transcript.rds", sep = "")
  } else if (data.org == "sce") {
    if (net.type == "string") {
      # only for yeast
      db.path <-
        paste(lib.path, data.org, "/entrez_embl_gene.rds", sep = "")  
    } else{
      db.path <- paste(lib.path, data.org, "/entrez_uniprot.rds", sep = "")     
    }
  } else if (net.type %in% c("innate", "irefinx", "rolland")) {
    db.path <- paste(lib.path, data.org, "/entrez_uniprot.rds", sep = "")   
  } else{
    if (data.org %in% c("bsu", "tbr", "cel", "dme", "eco", "pfa")) {
      db.path <- paste(lib.path, data.org, "/entrez_string.rds", sep = "")    
    } else if (data.org %in% c("mmu", "hsa") && net.type == "string") {
      db.path <- paste(lib.path, data.org, "/entrez.rds", sep = "")     
    } else {
      db.path <-
        paste(lib.path, data.org, "/entrez_embl_protein.rds", sep = "")
    }
  }
  db.map <-  readRDS(db.path);
  if (direction == "id") {
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
    ppi.mat <- db.map[hit.inx,];  
    # fix the factor col related to library issue
    i <- sapply(ppi.mat, is.factor)
    ppi.mat[i] <- lapply(ppi.mat[i], as.character)
    ppi.mat = as.vector(ppi.mat$accession)
  } else{
    hit.inx <- match(q.vec, db.map[, "accession"])
    ppi.mat <- db.map[hit.inx,]   
    # fix the factor col related to library issue
    i <- sapply(ppi.mat, is.factor)
    ppi.mat[i] <- lapply(ppi.mat[i], as.character)
    ppi.mat = as.vector(ppi.mat$gene_id)
  }
  return(ppi.mat)
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
    if(isKo){
      kegg.path <- paste(lib.path, "integ.rds", sep="");
    }else{
      kegg.path <- paste(lib.path, data.org, "/integ.rds", sep="");
    }
  }else if(type == "ginteg"){
    kegg.path <- paste(lib.path,"microbiome", "/integ.rds", sep="");
  }else if(type == "gkeggm"){
    kegg.path <- paste(lib.path,"microbiome",  "/keggm.rds", sep="");
  }else{
    if(isKo){
      kegg.path <- paste(lib.path,"microbiome", "/keggm.rds", sep="");
    }else{
      kegg.path <- paste(lib.path, data.org,"/keggm.rds", sep="");
    }
  }
  kegg.anot <- readRDS(kegg.path)
  current.setlink <- " "
  current.mset <- kegg.anot;
  current.setlink <<- current.setlink;
  current.setids <<- names(kegg.anot);
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.geneset));
}

PerformKeggEnrichment <- function(dataSetObj=NA, file.nm, fun.type, ids){
 dataSet <- .get.nSet(dataSetObj);  
 if(ids=="Not_applicable"){
   ora.vec <- keggp.allfeatures
   names(ora.vec) <- ora.vec;    
 } else {
   ora.vec <- unlist(strsplit(ids, "; "));
   names(ora.vec) <- ora.vec;
   ora.vecu <<- ora.vec
 }
  noGene <- F
  noMet <- F
  if(length(rownames(dataSet$seed[["ko"]]))== 0 && length(rownames(dataSet$seed[["gene"]]))== 0){
    noGene <- T
  }else if(length(rownames(dataSet$met.seed))== 0){
    noMet <- T
  }
  
  if(fun.type == "keggm" && noMet){
    return(2)
  }else if(fun.type != "keggm" && fun.type != "integ" && noGene){
    return(3)
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
  
  containsKo = length(rownames(dataSet$seed[["ko"]]))>0
  containsGene = length(rownames(dataSet$seed[["gene"]]))>0
  containsMet = length(rownames(dataSet$met.seed))>0
  if(containsMet && (containsGene || containsKo)){
    met.vec <- dataSet$met.seed
    if(isKo){
      gene.vec <- dataSet$seed[["ko"]]
    }else{
      gene.vec <- dataSet$seed[["gene"]]
    }
    res <- .pathwayMetGenePair(met.vec, gene.vec)
    nohits.vec = vector()
    
    for(i in 1:length(hits.query)){
      query.vec = hits.query[[i]]
      table.name = names(hits.query)[i]
      res1 = res[which(res[,2] %in% query.vec),]
      res2 = res[which(res[,5] %in% query.vec),]

      combinedres = rbind(res1, res2)
      combinedres = unique(combinedres)
      if(length(rownames(combinedres)) >0){
        combinedres$pathway = table.name
        if(i==1){
          alltables = combinedres
        }else{
          alltables = rbind(alltables, combinedres)
        }
      }
    }
    fast.write.csv(alltables, paste0("pairs_", file.nm, ".csv"), row.names = F)
  }
  return(.set.nSet(dataSetObj));
}

doEntrez2KoMapping <- function(entrez.vec){
  db.path <- paste(lib.path,data.org,"/koConv.rds", sep="");
  db.map <-  readRDS(db.path);
  hit.inx <- match(entrez.vec, db.map[, "entrez"]);
  entrez.vec[hit.inx] <- db.map[hit.inx, "ko"];
  return(entrez.vec);
}

PerformEnrichAnalysisKegg <- function(dataSetObj=NA, file.nm, fun.type, ora.vec){
  dataSet <- .get.nSet(dataSetObj); 
  # prepare lib
  if(tolower(fun.type) == 'kegg'){ 
    LoadKEGGLib();
  }else if(tolower(fun.type) == 'ko'){ 
    LoadKEGGKO_lib();
  }else if(tolower(fun.type) == 'reactome'){ 
    LoadREACTOMELib();
  }else if(tolower(fun.type) == 'motif'){ 
    LoadMotifLib();
  }else if(tolower(fun.type) %in% c("bp", "cc", "mf","panth")){ 
    LoadGOLib(fun.type);
  }else{
    LoadKEGGLibOther(fun.type)
  }
  
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
  require(RJSONIO);
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
  json.mat <- toJSON(json.res, .na='null');
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
  require(RJSONIO);
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
  json.mat <- toJSON(json.res, .na='null');
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
  if(isKo || data.org == "microbiome"){
    table.nm <- "ko_network"
    resMet <- QueryKoSQLiteNet(table.nm, rownames(met.vec), "inverse");
    resGene <- QueryKoSQLiteNet(table.nm, rownames(gene.vec), "direct");
    res = rbind(resMet, resGene)
    colnames(res) = c("entrez", "kegg", "met", "symbol")
  }else{
    table.nm <- paste(data.org, "kegg", sep="_");
    resMet <- QueryMetSQLiteNet(table.nm, rownames(met.vec), "direct");
    resGene <- QueryMetSQLiteNet(table.nm, rownames(gene.vec), "inverse");
    res <- rbind(resMet, resGene)
  }
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
    path <- paste0(sqlite.path, "ppi.sqlite")
    ppi.db <- dbConnect(SQLite(), path);
    query <- paste(shQuote(q.vec),collapse=",");
    if(grepl("string$", table.nm)){
        if(netbuild.opt == "pcsf"){
         statement <- paste("SELECT * FROM ",table.nm, " WHERE combined_score >=", min.score, " AND experimental > 0", sep="");
        }
        else if(requireExp){
            statement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))  AND combined_score >=", min.score, " AND experimental > 0", sep="");
        }else{
            statement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))  AND combined_score >=", min.score, sep="");
        }
    }else{
        if(netbuild.opt == "pcsf"){
         statement <- paste("SELECT * FROM ",table.nm, sep="");
        }else{
        statement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))", sep="");
        }
    }

    ppi.res <- .query.sqlite(ppi.db, statement);

    # remove dupliated edges
    # ppi.res <- ppi.res[!duplicated(ppi.res$row_id),]
    if(!data.org %in% c("mmu", "hsa") && net.type =="string" ){
        id1 <- doPpiIDMapping(as.vector(ppi.res$id1), "accession");
        id2 <- doPpiIDMapping(as.vector(ppi.res$id2), "accession");
        ppi.res$id1 <- id1;
        ppi.res$id2 <- id2;
        ppi.res <- ppi.res[!is.na(ppi.res),]
    }
    return(ppi.res);  
}

# table name is org code, id.type is column name
QueryMirSQLite <- function(table.nm, id.type, q.vec, inv, db.nm){
    require('RSQLite');
    path <- paste0(sqlite.path, "mir2gene.sqlite")
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
    mir.db <- dbConnect(SQLite(), path);
    table.nm = "human"
    query <- paste (shQuote(q.vec),collapse=",");

    statement <- paste("SELECT * FROM ", table.nm, " WHERE upid IN (",query,")", sep="");

    mir.dic <- .query.sqlite(mir.db, statement);
    mir.dic <- mir.dic[complete.cases(mir.dic),]; # get all records
    return(mir.dic);
}

QueryMicSQLite <- function(q.vec, table.nm, sql.nm, min.score, exclude.opt="currency"){
    print(exclude.opt);
    require('RSQLite');
    path <- paste0(sqlite.path, sql.nm)
    mir.db <- dbConnect(SQLite(), path);

    query <- paste (shQuote(q.vec),collapse=",");
    
    if(exclude.opt == "universal"){
    statement <- paste("SELECT * FROM ", table.nm, " WHERE ", table.nm, " IN (",query,")"," AND potential >=", min.score, " AND ExcluUni == 0", sep="");
    }else{
    statement <- paste("SELECT * FROM ", table.nm, " WHERE ", table.nm, " IN (",query,")"," AND potential >=", min.score, sep="");
    }

    mir.dic <- .query.sqlite(mir.db, statement);
    mir.dic <- mir.dic[complete.cases(mir.dic),]; # get all records

    met.ids <- mir.dic$metID
    met.path <- paste(lib.path, "microbiome", "/metInfo.rds", sep="");
    met.info <- readRDS(met.path);
    conv <- data.frame(metID=met.info$metID, KEGG=met.info$KEGG)
    mir.dic <- merge(mir.dic, conv, by="metID");

    mir.dic$KEGG[is.na(mir.dic$KEGG)] <- mir.dic$metID[is.na(mir.dic$KEGG)]
    return(mir.dic);
}

# table name is org code, id.type is column name
QueryMetSQLite <- function(table.nm, q.vec, type){
    require('RSQLite');
    path <- paste0(sqlite.path, "met.sqlite")
    lnc.db <- dbConnect(SQLite(), path);
    query <- paste (shQuote(q.vec),collapse=",");
    if(netbuild.opt == "pcsf"){
        statement <- paste("SELECT * FROM ", table.nm, sep="");
    }else if(net.inv != "inverse"){
        statement <- paste("SELECT * FROM ", table.nm, " WHERE entrez IN (",query,")", sep="");
    }else if(type == "meta"){  
        statement <- paste("SELECT * FROM ", table.nm, " WHERE met IN (",query,")", sep="");
    }else if(type == "kegg"){
        statement <- paste("SELECT * FROM ", table.nm, " WHERE kegg IN (",query,")", sep="");
    }else if(type == "name"){
        statement <- paste("SELECT * FROM ", table.nm, " WHERE name IN (",query,")", sep="");
    }else if(type == "bigg"){
        statement <- paste("SELECT * FROM ", table.nm, " WHERE bigg IN (",query,")", sep="");
    }else{
        statement <- paste("SELECT * FROM ", table.nm, " WHERE chebi IN (",query,")", sep="");
    }
    my.dic <- .query.sqlite(lnc.db, statement);
    return(my.dic);
}

QueryMetSQLiteNet <- function(table.nm, q.vec, inv){
    require('RSQLite');
    path <- paste0(sqlite.path, "met.sqlite")
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
    path <- paste0(sqlite.path, "met.sqlite")
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

QueryKoSQLiteNet <- function(table.nm, q.vec, inv){
    require('RSQLite');
    path <- paste0(sqlite.path, "ko.sqlite")
    lnc.db <- dbConnect(SQLite(), path);
    table.nm = "ko"
    query <- paste (shQuote(q.vec),collapse=",");
     
    if(netbuild.opt == "pcsf"){
    statement <- paste("SELECT * FROM ", table.nm, sep="");
    }else if(inv == "inverse"){
        statement <- paste("SELECT * FROM ", table.nm, " WHERE kegg IN (",query,")", sep="");
    }else{
        statement <- paste("SELECT * FROM ", table.nm, " WHERE ko IN (",query,")", sep="");
    } 

    my.dic <- .query.sqlite(lnc.db, statement);
    return(my.dic);
}

QueryTFSQLite <- function(table.nm, q.vec, inv){
    require('RSQLite');
    path <- paste0(sqlite.path, "tf2gene.sqlite")
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


Query.mGWASDB <- function(db.path, q.vec, table.nm, col.nm){
  require('RSQLite');
  db.path <<- db.path;
  
  db.path <- paste0(db.path, ".sqlite");
  mir.db <- dbConnect(SQLite(), db.path);

  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", table.nm, " WHERE ", col.nm," IN (", query, ")", sep="");
  mir.dic <- .query.sqlite(mir.db, statement);

  # remove duplicates

  dup.inx <- duplicated(mir.dic$mgwas);
  mir.dic <- mir.dic[!dup.inx, ];
  dataSet$snpTable <- mir.dic;
  return(mir.dic);
}