  
my.reg.enrich <- function(file.nm, fun.type, ora.vec, netInv, idType, sourceView="2d"){
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
  }  else if(idType == "entrez") {
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

  type = "regNetwork";
  # write csv
  csv.nm <- paste(file.nm, ".csv", sep="");
  fast.write.csv(resTable1, file=csv.nm, row.names=F);
  fast.write.csv(resTable1, file=paste0(type, "_enr_table.csv"), row.names=F);

  #record table for report
  dataSet$imgSet$enrTables[[type]] <- list()
  dataSet$imgSet$enrTables[[type]]$table <- resTable1;
  dataSet$imgSet$enrTables[[type]]$res.mat <- resTable1[,3, drop=F];
  dataSet$imgSet$enrTables[[type]]$sourceView <- sourceView;

  dataSet$imgSet$enrTables[[type]]$library <- fun.type
  dataSet$imgSet$enrTables[[type]]$algo <- "Overrepresentation Analysis"
  dataSet <<- dataSet

  return(1);
}