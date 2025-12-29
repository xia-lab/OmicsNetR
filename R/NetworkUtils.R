##################################################
## R scripts for OmicsNet
## Description: network analysis methods
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################


#' Create network from only input (seeds)
#'
#' @param dataSetObj Input the name of the created dataSetObj (see Init.Data)
#'
#' @export
#'
BuildSeedProteinNet <- function(dataSetObj=NA){
  dataSet <- .get.nSet(dataSetObj);
  nodes <- V(overall.graph)$name;
  hit.inx <- nodes %in% names(expr.vec);
  nodes2rm <- nodes[!hit.inx];

  g <- simplify(delete.vertices(overall.graph, nodes2rm));

  nodeList <- get.data.frame(g, "vertices");
  nodeList <- nodeList[,1:2];
  colnames(nodeList) <- c("Id", "Label");
  fast.write.csv(nodeList, file="orig_node_list.csv");
  nd.inx <- omics.net$node.data[,1] %in% nodeList[,1];

  edgeList <- get.data.frame(g, "edges");
  edgeList <- edgeList[,1:2];
  colnames(edgeList) <- c("Source", "Target");
  fast.write.csv(edgeList, file="orig_edge_list.csv");
  dataSet <- .decomposeGraph(dataSet, g);
  if(!is.null(dataSet$substats)){
    # update omics.net
    if(.on.public.web){
    .set.nSet(dataSet);
    return(c(length(seed.genes),length(dataSet$seeds.proteins), nrow(nodeList), nrow(edgeList), length(ppi.comps), dataSet$substats));
    }else{
    return(.set.nSet(dataSet));
    }
  }else{
    return(c(0));
  }
}


#' Create igraph object from edgelist created from database selection and decompose into connected subnetworks
#'
#' @export
#'
CreateGraph <- function(dataSetObj=NA){
  dataSet <- .get.nSet(dataSetObj);
  require('igraph');
  node.list <- omics.net$node.data;
  edge.list <- omics.net$edge.data;
  seed.proteins <- omics.net$node.data[,1];
  overall.graph <- simplify(graph.data.frame(edge.list, directed=FALSE)) #, vertices=node.list));

  # add node expression value
  newIDs <- dataSet$seeds.proteins;

  match.index <- match(V(overall.graph)$name, newIDs);
  expr.vals <- dataSet$seeds.expr[match.index];
  names(expr.vals) <- rownames(dataSet$seeds.expr)[match.index];
  expr.vals <- unlist(expr.vals);
  expr.vec <- vector();
  if(length(dataSet$seed) > 0){
    for(i in 1:length(dataSet$seed)){
      res <- unlist(dataSet$seed[[i]][,1])
      names(res) <- rownames(dataSet$seed[[i]])
      expr.vec <- c(expr.vec, res);
    }
  }

  match.index <- names(expr.vec) %in% V(overall.graph)$name;
  expr.vec <- expr.vec[match.index];
  expr.vec <<- expr.vec;
  current.overall.graph <- overall.graph;

  overall.graph <- suppressWarnings(set.vertex.attribute(overall.graph, "abundance", index = V(overall.graph), value = expr.vals));
  overall.graph <<- overall.graph;
  dataSet <- .decomposeGraph(dataSet, overall.graph);
  seed.proteins <<- seed.proteins
  if(!.on.public.web){
    return(.set.nSet(dataSet));
  }else if(!is.null(dataSet$substats)){
    .set.nSet(dataSet);
    message("Network Graph has been created!")
    return(c(length(seed.genes), length(seed.proteins), net.stats$Node[1], net.stats$Edge[1], length(ppi.comps), dataSet$substats));
  }else{
    message("Network Graph created!")
    overall.graph <<- current.overall.graph;
    return(0);
  }
}

#' Prepare the json file for network visualization
#'
#' @param dataSetObj Input the name of the created dataSetObj (see Init.Data)
#' @param net.nm Name of subnetwork created from network building (i.e subnetwork1)
#' @param json.nm Name of json file to be exported for network visualization
#'
#' @export
#'
PrepareNetwork <- function(dataSetObj=NA, net.nm, json.nm){
  require("igraph");
  dataSet <- .get.nSet(dataSetObj);

  if(is.null(net.nm)){
    net.nm <- names(ppi.comps)[1];
  }

  my.ppi <- ppi.comps[[net.nm]];
  nd.nms <- V(my.ppi)$name;
  current.net.nm <<- net.nm;

  if(uploadedGraph && dataSet$fileType != "jsonOmicsnet"){
    convertIgraph2JSONFromFile(net.nm, json.nm, 3);
    message("Conversion from Graph object into Json file completed successfully!")
    return(.set.nSet(dataSet));
  }else{
    GeneAnotDB <- doProteinIDMapping(nd.nms, "entrez");
    entrezIDs <- GeneAnotDB[,1];
    names(entrezIDs) <- nd.nms;
    current.anot <<- entrezIDs;
    convertIgraph2JSON(dataSet, net.nm, json.nm, FALSE);
    message("Conversion from Graph object into Json file completed successfully!")
    return(.set.nSet(dataSet));
  }
}

PrepareMinNetwork <- function(dataSetObj=NA, net.nm, json.nm){
  require("igraph");
  dataSet <- .get.nSet(dataSetObj);
  my.ppi <-SteinerTree_cons(seed.genes_entrez, ppi.comps[[net.nm]],2);
  ppi.comps$minimumNet <<- my.ppi;
  nd.nms <- V(my.ppi)$name;
  GeneAnotDB <- doProteinIDMapping(nd.nms, "entrez");
  entrezIDs <- GeneAnotDB[,1];
  names(entrezIDs) <- nd.nms;
  current.anot <<- entrezIDs;
  current.net.nm <<- "minimumNet";
  convertIgraph2JSON(dataSet, "minimumNet", json.nm, FALSE);
  return(.set.nSet(dataSet));
}

# from node ID (uniprot) return entrez IDs (one to many)
GetNodeEntrezIDs <- function(uniprotID){
  enIDs <- current.anot[uniprotID];
  enIDs <- paste(enIDs, collapse="||");
  enIDs;
}

GetNetsName <- function(){
  rownames(net.stats);
}

GetNetsNameString <- function(){
  paste(rownames(net.stats), collapse="||");
}

GetNetsTypes <- function(){
  paste(types_arr, collapse="??");
}

GetNetsEdgeNum <- function(){
  as.numeric(net.stats$Edge);
}

GetNetsNodeNum <- function(){
  as.numeric(net.stats$Node);
}

GetNetsQueryNum <- function(){
  as.numeric(net.stats$Query);
}

GetIndInputNames <- function(){
  as.vector(dataSet$ind.net.stats$Input);
}

GetIndInteractionValue <- function(){
  as.vector(dataSet$ind.net.stats$NetworkValue);
}

GetIndInteractionNames <- function(){
  as.vector(dataSet$ind.net.stats$Network);
}

GetIndNetsEdgeNum <- function(){
  as.integer(dataSet$ind.net.stats$Edge);
}

GetIndNetsNodeNum <- function(){
  as.integer(dataSet$ind.net.stats$Node);
}

GetIndNetsQueryNum <- function(){
  as.integer(dataSet$ind.net.stats$Query);
}

#' Compute the shortest path between two nodes
#'
#' @param from ID of from node
#' @param to ID of target node
#' @param intermediate if true check whether from and target node connects to each other
#'
GetShortestPaths <- function(from, to, intermediate="false"){
  current.net <- ppi.comps[[current.net.nm]];

  paths <- igraph::get.all.shortest.paths(current.net, from, to)$res;
  if(length(paths) == 0){
    return (paste("No connection between the two nodes!"));
  }

  path.vec <- vector(mode="character", length=length(paths));
  for(i in 1:length(paths)){
    path.inx <- paths[[i]];
    path.ids <- V(current.net)$name[path.inx];
    path.sybls <- lblsu[path.inx];
    pids <- paste(path.ids, collapse="->");
    psbls <- paste(path.sybls, collapse="->");
    path.vec[i] <- paste(c(pids, psbls), collapse=";")
  }

  if(intermediate == "true"){
    edge.list <- as_edgelist(current.net);
    edges.from <- edge.list[which(edge.list[,1] == from), ];
    edges.from <- rbind(edges.from, edge.list[which(edge.list[,2] == from), ]);
    edges.to <- edge.list[which(edge.list[,1] == to), ];
    edges.to  <- rbind(edges.to, edge.list[which(edge.list[,2] == to), ]);
    edges.from.universe <- unique(as.vector(edges.from));
    edges.to.universe <- unique(as.vector(edges.to));
    intersects <- intersect(edges.from.universe, edges.to.universe);
    if(length(intersects) != 0){
      for(i in 1:length(intersects)){
        element <- paste(from, intersects[i], to, sep="->")
        if(!element %in% path.vec){
          path.vec[length(path.vec) + i] <- paste(from, intersects[i], to, sep="->");
        }
      }
    }
  }

  if(length(path.vec) > 50){
    path.vec <- path.vec[1:50];
  }

  all.paths <- paste(path.vec, collapse="||");
  return(all.paths);
}

###################################
# Adapted from netweavers package
###################
weightCalc <- function(filt_net){
  all_gene_int <- merge(table(filt_net[,2]),table(filt_net[,1]),all=TRUE);
  neighbors <- aggregate(all_gene_int[,2],by=list(all_gene_int[,1]),sum);
  names(neighbors) <- c("Protein","n");
  data.frame(Protein=neighbors$Protein,weight=1/neighbors$n);
}

measureCalc <- function(protein_pvals,protein_wts,weightamt = 10){
  prot_data <- merge(protein_pvals,protein_wts,by='Protein',all=TRUE);
  prot_data$measure <- -log((prot_data$pvalue)*(prot_data$weight)^(1/weightamt));
  net_no_pval <- is.na(prot_data$measure) & !is.na(prot_data$weight);
  prot_data[net_no_pval,'measure'] <- 0;
  prot_data;
}

scoreClusters <- function(clust_2_protein,prot_meas){
  clust_score_total <- apply(clust_2_protein[,-c(1,2)],1,function(x){do.call(m,list(prot_meas[prot_meas$Protein%in%x[!is.na(x)],'measure']))});
  names(clust_score_total) <- clust_2_protein[,1];
  clust_score_total;
}

#' Compute minimum connected subnetwork from input nodes using shortest path based approach
#'
#' @param dataSetObj Input the name of the created dataSetObj (see Init.Data)
#' @param nodeids The ids of input nodes
#' @param dim dimension of network results to be exported
#'
#' @export
#'
ExtractModule<- function(dataSetObj=NA, nodeids, dim="3"){
  library(igraph);
  dim = as.numeric(dim);
  set.seed(8574);
  nodes <- strsplit(nodeids, ";")[[1]];

  g <- ppi.comps[[current.net.nm]];
  # try to see if the nodes themselves are already connected
  hit.inx <- V(g)$name %in% nodes;
  gObj <- induced.subgraph(g, V(g)$name[hit.inx]);

  # now find connected components
  comps <-decompose.graph(gObj, min.vertices=1);

  if(length(comps) == 1){ # nodes are all connected
    g <- comps[[1]];
  }else{
    # extract modules
    paths.list <-list();
    sd.len <- length(nodes);
    for(pos in 1:sd.len){
      paths.list[[pos]] <- igraph::get.shortest.paths(g, nodes[pos], nodes[-(1:pos)])$vpath;
    }
    nds.inxs <- unique(unlist(paths.list));
    nodes2rm <- V(g)$name[-nds.inxs];
    g <- simplify(delete.vertices(g, nodes2rm));
  }
  nodeList <- get.data.frame(g, "vertices");
  if(nrow(nodeList) < 3){
    return ("NA");
  }

  module.count <- module.count + 1;
  module.nm <- paste("module", module.count, sep="");
  colnames(nodeList) <- c("Id", "Label");
  ndFileNm = paste(module.nm, "_node_list.csv", sep="");
  fast.write.csv(nodeList, file=ndFileNm, row.names=F);

  edgeList <- get.data.frame(g, "edges");
  edgeList <- cbind(rownames(edgeList), edgeList);
  colnames(edgeList) <- c("Id", "Source", "Target");
  edgFileNm = paste(module.nm, "_edge_list.csv", sep="");
  fast.write.csv(edgeList, file=edgFileNm, row.names=F);

  filenm <- paste(module.nm, ".json", sep="");

  # record the module
  ppi.comps[[module.nm]] <<- g;
  UpdateSubnetStats();

  module.count <<- module.count
  if(uploadedGraph){
    convertIgraph2JSONFromFile(module.nm, filenm, dim);
  }else{
    convertIgraph2JSON(dataSet, module.nm, filenm, FALSE, dim);
  }

  if(.on.public.web){
    .set.nSet(dataSet)
    return (filenm);
  }else{
    return (.set.nSet(dataSet))
  }
}

SearchNetDB <- function(dataSetObj, protein.vec, orig.input, inputType, netw.type, db.type, require.exp=TRUE,
                        min.score = 900, netInv, zero=FALSE, snpRegion=FALSE){
  
  dataSet <- .get.nSet(dataSetObj);
  if(inputType == "peak"){
    netInv = "direct";
  }

  edge.res <-data.frame();

  require("RJSONIO");

  # now do the database search
  if(netw.type == "ppi" || netw.type == "gene" || netw.type == "protein"){
    table.nm = paste(data.org, db.type, sep="_");
    message("We are going to query the data table: ", table.nm, " from ppi database..");
    # Check if OmniPath database is selected
    if(grepl("omnipath$", table.nm)){
      res <- QueryOmnipathPpiSQLite(sqlite.path, data.org, protein.vec);
    }else{
      res <- QueryPpiSQLite(table.nm, protein.vec, require.exp, min.score);
    }

    if(dataSet$ppiZero){
      if(inputType == "gene"){
        universe.vec <- c(protein.vec,rownames(dataSet$exp.mat$protein))
      }else if(inputType == "protein"){
        universe.vec <- c(protein.vec,rownames(dataSet$exp.mat$gene))
      }else{
        universe.vec <- protein.vec;
      }

      hit.inx1 <- res[,1] %in% universe.vec;
      hit.inx2 <- res[,2] %in% universe.vec;
      res <- res[(hit.inx1 & hit.inx2),];
      n.ids <- c(res[,1], res[,2])
    }

    # no hits
    if(nrow(res)==0){ 
        return(c(0,0)); 
    }
    
    edge.res <- data.frame(Source=res$id1,Target=res$id2, stringsAsFactors=FALSE);
    row.names(edge.res) <- 1:nrow(res);
    fast.write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);

    node.ids <- c(res[,1], res[,2])
    node.nms <- c(res[,3], res[,4])

    if(netw.type == "gene" || netw.type == "ppi" ){
      net.info$protein.ids <- unique(c(net.info$protein.ids, unique(node.ids)));
    }else if(netw.type == "protein"){
      net.info$protein.ids <- c(net.info$protein.ids, unique(node.ids));
    }else{
      net.info$snp.ids <- c(net.info$snp.ids, unique(node.ids));
    }

  } else if(netw.type == "snp"){
  
    if(!exists("my.snp.query")){ # public web on same user dir
        compiler::loadcmp("../../rscripts/OmicsNetR/R/utils_snp.Rc");
    }
    res <- my.snp.query(dataSet, protein.vec, db.type, netInv, zero, snpRegion)
    edge.res <- res$edge.res;
    net.info <- res$net.info;
    node.ids <- res$node.ids;
    node.nms <- res$node.nms;

  } else if (netw.type == "tf") {

    table.nm <- paste(data.org, db.type, sep="_");
    # Check if OmniPath database is selected
    if(db.type == "omnipath"){
      res <- QueryOmnipathTfSQLite(sqlite.path, data.org, protein.vec);
    }else{
      res <- QueryTFSQLite(table.nm, protein.vec, netInv);
    }
    # no hits
    if(nrow(res)==0){ return(c(0,0)); }

    edge.res <- data.frame(Source=res[,"tfid"],Target=res[,"entrez"], stringsAsFactors=FALSE);
    if(nrow(res)!=0){
      row.names(edge.res) <- 1:nrow(res);
    }
    fast.write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    node.ids <- c(res[,"entrez"], res[,"tfid"])
    node.nms <- c(res[,"symbol"], res[,"tfname"]);

    if(netInv == "direct"){
      net.info$gene.ids <- unique(c(net.info$gene.ids, node.ids[!node.ids %in% protein.vec]));
      net.info$tf.ids <-  unique(res[,"tfid"]);
    }else{
      net.info$tf.ids <- unique(res[,"tfid"]);
      if(!zero){
        net.info$gene.ids <- c(unique(protein.vec));
      }
    }

  } else if(netw.type == "mir") { # in miRNA, table name is org code, colname is id type

    if(db.type =="targetscan"){
      table.nm <- data.org;
    }else{
      table.nm <- data.org;
    }
    # Check if OmniPath database is selected
    if(db.type == "omnipath"){
      res <- QueryOmnipathMirSQLite(sqlite.path, data.org, protein.vec);
    }else{
      res <- QueryMirSQLite(table.nm, "mir_id", protein.vec, netInv, db.type);
    }
    if(nrow(res)==0){ return(c(0,0)); }
    # no hits

    edge.res <- data.frame(Source=res[,"mir_id"],Target=res[,"entrez"],stringsAsFactors = FALSE)
    if(nrow(res)!=0){
      row.names(edge.res) <- 1:nrow(res);
    }
    fast.write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);

    node.ids <- c(res[,"mir_id"], res[,"entrez"]);
    node.nms <- c(res[,"mir_acc"], res[,"symbol"]);

    if(netInv == "direct"){
      net.info$gene.ids <- unique(c(net.info$gene.ids, node.ids[!node.ids %in% protein.vec]));
      net.info$mir.ids <- unique(res[,"mir_id"]);
    }else{
      net.info$mir.ids <- unique(res[,"mir_id"]);
      if(!zero){
        net.info$gene.ids <- c(net.info$gene.ids, unique(res[,"entrez"]));
      }
    }

  } else if(netw.type == "met") {
    ## metabolite-protein
    if(db.type == "keggp"){
      table.nm <- paste("keggp", sep="_");
    }else{
      table.nm <- paste(data.org, db.type, sep="_");
    }

    res <- QueryMetSQLiteNet(table.nm, protein.vec, netInv);
    # no hits
    if(nrow(res)==0){ return(c(0,0)); }

    if(db.type =="keggp"){ # project to kegg
      #res1 = res[which(res$entrez %in% c(dataSet$seeds.proteins)),];
      #res2 = res[which(res$kegg %in% c(dataSet$seeds.proteins)),];
      #res = rbind(res1, res2);
    }

    edge.res <- data.frame(Source=res[,"kegg"],
                           Target=res[,"entrez"],
                           stringsAsFactors=FALSE);

    if(nrow(res)!=0){
      row.names(edge.res) <- 1:nrow(res);
    }
    fast.write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    node.ids <- c(res[,"kegg"], res[,"entrez"]);
    node.nms <- c(res[,"met"], res[,"symbol"]);

    prot.ids <- unique( res[,"entrez"]);
    rm.inx <- prot.ids %in% net.info$met.ids;
    prot.ids <- prot.ids[!rm.inx];

    if(netInv == "direct"){
      if(!is.null(net.info$peak.ids)){
        not.gene.ids <- c(net.info$peak.ids,
                          net.info$met.ids, net.info$kncmpd.ids)
        net.info$met.ids <- c(net.info$met.ids, unique(res[,"kegg"][!res[,"kegg"] %in% not.gene.ids]));
        net.info$protein.ids <- c(net.info$protein.ids, unique( prot.ids));
      }else{
        net.info$met.ids <- c(net.info$met.ids, unique(res[,"kegg"]));
        net.info$protein.ids <- c(net.info$protein.ids, unique( prot.ids));
      }
    }else{
      net.info$met.ids <- c(net.info$met.ids, unique(res[,"kegg"]));
      if(!zero){
        net.info$protein.ids <- c(net.info$protein.ids, unique(prot.ids));
      }
    }

  } else if(netw.type == "peak") {
    PeakSet <- qs:::qread("PeakSet_net.qs");
    net.info$met.ids <- PeakSet$mets;
    net.info$peak.ids <- PeakSet$put.mets;
    seed.genes <<- unique(c(seed.genes, PeakSet$mets));    
    node.ids <- PeakSet$nodes.df[,1];
    node.nms <- PeakSet$nodes.df[,2];
    edge.res <- PeakSet$edges.df;

  } else if(netw.type == "mic") { # in mic
    # in mic
    sql.name <-paste0("omicsnet_", db.type, ".sqlite");
    table.nm <- mic.taxa;

    res <- QueryMicSQLite(protein.vec, table.nm, sql.name, dataSet$mic.thresh, dataSet$currExclude,dataSet$uniExclude,dataSet$orphExclude);
    if(nrow(res)==0){ return(c(0,0)); }
    # no hits

    colnames(res)[colnames(res) == table.nm] <- "taxa";
    edge.res.score <- data.frame(Source=res[,"taxa"],Target=res[,"KEGG"],Score=res[,"potential"],stringsAsFactors = FALSE);
    edge.res <- edge.res.score[,-3];
    if(nrow(res)!=0){
      row.names(edge.res) <- 1:nrow(res);
    }
    fast.write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    node.ids <- c(res[,"taxa"], res[,"KEGG"]);
    node.nms <- c(res[,"taxa"], res[,"metabolite"]);

    if(netInv == "direct"){
      net.info$met.ids <- unique(c(net.info$met.ids, res[,"KEGG"]));
      net.info$mic.ids <- unique(protein.vec);
    }else{
      net.info$mic.ids <- node.ids[!node.ids %in% protein.vec];
      if(!zero){
        net.info$met.ids <- c(net.info$met.ids, res[,"KEGG"]);
      }
    }
    require("igraph");
    g <- simplify(graph.data.frame(edge.res, directed=FALSE)) #, vertices=node.list));

    met.ids <- unique(net.info$met.ids);
    met.microbe.list <- list();
    met.microbe.score.list <- list();
    mic.table <- list();
    for( i in 1:length(met.ids)){
      met <- met.ids[i]
      if(length(which(edge.res$Target == met))>1){
        met.microbe.list[[met]] <- as.vector(edge.res[which(edge.res$Target == met), 1])
        met.microbe.score.list[[met]] <- as.vector(edge.res.score[which(edge.res.score$Target == met), 3])
      }
    }
    for(i in 1:nrow(res)){
      mic.table[[i]] <- unlist(res[i, ]);
      exp.mat <- dataSet$exp.mat$mic;
      exp <- exp.mat[rownames(exp.mat) == res[i,"taxa"]];
      if(exp != 0){
        mic.table[[i]]["AbundPred"] <- as.numeric(exp) * as.numeric(mic.table[[i]]["potential"]);
        mic.table[[i]]["potential"] <- as.numeric(mic.table[[i]]["potential"]);
      }
    }

    micSet <- list();
    micSet[["met.mic"]] <- met.microbe.list;
    micSet[["met.mic.score"]] <- met.microbe.score.list;
    micSet[["met.mic.table"]] <- mic.table;
    qs::qsave(micSet, "micSet.qs");

    dataSet <<- dataSet;
  } else if (netw.type == "m2m" && inputType == "peak") {
    if (db.type == "kegg"){
      table.nm <- paste(data.org, "KEGG_m2m", sep="_");
    } else if (db.type == "keggp") {
      table.nm <- "ko_KEGG_m2m";
    } else if (db.type == "agora") {
      table.nm <- "agora_kegg_m2m";
    } else if (db.type == "embl") {
      table.nm <- "embl_kegg_m2m";
    } else {
      ## for recon 3
      table.nm <- "hsa_recon3_m2m";
    }

    if (data.org == "microbiome" | data.org == "NA") {
      #TODO: to consider this option later
      table.nm <- db.type;
      data.org <<- "microbiome";
    }

    res <- extendMetPeakNetwork(table.nm);
    if(is.list(res)) {
      edge.res <- res$edge.res;
      #net.info <- res$net.info;
    } else {
      edge.res <- NULL;
    }

  } else if (netw.type == "m2m" && inputType != "peak"){

    if (db.type == "kegg"){
      table.nm <- paste(data.org, "KEGG_m2m", sep="_");
    } else if (db.type == "keggp") {
      table.nm <- "ko_KEGG_m2m";
    } else if (db.type == "agora") {
      table.nm <- "agora_kegg_m2m";
    } else if (db.type == "embl") {
      table.nm <- "embl_kegg_m2m";
    } else {
      ## for recon 3
      table.nm <- "hsa_recon3_m2m";
    }

    res <- QueryMicM2mSQLiteNet(table.nm, protein.vec);

    # no hits
    if(nrow(res)==0){ return(c(0,0)); }

    edge.res <- data.frame(Source=res[,"sourceID"],
                           Target=res[,"productID"],
                           stringsAsFactors=FALSE);

    if(nrow(res)!=0){
      row.names(edge.res) <- 1:nrow(res);
    }
    fast.write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    node.ids <- c(res[,"sourceID"], res[,"productID"])
    node.nms <- c(res[,"sourceNM"], res[,"productNM"]);

    if(netInv == "direct"){
      if("peak" %in% dataSet$type){
        #not.kncmpd.ids <- c(net.info$peak.ids, net.info$met.ids);
        #net.info$kncmpd.ids <- c(net.info$kncmpd.ids, node.ids[!node.ids %in% not.kncmpd.ids]);
      }else{
        #net.info$met.ids <- unique(protein.vec);
        #net.info$kncmpd.ids <- c(net.info$kncmpd.ids, node.ids[!node.ids %in% protein.vec]);
      }
    }else{
      #net.info$met.ids <- node.ids[!node.ids %in% protein.vec];
      if(!zero){
        #net.info$kncmpd.ids <- c(net.info$kncmpd.ids, unique(protein.vec));
      }
    }
  }

  if(length(edge.res)>0){
    netw_input_type <- paste0(netw.type,"_", inputType,"_", data.org, "_", orig.input);
    edgeu.res.list[[netw_input_type]]$table <- unique(edge.res[,c(1,2)]);
    edgeu.res.list[[netw_input_type]]$database <- makeReadable(db.type);
    edgeu.res.list <<- edgeu.res.list;
    nodeu.ids <<- c(nodeu.ids, node.ids);
    nodeu.nms <<- c(nodeu.nms, node.nms);
  }
  net.info <<- net.info;
  return(1);
}

SteinerTree_cons <- function(terminal_nodes, PPI_graph, run_times) {

  color = NULL
  terminal_nodes = na.omit(terminal_nodes)
  V(PPI_graph)$color = "yellow"
  hits.inx = V(PPI_graph)$name %in% terminal_nodes
  V(PPI_graph)[hits.inx]$color = "red"
  terminals = V(PPI_graph)[hits.inx]
  steinertmin = vector()
  steinertrees = list()
  for(runs in 1: run_times)
  {edges = c()
  prob = sample(1:length(terminals), 1, replace = FALSE)
  subtree = terminals$name[[prob]]
  nsubtree = setdiff(terminals$name, subtree)
  tparam = 1
  while(tparam <= length(terminals))
  {
    paths = igraph::get.all.shortest.paths(PPI_graph,subtree, nsubtree)
    if(length(paths$res)>1)
    {
      paths_length = sapply(paths$res, length)
      sp = paths$res[which(paths_length == min(paths_length))][[1]]
      subtree = igraph::union(subtree, V(PPI_graph)$name[sp])
      nsubtree = setdiff(nsubtree, V(PPI_graph)$name[sp])
    }else{
      subtree = subtree
      nsubtree = nsubtree
    }
    tparam = tparam+1
  }
  steinert =mst(induced.subgraph(PPI_graph, subtree))
  for(i in length(which(V(steinert)$color == "yellow"))+1)
  { degr = degree(steinert, v = V(steinert), mode = c("all"))
  todel = names(which(degr == 1))
  todel = todel[which(!todel %in% terminals$name)]
  if(length(todel) > 0)
  {steinert = delete.vertices(steinert, todel)}
  }
  steinertrees[[runs]] = steinert
  steinertmin[runs]  = length(V(steinert)$name)
  }
  return(steinertrees[[which(steinertmin == min(steinertmin))[1]]])
}

FilterBipartiNet <- function(dataSetObj=NA, nd.type, min.dgr, min.btw){
  library(igraph);
  dataSet <- .get.nSet(dataSetObj);
  all.nms <- V(overall.graph)$name;
  edge.mat <- get.edgelist(overall.graph);
  dgrs <- degree(overall.graph);
  nodes2rm.dgr <- nodes2rm.btw <- NULL;

  if(nd.type == "gene"){
    hit.inx <- all.nms %in% c(net.info$gene.ids, net.info$protein.ids);
  }else if(nd.type=="tf"){
    hit.inx <- all.nms %in% net.info$tf.ids;
  }else if(nd.type=="mir"){
    hit.inx <- all.nms %in% net.info$mir.ids;
  }else if(nd.type=="met"){
    hit.inx <- all.nms %in% net.info$met.ids;
  }else{ # all
    hit.inx <- rep(TRUE, length(all.nms));
  }

  if(min.dgr > 0){
    rm.inx <- dgrs <= min.dgr & hit.inx;
    nodes2rm.dgr <- V(overall.graph)$name[rm.inx];
  }
  if(min.btw > 0){
    btws <- betweenness(overall.graph);
    rm.inx <- btws <= min.btw & hit.inx;
    nodes2rm.btw <- V(overall.graph)$name[rm.inx];
  }

  nodes2rm <- unique(c(nodes2rm.dgr, nodes2rm.btw));
  overall.graph <- simplify(delete.vertices(overall.graph, nodes2rm));
  current.msg <<- paste("A total of", length(nodes2rm) , "was reduced.");
  dataSet <- .decomposeGraph(dataSet, overall.graph);
  if(!is.null(dataSet$substats)){
    overall.graph <<- overall.graph;
    if(.on.public.web){
    .set.nSet(dataSet);
    return(c(length(seed.genes),length(seed.proteins), vcount(overall.graph), ecount(overall.graph), length(ppi.comps), dataSet$substats));
    }else{
    return (.set.nSet(dataSet))
    }
  }else{
    return (.set.nSet(dataSet))
  }
}


#' Compute minimum connected network composed of seed nodes using shortest path based approach
#'
#' @param dataSetObj Input the name of the created dataSetObj (see Init.Data)
#' @param max.len Maximum number of seeds, if more than this number, take top 100 seed nodes based on degrees
#'
#' @return check overall.graph object for result
#' @export
#'
BuildMinConnectedGraphs <- function(dataSetObj=NA, max.len = 200){
  library(igraph);

  dataSet <- .get.nSet(dataSetObj);
  set.seed(8574);
  # first get shortest paths for all pair-wise seeds
  hit.inx <- unique(names(expr.vec)) %in% V(overall.graph)$name;
  my.seeds = names(expr.vec)[hit.inx];
  seed.proteins = my.seeds;
  sd.len <- length(my.seeds);
  paths.list <-list();

  # first trim overall.graph to remove no-seed nodes of degree 1
  dgrs <- degree(overall.graph);
  keep.inx <- dgrs > 1 | (names(dgrs) %in% my.seeds);
  nodes2rm <- V(overall.graph)$name[!keep.inx];
  overall.graph <-  simplify(delete.vertices(overall.graph, nodes2rm));

  # need to restrict the operation b/c get.shortest.paths is very time consuming
  # for top max.len highest degrees
  if(sd.len > max.len){
    hit.inx <- names(dgrs) %in% my.seeds;
    sd.dgrs <- dgrs[hit.inx];
    sd.dgrs <- rev(sort(sd.dgrs));
    # need to synchronize all (seed.proteins) and top seeds (my.seeds)
    seed.proteins <- names(sd.dgrs);
    my.seeds <- seed.proteins[1:max.len];
    sd.len <- max.len;
    current.msg <<- paste("The minimum connected network was computed using the top", sd.len, "seed proteins in the network based on their degrees.");
  }else{
    current.msg <<- paste("The minimum connected network was computed using all seed proteins in the network.");
  }
  # now calculate the shortest paths for
  # each seed vs. all other seeds (note, to remove pairs already calculated previously)
  for(pos in 1:sd.len){
    paths.list[[pos]] <- igraph::get.shortest.paths(overall.graph, my.seeds[pos], seed.proteins[-(1:pos)])$vpath;
  }
  nds.inxs <- unique(unlist(paths.list));
  nodes2rm <- V(overall.graph)$name[-nds.inxs];
  g <- simplify(delete.vertices(overall.graph, nodes2rm));

  nodeList <- get.data.frame(g, "vertices");
  colnames(nodeList) <- c("Id", "Label");
  fast.write.csv(nodeList, file="orig_node_list.csv", row.names=F);

  edgeList <- get.data.frame(g, "edges");
  edgeList <- cbind(rownames(edgeList), edgeList);
  colnames(edgeList) <- c("Id", "Source", "Target");
  fast.write.csv(edgeList, file="orig_edge_list.csv", row.names=F);

  path.list <- NULL;
  dataSet <- .decomposeGraph(dataSet, g);
  if(!is.null(dataSet$substats)){
    overall.graph <<- g;
    if(.on.public.web){
      .set.nSet(dataSet);
      return(c(length(seed.genes),length(seed.proteins), nrow(nodeList), nrow(edgeList), length(ppi.comps), dataSet$substats));
    }else{
      return(.set.nSet(dataSet));
    }
  }else{
    return(0);
  }
}

#' Compute minimum connect subnetwork based on input nodes using prize-collecting steiner forest approach
#'
#' @param dataSetObj Input the name of the created dataSetObj (see Init.Data)
#'
#' @return check overall.graph object for result
#' @export
#'
BuildPCSFNet <- function(dataSetObj=NA){
  library(igraph);
  dataSet <- .get.nSet(dataSetObj);
  edg <- as.data.frame(get.edgelist(overall.graph));
  edg$V3 <- rep(1, nrow(edg));
  colnames(edg) <- c("from", "to", "cost");

  node_names <- unique(c(as.character(edg[,1]),as.character(edg[,2])))
  ppi <- graph.data.frame(edg[,1:2],vertices=node_names,directed=F)
  E(ppi)$weight <- as.numeric(edg[,3])
  ppi <- simplify(ppi)

  #if(sum(expr.vec) == 0){ # make sure weights are not 0?!
    hit.inx <- unique(names(expr.vec)) %in% V(overall.graph)$name;
    my.seeds = names(expr.vec)[hit.inx];
    #print(my.seeds);
    expr.vec <- setNames(rep(5, length(my.seeds)), my.seeds)
  #}
  #expr.vec <- abs(expr.vec)

  g <- Compute.SteinerForest(ppi, expr.vec, w = 5, b = 100, mu = 0.0005);

  nodeList <- get.data.frame(g, "vertices");
  colnames(nodeList) <- c("Id", "Label");
  fast.write.csv(nodeList, file="orig_node_list.csv", row.names=F);

  edgeList <- get.data.frame(g, "edges");
  edgeList <- cbind(rownames(edgeList), edgeList);
  colnames(edgeList) <- c("Id", "Source", "Target");
  fast.write.csv(edgeList, file="orig_edge_list.csv", row.names=F);

  path.list <- NULL;
  dataSet <- .decomposeGraph(dataSet, g);
  if(!is.null(dataSet$substats)){
    overall.graph <<- g;
    current.msg<<- "Steiner Forest was completed successfully";
    if(.on.public.web){
      .set.nSet(dataSet)
      return(c(length(seed.genes),length(dataSet$seeds.proteins), nrow(nodeList), nrow(edgeList), length(ppi.comps), dataSet$substats));
    }else{
      return(.set.nSet(dataSet));
    }
  }else{
    return(-1);
  }
}

# for a given graph, obtain the smallest subgraphs that contain
# all the seed nodes. This is acheived by iteratively remove
# the marginal nodes (degree = 1) that are not in the seeds

.decomposeGraph <- function(dataSetObj=NA, gObj, minNodeNum = 2){
  dataSet <- .get.nSet(dataSetObj);
  # now decompose to individual connected subnetworks
  if(uploadedGraph == "false"){
    comps <-decompose.graph(gObj, min.vertices=minNodeNum);
  }else{
    if(gsize(gObj)>0){
      comps <-decompose.graph(gObj, min.vertices=minNodeNum);

    }else{
      comps = list()
      comps[[1]] <- overall.graph
    }
  }
  if(length(comps) == 0){
    current.msg <<- paste("No subnetwork was identified with at least", minNodeNum, "nodes!");
    return(NULL);
  }

  # first compute subnet stats
  dataSet <- .computeSubnetStats(dataSet, comps);
  net.stats <- dataSet$net.stats;
  ord.inx <- order(net.stats[,1], decreasing=TRUE);
  net.stats <- net.stats[ord.inx,];
  dataSet$type.nums <- dataSet$type.nums[ord.inx]
  dataSet$query.nums <- dataSet$query.nums[ord.inx]
  dataSet$query.nums[dataSet$query.nums == ""] <- net.stats$Query[dataSet$query.nums == ""]
  dataSet$type.nums[dataSet$type.nums == ""] <- net.stats$Node[dataSet$type.nums == ""]

  comps <- comps[ord.inx];
  names(comps) <- rownames(net.stats) <- paste("subnetwork", 1:length(comps), sep="");

  hit.inx <- net.stats$Node >= minNodeNum;

  comps <- comps[hit.inx];
  sub.stats <- unlist(lapply(comps, vcount));
  # now record
  ppi.comps <<- comps;
  net.stats <<- net.stats;
  dataSet$substats <- sub.stats
  return(dataSet);
}

.computeSubnetStats <- function(dataSetObj=NA, comps){
  dataSet <- .get.nSet(dataSetObj);
  type.nums <- vector();
  query.nums <- vector();
  net.stats <- as.data.frame(matrix(0, ncol = 3, nrow = length(comps)));
  colnames(net.stats) <- c("Node", "Edge", "Query");
  for(i in 1:length(comps)){
    g <- comps[[i]];
    my.stat <- GetNetStatByType(g);
    type.nums <- c(type.nums, my.stat$node.num)
    query.nums <- c(query.nums, my.stat$query.num);
    net.stats[i,] <- c(vcount(g),ecount(g),sum((unique(dataSet$seeds.proteins)) %in% V(g)$name));
  }
  dataSet$type.nums <- type.nums;
  dataSet$query.nums <- query.nums;
  dataSet$net.stats <- net.stats
  return(dataSet);
}

UpdateSubnetStats <- function(){
  old.nms <- names(ppi.comps);
  dataSet <- .computeSubnetStats(dataSet, ppi.comps);
  net.stats <- dataSet$net.stats;
  ord.inx <- order(net.stats[,1], decreasing=TRUE);
  net.stats <- net.stats[ord.inx,];
  dataSet$type.nums <- dataSet$type.nums[ord.inx]
  dataSet$query.nums <- dataSet$query.nums[ord.inx]
  rownames(net.stats) <- old.nms[ord.inx];
  net.stats <<- net.stats;
}


#' Detect graph communities using function from igraph
#'
#' @param method algorithm choice: walktrap, infomap, label propagation
#' @param use.weight Boolean consider edge weight: FALSE
#'
#' @export
FindCommunities <- function(method="infomap", use.weight=FALSE){
  library(igraph);
  set.seed(1);
  # make sure this is the connected
  current.net <- ppi.comps[[current.net.nm]];
  g <- current.net;
  if(!is.connected(g)){
    g <- decompose.graph(current.net, min.vertices=2)[[1]];
  }
  total.size <- length(V(g));

  if(use.weight){ # this is only tested for walktrap, should work for other method
    # now need to compute weights for edges
    egs <- get.edges(g, E(g)); #node inx
    nodes <- V(g)$name;
    # conver to node id
    negs <- cbind(nodes[egs[,1]],nodes[egs[,2]]);

    # get min FC change
    base.wt <- min(abs(seed.expr))/10;

    # check if user only give a gene list without logFC or all same fake value
    if(length(unique(seed.expr)) == 1){
      seed.expr <- rep(1, nrow(negs));
      base.wt <- 0.1; # weight cannot be 0 in walktrap
    }

    wts <- matrix(base.wt, ncol=2, nrow = nrow(negs));
    for(i in 1:ncol(negs)){
      nd.ids <- negs[,i];
      hit.inx <- match(names(seed.expr), nd.ids);
      pos.inx <- hit.inx[!is.na(hit.inx)];
      wts[pos.inx,i]<- seed.expr[!is.na(hit.inx)]+0.1;
    }
    nwt <- apply(abs(wts), 1, function(x){mean(x)^2})
  }

  if(method == "walktrap"){
    fc <- walktrap.community(g);
  }else if(method == "infomap"){
    fc <- infomap.community(g);
  }else if(method == "labelprop"){
    fc <- label.propagation.community(g);
  }else{
    return ("NA||Unknown method!");
  }

  if(length(fc) == 0){
    return ("NA||No communities were detected!");
  }

  # only get communities
  communities <- communities(fc);
  community.vec <- vector(mode="character", length=length(communities));
  com.vec <- list()
  gene.community <- NULL;
  qnum.vec <- NULL;
  pval.vec <- NULL;
  psize.vec <- vector();
  rowcount <- 0;
  nms <- V(g)$name;
  hit.inx <- match(nms, omics.net$node.data[,1]);
  sybls <- omics.net$node.data[hit.inx,2];
  names(sybls) <- V(g)$name;
  for(i in 1:length(communities)){
    # update for igraph 1.0.1
    path.ids <- communities[[i]];
    psize <- length(path.ids);
    if(psize < 10){
      next; # ignore very small community
    }
    hits <- seed.proteins %in% path.ids;
    qnums <- sum(hits);
    if(qnums == 0){
      next; # ignor community containing no queries
    }

    rowcount <- rowcount + 1;
    pids <- paste(path.ids, collapse="->");
    #path.sybls <- V(g)$Label[path.inx];
    path.sybls <- sybls[path.ids];
    com.mat <- cbind(path.ids, path.sybls, rep(i, length(path.ids)));
    gene.community <- rbind(gene.community, com.mat);
    qnum.vec <- c(qnum.vec, qnums);
    psize.vec <- c(psize.vec, psize);
    # calculate p values (comparing in- out- degrees)
    subgraph <- induced.subgraph(g, path.ids);
    in.degrees <- degree(subgraph);
    out.degrees <- degree(g, path.ids) - in.degrees;
    ppval <- wilcox.test(in.degrees, out.degrees)$p.value;
    ppval <- signif(ppval, 3);
    pval.vec <- c(pval.vec, ppval);

    # calculate community score
    community.vec[rowcount] <- paste(c(psize, qnums, ppval, pids), collapse=";");
    com.vec[[rowcount]] <- path.ids
  }

  #ord.inx <- order(pval.vec, decreasing=F);
  #community.vec <- community.vec[ord.inx];
  #qnum.vec <- qnum.vec[ord.inx];
  ord.inx <- order(psize.vec, decreasing=T);
  #print(ord.inx);
  community.vec <- community.vec[ord.inx];

  all.communities <- paste(community.vec, collapse="||");
  df <- convertModuleToDF(all.communities);
  df <- df[,-c(3,5)];


  #record table for report
  type = "module";
  dataSet$imgSet$enrTables[[type]] <- list()
  dataSet$imgSet$enrTables[[type]]$table <- df;

  df1<- data.frame(Size=df$Size,Pval=df$P.value)
  rownames(df1) <- df$Module;
  dataSet$imgSet$enrTables[[type]]$res.mat <- df1;
  dataSet$imgSet$enrTables[[type]]$library <- "";
  dataSet$imgSet$enrTables[[type]]$algo <- method;

  dataSet <<- dataSet;
  colnames(gene.community) <- c("Id", "Label", "Module");
  fast.write.csv(gene.community, file="module_table.csv", row.names=F);
  return(all.communities);
}

convertModuleToDF <- function(dataString) {
  # Split the string into lines based on the '||' separator
  lines <- strsplit(dataString, "\\|\\|")[[1]]
  
  # Initialize an empty list to store each row's data
  rows <- list()
  
  # Iterate over each line and split further by ';' to extract individual fields
    i <- 1;
  for (line in lines) {
    parts <- strsplit(line, ";")[[1]]
    if (length(parts) >= 4) {  # Ensure there are enough parts to form a complete row
      # Extract and store the data for this row
      rows[[length(rows) + 1]] <- list(
        Module = paste0("Module ", i),
        Size = parts[1],
        Name = as.numeric(parts[2]),
        `P-value` = as.numeric(parts[3]),
        IDs = parts[4]
      )
    i = i +1;
    }
  }
  
  # Convert the list of rows into a dataframe
  df <- do.call(rbind, lapply(rows, function(row) as.data.frame(row, stringsAsFactors = FALSE)))
  # Column names are already set via the list names, so this line is actually redundant
  return(df)
}


community.significance.test <- function(graph, vs, ...) {
  subgraph <- induced.subgraph(graph, vs)
  in.degrees <- degree(subgraph)
  out.degrees <- degree(graph, vs) - in.degrees
  wilcox.test(in.degrees, out.degrees, ...)
}

UpdateNetworkLayout3D <- function(algo, filenm){
  # get layers
  current.net <- ppi.comps[[current.net.nm]];

  pos.xyz <- PerformLayOut(current.net.nm, algo);
  nms <- V(current.net)$name;
  nodes <- vector(mode="list");
  for(i in 1:length(nms)){
    nodes[[i]] <- list(
      id=nms[i],
      x=pos.xyz[i,1],
      y=pos.xyz[i,2],
      z=pos.xyz[i,3]
    );
  }
  # now only save the node pos to json
  require("RJSONIO");
  netData <- list(nodes=nodes);
  sink(filenm);
  cat(RJSONIO::toJSON(netData));
  sink();
  return(filenm);
}

UpdateNetworkLayout <- function(algo, filenm, focus){
  # get layers
  current.net <- ppi.comps[[current.net.nm]];

  pos.xyz <- PerformLayOut(current.net.nm, algo, focus);
  nms <- V(current.net)$name;
  nodes <- vector(mode="list");
  for(i in 1:length(nms)){
    nodes[[i]] <- list(
      id=nms[i],
      x=pos.xyz[i,1],
      y=pos.xyz[i,2]
    );
  }
  # now only save the node pos to json
  require("RJSONIO");
  netData <- list(nodes=nodes);
  sink(filenm);
  cat(RJSONIO::toJSON(netData));
  sink();
  return(filenm);
}

PerformLayOut <- function(net.nm, algo, focus){
  require("igraph")

  g <- ppi.comps[[net.nm]];
  vc <- vcount(g);
  if(algo == "Default"){
    if(vc > 5000) {
      pos.xy <- layout_with_lgl(g);
    }else if(vc < 150){
      pos.xy <- layout_with_kk(g);
    }else{
      pos.xy <- layout_with_fr(g);
    }
  }else if(algo == "FrR"){
    pos.xy <- layout_with_fr(g, area=34*vc^2);
  }else if(algo == "random"){
    pos.xy <- layout_randomly (g);
  }else if(algo == "lgl"){
    pos.xy <- layout_with_lgl(g);
  }else if(algo == "gopt"){
    # this is a slow one
    if(vc > 3000) {
      maxiter <- 50;
    }else if(vc > 2000) {
      maxiter <- 100;
    }else if(vc > 1000) {
      maxiter <- 200;
    }else{
      maxiter <- 500;
    }
    pos.xy <- layout_with_graphopt(g, niter=maxiter);
  }else if(algo == "fr"){
    pos.xy <- layout_with_fr(g, dim=3, niter=500)
  }else if(algo == "kk"){
    pos.xy <- layout_with_kk(g, dim=3, maxiter=500)
  }else if(algo == "tree"){
    l <- layout_with_sugiyama(g, vgap=vc/4)
    pos.xy <- -l$layout
  }else if(algo == "circular_tripartite"){
    require("ggforce")
    l <- layout_with_sugiyama(g, layers = as.numeric(V(g)$layers)*(vc/3) +30)
    layout <- l$layout

    radial <- radial_trans(
      r.range = rev(range(layout[,2])),
      a.range = range(layout[,1]),
      offset = 0
    )
    coords <- radial$transform(layout[,2], layout[,1])
    layout[,1] <- coords$x
    layout[,2] <- coords$y
    pos.xy <-layout
  }else if(algo == "tripartite"){
    l <- layout_with_sugiyama(g, layers = as.numeric(V(g)$layers)*(vc/4))
    pos.xy <- -l$layout[,2:1]
  }else if(algo == "concentric"){
    require("graphlayouts")
    # the fist element in the list for concentric is the central node.
    if(focus==""){
      inx=1;
    }else{
      inx = which(V(g)$name == focus)
    }
    coords <- layout_with_focus(g,inx)
    pos.xy <- coords$xy
  }else if(algo == "backbone"){
    require("graphlayouts")
    if(length(V(g)$name)<2000){
      coords <- layout_with_stress(g)
      pos.xy <- coords
    }else{
      coords <-layout_with_sparse_stress(g,pivots=100)
      pos.xy <- coords
    }
  }
  pos.xy;
}

# Adapted from PCSF
# https://github.com/IOR-Bioinformatics/PCSF
Compute.SteinerForest <- function(ppi, terminals, w = 2, b = 1, mu = 0.0005, dummies){

  # Gather the terminal genes to be analyzed, and their scores
  terminal_names <- names(terminals)
  terminal_values <- as.numeric(terminals)

  # Incorporate the node prizes
  node_names <- V(ppi)$name
  node_prz <- vector(mode = "numeric", length = length(node_names))
  index <- match(terminal_names, node_names)
  percent <- signif((length(index) - sum(is.na(index)))/length(index)*100, 4)
  if (percent < 5){
    print("Less than 1% of your terminal nodes are matched in the interactome!");
    return(NULL);
  }
  paste0("  ", percent, "% of your terminal nodes are included in the interactome\n");
  terminal_names <- terminal_names[!is.na(index)]
  terminal_values <- terminal_values[!is.na(index)]
  index <- index[!is.na(index)]
  node_prz[index] <-  terminal_values

  if(missing(dummies)||is.null(dummies)||is.na(dummies)){
    dummies <- terminal_names #re-assign this to allow for input
  }

  ## Prepare input file for MST-PCSF implementation in C++

  # Calculate the hub penalization scores
  node_degrees <- igraph::degree(ppi)
  hub_penalization <- - mu*node_degrees

  # Update the node prizes
  node_prizes <- b*node_prz
  index <- which(node_prizes==0)
  node_prizes[index] <- hub_penalization[index]

  # Construct the list of edges
  edges <- ends(ppi,es = E(ppi))
  from <- c(rep("DUMMY", length(dummies)), edges[,1])
  to <- c(dummies, edges[,2])

  cost <- c(rep(w, length(dummies)), E(ppi)$weight)

  #PCSF will faill if there are NAs in weights, this will check and fail gracefully
  if(any(is.na(E(ppi)$weight))){
    print("NAs found in the weight vector!");
    return (NULL);
  }

  ## Feed the input into the PCSF algorithm
  if(.on.public.web){
    #XiaLabCppLib is required here, but load before
    output <- XiaLabCppLib::call_sr(from,to,cost,node_names,node_prizes)
  } else {
    #For R package, this function has already been included internally
    output <- XiaLabCppLib::call_sr(from,to,cost,node_names,node_prizes)
  }


  # Check the size of output subnetwork and print a warning if it is 0
  if(length(output[[1]]) != 0){

    # Contruct an igraph object from the MST-PCSF output
    e <- data.frame(output[[1]], output[[2]], output[[3]])
    e <- e[which(e[,2]!="DUMMY"), ]
    names(e) <- c("from", "to", "weight")

    # Differentiate the type of nodes
    type <- rep("Steiner", length(output[[4]]))
    index <- match(terminal_names, output[[4]])
    index <- index[!is.na(index)]
    type[index] <- "Terminal"

    v <- data.frame(output[[4]], output[[5]], type)
    names(v) <- c("terminals", "prize", "type")
    subnet <- graph.data.frame(e,vertices=v,directed=F)
    #E(subnet)$weight <- as.numeric(output[[3]])
    subnet <- delete_vertices(subnet, "DUMMY")
    subnet <- delete_vertices(subnet, names(which(degree(subnet)==0)));
    return(subnet);
  } else{
    print("Subnetwork can not be identified for a given parameter set")
    return(NULL);
  }
}

convertIgraph2JSON <- function(dataSetObj=NA, net.nm, filenm, thera="FALSE", dim=3){
  if(!exists("my.convert.igraph")){ # public web on same user dir
    compiler::loadcmp("../../rscripts/OmicsNetR/R/utils_convertIgraph2JSON.Rc");
  }
  return(my.convert.igraph(dataSet, net.nm, filenm, thera, dim));
}

# initialize or reset network objects
initNetwork <- function(){
  nodeu.type <<- vector()
  edgeu.res <<- data.frame()
  net.info <<- list();
  nodeu.ids <<- vector()
  nodeu.nms <<- vector()
  expr.vec <<- vector()
  seed.genes <<- vector();
  ppi.comps <<- list();
  overall.graph <<- "";
  net.stats <<- "";
  edgeu.res.list <<- list();
  if(exists("PeakSet",envir = .GlobalEnv)) {
    rm(PeakSet)
  }
}

#' Resetting individual omics networks generated from database selection step
#'
#' @param dataSetObj Input the name of the created dataSetObj (see Init.Data)
#'
#' @export
#'
resetNetwork <- function(dataSetObj=NA){
  dataSet <- .get.nSet(dataSetObj);
  ppi.comps <<- list();
  overall.graph <<- "";
  net.stats <<- "";
  return(.set.nSet(dataSet));
}

SetStringParam <-function(dataSetObj=NA, require.exp=T, min.score=900){
  dataSet <- .get.nSet(dataSetObj);
  dataSet$require.exp <- require.exp;
  dataSet$min.score <- min.score;
  return(.set.nSet(dataSet));
}

#' Query database using input list or previously computed network
#'
#' @param dataSetObj Input the name of the created dataSetObj (see Init.Data)
#' @param type Network type (gene, met, mir, tf, mic, peak, m2m, snp)
#' @param dbType Database name (i.e innatedb)
#' @param inputType Omics type of input features (gene, protein, met, tf, mic, peak, snp)
#'
#' @export
#'
QueryNet <- function(dataSetObj=NA, type="gene", dbType="default", inputType="gene"){
  dataSet <- .get.nSet(dataSetObj);

  orig.inputType <- inputType

  containMsg <- 0;

  require.exp <- dataSet$require.exp
  min.score <- dataSet$min.score

  snpRegion <- dataSet$snpRegion
  edge.res <- data.frame();


  #setting input type so it's displayed correctly on db selection page
  if(type ==  "m2m" ){
    inputType <- "met";
  }else if(inputType %in% c("mic") && type !=  inputType || inputType %in% c("peak") && type !=  "peak" ){
    inputType <- "met";
  }else if(inputType == "snp" && type != "snp"){
    inputType <- "gene";
  }else if (inputType %in% c("mir","tf","met") && type != inputType){
    inputType <- "gene";
  }

  nodes.query = "";

  if(!exists('edgeu.res.list')){
    edgeu.res.list <<- list();
  }

  if (length(edgeu.res.list) > 0) {
    nm <- paste0(type, "_", inputType, "_", data.org, "_", orig.inputType);

    for (i in 1:length(edgeu.res.list)) {
      curr.nm <- names(edgeu.res.list)[i];
      curr.orig.input <- strsplit(curr.nm, "_")[[1]][4];
      curr.type <- strsplit(curr.nm, "_")[[1]][1];
      # to make sure network expansion only uses nodes from primary network depending on input type
      if (curr.nm != nm && curr.orig.input == orig.inputType && .matchInputToPrimaryNet(orig.inputType, curr.type)) {
        edge.res <- rbind(edge.res, edgeu.res.list[[i]]$table);
      }
    }
    if(length(edge.res) > 0){
      nodes.query <- c(edge.res[, 1], edge.res[, 2]);
    }
  }

  old.edges <- edge.res;
  if(orig.inputType == type || orig.inputType %in% c("gene","protein")){
    query.mat <- dataSet$exp.mat[[inputType]];
  }else{
    query.mat <- data.frame();
  }

  if (length(query.mat) == 0) {

    if(.matchInputToPrimaryNet(orig.inputType, type)){
      query.vec <- unique(dataSet$seeds.proteins);
    }else{
      query.vec <- unique(nodes.query);
    }

    if(length(query.vec) == 0 && length(dataSet$exp.mat[["gene"]])>0){
      query.vec <- rownames(dataSet$exp.mat[["gene"]])
      seed.genes <<- c(seed.genes, query.vec);
    }else if(length(query.vec) == 0 && length(dataSet$exp.mat[["ko"]])>0){
      query.vec <- rownames(dataSet$exp.mat[["ko"]])
      seed.genes <<- c(seed.genes, query.vec);
    }

  } else {
    query.vec <- rownames(dataSet$exp.mat[[inputType]]);
    seed.genes <<- c(seed.genes, query.vec);
  }

  protein.vec <- unique(as.character(query.vec));

  if(length(protein.vec) == 0){
    protein.vec <- rownames(dataSet$exp.mat[[inputType]]);
  }


  inv <- .getQueryDir(inputType, type);

  cat(inputType, type, dbType, require.exp, min.score, inv, FALSE,snpRegion, "\n");
  SearchNetDB(dataSet, protein.vec, orig.inputType, inputType, type, dbType, require.exp, min.score, inv, FALSE, snpRegion);

  node.res <- data.frame(Id=nodeu.ids, Label=nodeu.nms);

  edge.res <- data.frame();
  if(length(edgeu.res.list) > 0){
    for(i in 1:length(edgeu.res.list)){
      edge.res <- rbind(edge.res, edgeu.res.list[[i]]$table);
    }
  }

  if(type == "m2m" && orig.inputType == "peak") {
    # to remove the useless expansion
    resclean <- redundancyClean(node.res, edge.res);
    node.res <- resclean$node.res;
    edge.res <- resclean$edge.res;
  }

  #if(length(edge.res)>0 && length(old.edges)>0){
  if(all(dim(edge.res) == dim(old.edges))){
    current.msg <<- paste("No interactions detected. Please make sure to start with a database containing the molecule type of the current queries!")
    containMsg <- 1;
    return(c(0, 0, containMsg));
  }
  #}

  omics.net <<- list(
    netw.type="abc",
    order=1,
    seeds=" ",
    table.nm=" ",
    node.data = node.res,
    edge.data = edge.res,
    require.exp = F,
    min.score = 900
  );

  dataSet$seeds.proteins <- unique(c(dataSet$seeds.proteins, seed.genes));

  dataSet <<- dataSet;

  ComputeIndSubnetStats(dataSet);

  if(.on.public.web){
    .set.nSet(dataSet);
    return(c(nrow(node.res), nrow(edge.res), containMsg));
  }else{
    message("Query Net Finished.")
    return(.set.nSet(dataSet));
  }
}

GetNetworkTopology <- function(netnm){
  g <- ppi.comps[[netnm]];
  globalProperties <-list();
  globalProperties[["Diameter"]] <-diameter(g);
  globalProperties[["Radius"]] <-radius(g);
  globalProperties[["Average path length"]] <-signif(mean_distance(g), 3);
  globalProperties[["Clustering coefficient"]] <- signif(transitivity(g, type="global"), 3);
  propertiesVector <- c(globalProperties[[1]], globalProperties[[2]], globalProperties[[3]], globalProperties[[4]]);
  return(propertiesVector);
}

#' PlotDegreeHistogram
#'
#' @param imgNm image name
#' @param netNm network name
#' @param dpi dpi value
#' @param format format
#' @export
#'
PlotDegreeHistogram <- function(imgNm, netNm = "NA", dpi=72, format="png"){
  library(igraph);
  require("Cairo")
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  Cairo(file=imgNm, width=400, height=400, type="png", bg="white");
  require("ggplot2")
  if(netNm != "NA"){
    overall.graph <- ppi.comps[[netNm]];
  }
  G.degrees <- degree(overall.graph)
  G.degree.histogram <- as.data.frame(table(G.degrees))
  G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])

  p <- ggplot(G.degree.histogram, aes(x = G.degrees, y = Freq)) +
    geom_point() +
    scale_x_continuous("Degree\n(nodes containing that amount of connections)",
                       breaks = c(1, 3, 10, 30, 100, 300),
                       trans = "log10") +
    scale_y_continuous("Frequency\n(number of nodes)",
                       breaks = c(1, 3, 10, 30, 100, 300, 1000),
                       trans = "log10") +
    ggtitle("Degree Distribution (log-log)") +
    theme_bw()  +
    theme(plot.title = element_text(hjust = 0.5));
  print(p);
  dev.off();

  dataSet$imgSet$degreeHistogram <- list();
  dataSet$imgSet$degreeHistogram$netName <- netNm;
  dataSet <<- dataSet;
  return(1);
}

#' PlotBetweennessHistogram
#'
#' @param imgNm image name
#' @param netNm network name
#' @param dpi dpi value
#' @param format format
#' @export
#' @import Cairo
PlotBetweennessHistogram <- function(imgNm, netNm = "NA",dpi=72, format="png"){
  require("Cairo")
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  Cairo(file=imgNm, width=400, height=400, type="png", bg="white");
  require("ggplot2")
  if(netNm != "NA"){
    overall.graph <- ppi.comps[[netNm]];
  }
  G.degrees <- betweenness(overall.graph);
  G.degree.histogram <- as.data.frame(table(G.degrees));
  G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1]);

  p <- ggplot(G.degree.histogram, aes(x = G.degrees, y = Freq)) +
    geom_point() +
    scale_x_continuous("Betweenness\n(nodes with that amount of betweenness)",
                       breaks = c(1, 3, 10, 30, 100, 300,1000,3000,10000,30000),
                       trans = "log10") +
    scale_y_continuous("Frequency\n(number of nodes)",
                       breaks = c(1, 3, 10, 30, 100, 300, 1000),
                       trans = "log10") +
    ggtitle("Betweenness Distribution (log-log)") +
    theme_bw()  +
    theme(plot.title = element_text(hjust = 0.5));
  print(p);
  dev.off();


  dataSet$imgSet$betwennessHistogram <- list();
  dataSet$imgSet$betwennessHistogram$netName <- netNm;
  dataSet <<- dataSet;
  return(1);
}

#' PreparePeaksNetwork
#'
#' @param dataSetObj dataSetObj
#' @export

PreparePeaksNetwork <- function(dataSetObj=NA){
    if(!exists("my.peak.net")){ # public web on same user dir
        compiler::loadcmp("../../rscripts/OmicsNetR/R/utils_peak_net.Rc");
    }
    res <- my.peak.net(dataSetObj);
    return(res);
}


#' Guilt-by-association analysis using Random Walk With Restart
#'
#' @param fileNm Input file name for exporting result table as .csv
#' @param method Method name, only "rwr" is supported
#' @param queryType seeds or highlighted nodes for report saving purposes
#' @param nodeids The IDs of nodes set as seed nodes "; " separated.
#'
#' @export
DoGba <- function(fileNm="NA", method="rwr", queryType="seed", nodeids){
  require("RandomWalkRestartMH")
  require("igraph")

  if(queryType == "seed"){

  nodes <- strsplit(nodeids, ",")[[1]];
  }else{
  nodes <- strsplit(nodeids, ",")[[1]];
  }
  g <- ppi.comps[[current.net.nm]];
  PPI_MultiplexObject <- RandomWalkRestartMH::create.multiplex(list(PPI=g))

  AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
  AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)
  seeds <- nodes;

  RWR_PPI_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI,
                                                   PPI_MultiplexObject,seeds);
  resObj <- RWR_PPI_Results$RWRM_Results;
  if(nrow(resObj) > 100){
    resObj <- resObj[c(1:100),];
  }
  require("RJSONIO");
  sink(fileNm);
  cat(toJSON(resObj));
  sink();

  csvNm <- paste0(gsub(".json", "", fileNm), ".csv");
  df <- RWR_PPI_Results$RWRM_Results;
  colnames(df) <- c("IDs", "Score");
  sym.vec <- doEntrez2SymbolMapping(df$IDs);

  df <- cbind(Name=sym.vec, df);
  fast.write.csv(df, file=csvNm);
  fast.write.csv(df, file="gba_table.csv");

  #record table for report
  type = "gba";
  dataSet$imgSet$enrTables[[type]] <- list();
  dataSet$imgSet$enrTables[[type]]$table <- df;
  dataSet$imgSet$enrTables[[type]]$res.mat<- df[,3, drop=F];

  dataSet$imgSet$enrTables[[type]]$library <- "";
  dataSet$imgSet$enrTables[[type]]$query <- queryType; 
  dataSet$imgSet$enrTables[[type]]$algo <- "Random walk with restart";
  dataSet <<- dataSet;

  if(.on.public.web){
    return(1);
  }else{
    return(RWR_PPI_Results$RWRM_Results);
  }
}

#' FilterByTissue
#' @description Filter gene/protein by tissue
#'
#' @param dataSetObj dataSetObj
#' @param type Input database name
#' @param tissue Input the tissue name
#'
#' @export
#'
FilterByTissue <- function(dataSetObj=NA, type, tissue){
  dataSet <- .get.nSet(dataSetObj);
  gene.nms <- V(overall.graph)$name;

  tissue.mat <- queryFilterDB(type, data.org);
  hit.inx <- tissue.mat$Tissue %in% c(tissue, "All");

  tissue.genes <- tissue.mat[,1][!hit.inx];
  nodes2rm <- gene.nms[which(gene.nms %in% tissue.genes)]

  g <- simplify(delete.vertices(overall.graph, nodes2rm));

  nodeList <- get.data.frame(g, "vertices");
  nodeList <- nodeList[,1:2];
  colnames(nodeList) <- c("Id", "Label");

  fast.write.csv(nodeList, file="orig_node_list.csv");
  nd.inx <- omics.net$node.data[,1] %in% nodeList[,1];

  edgeList <- get.data.frame(g, "edges");
  edgeList <- edgeList[,1:2];
  colnames(edgeList) <- c("Source", "Target");
  fast.write.csv(edgeList, file="orig_edge_list.csv");

  overall.graph <- g
  dataSet <- .decomposeGraph(dataSet, g);
  if(!is.null(dataSet$substats)){
    overall.graph <<- overall.graph;
    if(.on.public.web){
      .set.nSet(dataSet);
      return(c(length(seed.genes),length(seed.proteins), vcount(overall.graph), ecount(overall.graph), length(ppi.comps), dataSet$substats));
    }else{
      return(.set.nSet(dataSet));
    }
  }else{
    return(0);
  }
}

#' Filter peaks by p-value cutoff for peak processing
#'
#' @param pvaluecutoff Input the p-value cutoff
#'
#' @export
#'
FilterByPvalue <- function(pvaluecutoff){
  data <- qs::qread("PeakSet_data.qs")

  if(exists("dataSet",envir = .GlobalEnv)) {
    dataSetObj <- get("dataSet", envir = .GlobalEnv)
  } else {
    dataSetObj <- NA;
  }

  dataSet <- .get.nSet(dataSetObj);
  all.nms <- V(overall.graph)$name;

  nms_ori <- PeakSet[["NetID_output"]][["formula"]];
  nms_ori[PeakSet[["NetID_output"]][["KEGG"]] !=''] <- PeakSet[["NetID_output"]][["KEGG"]][PeakSet[["NetID_output"]][["KEGG"]] !=''];
  nonsigPeaks <- data$id[data$pvalue > pvaluecutoff];
  nodes2rm <- nms_ori[nonsigPeaks];
  nodes2rm <- nodes2rm[nodes2rm != "Unknown"];
  res <- vapply(all.nms, function(x){x %in% nodes2rm}, FUN.VALUE = logical(length = 1));
  nodes2rm <- all.nms[res];

  g <- simplify(delete.vertices(overall.graph, nodes2rm));

  nodeList <- get.data.frame(g, "vertices");
  lbss <- doKegg2NameMapping(nodeList[,1]);
  nodeList <- cbind(nodeList[,1], label = lbss);
  colnames(nodeList) <- c("Id", "Label");

  #fast.write.csv(nodeList, file="orig_node_list.csv");
  nd.inx <- omics.net$node.data[,1] %in% nodeList[,1];

  edgeList <- get.data.frame(g, "edges");
  edgeList <- edgeList[,1:2];
  colnames(edgeList) <- c("Source", "Target");
  #fast.write.csv(edgeList, file="orig_edge_list.csv");

  # update omics.net
  overall.graph <- g;
  dataSet <- .decomposeGraph(dataSet, g);
  if(!is.null(dataSet$substats)){
    overall.graph <<- overall.graph;
    if(.on.public.web){
      .set.nSet(dataSet);
      return(c(length(seed.genes),length(seed.proteins), vcount(overall.graph), ecount(overall.graph), length(ppi.comps), dataSet$substats));
    }else{
      return(.set.nSet(dataSet));
    }
  }else{
    return(0);
  }
}

#' queryFilterDB
#'
#' @param type database type
#' @param org organism key
#' @export
queryFilterDB <- function(type, org){
  require('RSQLite');
  if(!PrepareSqliteDB(paste(sqlite.path, "tissue_filter.sqlite", sep=""), .on.public.web)){
    stop("Sqlite database is missing, please check your internet connection!");
  }
  conv.db <- dbConnect(SQLite(), paste(sqlite.path, "tissue_filter.sqlite", sep=""));
  db.map <- dbReadTable(conv.db, paste0(data.org,"_",type));
  dbDisconnect(conv.db); cleanMem();

  return(db.map)
}

#' SetSnp2GeneOpts
#' @param snpDBType snpDBType
#' @param phescOpt phescOpt
#' @param vepOpt vepOpt
#' @param vepNum vepNum
#' @param vepDis vepDis
#' @export
SetSnp2GeneOpts <- function(snpDBType, phescOpt, vepOpt, vepNum,vepDis){
    if(snpDBType=="PhenoScanner"){
        dataSet$snp2gene$phesc.opt <<- phescOpt;
    }else if(snpDBType=="vep"){
        dataSet$snp2gene$vep.opt <<- vepOpt;
        dataSet$snp2gene$vep.num <<- vepNum;
        dataSet$snp2gene$vep.dis <<- vepDis;
    }
}

#' SetMetPotentialOpts
#' @param metPotentialThresh metPotentialThresh
#' @param currExclude currExclude
#' @param uniExclude uniExclude
#' @param orphExclude orphExclude
#' @export
SetMetPotentialOpts <- function(metPotentialThresh, currExclude, uniExclude,orphExclude){
    dataSet$mic.thresh <<- metPotentialThresh;
    dataSet$currExclude <<- currExclude;
    dataSet$uniExclude <<- uniExclude;
    dataSet$orphExclude <<- orphExclude;
}

#for computing individual networks before merging
ComputeIndSubnetStats <- function(dataSetObj=NA){
  dataSet <- .get.nSet(dataSetObj);
  if(length(edgeu.res.list) == 0){
    return(0)
  }
  net.stats <- as.data.frame(matrix(0, ncol = 6, nrow = length(edgeu.res.list)));
  colnames(net.stats) <- c("Input", "Network","Node", "Edge", "Query", "NetworkValue");
  for(i in 1:length(edgeu.res.list)){
    seed.type <- strsplit(names(edgeu.res.list)[i], "_")[[1]][2];
    seeds <- unique(rownames(dataSet$exp.mat[[seed.type]]));
    edge.df <- edgeu.res.list[[i]]$table;
    edge.num <- dim(edge.df)[1];
    nodes <- unique(unname(unlist(edge.df)));
    node.num <- length(nodes);
    query.num <- sum(seeds %in% nodes);
    netw.type <- strsplit(names(edgeu.res.list[i]), "_")[[1]][1];
    input.type <- strsplit(names(edgeu.res.list[i]), "_")[[1]][2];
    inputName <- convertInputTypes2Names(input.type);
    interactionName <- paste0(convertInteraction2Names(netw.type), " (", edgeu.res.list[[i]]$database ,")");
    net.stats[i,] <- c(inputName, interactionName, node.num,edge.num,query.num, names(edgeu.res.list)[i]);
  }
  saveRDS(dataSet$seeds.proteins, "seeds.rds");

  dataSet$ind.net.stats <- net.stats;
  return(.set.nSet(dataSet));
}

DeleteIndNet <- function(netNm){
  inx <- which(names(edgeu.res.list) == netNm);

  edgeu.res.list <<- edgeu.res.list[-inx];
  resetNetwork();

  edge.res <- data.frame();
  if(length(edgeu.res.list) > 0){
    for(i in 1:length(edgeu.res.list)){
      edge.res <- rbind(edge.res, edgeu.res.list[[i]]$table);
    }
  }
  omics.net$edge.data <- edge.res;
  omics.net <<- omics.net;
  return(1)
}

#' PrepareGraph
#' @param net.nm net.nm
#' @param type type
#' @param export export
#' @export
PrepareGraph <- function(net.nm, type="", export=T){
  #type: subnetwork (networkbuilder page) or ind (database selection page)
  require("igraph");
  if(type == "subnetwork"){
    g <- ppi.comps[[net.nm]];
  }else{
    g <- simplify(graph.data.frame(edgeu.res.list[[net.nm]]$table, directed=FALSE));
    net.value <- net.nm;
    net.nm <- paste0("network_", net.nm);
  }
  file.nm <- paste(net.nm, ".csv", sep="");
  edge.res <- as.data.frame(get.edgelist(g))

  for(i in 1:ncol(edge.res)){
    nms <- edge.res[,i];
    hit.inx <- match(nms, omics.net$node.data[,1]);

    if(length(which(is.na(hit.inx)))>0){
      hit.inx <- match(nms, omics.net$node.data[,2]);
      lbls <- omics.net$node.data[hit.inx, 1];
    }else{
      lbls <- omics.net$node.data[hit.inx, 2];
    }
    edge.res <- cbind(edge.res, lbls);

  }


  colnames(edge.res) <- c("Id1", "Id2", "Name1", "Name2");
  if(export){
    fast.write.csv(edge.res, file=file.nm,row.names=FALSE);
  }else{
    edgeu.res.list.init <<- edgeu.res.list;
    dataSet$viewTable[[net.value]] <- edge.res;
    dataSet <<- dataSet;
  }

  return(file.nm);
}

#' SetPpiZero
#' @param ppiZero ppiZero
#' @export
SetPpiZero <- function(ppiZero){
  dataSet$ppiZero <<- ppiZero;
}


#' Save graph file in OmicsNet json format
#'
#' @param fileNm The file name of output json file
#'
#' @export
#'
SaveNetworkJson <- function(fileNm){
  dataSet <- .get.nSet(dataSet);
  obj <- list();
  obj$dataSet <- list();
  seed.list <- dataSet$seed;

  for ( i in 1:length(seed.list)){
    df <- seed.list[[i]];
    df1 <- cbind(df, rownames(df))
    colnames(df1) <- c("expr", "id");
    df1 <- as.data.frame(df1);
    seed.list[[i]] <- df1;
  }
  obj$allOrgs <- dataSet$allOrgs;
  obj$snpRegion <- dataSet$snpRegion;
  obj$seed.genes <- seed.genes
  obj$seed.proteins <- seed.proteins
  obj$dataSet$seed <- seed.list;
  obj$dataSet$seeds.expr <- dataSet$seeds.expr;
  obj$dataSet$seeds.proteins <- dataSet$seeds.proteins;
  obj$dataSet$type.nums <- dataSet$type.nums;
  obj$dataSet$query.nums <- dataSet$query.nums;
  obj$net.info <- net.info;
  obj$omics.net <- omics.net;
  obj$data.org <- data.org;

  require("rjson");
  sink(fileNm);
  cat(toJSON(obj));
  sink();
}

#gene, met, mir, tf, mic, peak, m2m, snp
.matchInputToPrimaryNet <- function(inputType, netType){
  bool <- F;
  if(inputType == "protein"){
    bool <- netType == "gene";
  }else{
    bool <- netType == inputType;
  }
  return(bool);
}

.getQueryDir <- function(inputType, netType){
  if(inputType == netType){
    direction <- "direct";
  }else{
    direction <- "inverse";
  }
  return(direction);
}


#' SanityCheckSelection
#'
#' @param dataSetObj dataSetObj
#' @export
SanityCheckSelection <- function(dataSetObj=NA){
  dataSet <- .get.nSet(dataSet);
  type.vec <- names(dataSet$exp.mat);
  computed.ok <- 1; #check if every input type has computed for network at least once
  overlap.ok <- 0; #check if there is any overlap between omics networks

  #check compute
  input.from.net.vec <- vector();
  for(i in 1:length(edgeu.res.list)){
    curr.nm <- names(edgeu.res.list)[i];
    curr.orig.input <- strsplit(curr.nm, "_")[[1]][4];
    input.from.net.vec <- c(input.from.net.vec, curr.orig.input);
  }

  input.from.net.vec <- unique(input.from.net.vec);

  for(i in 1:length(type.vec)){
    curr.type <- type.vec[i];
    if(!curr.type %in% input.from.net.vec){
      computed.ok <- 0;
    }
  }

  #check overlap
  e.list.len <- length(edgeu.res.list);
  if(length(type.vec) > 1 && e.list.len>1){
    node.res.list <- lapply(edgeu.res.list, function(x){
      unique(unlist(x$table))
    })
    num.vec <- seq.int(e.list.len);
    cat(e.list.len, "==elen");
    if(e.list.len==2){
      res <- intersect(node.res.list[[1]], node.res.list[[2]])
      cat("length=====", length(res));

      if(length(res)>0){
        overlap.ok <- 1;
      }
    }else{
      pairs <- t(combn(num.vec,2));
      for(i in 1:nrow(pairs)){
        pair <- pairs[i,];
        res <- intersect(node.res.list[[ pair[1] ]], node.res.list[[ pair[2] ]]);
        cat("length=====", length(res));
        if(length(res)>0){
          overlap.ok <- 1;
        }
      }
    }
  }else{
    overlap.ok <- 1;
  }
  if(!.on.public.web){
    return(.set.nSet(dataSet));
  }else{
    return(c(computed.ok,overlap.ok));
  }
}

getNodeTypesByNetType <- function(type){
    if(type == "gene"){
        ndTypes <- c("gene")
    }else if(type == "tf"){
        ndTypes <- c("gene", "tf")
    }else if(type == "mir"){
        ndTypes <- c("gene", "mir")
    }else if(type == "met"){
        ndTypes <- c("gene", "met")
    }else if(type == "mic"){
        ndTypes <- c("mic", "met")
    }else if(type == "peak"){
        ndTypes <- c("putative_met", "met")
    }else if(type == "snp"){
        ndTypes <- c("gene", "snp")
    }
}

# generate subnetwork stats with node type information
GetNetStatByType <- function(g){
  nd.queries <- V(g)$name;
  uniq.ins <- unique(dataSet$seeds.proteins);
  sd.queries <- uniq.ins[uniq.ins %in% nd.queries];

  #check compute
  net.type.vec <- vector();

  if(length(edgeu.res.list) == 0){
    my.stat <- list(
      node.num = vcount(g),
      query.num = sum((unique(dataSet$seeds.proteins)) %in% V(g)$name)
    );
    return(my.stat);
  }

  for(i in 1:length(edgeu.res.list)){
    curr.nm <- names(edgeu.res.list)[i];
    curr.net.type <- strsplit(curr.nm, "_")[[1]][1];
    curr.node.types <- getNodeTypesByNetType(curr.net.type);
    net.type.vec <- c(net.type.vec, curr.node.types);
  }

    vec <- net.type.vec
    vec <- unique(unlist(vec))
    sd.res <- "";
    nd.res <- "";
    for( i in 1:length(vec)){
      if(vec[i] == "mir"){
        nms = net.info$mir.ids
        lbl = "miRNA"
      }else if(vec[i] == "gene"){
        nms = c(net.info$gene.ids, net.info$protein.ids)
        lbl = "mRNA/protein";
      }else if(vec[i] == "tf"){
        nms = net.info$tf.ids
        lbl = "TF";
      }else if(vec[i] == "mic"){
        nms = net.info$mic.ids
        lbl = "Taxon";
      }else if(vec[i] == "met"){
        nms = net.info$met.ids
        lbl = "Metabolite";
      }else if(vec[i] == "snp"){
        nms = net.info$snp.ids
        lbl = "SNP";
      }else if(vec[i] == "putative_met"){
        nms = net.info$peak.ids;
        lbl = "Putative metabolite";
      }

      nms <- unique(nms);
      if(sum(nms %in% nd.queries)>0 && !grepl(lbl, nd.res)){
        nd.res <- paste0(lbl,": ", sum(nms %in% nd.queries), "; ", nd.res)
      }
      if(sum(nms %in% sd.queries)>0){
        sd.res <- paste0(lbl,": ", sum(nms %in% sd.queries), "; ", sd.res)
      }
    }
    my.stat <- list(
      node.num = nd.res,
      query.num = sd.res
    );

  return(my.stat);
}

## subnetwork seed node type and number
GetQueryNum <-function(){
  return(dataSet$query.nums)
}

## subnetwork overall node type and number
GetTypeNum <-function(){
  return(dataSet$type.nums)
}
