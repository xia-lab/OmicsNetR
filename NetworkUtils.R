##################################################
## R scripts for OmicsNet 
## Description: network analysis methods
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################


SetGraphFileName <- function(fileName){
  graph_name <<- fileName;
}

.prepareSigProteinJSON <- function(type){
  result.list <- .prepareListSeeds(type);
  return(result.list);
}

.prepareListSeeds <- function(type){
  type <<- type;
  protein.list <- list();
  gene.list <- list();
  # return a json array object
  # each object for a single dataset its sig proteins
  meta.vec <- meta.gene.vec <- meta.seed.expr <- NULL;
  
  GeneAnotDB <- NULL;
  if(type %in% c("gene", "geneonly", "protein")){
    gene.mat <- dataSet$mat[["gene"]];
    prot.mat <- dataSet$exp.mat[["gene"]];
  }else{
    gene.mat<-dataSet$mat[[type]];
    prot.mat<-dataSet$exp.mat[[type]];
  }
  if(length(prot.mat) == 0){
    gene.mat <- dataSet$mat[["gene"]];
    prot.mat <- dataSet$exp.mat[["gene"]];
  }
  
  meta.gene.vec <- c(meta.gene.vec, rownames(gene.mat));
  gene.list[[dataSet$name]] <- list(gene=rownames(gene.mat),logFC=unname(gene.mat[,1]));
  GeneAnotDB <- rbind(GeneAnotDB, dataSet$GeneAnotDB);
  meta.seed.expr <- c(meta.seed.expr, prot.mat[,1]);
  protein.vec <- prot.mat[,1];
  meta.vec <- c(meta.vec, names(protein.vec));
  if(length(protein.vec) == 1){
    protein.vec <- as.matrix(protein.vec)
  }   
  protein.list[[dataSet$name]] <-protein.vec
  
  gene.list$name <- dataSet$name;
  seed.genes <<- unique(meta.gene.vec);
  
  meta.seed.df <- as.matrix(meta.seed.expr);
  rownames(meta.seed.df) <- names(meta.seed.expr);
  
  seed.expr <- RemoveDuplicates(meta.seed.df, "max", quiet=F); 
  seed.expr <<- seed.expr[,1];
  protein.vec <- unique(meta.vec);
  
  list(
    gene.list = gene.list,
    protein.list = protein.list,
    protein.vec = protein.vec,
    type = type
  );
}

# zero-order network - create ppi nets from only input (seeds)
BuildSeedProteinNet <- function(){
  nodes = V(overall.graph)$name;
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
  substats <- DecomposeGraph(g);    
  if(is.null(substats)){
    return(0) 
  }else{
    # update omics.net 
    omics.net$order = 0;
    omics.net$node.data <- nodeList;
    omics.net$edge.data <- edgeList;
    omics.net <<- omics.net;
    return(1)
  }
}

# create igraph from the edgelist saved from graph DB
# and decompose into subnets
CreateGraph <- function(){
  require('igraph');
  node.list <- omics.net$node.data;
  edge.list <- omics.net$edge.data;
  seed.proteins <- omics.net$node.data[,1];
  overall.graph <- simplify(graph.data.frame(edge.list, directed=FALSE)) #, vertices=node.list));
  
  # add node expression value
  if(omics.net$db.type == "ppi"){
    newIDs <- doPpiIDMapping(names(seed.expr), "id")$accession;
  }else{# all entrez in mirNet
    newIDs <- dataSet$seeds.proteins;
  }
  match.index <- match(V(overall.graph)$name, newIDs);
  expr.vals <- dataSet$seeds.expr[match.index];
  names(expr.vals) <- rownames(dataSet$seeds.expr)[match.index];
  expr.vec <- vector();
  for(i in 1:length(dataSet$seed)){
    expr.vec <- c(expr.vec, unlist(dataSet$seed[[i]][,1]));
  }
  match.index <- names(expr.vec) %in% V(overall.graph)$name;
  expr.vec <- expr.vec[match.index];
  expr.vec <<- expr.vec;
  current.overall.graph <- overall.graph;
  overall.graph <- set.vertex.attribute(overall.graph, "abundance", index = V(overall.graph), value = expr.vals);
  overall.graph <<- overall.graph;
  substats <- DecomposeGraph(overall.graph);      
  g <<- overall.graph
  seed.proteins <<- seed.proteins
  if(!is.null(substats)){
    return(c(length(seed.genes), length(seed.proteins), net.stats$Node[1], net.stats$Edge[1], length(ppi.comps), substats));        
  }else{
    overall.graph <<- current.overall.graph;
    return(0);
  }
}

PrepareTheraNet <- function(dataSetObj=NA, net.nm, json.nm, theraids, regids, regnms){
  dataSet <- .get.nSet(dataSetObj);  
  net.nm <- gsub(".json", "", net.nm);
  
  if(is.null(net.nm)){
    net.nm <- names(ppi.comps)[1];
  }
  
  reg.count <- reg.count + 1;
  reg.nm <- paste("addedReg", reg.count, sep="");
  
  my.ppi <- ppi.comps[[net.nm]];
  net.nm <<- net.nm;
  json.nm <<- json.nm;
  theraids <<- theraids;
  regids <<- regids;
  regnms <<- regnms;
  
  theraids <- strsplit(theraids, ";");
  theraids <- theraids[[1]];
  theraids <- strsplit(theraids, ",");
  
  regids <- strsplit(regids, ",");
  regids <- regids[[1]]
  regnms <- strsplit(regnms, ",");
  regnms <- regnms[[1]]
  regnms <<- regnms;
  regids <<- regids;
  my.ppi <- ppi.comps[[net.nm]];
  
  my.ppi <- my.ppi + vertices(regids);
  for(i in 1:length(regids)){
    for(j in 1:length(theraids[[i]])){
      my.ppi <- my.ppi + edge(regids[[i]], theraids[[i]][j]);
    }
  }
  my.ppi <<- my.ppi;
  edge.data <- data.frame(regids, regnms, stringsAsFactors=FALSE);
  
  colnames(edge.data) <- c("Id", "Label");
  omics.net$node.data <<- rbind(omics.net$node.data, edge.data);
  filenm <- paste(reg.nm, ".json", sep="");
  reg.count <<- reg.count
  ppi.comps[[reg.nm]] <<- my.ppi;
  UpdateSubnetStats();
  net.info$reg.ids <<- regids;
  convertIgraph2JSON(dataSet, reg.nm, filenm, TRUE);
  if(.on.public.web){
  .set.nSet(dataSetObj);  
  return(filenm);
  }else{
  return(.set.nSet(dataSetObj));  
  }
 
}

PrepareNetwork <- function(dataSetObj=NA, net.nm, json.nm){
  dataSet <- .get.nSet(dataSetObj); 
  if(is.null(net.nm)){
    net.nm <- names(ppi.comps)[1];
  }
  
  my.ppi <- ppi.comps[[net.nm]];
  nd.nms <- V(my.ppi)$name;
  
  GeneAnotDB <- doProteinIDMapping(nd.nms, "entrez");
  
  entrezIDs <- GeneAnotDB[,1];
  names(entrezIDs) <- nd.nms;
  current.anot <<- entrezIDs;
  current.net.nm <<- net.nm;
  if(uploadedGraph){
    convertIgraph2JSONFromFile(net.nm, json.nm, 3);
  }else{
    convertIgraph2JSON(dataSet, net.nm, json.nm, FALSE);
  }

  return(.set.nSet(dataSetObj));
}

PrepareMinNetwork <- function(dataSetObj=NA, net.nm, json.nm){
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

GetShortestPaths <- function(from, to, intermediate="false"){
  current.net <- ppi.comps[[current.net.nm]];
  
  paths <- get.all.shortest.paths(current.net, from, to)$res;
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
    edge.list <- as_edgelist(current.net)
    edges.from <- edge.list[which(edge.list[,1] == from), ]
    edges.from <- rbind(edges.from, edge.list[which(edge.list[,2] == from), ])
    edges.to <- edge.list[which(edge.list[,1] == to), ]
    edges.to  <- rbind(edges.to, edge.list[which(edge.list[,2] == to), ])
    edges.from.universe <- unique(as.vector(edges.from))
    edges.to.universe <- unique(as.vector(edges.to))
    intersects <- intersect(edges.from.universe, edges.to.universe)
    if(length(intersects) != 0){
      for(i in 1:length(intersects)){
        element <- paste(from, intersects[i], to, sep="->")
        if(!element %in% path.vec){
          path.vec[length(path.vec) + i] <- paste(from, intersects[i], to, sep="->")
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
  all_gene_int <- merge(table(filt_net[,2]),table(filt_net[,1]),all=TRUE)
  neighbors <- aggregate(all_gene_int[,2],by=list(all_gene_int[,1]),sum)
  names(neighbors) <- c("Protein","n")
  data.frame(Protein=neighbors$Protein,weight=1/neighbors$n)
}

measureCalc <- function(protein_pvals,protein_wts,weightamt = 10){
  prot_data <- merge(protein_pvals,protein_wts,by='Protein',all=TRUE)
  prot_data$measure <- -log((prot_data$pvalue)*(prot_data$weight)^(1/weightamt))
  net_no_pval <- is.na(prot_data$measure) & !is.na(prot_data$weight)
  prot_data[net_no_pval,'measure'] <- 0
  prot_data
}

scoreClusters <- function(clust_2_protein,prot_meas){
  clust_score_total <- apply(clust_2_protein[,-c(1,2)],1,function(x){do.call(m,list(prot_meas[prot_meas$Protein%in%x[!is.na(x)],'measure']))})
  names(clust_score_total) <- clust_2_protein[,1]
  clust_score_total
}

ExtractModule<- function(dataSetObj=NA, nodeids, dim="3"){
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
      paths.list[[pos]] <- get.shortest.paths(g, nodes[pos], nodes[-(1:pos)])$vpath;
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
  if(uploadedGraph == "false"){
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

SearchNetDBX <- function(db.type, table.nm, require.exp=TRUE, 
                         min.score = 900, netInv, zero=FALSE){
  
  if("peak" %in% dataSet$type){
    netInv = "direct";
  }
  
  edge.res <-data.frame();
  result.list <- result.listu;
  protein.vec <- result.list$protein.vec;
  cat(db.type, netInv, "\n")
  
  require(RJSONIO);
  
  # now do the database search
  if(db.type == "ppi" || db.type == "gene" || db.type == "protein"){
    table.nm = paste(data.org, net.type, sep="_");
    if(!data.org %in% c("mmu", "hsa") && net.type == "string"){
      protein.vec <- doPpiIDMapping(protein.vec, "id");
    }
    
    res <- QueryPpiSQLite(table.nm, protein.vec, require.exp, min.score);
    if(net.order == "zero" || zero){
      hit.inx1 <- res[,1] %in% protein.vec;
      hit.inx2 <- res[,2] %in% protein.vec;
      res <- res[(hit.inx1 & hit.inx2),];
      n.ids <- c(res[,1], res[,2])
    }
    
    edge.res <- data.frame();
    if(length(edgeu.res.list)>0){
      for(i in 1:length(edgeu.res.list)){
        edge.res <- rbind(edge.res, edgeu.res.list[[i]]);
      }
      nodes.current <- unique(c(edge.res[,1], edge.res[,2]));
    }else{
      nodes.current <- vector();
    }
    #print(length(nodes.current) > 0 && !anchor_type %in% c("gene","protein") && !"snp" %in% names(edgeu.res.list))
    
    if(length(nodes.current) > 0 && !anchor_type %in% c("gene","protein") && !"snp" %in% names(edgeu.res.list) ){
      inx = which(res$id1 %in% nodes.current & res$id2 %in% nodes.current)
      res <- res[inx,]
    }
    res = na.omit(res);
    
    # no hits
    if(nrow(res)==0){ return(c(0,0)); }
    edge.res <- data.frame(Source=res$id1,Target=res$id2, stringsAsFactors=FALSE);
    if(nrow(res)!=0){
      row.names(edge.res) <- 1:nrow(res);
    }
    fast.write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    node.ids <- c(res[,1], res[,2])
    node.nms <- c(res[,3], res[,4])
    
    if(anchor_type == "gene" || anchor_type == "protein"){
      net.info$int.ids <- unique(node.ids);
      seed.genes_entrez <<- net.info$int.ids;
    }
    
    edgeu.res <<- rbind(edgeu.res, edge.res)
    nodeu.ids <<- c(nodeu.ids, node.ids)   
    nodeu.nms <<- c(nodeu.nms, node.nms)
    if(db.type == "gene" || db.type == "ppi" ){
      net.info$protein.ids <- unique(c(net.info$protein.ids, unique(node.ids)));
    }else if(db.type == "protein"){
      net.info$protein.ids <- c(net.info$protein.ids, unique(node.ids));
    }else{
      net.info$snp.ids <- c(net.info$snp.ids, unique(node.ids));
    }
    net.info$int.ids <- c(net.info$int.ids, node.ids[!(node.ids %in% unique(protein.vec))]);
    if(!zero){
      net.info$ppi.ids <- unique(node.ids);
    }
  } else if(db.type == "snp"){
    require('RSQLite');
    db.path <- paste(sqlite.path, "snp_annot", sep="");
    table.nm <- "snp_annot";
    col.nm <- "rsid";
    res <- Query.mGWASDB(db.path, protein.vec, table.nm, col.nm);
    
    snp.list <- list()
    snp.type.list <- list()
    
    for( i in 1:length(unique(res$entrez))){
      geneChr <- unique(res$entrez)[i];
      snp.list[[geneChr]] <- as.vector(res[which(res$entrez == geneChr), "rsid"])
      snp.type.list[[geneChr]] <- as.vector(res[which(res$entrez == geneChr), "consequence"])
    }
    snpTable <<- snp.list;
    
    edge.res <- data.frame(Source=res[,"rsid"],Target=res[,"entrez"], stringsAsFactors=FALSE);              
    if(nrow(res)!=0){
      row.names(edge.res) <- 1:nrow(res);
    }
    fast.write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    node.ids <- c(res[,"rsid"], res[,"entrez"])
    symb <- doEntrez2SymbolMapping(res[,"entrez"]);
    node.nms <- c(res[,"rsid"], symb);
    
    edgeu.res <<- rbind(edgeu.res, edge.res)
    edgeu.res <<- edgeu.res[!duplicated(edgeu.res), ]
    nodeu.ids <<- c(nodeu.ids, node.ids)
    nodeu.nms <<- c(nodeu.nms, node.nms)
    if(netInv == "direct"){
      net.info$gene.ids <- unique(c(net.info$gene.ids, node.ids[!node.ids %in% protein.vec]))
      net.info$snp.ids <- unique(protein.vec)
    }else{
      net.info$snp.ids <- unique(res[,"rsid"])
      if(!zero){
        net.info$gene.ids <- c(unique(protein.vec))
      }
    }
    
    if(anchor_type == "snp"){      
      net.info$int.ids <- unique(node.ids);
      seed.genes_entrez <<- net.info$int.ids;
    }
    
    net.info$snpi.ids <- unique(node.ids);
    
  } else if (db.type == "tf") {
    
    table.nm <- paste(data.org, tf.type, sep="_");
    res <- QueryTFSQLite(table.nm, protein.vec, netInv);
    # no hits
    if(nrow(res)==0){ return(c(0,0)); }
    
    if(length(dataSet$type) >1 && anchor_type != "tf" && length(net.info$int.ids) != 0){
      if(build.opt == "seed"){
        res = res[which(res$entrez %in% seed.genes),];
      }else{
        res = res[which(res$entrez %in% nodeu.ids),];
      }
    }
    
    edge.res <- data.frame(Source=res[,"tfid"],Target=res[,"entrez"], stringsAsFactors=FALSE);              
    if(nrow(res)!=0){
      row.names(edge.res) <- 1:nrow(res);
    }
    fast.write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    node.ids <- c(res[,"entrez"], res[,"tfid"])
    node.nms <- c(res[,"symbol"], res[,"tfname"]);
    
    edgeu.res <<- rbind(edgeu.res, edge.res)
    edgeu.res <<- edgeu.res[!duplicated(edgeu.res), ]
    nodeu.ids <<- c(nodeu.ids, node.ids)
    nodeu.nms <<- c(nodeu.nms, node.nms)
    if(netInv == "direct"){
      net.info$gene.ids <- unique(c(net.info$gene.ids, node.ids[!node.ids %in% protein.vec]))
      net.info$tf.ids <-  unique(res[,"tfid"])
    }else{
      net.info$tf.ids <- unique(res[,"tfid"])
      if(!zero){
        net.info$gene.ids <- c(unique(protein.vec))
      }
    }
    if(anchor_type == "tf"){      
      net.info$int.ids <- unique(node.ids);
      seed.genes_entrez <<- net.info$int.ids;
    }
        
  } else if(db.type == "mir") { # in miRNA, table name is org code, colname is id type   
    
    if(mir.type =="targetscan"){
      table.nm <- data.org
    }else{
      table.nm <- data.org
    }    
    res <- QueryMirSQLite(table.nm, "mir_id", protein.vec, netInv, mir.type);
    if(nrow(res)==0){ return(c(0,0)); }
    # no hits
    
    if(length(dataSet$type) >1 && anchor_type != "mir" && length(net.info$int.ids) != 0) {
      if(build.opt == "seed"){
        res = res[which(res$entrez %in% seed.genes),];
      }else{
        res = res[which(res$entrez %in% nodeu.ids),];
      }
    }
    
    edge.res <- data.frame(Source=res[,"mir_id"],Target=res[,"entrez"],stringsAsFactors = FALSE)        
    if(nrow(res)!=0){
      row.names(edge.res) <- 1:nrow(res);
    }
    fast.write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    
    node.ids <- c(res[,"mir_id"], res[,"entrez"]);
    node.nms <- c(res[,"mir_acc"], res[,"symbol"]);
    
    
    edgeu.res <<- rbind(edgeu.res, edge.res)
    edgeu.res <<- edgeu.res[!duplicated(edgeu.res), ]
    nodeu.ids <<- c(nodeu.ids, node.ids)
    nodeu.nms <<- c(nodeu.nms, node.nms)
    if(netInv == "direct"){
      net.info$gene.ids <- unique(c(net.info$gene.ids, node.ids[!node.ids %in% protein.vec]))
      net.info$mir.ids <- unique(res[,"mir_id"])
    }else{
      net.info$mir.ids <- unique(res[,"mir_id"])
      if(!zero){
        net.info$gene.ids <- c(net.info$gene.ids, unique(res[,"entrez"]))
      }
    }
    
    if(anchor_type == "mir"){      
      net.info$int.ids <- unique(node.ids);
      seed.genes_entrez <<- net.info$int.ids;
    }
    
  } else if (db.type == "ko"){
    table.nm <- "ko"
    #print(table.nm)
    if(netInv == "direct"){
      protein.vec <- rownames(dataSet$exp.mat[["ko"]]);
      seed.genes <<- c(seed.genes, protein.vec); 
    }
    res <- QueryKoSQLiteNet(table.nm, protein.vec, netInv);
    
    # no hits
    if(nrow(res)==0){ return(c(0,0)); }
    
    edge.res <- data.frame(Source=res[,"ko"],Target=res[,"kegg"], stringsAsFactors=FALSE);
    
    if(nrow(res)!=0){
      row.names(edge.res) <- 1:nrow(res);
    }
    fast.write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    node.ids <- c(res[,"ko"], res[,"kegg"])
    node.nms <- c(res[,"enzyme"], res[,"compound"]);
    
    edgeu.res <<- rbind(edgeu.res, edge.res)
    nodeu.ids <<- c(nodeu.ids, node.ids)
    nodeu.nms <<- c(nodeu.nms, node.nms)
    if(netInv == "direct"){
      net.info$protein.ids <- unique(protein.vec)
      net.info$met.ids <- node.ids[!node.ids %in% protein.vec]
    }else{
      net.info$protein.ids <- node.ids[!node.ids %in% protein.vec]
      if(!zero){
        net.info$protein.ids <- c(net.info$gene.ids, unique(protein.vec))
      }
    }
    
    if(anchor_type == "ko"){      
      net.info$int.ids <- unique(node.ids);
      seed.genes_entrez <<- net.info$int.ids;
    }
  } else if(db.type == "met") {
    ## metabolite-protein
    if(met.type == "keggp"){
      table.nm <- paste("keggp", sep="_");
    }else{
      table.nm <- paste(data.org, met.type, sep="_");
    }
    cat(table.nm, "===========", "\n");
    res <- QueryMetSQLiteNet(table.nm, protein.vec, netInv);
    # no hits
    if(nrow(res)==0){ return(c(0,0)); }
    if(!"peak" %in% dataSet$type && !"mic" %in% dataSet$type){
      if(length(dataSet$type) >1 && anchor_type != "met" && length(net.info$int.ids) != 0){
        if(build.opt == "seed"){
          res = res[which(res$entrez %in% seed.genes),];
        }else{
          res = res[which(res$entrez %in% nodeu.ids),];
        }
      }
    }
    
    if(met.type =="keggp"){ # project to kegg
      res1 = res[which(res$entrez %in% c(dataSet$seeds.proteins)),];
      res2 = res[which(res$kegg %in% c(dataSet$seeds.proteins)),];
      res = rbind(res1, res2);
    }
    
    edge.res <- data.frame(Source=res[,"kegg"],
                           Target=res[,"entrez"], 
                           stringsAsFactors=FALSE);
    
    if(nrow(res)!=0){
      row.names(edge.res) <- 1:nrow(res);
    }
    fast.write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    node.ids <- c(res[,"kegg"], res[,"entrez"])
    node.nms <- c(res[,"met"], res[,"symbol"]);
    
    if("mic" %in% names(dataSet$exp.mat)){
      seed.genes <<- c(seed.genes, res[,"kegg"]);
    }
    
    edgeu.res <<- rbind(edgeu.res, edge.res)
    nodeu.ids <<- c(nodeu.ids, node.ids)
    nodeu.nms <<- c(nodeu.nms, node.nms)
    if(netInv == "direct"){
      cat(curr.types, "=====", "\n")
      if("peak" %in% curr.types){
        not.gene.ids <- c(net.info$peak.ids, 
                          net.info$met.ids, net.info$kncmpd.ids)
        net.info$met.ids <- c(net.info$met.ids, unique(res[,"kegg"][!res[,"kegg"] %in% not.gene.ids]))
        net.info$protein.ids <- c(net.info$protein.ids, 
                                  unique( res[,"entrez"]))
        net.info$int.ids <- net.info$gene.ids
      }else{
        net.info$met.ids <- c(net.info$met.ids, unique(res[,"kegg"]))
        net.info$protein.ids <- c(net.info$protein.ids, 
                                  unique( res[,"entrez"]))
      }
    }else{
      net.info$met.ids <- node.ids[!node.ids %in% protein.vec]
      if(!zero){
        net.info$protein.ids <- c(net.info$protein.ids, unique(protein.vec))
      }
    }
    
    if(anchor_type == "met"){      
      net.info$int.ids <- unique(node.ids);
      seed.genes_entrez <<- net.info$int.ids;
    }
  } else if(db.type == "peak") {
    
    net.info$met.ids <- PeakSet$mets
    net.info$peak.ids <- PeakSet$put.mets 
    seed.genes <<- unique(c(seed.genes, PeakSet$mets));
    node.ids <- PeakSet$nodes.df[,1]
    node.nms <- PeakSet$nodes.df[,2]
    edge.res <- PeakSet$edges.df
    edgeu.res <<- rbind(edgeu.res, edge.res)
    nodeu.ids <<- c(nodeu.ids, node.ids)
    nodeu.nms <<- c(nodeu.nms, node.nms)
    
    net.info$meti.ids <- unique(node.ids);
    
    if(anchor_type == "peak"){      
      net.info$int.ids <- unique(node.ids);
      seed.genes_entrez <<- net.info$int.ids;
    }
  } else if(db.type == "mic") { # in mic
    # in mic
    sql.name <-paste0("omicsnet_", dataSet$mic.type, ".sqlite");
    table.nm <- mic.taxa;
    res <- QueryMicSQLite(protein.vec, table.nm, sql.name, dataSet$mic.thresh, dataSet$mic.exclude.opt);
    if(nrow(res)==0){ return(c(0,0)); }
    # no hits
    if(length(dataSet$type) >1 && anchor_type != "mic") {
      if(build.opt == "seed"){
        res = res[which(res$entrez %in% seed.genes),];
      }else{
        res = res[which(res$entrez %in% nodeu.ids),];
      }
    }
    colnames(res)[colnames(res) == table.nm] <- "taxa"
    edge.res <- data.frame(Source=res[,"taxa"],Target=res[,"KEGG"],Score=res[,"potential"],stringsAsFactors = FALSE)        
    if(nrow(res)!=0){
      row.names(edge.res) <- 1:nrow(res);
    }
    fast.write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    
    node.ids <- c(res[,"KEGG"]);
    node.nms <- c(res[,"metabolite"]);
    
    nodeu.ids <<- c(nodeu.ids, node.ids)
    nodeu.nms <<- c(nodeu.nms, node.nms)
    if(netInv == "direct"){
      net.info$met.ids <- unique(c(net.info$met.ids, node.ids[!node.ids %in% protein.vec]))
      net.info$mic.ids <- unique(protein.vec)
    }else{
      net.info$mic.ids <- node.ids[!node.ids %in% protein.vec]
      if(!zero){
        net.info$met.ids <- c(net.info$met.ids, unique(protein.vec))
      }
    }
    net.info$mici.ids <- unique(node.ids);
    
    library(igraph);
    g <- simplify(graph.data.frame(edge.res, directed=FALSE)) #, vertices=node.list));
    dgrs <- degree(g);
    dataSet$mic.dgrs <<- dgrs;
    met.ids <- net.info$met.ids
    met.microbe.list <- list()
    met.microbe.score.list <- list()
    mic.table <- list()
    for( i in 1:length(met.ids)){
      met <- met.ids[i]
      met.microbe.list[[met]] <- as.vector(edge.res[which(edge.res$Target == met), 1])
      met.microbe.score.list[[met]] <- as.vector(edge.res[which(edge.res$Target == met), 3])
    }
    for(i in 1:nrow(res)){
      mic.table[[i]] <- unlist(res[i, ]);
      exp.mat <- dataSet$exp.mat$mic;
      exp <- exp.mat[rownames(exp.mat) == res[i,"taxa"]];
      if(exp != 0){
        mic.table[[i]]["AbundPred"] <- as.numeric(exp) * as.numeric(mic.table[[i]]["potential"])
        mic.table[[i]]["potential"] <- as.numeric(mic.table[[i]]["potential"])
      }
    }
    
    
    dataSet$met.mic <- met.microbe.list;
    dataSet$met.mic.score <- met.microbe.score.list;
    dataSet$met.mic.table <- mic.table;
    
    dataSet <<- dataSet;
  } else if (db.type == "m2m") {
    if (m2m.type == "kegg"){
      table.nm <- paste(data.org, "KEGG_m2m", sep="_");
    } else if (m2m.type == "keggp") {
      table.nm <- "ko_KEGG_m2m";
    } else {
      ## for recon 3
      table.nm <- "hsa_recon3_m2m";
    }
    
    if (data.org == "microbiome" | data.org == "NA") {
      #TODO: to consider this option later
      table.nm <-"ko_KEGG_m2m";
      data.org <<- "microbiome";
    }
    
    res <- extendMetPeakNetwork(table.nm)
    if(is.list(res)) {
      edge.res <- res$edge.res;
      net.info <- res$net.info;
    } else {
      edge.res <- NULL;
    }
    
  } else if (db.type == "m2m.mic"){
    
    if(m2m.type == "keggp" | m2m.type == "kegg"){
      table.nm <- paste(data.org, "KEGG_m2m", sep="_");
    } else {
      table.nm <- paste(data.org, "m2m", sep="_");
    }
    cat("data.org is", data.org, "\n")
    if (data.org == "microbiome" | data.org == "NA") {
      #TODO: to consider this option later
      #table.nm <-"ko_KEGG_m2m";
      #data.org <<- "microbiome";
    }
    res <- QueryM2mSQLiteNet(table.nm, protein.vec, netInv);
    
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
    node.nms <- c(res[,"productNM"], res[,"productNM"]);
    
    edgeu.res <<- rbind(edgeu.res, edge.res)
    nodeu.ids <<- c(nodeu.ids, node.ids)
    nodeu.nms <<- c(nodeu.nms, node.nms)
    
    if(netInv == "direct"){
      if("peak" %in% dataSet$type){
        not.kncmpd.ids <- c(net.info$peak.ids, 
                            net.info$met.ids)
        net.info$kncmpd.ids <- c(net.info$kncmpd.ids, 
                                 node.ids[!node.ids %in% not.kncmpd.ids])
        net.info$int.ids <- net.info$kncmpd.ids
      }else{
        net.info$met.ids <- unique(protein.vec)
        net.info$kncmpd.ids <- c(net.info$kncmpd.ids, 
                                 node.ids[!node.ids %in% protein.vec])
        
      }
    }else{
      net.info$met.ids <- node.ids[!node.ids %in% protein.vec]
      if(!zero){
        net.info$kncmpd.ids <- c(net.info$kncmpd.ids, unique(protein.vec))
      }
    }
    net.info$meti.ids <- unique(node.ids);
    cat("m2m.mic running reaches here \n")
  }
  
  if(length(edge.res)>0 && db.type != "mic"){
    edgeu.res.list[[db.type]] <<- unique(edge.res[,c(1,2)]);
  }
  net.info <<- net.info;
  return(1);
}

SetTfType <- function(tfType){
  tf.type <<- tfType;
}

SetMetType <- function(metType){
  met.type <<- metType;
}

SetM2mType <- function(metType){
  m2m.type <<- metType;
}

SetMirType <- function(mirType){
  mir.type <<- mirType;
}

SetNetOrder <- function(netOrder){
  net.order <<- netOrder;
}

SetBuildOpt<- function(opt){
  build.opt <<- opt;
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
    paths = get.all.shortest.paths(PPI_graph,subtree, nsubtree)
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

FilterBipartiNet <- function(nd.type, min.dgr, min.btw){
  
  all.nms <- V(overall.graph)$name;
  edge.mat <- get.edgelist(overall.graph);
  dgrs <- degree(overall.graph);
  nodes2rm.dgr <- nodes2rm.btw <- NULL;
  
  if(nd.type == "gene"){
    hit.inx <- all.nms %in% net.info$int.ids;
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
  substats <- DecomposeGraph(overall.graph);
  if(!is.null(substats)){
    overall.graph <<- overall.graph;
    return(c(length(seed.genes),length(seed.proteins), vcount(overall.graph), ecount(overall.graph), length(ppi.comps), substats));
  }else{
    return(0);
  }
}

SetAnchorType <- function(type){
  anchor_type <<- type;
}

BuildMinConnectedGraphs <- function(dataSetObj=NA, max.len = 200){
  dataSet <- .get.nSet(dataSetObj);  
  set.seed(8574);
  # first get shortest paths for all pair-wise seeds
  hit.inx <- unique(dataSet$seeds.proteins) %in% V(overall.graph)$name;
  my.seeds = dataSet$seeds.proteins[hit.inx];
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
    paths.list[[pos]] <- get.shortest.paths(overall.graph, my.seeds[pos], seed.proteins[-(1:pos)])$vpath;
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
  substats <<- DecomposeGraph(g);
  if(!is.null(substats)){
    overall.graph <<- g;
    if(.on.public.web){
      .set.nSet(dataSet)
      return(c(length(seed.genes),length(seed.proteins), nrow(nodeList), nrow(edgeList), length(ppi.comps), substats));
    }else{
      return(.set.nSet(dataSet));
    }
  }else{
    return(0);
  }
}

BuildPCSFNet <- function(dataSetObj=NA){
  dataSet <- .get.nSet(dataSetObj);  
  edg <- as.data.frame(get.edgelist(overall.graph));
  edg$V3 <- rep(1, nrow(edg));
  colnames(edg) <- c("from", "to", "cost");
  
  node_names <- unique(c(as.character(edg[,1]),as.character(edg[,2])))
  ppi <- graph.data.frame(edg[,1:2],vertices=node_names,directed=F)
  E(ppi)$weight <- as.numeric(edg[,3])
  ppi <- simplify(ppi)
  
  if(sum(expr.vec) == 0){ # make sure weights are not 0?!
    expr.vec <- expr.vec +1
  }
  expr.vec <- abs(expr.vec)
  
  g <- Compute.SteinerForest(ppi, expr.vec, w = 5, b = 100, mu = 0.0005);
  
  nodeList <- get.data.frame(g, "vertices");
  colnames(nodeList) <- c("Id", "Label");
  fast.write.csv(nodeList, file="orig_node_list.csv", row.names=F);
  
  edgeList <- get.data.frame(g, "edges");
  edgeList <- cbind(rownames(edgeList), edgeList);
  colnames(edgeList) <- c("Id", "Source", "Target");
  fast.write.csv(edgeList, file="orig_edge_list.csv", row.names=F);
  
  path.list <- NULL;
  substats <- DecomposeGraph(g);
  if(!is.null(substats)){
    overall.graph <<- g;
    current.msg<<- "Steiner Forest was completed successfully";
    if(.on.public.web){
      .set.nSet(dataSet)
      return(c(length(seed.genes),length(dataSet$seed.proteins), nrow(nodeList), nrow(edgeList), length(ppi.comps), substats));
    }else{
      return(.set.nSet(dataSet));
    }
  }else{
    return(-1);
  }
}

SetNetBuildOpt <- function(opt){
  netbuild.opt <<- opt
}


# for a given graph, obtain the smallest subgraphs that contain
# all the seed nodes. This is acheived by iteratively remove 
# the marginal nodes (degree = 1) that are not in the seeds

DecomposeGraph <- function(gObj, minNodeNum = 3){
  # now decompose to individual connected subnetworks
  if(netbuild.opt == "zero"){
    minNodeNum =2;
  }
  if(uploadedGraph == "false"){
    comps <-decompose.graph(gObj, min.vertices=minNodeNum);
  }else{
    if(gsize(gObj)>0 || met.type != "keggp"){ # do not decompose if select kegg projection
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
  net.stats <- ComputeSubnetStats(comps);
  ord.inx <- order(net.stats[,1], decreasing=TRUE);
  net.stats <- net.stats[ord.inx,];
  comps <- comps[ord.inx];
  names(comps) <- rownames(net.stats) <- paste("subnetwork", 1:length(comps), sep="");
  
  hit.inx <- net.stats$Node >= minNodeNum;
  
  comps <- comps[hit.inx];
  sub.stats <- unlist(lapply(comps, vcount));
  # now record
  ppi.comps <<- comps;
  net.stats <<- net.stats;
  
  return(sub.stats);
}

ComputeSubnetStats <- function(comps){
  net.stats <- as.data.frame(matrix(0, ncol = 3, nrow = length(comps)));
  colnames(net.stats) <- c("Node", "Edge", "Query");
  for(i in 1:length(comps)){
    g <- comps[[i]];
    if("mic" %in% dataSet$type){
      net.stats[i,] <- c(vcount(g),ecount(g),sum((unique(net.info$met.ids)) %in% V(g)$name));
    }else{
      net.stats[i,] <- c(vcount(g),ecount(g),sum((unique(dataSet$seeds.proteins)) %in% V(g)$name));
    }
  }
  return(net.stats);
}

UpdateSubnetStats <- function(){
  old.nms <- names(ppi.comps);
  net.stats <- ComputeSubnetStats(ppi.comps);
  ord.inx <- order(net.stats[,1], decreasing=TRUE);
  net.stats <- net.stats[ord.inx,];
  rownames(net.stats) <- old.nms[ord.inx];
  net.stats <<- net.stats;
}


# support walktrap, infomap and lab propagation
FindCommunities <- function(method="infomap", use.weight=FALSE){
  
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
  
  if(length(fc) == 0 || modularity(fc) == 0){
    return ("NA||No communities were detected!");
  }
  
  # only get communities
  communities <- communities(fc);
  community.vec <- vector(mode="character", length=length(communities));
  com.vec <- list()
  gene.community <- NULL;
  qnum.vec <- NULL;
  pval.vec <- NULL;
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
    
    # calculate p values (comparing in- out- degrees)
    #subgraph <- induced.subgraph(g, path.inx);
    subgraph <- induced.subgraph(g, path.ids);
    in.degrees <- degree(subgraph);
    #out.degrees <- degree(g, path.inx) - in.degrees;
    out.degrees <- degree(g, path.ids) - in.degrees;
    ppval <- wilcox.test(in.degrees, out.degrees)$p.value;
    ppval <- signif(ppval, 3);
    pval.vec <- c(pval.vec, ppval);
    
    # calculate community score
    community.vec[rowcount] <- paste(c(psize, qnums, ppval, pids), collapse=";");
    com.vec[[rowcount]] <- path.ids
  }
  
  ord.inx <- order(pval.vec, decreasing=F);
  community.vec <- community.vec[ord.inx];
  qnum.vec <- qnum.vec[ord.inx];
  ord.inx <- order(qnum.vec, decreasing=T);
  community.vec <- community.vec[ord.inx];
  if(1==2){
    metList <- rownames(dataSet$exp.mat[["met"]]);
    hits.query <- lapply(com.vec, 
                         function(x) {
                           x[x %in% metList];
                         }
    );
    toRemove <- vector()
    inxs <- vector();
    for(i in 1:length(hits.query)){
      if(length(hits.query[[i]]) < 10){
        toRemove = c(toRemove, com.vec[[i]])
        inxs = c(inxs, i)
      }
    }
    com.vec = com.vec[inxs]
    community.vec= community.vec[inxs]
    #g <- ppi.comps[[current.net.nm]];
    #g <-delete_vertices(g, unlist(toRemove))
    #ppi.comps[[current.net.nm]] <<- g
  }
  
  all.communities <- paste(community.vec, collapse="||");
  colnames(gene.community) <- c("Id", "Label", "Module");
  fast.write.csv(gene.community, file="module_table.csv", row.names=F);
  return(all.communities);
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
  library(RJSONIO);
  netData <- list(nodes=nodes);
  sink(filenm);
  cat(toJSON(netData));
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
  library(RJSONIO);
  netData <- list(nodes=nodes);
  sink(filenm);
  cat(toJSON(netData));
  sink();
  return(filenm);
}

PerformLayOut <- function(net.nm, algo, focus){
  require(igraph)
  
  g <- ppi.comps[[net.nm]];
  vc <- vcount(g);
  if(algo == "Default"){
    if(vc > 3000) {
      pos.xy <- layout.lgl(g, maxiter = 100);
    }else if(vc > 2000) {
      pos.xy <- layout.lgl(g, maxiter = 150);
    }else if(vc > 1000) {
      pos.xy <- layout.lgl(g, maxiter = 200);
    }else if(vc < 150){
      pos.xy <- layout.kamada.kawai(g);
    }else{
      pos.xy <- layout.fruchterman.reingold(g);
    }
  }else if(algo == "FrR"){
    pos.xy <- layout.fruchterman.reingold(g);
  }else if(algo == "random"){
    pos.xy <- layout.random(g);
  }else if(algo == "lgl"){
    if(vc > 3000) {
      pos.xy <- layout.lgl(g, maxiter = 100);
    }else if(vc > 2000) {
      pos.xy <- layout.lgl(g, maxiter = 150);
    }else {
      pos.xy <- layout.lgl(g, maxiter = 200);
    }
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
    pos.xy <- layout.graphopt(g, niter=maxiter);
  }else if(algo == "fr"){
    pos.xy <- layout_with_fr(g, dim=3, niter=500)
  }else if(algo == "kk"){
    pos.xy <- layout_with_kk(g, dim=3, maxiter=500)
  }else if(algo == "tree"){
    l <- layout_with_sugiyama(g, vgap=vc/4)
    pos.xy <- -l$layout
  }else if(algo == "circular_tripartite"){
    library(ggforce)
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
    library(graphlayouts)
    # the fist element in the list for concentric is the central node.
    if(focus==""){
      inx=1;
    }else{
      inx = which(V(g)$name == focus)
    }
    coords <- layout_with_focus(g,inx)
    pos.xy <- coords$xy
  }else if(algo == "backbone"){
    library(graphlayouts)
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
  output <- XiaLabCppLib::call_sr(from,to,cost,node_names,node_prizes)
  
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
    E(subnet)$weight <- as.numeric(output[[3]])
    subnet <- delete_vertices(subnet, "DUMMY")
    subnet <- delete_vertices(subnet, names(which(degree(subnet)==0)));
    return(subnet);
  } else{
    print("Subnetwork can not be identified for a given parameter set")
    return(NULL);
  }
}

convertIgraph2JSON <- function(dataSetObj=NA, net.nm, filenm, thera="FALSE", dim=3){
  dataSet <- .get.nSet(dataSetObj);
  if(!exists("my.convert.igraph")){ # public web on same user dir
    compiler::loadcmp("../../rscripts/omicsnetr/_utils_convertIgraph2JSON.Rc");    
  }
  return(my.convert.igraph(dataSet, net.nm, filenm, thera, dim));
}

initNet <- function(){
  nodeu.type <<- vector()
  edgeu.res <<- data.frame()
  net.info <<- list();
  nodeu.ids <<- vector()
  nodeu.nms <<- vector()
  expr.vec <<- vector()
  expand.types <- vector();
  result.listu <<- list();
  seed.genes <<- vector();
  ppi.comps <<- list();
  overall.graph <<- "";
  net.stats <<- "";
  curr.types <<- vector()
  queried.types <<- vector();
  edgeu.res.list <<- list()
}

QueryNet <- function(dataSetObj=NA, type="gene", seedOnly=T, require.exp, min.score){
  dataSet <- .get.nSet(dataSetObj);  
  containMsg <- 0;
  snpPeakMicBool <- F;
  build.opt <<- "all";
  
  edge.res <- data.frame();
  
  if("mic" %in% curr.types && type == "mic"){
    initNet();
  }
  
  if(length(edgeu.res.list) > 0){
    for(i in 1:length(edgeu.res.list)){
      edge.res <- rbind(edge.res, edgeu.res.list[[i]]);
    }
    nodes.query <- c(edge.res[,1], edge.res[,2])
  }
  
  old.edges <- edge.res;
  query.mat <- dataSet$exp.mat[[type]];
  
  if (length(query.mat) == 0) {
    types.vec <- c("mic", "peak", "snp");
    for(i in 1:length(types.vec)){
      typeVal <- types.vec[i];
      if(length(dataSet$exp.mat[[typeVal]]) > 0){
        snpPeakMicBool <- T;
      }
    }
    inv = "inverse"
    if(seedOnly && !snpPeakMicBool && length(dataSet$exp.mat[["gene"]])>0){
      query.vec <- unique(seed.genes)
    }else{
      query.vec <- unique(nodeu.ids)
    }
    
    if(any(types.vec %in% names(dataSet$exp.mat))){
      if(type %in% c("mir", "tf", "met")){
        inv = "inverse"
      }else{
        inv = "direct"
      }
    }
    
    if(length(query.vec) == 0 && length(dataSet$exp.mat[["gene"]])>0){
      query.vec <- rownames(dataSet$exp.mat[["gene"]])
      inv = "inverse";
      seed.genes <<- c(seed.genes, query.vec); 
    }else if(length(query.vec) == 0 && length(dataSet$exp.mat[["ko"]])>0){
      query.vec <- rownames(dataSet$exp.mat[["ko"]])
      inv = "inverse";
      seed.genes <<- c(seed.genes, query.vec); 
    }
    
  } else {
    
    query.vec <- rownames(dataSet$exp.mat[[type]]);
    inv <- "direct";
    seed.genes <<- c(seed.genes, query.vec); 
  }
  
  result.listu$protein.vec <- unique(as.character(query.vec));
  result.listu$type <- type;
  result.listu <<- result.listu;
  
  if(!CheckQueryTypeMatch(result.listu$protein.vec, type)){
    current.msg<<- paste("Please make sure correct interaction type is selected!")
    containMsg <- 1;
    print("Please make sure correct interaction type is selected!");
  }
  
  if(type == "m2m" && exists("PeakSet")) {
    result.listu$protein.vec <- PeakSet[["nodes.df"]][["Id"]];
    result.listu$type <- "m2m";
    result.listu <<- result.listu;
  }
  
  if("mic" %in% curr.types && type %in% c("met")){ #m2p after mic is direct
    inv = "direct";
  }
  
  SearchNetDBX(type, "ppi", require.exp, min.score, inv, FALSE);
  
  if(type == "mic"){
    result.listu$protein.vec <- unique(as.character(net.info$met.ids));   
    result.listu$type <- "met";
    result.listu <<- result.listu;
    data.org <<- "microbiome";
    m2m.type <<- "mic";
    SearchNetDBX("m2m.mic", "ppi", TRUE, 900, "direct", FALSE);
  }else if (type == "met"){ ## if both KO and met uploaded
    if(length(dataSet$exp.mat[["ko"]])>0){
      result.listu$protein.vec <- rownames(dataSet$exp.mat[["ko"]]);   
      result.listu$type <- "met"
      result.listu <<- result.listu;
      seed.genes <<- c(seed.genes, result.listu$protein.vec); 
      SearchNetDBX("met", "ppi", TRUE, 900, "inverse", FALSE);
    }else if(length(dataSet$exp.mat[["gene"]])>0){
      result.listu$protein.vec <- rownames(dataSet$exp.mat[["gene"]]);   
      result.listu$type <- "met"
      result.listu <<- result.listu;
      seed.genes <<- c(seed.genes, result.listu$protein.vec); 
      SearchNetDBX("met", "ppi", TRUE, 900, "inverse", FALSE);
    }
  }
  
  node.res <- data.frame(Id=nodeu.ids, Label=nodeu.nms);
  
  edge.res <- data.frame();
  if(length(edgeu.res.list) > 0){
    for(i in 1:length(edgeu.res.list)){
      edge.res <- rbind(edge.res, edgeu.res.list[[i]]);
    }
  }
  
  if(type == "m2m") {
    # to remove the useless expansion
    resclean <- redundancyClean(node.res, edge.res)
    node.res <- resclean$node.res
    edge.res <- resclean$edge.res
  }
  
  if(length(edge.res)>0 && length(old.edges)>0){
    if(dim(edge.res) == dim(old.edges)){
      current.msg<<- paste("No interactions have been detected!")
      containMsg <- 1;
      return(c(0, 0, containMsg));
    }else{
      curr.types <<- unique(c(curr.types, type));
      containMsg <- 0;
    }
  }else{
    curr.types <<- unique(c(curr.types, type));
    containMsg <- 0;
  }
  
  omics.net <<- list(
    db.type="abc",
    order=1, 
    seeds=" ", 
    table.nm=" ", 
    node.data = node.res, 
    edge.data = edge.res,
    require.exp = F,
    min.score = 900
  );
  #print(curr.types);
  #print(c(nrow(node.res), nrow(edgeu.res), containMsg));
  dataSet$seeds.proteins <- unique(seed.genes);
  dataSet <<- dataSet
  
  if(.on.public.web){
    .set.nSet(dataSet);
    return(c(nrow(node.res), nrow(edge.res), containMsg));
    }else{
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
  
  PlotDegreeHistogram <- function(imgNm, netNm = "NA", dpi=72, format="png"){
    library(Cairo)
    dpi<-as.numeric(dpi)
    imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
    Cairo(file=imgNm, width=400, height=400, type="png", bg="white");
    library(ggplot2)
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
      theme(plot.title = element_text(hjust = 0.5))
    print(p)
    dev.off();
  }
  
  PlotBetweennessHistogram <- function(imgNm, netNm = "NA",dpi=72, format="png"){
    library(Cairo)
    dpi<-as.numeric(dpi)
    imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
    Cairo(file=imgNm, width=400, height=400, type="png", bg="white");
    library(ggplot2)
    if(netNm != "NA"){
      overall.graph <- ppi.comps[[netNm]];
    }
    G.degrees <- betweenness(overall.graph)
    
    G.degree.histogram <- as.data.frame(table(G.degrees))
    G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])
    
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
      theme(plot.title = element_text(hjust = 0.5))
    print(p)
    dev.off();
  }
  
  PreparePeaksNetwork <- function(){
    dataSet <- .get.nSet(dataSetObj);
    if(exists("PeakSet",envir = .GlobalEnv)) {
      PeakSet <- get("PeakSet", envir = .GlobalEnv)
    } else {
      return(0)
    }
    
    nodes <- PeakSet$nodes
    edges <- PeakSet$edges
    met.ids <- nodes$KEGGID
    met.ids[met.ids == ""] <- nodes$formula[met.ids == ""]
    PeakSet$nodes$KEGGID <- met.ids;
    art.ids <- PeakSet$nodes$name[which(PeakSet$nodes$class == "Artifact")]
    
    PeakSet$mets <- unique(PeakSet$nodes$KEGGID[which(PeakSet$nodes$class == "Metabolite")]);
    PeakSet$put.mets <- PeakSet$nodes$KEGGID[which(PeakSet$nodes$class == "Putative metabolite")];
    PeakSet$mets.ids <- PeakSet$nodes$name[which(PeakSet$nodes$class == "Metabolite")];
    PeakSet$put.mets.ids <- PeakSet$nodes$name[which(PeakSet$nodes$class == "Putative metabolite")];
    dataSet$seed[["peak"]] <- data.frame(rep(0, length(PeakSet$mets)));
    
    rownames(dataSet$seed[["peak"]]) <- PeakSet$mets;
    
    dataSet$seed[["peakRaw"]] <- data.frame(rep(0, length(PeakSet$NetID_output$peak_id)))
    
    rownames(dataSet$seed[["peakRaw"]]) <- PeakSet$NetID_output$peak_id
    
    cmp.nms <- PrepareInputList(dataSet, dataSet$seed[["peak"]], data.org, "peak", "kegg", "direct");
    
    PeakSet$nodes.df <- data.frame(Id=met.ids, Label=doKegg2NameMapping(met.ids))
    edges.df0 <- data.frame(from=edges[,1], to=edges[,2])
    edges.df01 <- edges.df0[!edges.df0[,1] %in% art.ids,];
    edges.df02 <- edges.df01[!edges.df01[,2] %in% art.ids,];
    conv1 <- data.frame(from=PeakSet$nodes$name, Source=PeakSet$nodes$KEGGID);
    conv2 <- data.frame(to=PeakSet$nodes$name, Target=PeakSet$nodes$KEGGID);
    edges.df1 <- merge(edges.df02, conv1, by = "from");
    edges.df2 <- merge(edges.df1, conv2, by = "to");
    PeakSet$edges.df <- data.frame(Source=edges.df2$Source, Target=edges.df2$Target);
    PeakSet <<- PeakSet;
    anchor_type <<- "peak";
    met.type <<- "kegg";
    
    dataSet <<- dataSet;
    return(cmp.nms)
  }
  
  SetMicType <- function(name){
    dataSet$mic.type <<- name;
  }
  
  SetMicExcludeOpt <- function(opt){
    dataSet$mic.exclude.opt <<- opt;
  }
  
  SetMicThresh <- function(thresh){
    dataSet$mic.thresh <<- thresh
  }
  
  # support rwr right now
  DoGba <- function(fileNm="NA", method="rwr", nodeids){
    library(RandomWalkRestartMH)
    library(igraph)
    
    nodes <- strsplit(nodeids, ",")[[1]];
    g <- ppi.comps[[current.net.nm]];
    PPI_MultiplexObject <- create.multiplex(list(PPI=g))
    
    AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
    AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)
    seeds <- nodes;
    
    RWR_PPI_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI,
                                                     PPI_MultiplexObject,seeds);
    resObj <- RWR_PPI_Results$RWRM_Results
    if(nrow(resObj) > 100){
      resObj <- resObj[c(1:100),]
    }
    require(RJSONIO);
    sink(fileNm);
    cat(toJSON(resObj));
    sink();
    
    csvNm <- paste0(gsub(".json", "", fileNm), ".csv");
    fast.write.csv(RWR_PPI_Results$RWRM_Results, file=csvNm);
    return(1)
  }
  
  
  FilterByTissue <- function(type, tissue){
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
    
    # update omics.net
    omics.net$order <- 1;
    omics.net$node.data <- nodeList;
    omics.net$edge.data <- edgeList;
    omics.net <<- omics.net;
    overall.graph <- g
    substats <- DecomposeGraph(overall.graph);
    if(!is.null(substats)){
      overall.graph <<- overall.graph;
      return(c(length(seed.genes),length(seed.proteins), vcount(overall.graph), ecount(overall.graph), length(ppi.comps), substats));
    }else{
      return(0);
    }
  }
  
  queryFilterDB <- function(type, org){
    require('RSQLite');
    conv.db <- dbConnect(SQLite(), paste(sqlite.path, "tissue_filter.sqlite", sep="")); 
    db.map <- dbReadTable(conv.db, paste0(data.org,"_",type));
    dbDisconnect(conv.db); cleanMem();
    
    return(db.map)
  }
  
