##################################################
## R scripts for OmicsNet
## Description: Data IO functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

ProcessOmicsNetJson <- function(dataSet, fileName) {
  require("RJSONIO");
  obj <- RJSONIO::fromJSON(fileName)
  seed <- obj$dataSet$seed
  for(i in 1:length(seed)){
    df <- seed[[i]];
    df1 <- data.frame(expr = df$expr, id = df$id)
    rownames(df1) <- df1$id;
    seed[[i]] <- df1;
  }
  dataSet$seed <- seed;
  dataSet$seeds.expr <- obj$dataSet$seeds.expr
  dataSet$seeds.proteins <- obj$dataSet$seeds.proteins
  seed.expr <<- obj$seed.expr
  omics.net <- obj$omics.net
  if(is.null(dim(obj$omics.net$node.data))){
  node.data.df <- data.frame(Id=omics.net$node.data, Label=omics.net$node.data);
  }else{
  node.data.df <- data.frame(Id=omics.net$node.data$Id, Label=omics.net$node.data$Label);
  }
  node.data.df <- unique(node.data.df);
  edge.data.df <- data.frame(Source=omics.net$edge.data$Source, Target=omics.net$edge.data$Target);
  omics.net$node.data <- node.data.df
  omics.net$edge.data <- edge.data.df
  seed.proteins <<- node.data.df[,1];
  seed.genes <<- node.data.df[,1];
  omics.net <<- omics.net;
  net.info <<- obj$net.info
  data.org <<- obj$data.org

  if(is.null(obj$snpRegion)){
    dataSet$snpRegion <- F;
  }else{
    dataSet$snpRegion <- obj$snpRegion;
  }

  if(is.null(obj$allOrgs)){
    dataSet$allOrgs <- data.org;
  }else{
    dataSet$allOrgs <- obj$allOrgs;
  }

  dataSet <<- dataSet;
  CreateGraph(dataSet);
  net.nm <- names(ppi.comps)[1];
  net.nmu <<- net.nm;
  current.net.nm <<- net.nm;
  g <- ppi.comps[[net.nm]];

  dataSet <- .computeSubnetStats(dataSet, ppi.comps);
  return(.set.nSet(dataSet));
}

#' ReadGraphFile
#'
#' @param dataSetObj dataSetObj
#' @param fileName fileName, Input name of graph file
#' @param fileType fileType, File type of graph file (".json", ".graphml", ".txt", ".sif")
#' @export
ReadGraphFile <- function(dataSetObj=NA, fileName, fileType) {
  if(!exists("my.read.graphfile")){ # public web on same user dir
    compiler::loadcmp("../../rscripts/OmicsNetR/R/utils_graph_io.Rc");
  }
  return(my.read.graphfile(dataSetObj, fileName, fileType));
}

convertIgraph2JSONFromFile <- function(net.nm, filenm, dim=3){

  if(!exists("my.convert.Igraph2JSONFromFile")){ # public web on same user dir
    compiler::loadcmp("../../rscripts/OmicsNetR/R/utils_graph_io.Rc");
  }
  return(my.convert.Igraph2JSONFromFile(net.nm, filenm, dim));
}

getGraphStatsFromFile <- function(){
  g <- ppi.comps[[net.nmu]];
  nms <- V(g)$name;
  edge.mat <- get.edgelist(g);
  return(c(length(nms), nrow(edge.mat)));
}
