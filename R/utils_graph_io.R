my.read.graphfile <- function(dataSetObj=NA, fileName, fileType) {
  dataSet <- .get.nSet(dataSetObj);
  require("igraph");
  types_arr <<- "";
  lbls <- "";
  dataSet$fileType <- fileType;
  current.msg <<- NULL;

  if(grepl("scatter", fileName)){
    require("RJSONIO")
    j = fromJSON(fileName)
    nm = "scatter3D.json"
    sink(nm);
    cat(toJSON(j));
    sink();
    return(1);
  }

  if(fileType == "jsonOmicsnet"){
    res <- ProcessOmicsNetJson(dataSet, fileName);
    return(res);
  }else if(fileType == "jsonGraph"){
    return(.set.nSet(dataSet));
  }else if(fileType == "graphml"){
    graphX = tryCatch({
      read_graph(fileName, format = "graphml")
    }, warning = function(w) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0);
    }, error = function(e) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0);
    }, finally = {

    })
  }else if(fileType == "sif"){
    graphX = tryCatch({
      read.sif(fileName, format="igraph", directed = FALSE)
    }, warning = function(w) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0);
    }, error = function(e) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0);
    }, finally = {

    })
  }else if(fileType == "txt"){
    df <- read.table(fileName, header=FALSE, stringsAsFactors = FALSE)
    df <- as.matrix(df);
    graphX = tryCatch({
      graph_from_edgelist(df)
    }, warning = function(w) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0);
    }, error = function(e) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0);
    }, finally = {

    })


  }else if(fileType == "json"){
    require("RJSONIO");
    dat <- fromJSON(fileName);
    if(!is.null(dat$data[["__Org"]])){
      data.org <<- dat$data[["__Org"]];
    }
    dfn <- unlist(dat$elements$nodes);
    conv <- data.frame(id1=dfn[which(names(dfn)=='data.id')], name1=dfn[which(names(dfn)=='data.name')]);
    dfe <- unlist(dat$elements$edges);
    dffe <- data.frame(id1=dfe[which(names(dfe) == "data.source")], id2=dfe[which(names(dfe) == "data.target")]);
    dfint <- merge(conv, dffe, by="id1");
    colnames(conv) <- c("id2", "name2");
    df <- merge(conv, dfint, by="id2");
    df <- df[,c("id1", "id2")];
    df <- as.matrix(df);

    graphX = tryCatch({
      graph_from_edgelist(df, directed=FALSE);
    }, warning = function(w) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0);
    }, error = function(e) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0);
    }, finally = {

    })
    nms <- V(graphX)$name
    hit.inx <- match(nms, conv$id2);
    lbls <- conv[hit.inx, "name2"];
    mode(lbls) <- "character";
    na.inx <- is.na(lbls);
    lbls[na.inx] <- nms[na.inx];

  }else {
    current.msg <<- "Unknown format, please make sure that the file is saved in the supported formats!";
    return(0);
  }

  if(!is_igraph(graphX)){
    current.msg <<- "Failed to parse your file, please make sure that the file is formatted correctly";
    return(0);
  }
  current.msg <<- "Sucessfully parsed your graph file!";

  nms <- V(graphX)$name;
  if(length(nms)<1){
    nms <- V(graphX)$id;
    graphX <- set_vertex_attr(graphX, "name", value=nms);
  }

  if(length(lbls) == 0){
    lbls <- nms;
  }

  node.data <- data.frame(nms, lbls);
  seed.proteins <<- nms;
  seed.genes <<- seed.proteins;
  e <- get.edgelist(graphX)
  edge.data <- data.frame(Source=e[,1], Target=e[,2])
  seed.expr <<- rep(0, length(node.data));
  dataSet$seed <- list();
  dataSet$seeds.proteins <- c();
  dataSet <<- dataSet;
  omics.net <<- list(
    netw.type="graph",
    order=1,
    seeds=nms,
    table.nm=" ",
    node.data = node.data,
    edge.data = edge.data
  );

  CreateGraph(dataSet);
  net.nm <- names(ppi.comps)[1];
  net.nmu <<- net.nm;
  current.net.nm <<- net.nm;
  g <- ppi.comps[[net.nm]];

  if(fileType == "json"){
    dataSet$nodeLabels <- conv
  }
  dataSet <- .computeSubnetStats(dataSet, ppi.comps);
  return(.set.nSet(dataSet));
}


read.sif <- function (sif.file, format = "graphNEL", directed = FALSE, header = TRUE, sep = "\t", ...) {

  net <- read.csv(file = sif.file, sep = sep, colClasses = "character", header = header, ...)

  # Assume form: node1 linktype node2 side.info..
  if ( ncol(net) > 2 ) {

    # remove NA nodes
    nas <- apply(net, 1, function (x) {any(is.na(x[c(1,3)]))})
    if (any(nas)) {
      net <- net[!nas, ];
      warning("NAs removed from network node list, ", sum(nas), " edges removed.");
    }

    net <- graph.edgelist(as.matrix(net[, -2]), directed = directed);

  } else if ( ncol(net) == 2 ) { # assume form: node1 node2

    # remove NA nodes
    nas <- apply(net, 1, function (x) {any(is.na(x))});
    if (any(nas)) {
      net <- net[!nas, ];
      warning("NAs removed from network node list, ", sum(nas), " edges removed.");
    }

    net <- graph.edgelist(cbind(net[,1],net[,2]), directed = directed);
  }

  if (format == "graphNEL") { net <- igraph.to.graphNEL(net) }
  # if (format == "igraph") { net <- igraph.from.graphNEL(igraph.to.graphNEL(net)) }

  net
}

