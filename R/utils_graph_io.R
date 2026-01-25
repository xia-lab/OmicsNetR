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

my.convert.Igraph2JSONFromFile <- function(net.nm, filenm, dim=3){
  dim = as.numeric(dim);
  g <- ppi.comps[[net.nm]];

  # annotation
  nms <- V(g)$name;
  lbls = nms;
  if(length(V(g)$Label) >0){
    lbls <- V(g)$Label;
  }else if(length(dataSet$nodeLabels) > 0){
    conv <- dataSet$nodeLabels;
    hit.inx <- match(nms, conv$id2);
    lbls <- conv[hit.inx, "name2"];
    mode(lbls) <- "character";
    na.inx <- is.na(lbls);
    lbls[na.inx] <- nms[na.inx];
  }
  #print(lbls);
  # setup shape (gene circle, other squares)
  shapes <- rep("circle", length(nms));

  # get edge data
  edge.mat <- get.edgelist(g);

  edge.mat = data.frame(edge.mat);
  edge.mat$color = "target";
  # edge.mat1 = as.matrix(edge.mat1);

  edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2], true_col = edge.mat[,3]);
  edge.exp <- get.edge.attribute(g, name="Expression", index = E(g));
  edge.pval <- get.edge.attribute(g, name="Pval", index = E(g));


  if(length(edge.pval) > 0){
    edge.sizes <- as.numeric(rescale2NewRange(abs(edge.pval), 0.3, 6));
    edge.mat <- cbind(edge.mat, size=edge.sizes);
    edge.label <-get.edge.attribute(g, name="Label", index = E(g));
    edge.mat <- cbind(edge.mat, label=edge.label);
  }

  if(length(edge.exp) > 0){
    bad.inx <- is.na(edge.exp) | edge.exp==0;
    edge.col <- ComputeColorGradient(edge.exp, "black", T);
    edge.col[bad.inx] <- "#d3d3d3";
    edge.mat[,4] <- edge.col;
  }
  # get the node data
  if(length(E(g)$weight)>0){
    E(g)$weight = abs(E(g)$weight);
    weight = E(g)$weight;
    ppi.comps[[net.nm]] <<- g;
  }
  node.btw <- as.numeric(betweenness(g));
  node.clo <- as.numeric(closeness(g));
  #node.adh <- as.numeric(adhesion(g));
  node.eig <- eigen_centrality(g);
  node.eig = as.numeric(node.eig$vector);

  node.dgr <- as.numeric(degree(g));
  node.exp <- as.numeric(get.vertex.attribute(g, name="Expression", index = V(g)));

  if(length(node.exp) == 0){
    node.exp <- rep(0,length(node.dgr));
  }

  node.type <- get.vertex.attribute(g, name="Type", index = V(g));

  if(length(node.type) == 0){
    node.type <- rep("Nodes",length(node.dgr));
  }

  # node size to degree values
  if(vcount(g)>500){
    min.size = 1;
  }else if(vcount(g)>200){
    min.size = 2;
  }else{
    min.size = 3;
  }

  minval = min(node.dgr, na.rm=T);
  maxval = max(node.dgr, na.rm=T);
  result = maxval-minval;

  if(result == 0){
    node.sizes <- rep((log(node.dgr))^2, length(nms));
    node.sizes2d <- node.sizes;
  }else{
    node.sizes <- as.numeric(rescale2NewRange((log(node.dgr))^2, min.size, 9))*3 +5;
    node.sizes2d <- as.numeric(rescale2NewRange((log(node.btw+1))^2, min.size, 12));
  }

  nsize <- get.vertex.attribute(g, name="Pval", index = V(g));
  if(length(nsize) > 0){
    nsize[is.na(nsize)] <- 0;
    node.sizes <- as.numeric(rescale2NewRange(abs(nsize), min.size, 9))*3 +5;
    node.sizes2d <- as.numeric(rescale2NewRange(abs(nsize), min.size, 12));
  }

  centered = T;
  notcentered = F;
  # update node color based on betweenness
  require("RColorBrewer");
  topo.val <- log(node.btw+1);
  topo.colsb <- ComputeColorGradient(topo.val, "black", notcentered);
  topo.colsw <-  ComputeColorGradient(topo.val, "white", notcentered);

  # color based on expression
  bad.inx <- is.na(node.exp) | node.exp==0;
  if(!all(bad.inx)){
    exp.val <- node.exp;
    node.colsb.exp <- ComputeColorGradient(exp.val, "black", centered);
    node.colsw.exp <- ComputeColorGradient(exp.val, "white", centered);
    node.colsb.exp[bad.inx] <- "#d3d3d3";
    node.colsw.exp[bad.inx] <- "#c6c6c6";
  }else{
    node.colsb.exp <- rep("#d3d3d3",length(node.exp));
    node.colsw.exp <- rep("#c6c6c6",length(node.exp));
  }

  node_attr = list.vertex.attributes(g);

  attr=list();
  for(j in 1:length(node_attr)){
    attr[[node_attr[j]]] = vertex_attr(g, node_attr[j]);
  }
  attr_names <- names(attr);
  attr_nd <- list();
  arr <- list();
  for(i in 1:length(node.sizes)){
    for(j in 1:length(attr)){
      attr_nd[node_attr[j]] = as.character(unlist(attr[node_attr[j]])[i])
    }
    arr[[i]] = attr_nd;
  }
  network_prop <- list();
  for(i in 1:length(node.sizes)){
    network_prop[[i]]  <- list(
      closeness = node.clo[i],
      eigen = node.eig[i]
    )
  }

  lblsu <<- lbls;
  pos.xy <- PerformLayOut(net.nm, "Default");

  if(length(V(g)$x) >0){
    pos.xy[,1] <-V(g)$x;
    pos.xy[,2] <-V(g)$y;
  }

  # now create the json object
  nodes <- vector(mode="list");

  require("stringr")
  displayedLabel<-lbls;
  long.inx <- which(stringr:::str_length(lbls) > 32);
  displayedLabel[long.inx] <- paste0(strtrim(lbls[long.inx],  rep(32, length(lbls[long.inx]))), "..." )

  for(i in 1:length(node.sizes)){
    nodes[[i]] <- list(
      id=nms[i],
      label=lbls[i],
      displayedLabel=displayedLabel[i],
      size=node.sizes[i],
      size2d=node.sizes2d[i],
      type="circle",
      types=node.type[i],
      molType=node.type[i],
      x2d = pos.xy[i,1],
      y2d = pos.xy[i,2],
      colorb=topo.colsb[i],
      colorw=topo.colsw[i],
      topocolb=topo.colsb[i],
      topocolw=topo.colsw[i],
      expcolb=node.colsb.exp[i],
      expcolw=node.colsw.exp[i],
      user=network_prop[[i]],
      attributes=list(
        expr = node.exp[i],
        degree=node.dgr[i],
        between=node.btw[i]
      )
    );
  }

  if(length(edge.exp) == 0){
    edge.exp <- rep("Nodes",length(E(g)$name));
  }

  # save node table
  nd.tbl <- data.frame(Id=nms, Degree=node.dgr, Betweenness=round(node.btw,2));
  # order
  ord.inx <- order(nd.tbl[,2], nd.tbl[,3], decreasing =TRUE);
  nd.tbl <- nd.tbl[ord.inx, ];
  fast.write.csv(nd.tbl, file="node_table.csv", row.names=FALSE);

  # covert to json
  require("RJSONIO");
  dg <- decompose(g)
  if(length(dg)>1){
    modules <- "NA";
  }else{
    modules <- FindCommunities("infomap", FALSE);
  }
  netData <- list(nodes=nodes, edges=edge.mat, modules=modules, conv="NA", gene.nms=c("NA"), prot.nms=c("NA"));
  if(exists("netPropU")){
    netData[["netProp"]] <- netPropU;
  }
  partialToBeSaved <<- c(partialToBeSaved, c(filenm));
  sink(filenm);
  cat(RJSONIO::toJSON(netData));
  sink();
}
