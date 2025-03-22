
#color shape layout important
my.convert.igraph <- function(dataSetObj=NA, net.nm, filenm, thera=FALSE, dim=3){
  #save.image("igraph.RData");
  if(exists("lbls",envir = .GlobalEnv)) {
    lbls <- get("lbls", envir = .GlobalEnv)
  } else {
    lbls <<- NA;
  }

  dataSet <- .get.nSet(dataSetObj);
  dim <- as.numeric(dim)
  types_arr <<- dataSet$type;

  if(thera){
    g <- my.ppi;
  }else{
    g <- ppi.comps[[net.nm]];
  }
  current.net.nm <<- net.nm;
  modules = "NA"

  nms <- V(g)$name;
  #save(omics.net,g, file = "omics.net___checking.rda")

    hit.inx <- match(nms, omics.net$node.data[,1]);
  lbls <- omics.net$node.data[hit.inx, 2];

  # setup shape (gene circle, other squares)
  if("peak" %in% dataSet$type){
    shapes <- rep("putative_met", length(nms));
  }else{
    shapes <- rep("gene", length(nms));
  }
  itypes <- rep("circle", length(nms));

  # get edge data
  edge.mat <- get.edgelist(g);
  edge.mat1 = data.frame(edge.mat)
  edge.mat1$color = "target"
  edge.mat1 = as.matrix(edge.mat1)

  edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2], color = edge.mat1[,3]);

  # now get coords
  pos.xy <- PerformLayOut(net.nm, "Default");
  # get the note data
  node.btw <- as.numeric(betweenness(g));
  node.clo <- as.numeric(closeness(g));
  node.eig <- eigen_centrality(g);
  node.eig = as.numeric(node.eig$vector);

  node.dgr <- as.numeric(degree(g));
  exp <-as.vector(get.vertex.attribute(g, name="abundance", index = V(g)));
  exp[is.na(exp)] <- 0;
  node.exp <- as.numeric(exp);
  exp2 <- parse_expression_data(dataSet, exp,nms);

  # node size to degree values
  if(vcount(g)>500){
    min.size = 2;
  }else if(vcount(g)>200){
    min.size = 3;
  }else{
    min.size = 3;
  }

  node.sizes <- as.numeric(rescale2NewRange((log(node.btw+1))^2, min.size, 9))*3 +5;
  node.sizes2d <- as.numeric(rescale2NewRange((log(node.btw+1))^2, min.size, 12));
  centered <- T;
  notcentered <- F;

  require("RColorBrewer");
  topo.val <- log(node.btw+1);
  topo.colsb <- topo.colsb1 <- rep("#D3D3D3", length(nms));
  topo.colsw <- topo.colsw1 <- rep("#D3D3D3", length(nms));

  # color based on expression
  bad.inx <- is.na(node.exp) | node.exp==0;
  if(!all(bad.inx)){
    exp.val <- node.exp;
    node.colsb.exp <- ComputeColorGradient(exp.val, "black", centered);
    node.colsw.exp <- ComputeColorGradient(exp.val, "white", centered);
    node.cols.pie <- ComputeColorGradient(exp2, "black", centered);
    node.colsb.exp[bad.inx] <- "#d3d3d3";
    node.colsw.exp[bad.inx] <- "#c6c6c6";
  }else{
    node.colsb.exp <- rep("#D3D3D3",length(node.dgr));
    node.colsw.exp <- rep("#C6C6C6",length(node.dgr));
    node.cols.pie <- rep("#D3D3D3",length(node.dgr));

  }

  gene.nms <- rownames(dataSet$seed[["gene"]] )
  prot.nms <- rownames(dataSet$seed[["protein"]] )

  # now update for bipartite network

  #edge.mat1 <<- edge.mat
  gene.inx <- nms %in% c(net.info$gene.ids, net.info$ko.ids);

  if(length(net.info$gene.ids) == 0 && length(net.info$protein.ids) == 0 ){
    #gene.inx = predicted.inx;
    #predicted.inx = rep(FALSE,length(node.dgr));
  }

  kncmpd.inx <- nms %in% net.info$kncmpd.ids
  prot.inx <- nms %in% net.info$protein.ids;
  tf.inx <- nms %in% net.info$tf.ids;
  mir.inx <- nms %in% net.info$mir.ids;
  peak.inx <- nms %in% net.info$peak.ids;
  met.inx <- nms %in% net.info$met.ids;
  reg.inx <- nms %in% net.info$reg.ids;
  mic.inx <- nms %in% net.info$mic.ids;
  snp.inx <- nms %in% rownames(dataSet$seed[["snp"]]);


  containsGP <- any(gene.nms %in% prot.nms); ## check if any overlap between gene and protein


  topo.colsw[gene.inx] <- "#FF8484";
  topo.colsw[kncmpd.inx] <- "#00ffff";
  topo.colsw[prot.inx] <- "#FF8484";
  topo.colsw[tf.inx] <- "#39FF14";
  topo.colsw[mir.inx] <- "#00f6ff";
  topo.colsw[peak.inx] <- "#D3D3D3";

  topo.colsw[met.inx] <- "#ffff00";
  topo.colsw[reg.inx] <- "#ff9900";
  topo.colsw[mic.inx] <- "#39FF14";
  topo.colsw[snp.inx] <- "#4B0082";
  topo.colsw[nms %in% gene.nms] <- "#d3d3d3";

  topo.colsb <- topo.colsw;

  if(containsGP){
    inx <- nms %in% gene.nms
    topo.colsb[inx] <- "#D3D3D3";
    topo.colsw[inx] <- "#D3D3D3";
    color.vec <- c("#D3D3D3", "#FF8484", "#39FF14","#00f6ff", "#D3D3D3", "#00ffff", "#ffff00", "#4B0082", "#39FF14");
  }else{
    color.vec <- c("#FF8484", "#FF8484", "#39FF14","#00f6ff", "#D3D3D3", "#00ffff", "#ffff00", "#4B0082", "#39FF14");
  }
  names(color.vec) <- c("gene", "protein", "tf", "mir", "peak", "kncmpd", "met", "reg", "mic");

  #if(any(predicted.inx)){
  #  types_arr <<- c(types_arr, "interactor");
  #}
  if(any(gene.inx) && !"gene" %in% types_arr){
    types_arr <<- c(types_arr, "gene");
  }

  shapes[gene.inx] <- "gene";
  shapes[prot.inx] <- "protein";
  shapes[tf.inx] <- "tf";
  shapes[mir.inx] <- "mir";
  shapes[kncmpd.inx] <- "kncmpd";

  shapes[peak.inx] <- "putative_met";
  shapes[met.inx] <- "met";
  shapes[reg.inx] <- "reg"
  shapes[mic.inx] <- "microbe"
  shapes[snp.inx] <- "snp";

  types <- rep("", length(shapes))
  types[gene.inx] <- paste(types[gene.inx],"gene");
  types[peak.inx] <- "putative_met"
  types[prot.inx] <- paste(types[prot.inx],"protein");
  types[tf.inx] <- "tf"
  types[mir.inx] <- "mir";
  types[met.inx] <- "met";
  types[snp.inx] <- "snp";
  types[kncmpd.inx] <- "kncmpd";

  types <- trimws(types);

  network_prop = list();
  for(i in 1:length(node.sizes)){
    network_prop[[i]]  <- list(
      closeness = node.clo[i],
      eigen = node.eig[i]
    )
  }

  seed.inx <- nms %in% names(expr.vec);

  seed_arr <- rep("notSeed",length(node.dgr));
  seed_arr[seed.inx] <- "seed";

  type <- rep(FALSE,length(node.dgr))
  type[seed.inx] = TRUE

  lblsu <<- lbls;
  node_attr = list.vertex.attributes(g);
  node_attr = node_attr[which(node_attr!="names")]

  # now create the json object
  nodes <- vector(mode="list");


  require("stringr")
  displayedLabel<-lbls;
  long.inx <- which(stringr:::str_length(lbls) > 32);
  displayedLabel[long.inx] <- paste0(strtrim(lbls[long.inx],  rep(32, length(lbls[long.inx]))), "..." )

  moltypes <- shapes;
  moltypes[nms %in% gene.nms] <- "gene";
  moltypes[nms %in% prot.nms] <- "protein";
  #moltypes[(nms %in% prot.nms) & (nms %in% gene.nms)] <- "gene_protein";
  for(i in 1:length(node.sizes)){
    nodes[[i]] <- list(
      id=nms[i],
      idx=i,
      label=lbls[i],
      displayedLabel = displayedLabel[i],
      size=node.sizes[i],
      size2d=node.sizes2d[i],
      type=shapes[i],
      molType=moltypes[i],
      types=types[i],
      seedArr =seed_arr[i],
      color=topo.colsb[i],
      colorb=topo.colsb[i],
      colorw=topo.colsw[i],
      topocolb=topo.colsb[i],
      topocolw=topo.colsw[i],
      expcolb=node.colsb.exp[i],
      expcolw=node.colsw.exp[i],
      x2d = pos.xy[i,1],
      y2d = pos.xy[i,2],
      user =network_prop[[i]],
      attributes=list(
        degree=node.dgr[i],
        between=node.btw[i],
        expr = node.exp[i],
        closeness = node.clo[i],
        eigen = node.eig[i]
      )
    );
  }

    if(containsGP){
      nodes[[i]]$exp2 = exp2[[i]];
      nodes[[i]]$expcolpie = node.cols.pie[[i]];
    }

  middleLayerType <- unique(shapes)[1]

  for(i in 1:length(unique(shapes))){
    inx <- shapes == unique(shapes)[i]
    if(i == 1){
    sumDgrs <- sum(node.dgr[inx])
    }else{
    if(sumDgrs < sum(node.dgr[inx])){
    middleLayerType <-unique(shapes)[i]
    }
    }
  }

  V(g)$layers <- as.numeric(as.factor(moltypes));
  ppi.comps[[net.nm]] <<- g;

  summary = vector(mode="list");
  summary[[1]] = list(nodes=vcount(g));
  summary[[2]] = list(edges=ecount(g));

  # save node table


  nd.tbl <- data.frame(Id=nms, Label=lbls, Degree=node.dgr, Betweenness=round(node.btw,2), Expression = node.exp, Closeness=node.clo, Eigenvalue=node.eig);
  # order
  ord.inx <- order(nd.tbl[,3], nd.tbl[,4], decreasing = TRUE)
  nd.tbl <- nd.tbl[ord.inx, ];

  netData <- list(nodes=nodes,typeVec=dataSet$gene_type_vec, edges=edge.mat, modules=modules, prot.nms=prot.nms,gene.nms=gene.nms, containsGP=containsGP, middleLayerType = middleLayerType);

  dataSet$imgSet$node_table <- nd.tbl;

  if(exists("PeakSet",envir = .GlobalEnv)) {
    PeakSet <- get("PeakSet", envir = .GlobalEnv)
    netData[["peakNodes"]] <- PeakSet$nodes;
    netData[["peakEdges"]] <- PeakSet$edges;
  }

  if(file.exists("micSet.qs")){
    micSet <- qs::qread("micSet.qs");
    inx = names(micSet$met.mic) %in% nms;
    netData[["metMic"]] <- micSet$met.mic[inx];
    netData[["metMicScore"]] <- micSet$met.mic.score[inx];
    micTbl <- list();
    j = 1;
    for(i in 1:length(micSet$met.mic.table)){
      if(micSet$met.mic.table[[i]]["KEGG"] %in% nms){
        micTbl[[j]] <- micSet$met.mic.table[[i]]
        j = j + 1;
      }
    }
    netData[["metMicTable"]] <- micTbl
    netData[["micThresh"]] <- dataSet$mic.thresh;
  }

  if(any(dataSet$gene_type_vec == 4)){
    netData[["snpTable"]] <- snpTable;
    netData[["snpTableConsequence"]] <- snpTableConsequence;
  }

  netData[["colorVec"]] <- color.vec;
  dataSet$jsonNms$network <- filenm
  partialToBeSaved <<- c(partialToBeSaved, c(filenm))
  sink(filenm);
  cat(RJSONIO::toJSON(netData));
  sink();

  dataSet <<- dataSet;
  return(.set.nSet(dataSet));
}

# Define the function
parse_expression_data <- function(dataSet, exp, nms) {
  # Create an empty list to store the parsed data
  parsed_data <- list()
  # Loop through each name in nms
  for (i in seq_along(nms)) {
    nm <- nms[i]
    
    # Initialize a vector to store the values
    values <- numeric(2)
    
    # Get the expression value
    values[1] <- exp[i]
    
    # Find the corresponding value in the dataSet$exp.mat matrices
    for (name in names(dataSet$exp.mat)) {
      mat <- dataSet$exp.mat[[name]]
      if (nm %in% rownames(mat)) {
        values[2] <- mat[nm, 1]
      } else {
        values[2] <- 0
      }
    }
    
    # Add the vector to the parsed_data list
    parsed_data[[i]] <- values
  }
  
  # Return the parsed data
  return(parsed_data)
}