
#color shape layout important
my.convert.igraph <- function(dataSetObj=NA, net.nm, filenm, thera=FALSE, dim=3){
  dataSet <- .get.nSet(dataSetObj);
  dim <- as.numeric(dim)
  types_arr <<- dataSet$type;
  
  if(thera){
    g <- my.ppi;
  }else{
    g <- ppi.comps[[net.nm]];
  }
  
  if(length(V(g)$name)>100){
    #modules = FindCommunities("infomap", FALSE);
    modules = "NA"
  }else{
    modules = "NA"
  }
  
  if(met.type == "keggp"){
    netData = list()
    netData[["ids"]] = c(rownames(dataSet$seed[["ko"]]),rownames(dataSet$seed[["gene"]]), rownames(dataSet$seed[["met"]]));
    netData[["expr"]] = c(as.character(dataSet$seed[["ko"]]),as.character(dataSet$seed[["gene"]]), as.character(dataSet$seed[["met"]]));
    if(sum(dataSet$seed[["ko"]])!=0){
      ko.colsb.exp <- ComputeColorGradient(as.vector(dataSet$seed[["ko"]]), "black", centered); 
    }else{
      ko.colsb.exp <- rep("#ffff00", length(as.vector(dataSet$seed[["ko"]])));
    }
    if(sum(dataSet$seed[["gene"]])!=0){
      gene.colsb.exp <- ComputeColorGradient(as.vector(dataSet$seed[["gene"]]), "black", centered); 
    }else{
      gene.colsb.exp <- rep("#ffff00", length(as.vector(dataSet$seed[["gene"]])));
    }
    if(sum(dataSet$met.seed)!=0){
      met.colsb.exp <- ComputeColorGradient(as.vector(dataSet$met.seed), "black", centered); 
    }else{
      met.colsb.exp <- rep("#ffff00", length(as.vector(dataSet$met.seed)));
    }
    netData[["col"]] = c(ko.colsb.exp,gene.colsb.exp, met.colsb.exp);
    ord.inx<-order(-abs(as.numeric(netData[["expr"]])));
    netData[["expr"]] = netData[["expr"]][ord.inx]
    netData[["ids"]] = netData[["ids"]][ord.inx]
    keggp.allfeatures <<- netData[["ids"]][ord.inx]
    #netData[["ids"]] = doEntrez2KoMapping(netData[["ids"]])
    containsKo = length(rownames(dataSet$seed[["ko"]]))>0
    containsGene = length(rownames(dataSet$seed[["gene"]]))>0
    containsMet = length(rownames(dataSet$met.seed))>0
    
    enrType = "keggm";
    if(length(dataSet$seed[["ko"]])>0 && length(dataSet$met.seed)>0 ){
      enrType = "integ";
    }else if(length(dataSet$seed[["ko"]])>0){
      enrType = "ko";
    }else{
      enrType = "keggm";
    }
    ora.vec = netData[["ids"]]
    names(ora.vec) <- ora.vec;
    PerformKeggEnrichment(dataSet, "omicsnet_enrichment_0", "keggm", ora.vec);
    
    #Create met/gene pairs table Gene Name, ID, p value, FC, Metabolite Name, ID, p value, FC
    
    
    lengthMet=length(rownames(dataSet$met.seed))
    gDat = vector()
    metDat = vector()
    if(containsKo){
      dataSet$seed[["ko"]] = dataSet$seed[["ko"]][order(-dataSet$seed[["ko"]]),] 
      gDat = dataSet$seed[["ko"]]
    }else if(containsGene){
      gDat = dataSet$seed[["gene"]][order(-abs(dataSet$seed[["gene"]])),]
    }
    if(containsMet){
      metDat= dataSet$seed[["met"]][order(-abs(dataSet$seed[["met"]])),]
    }
    lengthgDat = length(gDat)
    replacingMet = c(1:lengthMet)
    replacinggDat = c(1:lengthgDat)
    
    if(containsMet && (containsGene || containsKo)){
      combinedTable = matrix(NA, nrow=max(c(lengthgDat, lengthMet)), ncol=6);
      combinedTable[replacinggDat,1]=names(gDat)
      combinedTable[replacinggDat,2]=doEntrez2SymbolMapping(names(gDat))
      combinedTable[replacinggDat,3]=unname(gDat)
      combinedTable[replacingMet,4]=names(metDat)
      combinedTable[replacingMet,5]=doKegg2NameMapping(names(metDat))
      combinedTable[replacingMet,6]=unname(metDat)
      colnames(combinedTable) = c("Gene ID", "Gene Symbol", "Gene Expr.", "KEGG ID", "Name", "Met Expr.")
      netData[["containsBoth"]] = 1
    }else if(containsMet){
      combinedTable = matrix(NA, nrow=lengthMet, ncol=3);
      combinedTable[,1]=names(metDat)
      combinedTable[,2]=doKegg2NameMapping(names(metDat))
      combinedTable[,3]=unname(metDat)
      colnames(combinedTable) = c("ID", "Name", "Expr.")
      netData[["containsBoth"]] = 0
    }else{
      combinedTable = matrix(NA, nrow=lengthgDat, ncol=3);
      combinedTable[,1]=names(gDat)
      combinedTable[,2]=doEntrez2SymbolMapping(names(gDat))
      combinedTable[,3]=unname(gDat)
      colnames(combinedTable) = c("ID", "Name", "Expr.")
      netData[["containsBoth"]] = 0
    }
    fast.write.csv(combinedTable, file="feature_table.csv", row.names=FALSE);
    sink(filenm);
    cat(toJSON(netData));
    sink();
    
    return(1)
  } else {
    netData = list()
    netData[["ids"]] = c(rownames(dataSet$seed[["ko"]]),rownames(dataSet$seed[["gene"]]), rownames(dataSet$seed[["met"]]), rownames(dataSet$seed[["tf"]]));
    netData[["expr"]] = c(as.character(dataSet$seed[["ko"]]),as.character(dataSet$seed[["gene"]]), as.character(dataSet$seed[["met"]]), as.character(dataSet$seed[["tf"]]));
    ord.inx<-order(-abs(as.numeric(netData[["expr"]])));
    netData[["expr"]] = netData[["expr"]][ord.inx]
    netData[["ids"]] = netData[["ids"]][ord.inx]
    keggp.allfeatures <<- netData[["ids"]][ord.inx]
  }
  
  nms <- V(g)$name;
  
  hit.inx <- match(nms, omics.net$node.data[,1]);
  lbls <- omics.net$node.data[hit.inx, 2];
  
  
  # setup shape (gene circle, other squares)
  if("peak" %in% dataSet$type){
    shapes <- rep("putative_met", length(nms));
  }else{
    shapes <- rep("gene", length(nms));
  }
  itypes <- rep("circle", length(nms));
  seeds <- rep("circle", length(nms));
  
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
  
  if("mic" %in% dataSet$type){
    node.dgr <- degree(g);
    mic.dgrs <- dataSet$mic.dgrs
    mic.dgrs <- mic.dgrs[names(mic.dgrs) %in% names(node.dgr)]
    node.dgr[names(mic.dgrs)] <- mic.dgrs
    node.dgr <- as.numeric(node.dgr);
  }else{
    node.dgr <- as.numeric(degree(g));
  }
  node.exp <- as.numeric(get.vertex.attribute(g, name="abundance", index = V(g)));
  
  # node size to degree values
  if(vcount(g)>500){
    min.size = 1;
  }else if(vcount(g)>200){
    min.size = 2;
  }else{
    min.size = 3;
  }
  
  node.sizes <- as.numeric(rescale2NewRange((log(node.dgr))^2, min.size, 9))*3 +5;
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
    node.colsb.exp[bad.inx] <- "#d3d3d3"; 
    node.colsw.exp[bad.inx] <- "#c6c6c6"; 
  }else{
    node.colsb.exp <- rep("#D3D3D3",length(node.dgr)); 
    node.colsw.exp <- rep("#C6C6C6",length(node.dgr)); 
  }
  
  gene.nms <- rownames(dataSet$exp.mat[["gene"]] )[dataSet$gene_type_vec %in% c(1,2)];
  prot.nms <- rownames(dataSet$exp.mat[["gene"]] )[dataSet$gene_type_vec %in% c(2,3)]; 
  snp.nms <- rownames(dataSet$exp.mat[["snp"]] ) 
  snp.nms <- snp.nms[!is.na(snp.nms)]
  
  # now update for bipartite network
  
  #edge.mat1 <<- edge.mat
  gene.inx <- nms %in% c(net.info$gene.ids, net.info$ko.ids);
  
  #predicted.inx <- !(nms %in% c(net.info$gene.ids, net.info$tf.ids, net.info$mir.ids, net.info$met.ids, net.info$ko.ids, net.info$protein.ids));
  if(length(net.info$gene.ids) == 0 && length(net.info$protein.ids) == 0 ){
    #gene.inx = predicted.inx;
    #predicted.inx = rep(FALSE,length(node.dgr));
  }
  
  kncmpd.inx <- nms %in% net.info$kncmpd.ids
  #predicted.inx <- nms %in% net.info$int.ids;
  prot.inx <- nms %in% net.info$protein.ids;
  tf.inx <- nms %in% net.info$tf.ids;
  mir.inx <- nms %in% net.info$mir.ids;
  peak.inx <- nms %in% net.info$peak.ids;
  met.inx <- nms %in% net.info$met.ids;
  reg.inx <- nms %in% net.info$reg.ids;
  mic.inx <- nms %in% net.info$mic.ids
  snp.inx <- nms %in% rownames(dataSet$mat[["snp"]]);
  
  containsGP <- any(gene.nms %in% prot.nms) ## check if any overlap between gene and protein
  cat("containsGP = ", containsGP, "\n")

  genes.inx <- nms %in% rownames(dataSet$seed[["gene"]])[dataSet$gene_type_vec %in% c(1,2)];
  tfs.inx <- nms %in% rownames(dataSet$seed[["tf"]]);
  mirs.inx <- nms %in% rownames(dataSet$seed[["mir"]]);
  mets.inx <- nms %in% rownames(dataSet$seed[["met"]]);

  topo.colsw[gene.inx] <- "#FF8484";
  topo.colsw[kncmpd.inx] <- "#00ffff";
  topo.colsw[prot.inx] <- "#FF8484";
  topo.colsw[tf.inx] <- "#39FF14"; 
  topo.colsw[mir.inx] <- "#00f6ff"; 
  topo.colsw[peak.inx] <- "#D3D3D3";
  
  topo.colsw[met.inx] <- "#ffff00";
  topo.colsw[reg.inx] <- "#ff9900";
  topo.colsw[mic.inx] <- "#39FF14";
  
  
  topo.colsb <- topo.colsw;
  
  
  if(containsGP){
    topo.colsb[genes.inx] <- "#D3D3D3";
    topo.colsw[genes.inx] <- "#D3D3D3";
    color.vec <- c("#D3D3D3", "#FF8484", "#39FF14","#00f6ff", "#D3D3D3", "#00ffff", "#ffff00", "#ff9900", "#39FF14");
  }else{
    color.vec <- c("#FF8484", "#FF8484", "#39FF14","#00f6ff", "#D3D3D3", "#00ffff", "#ffff00", "#ff9900", "#39FF14");
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
  
  
  seeds[genes.inx] <- "gene";
  seeds[tfs.inx] <- c("tf",itypes[tf.inx]);
  seeds[mirs.inx] <- c("mir",itypes[mir.inx]);
  seeds[mets.inx] <- c("met", itypes[met.inx]);
  
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
  
  seed.inx <- nms %in% unique(dataSet$seeds.proteins);
  
  seed_arr <- rep("notSeed",length(node.dgr));
  seed_arr[seed.inx] <- "seed";
  
  type <- rep(FALSE,length(node.dgr))
  type[seed.inx] = TRUE
  
  lblsu <<- lbls;
  node_attr = list.vertex.attributes(g);
  node_attr = node_attr[which(node_attr!="names")] 
  
  # now create the json object
  nodes <- vector(mode="list");
  
  
  library(stringr)
  displayedLabel<-lbls;
  long.inx <- which(str_length(lbls) > 32);
  displayedLabel[long.inx] <- paste0(strtrim(lbls[long.inx],  rep(32, length(lbls[long.inx]))), "..." )
  
  for(i in 1:length(node.sizes)){
    nodes[[i]] <- list(
      id=nms[i], 
      idx=i,
      label=lbls[i],
      displayedLabel = displayedLabel[i],
      size=node.sizes[i], 
      size2d=node.sizes2d[i],
      type=shapes[i],
      molType=shapes[i],
      types=types[i],
      seed=seeds[i],
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
  
  summary = vector(mode="list");
  summary[[1]] = list(nodes=vcount(g));
  summary[[2]] = list(edges=ecount(g));
  
  # save node table
  nd.tbl <- data.frame(Id=nms, Label=lbls, Degree=node.dgr, Betweenness=round(node.btw,2));
  # order 
  ord.inx <- order(nd.tbl[,3], nd.tbl[,4], decreasing = TRUE)
  nd.tbl <- nd.tbl[ord.inx, ];
  fast.write.csv(nd.tbl, file="node_table.csv", row.names=FALSE);
  
  # covert to json
  require(RJSONIO);
  db.path <- paste(lib.path,data.org,"/entrez.rds", sep="");
  conv <- readRDS(db.path)
  netData <- list(nodes=nodes,typeVec=dataSet$gene_type_vec, edges=edge.mat, modules=modules, conv = conv, prot.nms=prot.nms,gene.nms=gene.nms, snp.nms=snp.nms);

  
  if(exists("PeakSet",envir = .GlobalEnv)) {
    PeakSet <- get("PeakSet", envir = .GlobalEnv)
    netData[["peakNodes"]] <- PeakSet$nodes;
    netData[["peakEdges"]] <- PeakSet$edges;
  } 

  if("mic" %in% dataSet$type ){
    inx = names(dataSet$met.mic) %in% nms;
    netData[["metMic"]] <- dataSet$met.mic[inx];
    netData[["metMicScore"]] <- dataSet$met.mic.score[inx];
    micTbl <- list();
    j = 1;
    for(i in 1:length(dataSet$met.mic.table)){
      if(dataSet$met.mic.table[[i]]["KEGG"] %in% nms){
        micTbl[[j]] <- dataSet$met.mic.table[[i]]
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
  dataSet$jsonNms$network <<- filenm
  partialToBeSaved <<- c(partialToBeSaved, c(filenm))
  sink(filenm);
  cat(toJSON(netData));
  sink();
  return(.set.nSet(dataSetObj));
}