  
my.peak.net <- function(dataSetObj=NA){
  #save.image("peaknet.RData");
  #print("peaknet1");

  dataSet <- .get.nSet(dataSetObj);
  if(exists("PeakSetDone",envir = .GlobalEnv)) {
    PeakSet <- get("PeakSetDone", envir = .GlobalEnv);
  } else if(exists("PeakSet",envir = .GlobalEnv)) {
    PeakSet <- get("PeakSet", envir = .GlobalEnv);
  } else {
    return(0L);
  }
  nodes <- PeakSet$nodes;
  edges <- PeakSet$edges;
  met.kg.ids <- nodes$KEGGID;

  tmpsIDs <- nodes$HMDBID;
  tmpsIDs <- tmpsIDs[tmpsIDs != ""]
  if(all(grepl("^C[0-9]+", tmpsIDs))){
    met.ids <- nodes$KEGGID;
    met.ids[met.ids == ""] <- nodes$formula[met.ids == ""];
    lbls <- doKegg2NameMapping(met.ids);
  } else if (all(grepl("^HMDB[0-9]+", tmpsIDs))) {
    met.ids <- nodes$HMDBID;
    met.ids[met.ids == ""] <- nodes$formula[met.ids == ""];
    lbls <- doHMDB2NameMapping(met.ids);
  } else if (all(grepl("^PBCM[0-9]+", tmpsIDs))) {
    met.ids <- nodes$HMDBID;
    met.ids[met.ids == ""] <- nodes$formula[met.ids == ""];
    lbls <- doPubchem2NameMapping(met.ids);
  }

  met.kg.ids[met.kg.ids == ""] <- nodes$formula[met.kg.ids == ""];
  PeakSet$nodes$KEGGID_realID <- met.kg.ids;
  PeakSet$nodes$KEGGID <- met.ids;
  art.ids <- PeakSet$nodes$name[which(PeakSet$nodes$class == "Artifact")]

  PeakSet$mets <- unique(PeakSet$nodes$KEGGID[which(PeakSet$nodes$class == "Metabolite")]);
  PeakSet$put.mets <- PeakSet$nodes$KEGGID[which(PeakSet$nodes$class == "Putative metabolite")];
  PeakSet$mets.ids <- PeakSet$nodes$name[which(PeakSet$nodes$class == "Metabolite")];
  PeakSet$put.mets.ids <- PeakSet$nodes$name[which(PeakSet$nodes$class == "Putative metabolite")];


  PeakSet$nodes.df <- data.frame(Id=met.ids, Label=lbls)
  edges.df0 <- data.frame(from=edges[,1], to=edges[,2])
  edges.df01 <- edges.df0[!edges.df0[,1] %in% art.ids,];
  edges.df02 <- edges.df01[!edges.df01[,2] %in% art.ids,];
  conv1 <- data.frame(from=PeakSet$nodes$name, Source=PeakSet$nodes$KEGGID);
  conv2 <- data.frame(to=PeakSet$nodes$name, Target=PeakSet$nodes$KEGGID);
  edges.df1 <- merge(edges.df02, conv1, by = "from");
  edges.df2 <- merge(edges.df1, conv2, by = "to");
  PeakSet$edges.df <- data.frame(Source=edges.df2$Source, Target=edges.df2$Target);

    if(length(which(grepl("HMDB", PeakSet$mets)))/length(PeakSet$mets) > 0){
      PeakSet$mets <- doHMDB2KEGGMapping(PeakSet$mets);
      PeakSet$nodes.df[,1] <- doHMDB2KEGGMapping(PeakSet$nodes.df[,1]);
      PeakSet$edges.df[,1] <- doHMDB2KEGGMapping(PeakSet$edges.df[,1]);
      PeakSet$edges.df[,2] <- doHMDB2KEGGMapping(PeakSet$edges.df[,2]);

    }

    if(length(which(grepl("PBCM", PeakSet$mets)))/length(PeakSet$mets) > 0){
      PeakSet$mets <- doPubchem2KEGGMapping(PeakSet$mets);
      PeakSet$nodes.df[,1] <- doPubchem2KEGGMapping(PeakSet$nodes.df[,1]);
      PeakSet$edges.df[,1] <- doPubchem2KEGGMapping(PeakSet$edges.df[,1]);
      PeakSet$edges.df[,2] <- doPubchem2KEGGMapping(PeakSet$edges.df[,2]);

    }

  dataSet$seed[["peak"]] <- data.frame(rep(0, length(unique(PeakSet$mets))));

  rownames(dataSet$seed[["peak"]]) <-unique(PeakSet$mets);
  dataSet$exp.mat[["peak"]] <- data.frame(rep(0, length(PeakSet$NetID_output$peak_id)));
  rownames(dataSet$exp.mat[["peak"]]) <- paste0(PeakSet$NetID_output$peak_id,"_", PeakSet$NetID_output$medMz);

  cmp.nms <- PrepareInputList(dataSet, rownames(dataSet$seed[["peak"]]), data.org, "peak", "kegg");
  
  qs:::qsave(PeakSet, "PeakSet_net.qs");
  dataSet <<- dataSet;
  if(.on.public.web){
    .set.nSet(dataSet);
    return(1L);
  }else{
    return(.set.nSet(dataSet));
  }
}
