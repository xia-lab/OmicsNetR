##################################################
## R scripts for OmicsNet
## Description: MS Peak Analysis Utilities
## Author: Zhiqiang Pang, zhiqiang.pang@mcgill.ca
###################################################

### ------ Content ------- ###
# Section I: MS peak data uploading and preprocessing
# Section II: Prepare MS nodes
# Section III: Prepare MS edges
# Section IV: Prepare & Optimize Network
# Section V: Output Network
# Section VI: Other Utilities
### ---- End of Content ----###

#' InitiatePeakSet
#'
#' @param ionMode ion mode, postive is 1 and negative is -1
#' @param ppm ppm value of mass
#' @param database database, can be "KEGG", "HMDB" or "Pubchem"
#' @param p_cutoff p value cutoff for filtration, default is 0.1
#' @export
#'

InitiatePeakSet <- function(ionMode = -1,
                            ppm = 10,
                            database = "KEGG",
                            p_cutoff = 0.1) {

    options(scipen=999);
    PeakSet <- list();
    PeakSet$ParamSet <- list(DB = database,
                             polarity = ionMode, # "1" is pos, "-1" is neg
                             ppm = ppm,
                             p_cutoff = p_cutoff,
                             data.org = data.org);
    PeakSet <<- PeakSet;
    #data.org <<- "ko";
    return(1);
}

#' ImportMSPeaks
#'
#' @param PeakFile PeakFile path
#' @export
#'
ImportMSPeaks <- function(PeakFile = NA) {
  require("dplyr")
  require("tidyr")
  t0 <- Sys.time();
  #RequiredPackageCheck(); # dangerous on a production server!
  cat("Running into ImportMSPeaks! \n");
  if(exists("PeakSet",envir = .GlobalEnv)) {
    PeakSet <- get("PeakSet", envir = .GlobalEnv)
  } else {
    InitiatePeakSet();
    PeakSet <- get("PeakSet", envir = .GlobalEnv)
  }
  if(is.na(PeakFile)) {return(0)}

  raw_data <- read.csv(PeakFile);
  if(ncol(raw_data) == 3) {
    first_sample_col_num <- 4;
  } else if(ncol(raw_data) == 4) {
    first_sample_col_num <- 5;
  } else {
    stop("Data format is wrong!")
  }

  raw_data <- cbind(id = c(1:nrow(raw_data)), raw_data);
  colnames(raw_data) <- c("id", "medMz", "medRt", "Sample",
                          "pvalue", "Unknown1", "Unknown2")[1:ncol(raw_data)];
  raw_data[,3] <- round(raw_data[,3], 3);

  PeakSet[["raw_data"]] <- raw_data;
  dataSet$exp.mat[["peak"]] <- data.frame(rep(0, nrow(raw_data)));
  rownames(dataSet$exp.mat[["peak"]]) <- paste0(raw_data$id, "_", raw_data$medMz);

  PeakSet$ParamSet$RTTol <- (max(raw_data[,3]) -
                               min(raw_data[,3]))*0.005; #0.5% of whole rt range as rt tol

  PeakSet$Data <- Peak_cleanup(Mset = PeakSet,
                               first_sample_col_num = first_sample_col_num);

  PeakSet$HMDB_library <- .importCmpdLib(PeakSet$ParamSet$DB);
  PeakSet$empirical_rules <- .importEmpiricalRule();
  PeakSet$propagation_rule <- .importPropagationRule();
  PeakSet$MS2_library <- .importMS2Lib();
  cat("Peak Data import done \n Time consumed this step: ", Sys.time() - t0, "\n")
  PeakSet <<- PeakSet;
  qs::qsave(PeakSet, file = "PeakSet_initial.qs")
  Sys.sleep(0.15);  # CRITICAL: Prevent race condition
  dataSet <<- dataSet;
  return(1);
}

#' UpdateCmpdDB
#'
#' @param dbNM database name, KEGG, HMDB or PubChem
#' @export
UpdateCmpdDB <- function(dbNM){
  PeakSet <<- qs::qread(file = "PeakSet_initial.qs");
  cat("Now the db is being changed as ", dbNM, "\n")
  if(dbNM == "kegg"){
    PeakSet$HMDB_library <- .importCmpdLib("KEGG");
  } else if (dbNM == "hmdb"){
    PeakSet$HMDB_library <- .importCmpdLib("HMDB");
  } else if(dbNM == "pubchem"){
    PeakSet$HMDB_library <- .importCmpdLib("Pubchem");
  }
  PeakSet <<- PeakSet;
  return(1)
}

#' PerformDataProcessing
#' @export
#'

PerformDataProcessing <- function() {

  cat("Running into PerformDataProcessing! \n");
  if(.on.public.web){
    ## This package is no longer required for local R package, because the source code has been included
    #requireNamespace("lc8")
  }

  require("stringr")

  t0 <- Sys.time();
  if(exists("PeakSet",envir = .GlobalEnv)) {
    PeakSet <- get("PeakSet", envir = .GlobalEnv)
  } else {
    return(0L)
  }
  cat("Start preparing node set...")
  NodeSet <- Initiate_nodeset(PeakSet)
  LibrarySet <- Initiate_libraryset(PeakSet)
  LibrarySet_expand <- Expand_libraryset(LibrarySet, Mset = PeakSet)

  #StructureSet = Initilize_empty_structureset(NodeSet)
  StructureSet <- initialize_structureset(NodeSet)
  StructureSet <- Match_library_structureset(LibrarySet,
                                             LibrarySet_expand,
                                             StructureSet,
                                             NodeSet,
                                             ppm_tol =
                                               PeakSet$ParamSet$ppm*1e-6)

  cat("done!\nStart preparing edge set...\n")
  # return(testtingg(Mset = PeakSet,
  #                  NodeSet,
  #                  mz_tol_abs = 2e-4,
  #                  mz_tol_ppm = PeakSet$ParamSet$ppm,
  #                  rt_tol_bio = Inf,
  #                  rt_tol_nonbio = PeakSet$ParamSet$RTTol))

  EdgeSet <- Initiate_edgeset(Mset = PeakSet,
                              NodeSet,
                              mz_tol_abs = 2e-4,
                              mz_tol_ppm = PeakSet$ParamSet$ppm,
                              rt_tol_bio = Inf,
                              rt_tol_nonbio = PeakSet$ParamSet$RTTol)
  ## errr in the function above ---
  EdgeSet_expand <- Expand_edgeset(EdgeSet,
                                   NodeSet,
                                   StructureSet,
                                   RT_cutoff = PeakSet$ParamSet$RTTol,
                                   inten_cutoff = 1e4,
                                   types = c("ring_artifact",
                                             "oligomer_multicharge",
                                             "heterodimer"),
                                   Mset = PeakSet,
                                   LibrarySet = LibrarySet)
  EdgeSet_all = merge_edgeset(EdgeSet, EdgeSet_expand)
  cat("done!\nTime consumed this step: ", Sys.time() - t0, "\n");
  PeakSet$LibrarySet <- LibrarySet;
  PeakSet$NodeSet <- NodeSet;
  PeakSet$EdgeSet <- EdgeSet;
  PeakSet$StructureSet <- StructureSet;
  PeakSet$EdgeSet_all <- EdgeSet_all;
  PeakSet <<- PeakSet;
  print("Peak data processing completed!");
  return(1L)
}

#' PerformPropagation
#' @export
PerformPropagation <- function() {
  cat("Running into PerformPropagation! \n");
  t0 <- Sys.time();
  if(exists("PeakSet",envir = .GlobalEnv)) {
    PeakSet <- get("PeakSet", envir = .GlobalEnv)
  } else {
    return(0L)
  }
  cat("Start Propagating from seeds...\n")
  StructureSet <- Propagate_structureset(PeakSet,
                                         biotransform_step = 2,
                                         artifact_step = 3,
                                         propagation_ppm_threshold = PeakSet$ParamSet$ppm*1e-6,
                                         propagation_abs_threshold = 0e-4,
                                         record_RT_tol = PeakSet$ParamSet$RTTol,
                                         record_ppm_tol = PeakSet$ParamSet$ppm*1e-6/2);
  PeakSet$StructureSet <- StructureSet;
  cat("done!\nStart scoring structure ...");
  StructureSet_df <- Score_structureset(PeakSet);
  StructureSet_df <- Score_structure_propagation(StructureSet_df,
                                                 artifact_decay = -0.5)
  PeakSet$StructureSet_df <- StructureSet_df;
  cat("done!\n Time consumed this step: ", Sys.time() - t0, "\n")
  PeakSet <<- PeakSet;
  print("Peak data constructing propagation completed!");
  return(1L);
}

#' PerformGlobalOptimization
#' @export
PerformGlobalOptimization <- function() {
  cat("Running into PerformGlobalOptimization! \n");
  t0 <- Sys.time();
  if(exists("PeakSet",envir = .GlobalEnv)) {
    PeakSet <- get("PeakSet", envir = .GlobalEnv);
  } else {
    return(0L);
  }
  cat("Start Preparing ILP Set for optimization...\n")
  PeakSet$ILPSet <- .prepareILPSet(PeakSet);
  PeakSet <- .networkConstruct(PeakSet);
  cat("Global done!\n Time consumed this step: ", Sys.time() - t0, "\n")
  PeakSet <<- PeakSet;
  print("Global optimization completed!");
  return(1L);
}

#' PerformAnnotation
#' @export
PerformAnnotation <- function() {
  cat("Running into PerformAnnotation! \n");
    t0 <- Sys.time();
    if(exists("PeakSet",envir = .GlobalEnv)) {
        PeakSet <- get("PeakSet", envir = .GlobalEnv)
    } else {
        return(0L)
    }
    cat("Start Perform Network Annotation...\n")
    ResSet <- path_annotation(PeakSet$ILPSet,
                              PeakSet$NetworkSet,
                              solution = "ilp_solution")

    PeakSet$NetID_output =  get_NetID_output(ResSet$ilp_nodes, simplified = T)
    PeakSet$NetID_output <- cbind(PeakSet$NetID_output,
                                  HMDB = matchHMDB(ResSet$core_annoted_met,
                                                   Restable = PeakSet$NetID_output, 0, 5))
    PeakSet$NetID_output <- cbind(PeakSet$NetID_output,
                                  KEGG = convert2KEGG(PeakSet[["NetID_output"]][["HMDB"]],
                                                      PeakSet$HMDB_library))
    fast.write.csv(PeakSet$NetID_output, file = "peak_annotation.csv", row.names = FALSE)

    vs <- vertex.attributes(PeakSet$NetworkSet$cleanNetwork);
    mhdbvalues <- matchHMDB(ResSet$core_annoted_met,
                            Restable = as.data.frame(vs),
                            1,2)
    PeakSet$NetworkSet$cleanNetwork <-
      set_vertex_attr(PeakSet$NetworkSet$cleanNetwork, "HMDBID",
                      value = mhdbvalues)
    PeakSet$NetworkSet$cleanNetwork <-
      set_vertex_attr(PeakSet$NetworkSet$cleanNetwork, "KEGGID",
                      value = convert2KEGG(mhdbvalues,
                                           PeakSet$HMDB_library))

    nods <- igraph::as_data_frame(PeakSet[["NetworkSet"]][["cleanNetwork"]], what = "vertices");
    nods <- nods[,-c(12:16)]
    egs <- igraph::as_data_frame(PeakSet[["NetworkSet"]][["cleanNetwork"]], what = "edges");
    egs <- egs[,-c(13:23)]
    nods <- .filterCurrency(nods);
    PeakSet$nodes <- nods;
    PeakSet$edges <- egs;
    PeakSet <- .filterNonSigCMPDs(PeakSet);
    PeakSet <- .enforceNodeAdduction(PeakSet);
    qs::qsave(PeakSet$Data, file = "PeakSet_data.qs")
    PeakSet <<- PeakSet <- .cleanPeakSet(PeakSet);
    qs::qsave(PeakSet, file = "PeakSet_done.qs");
    Sys.sleep(0.15);  # CRITICAL: Prevent race condition
    cat("Annotation done!\n Time consumed this step: ", Sys.time() - t0, "\n")
    print("Peak annotation completed!");
    return(1L)
}

GetFastPeak <- function(){
    ImportMSPeaks("../../data/test/ibd_peaks_p.csv");
    PeakSetDone <<- qs::qread("../../data/test/PeakSet_done.qs");
    PeakSet_data <- qs::qread("../../data/test/PeakSet_data.qs");
    qs::qsave(PeakSet_data, "PeakSet_data.qs");
}

#### Other internal functions

.importCmpdLib <- function(DB) {

  if(exists(".on.public.web",envir = .GlobalEnv)) {
    .on.public.web <- get(".on.public.web", envir = .GlobalEnv)
  } else {
    .on.public.web <- FALSE;
  }

  if(!.on.public.web){
    if(DB == "HMDB"){
      file_path <- system.file('db/hmdb_lib.qs', package = "OmicsNetR")
    } else if(DB == "KEGG") {
      file_path <- system.file('db/kegg_lib.qs', package = "OmicsNetR")
    } else if(DB == "Pubchem"){
      file_path <- system.file('db/pubchem_lib.qs', package = "OmicsNetR")
    }
    qs::qread(file_path)
  } else {
    if(DB == "HMDB"){
      qs::qread("../../data/lib/hmdb_lib.qs")
    } else if(DB == "KEGG") {
      qs::qread("../../data/lib/kegg_lib.qs")
    } else if(DB == "Pubchem"){
      qs::qread("../../data/lib/pubchem_lib.qs")
    }
  }
}
.importEmpiricalRule <- function() {
  if(exists(".on.public.web",envir = .GlobalEnv)) {
    .on.public.web <- get(".on.public.web", envir = .GlobalEnv)
  } else {
    .on.public.web <- FALSE;
  }
  if(!.on.public.web){
    file_path <- system.file('db/empirical_rule.qs', package = "OmicsNetR")
    qs::qread(file_path)
  } else {
    qs::qread("../../data/lib/empirical_rule.qs")
  }
}
.importPropagationRule <- function() {
  if(exists(".on.public.web",envir = .GlobalEnv)) {
    .on.public.web <- get(".on.public.web", envir = .GlobalEnv)
  } else {
    .on.public.web <- FALSE;
  }
  if(!.on.public.web){
    file_path <- system.file('db/propagation_rule.qs', package = "OmicsNetR")
    qs::qread(file_path)
  } else {
    qs::qread("../../data/lib/propagation_rule.qs")
  }
}
.importMS2Lib <- function() {
  if(exists(".on.public.web",envir = .GlobalEnv)) {
    .on.public.web <- get(".on.public.web", envir = .GlobalEnv)
  } else {
    .on.public.web <- FALSE;
  }
  if(!.on.public.web){
    file_path <- system.file('db/ms2_lib.qs', package = "OmicsNetR")
    qs::qread(file_path)
  } else {
    qs::qread("../../data/lib/ms2_lib.qs")
  }
}

.prepareILPSet <- function(PeakSet) {
  
  ILPSet = list()
  
  print("Start initiate_ilp_nodes")
  ILPSet[["ilp_nodes"]] <-
    initiate_ilp_nodes(PeakSet);
  print("Start initiate_ilp_edges")
  
  ILPSet[["ilp_edges"]] <-
    initiate_ilp_edges(PeakSet$EdgeSet_all,
                       ILPSet,
                       Exclude = ""); #can be faster
  
  print("initiate_heterodimer_ilp_edges")
  ILPSet[["heterodimer_ilp_edges"]] <-
    initiate_heterodimer_ilp_edges(PeakSet$EdgeSet_all,
                                   ILPSet,
                                   PeakSet$NodeSet);
  
  print("Finish ILPSet initialization")
  
  # Score
  ILPSet[["ilp_nodes"]] <-
    score_ilp_nodes(ILPSet,
                    metabolite_score = 0.1,
                    putative_metabolite_score = 0,
                    artifact_score = 0,
                    unknown_score = -0.5)
  
  ILPSet[["ilp_edges"]] <- score_ilp_edges(ILPSet, PeakSet$NodeSet)
  
  ILPSet[["heterodimer_ilp_edges"]] <-
    score_heterodimer_ilp_edges(ILPSet,
                                type_score_heterodimer = 0,
                                MS2_score_experiment_fragment = 0.5)
  print("Finish ILPSet scoring")
  
  ILPSet[["para"]] <- Initiate_cplexset(ILPSet)
  
  cat("Done!\n")
  return(ILPSet)
}

.networkConstruct <- function(PeakSet) {
    require("lpsymphony");
    CplexSet <- PeakSet$ILPSet;
    mat <- CplexSet$para$mat;
    max <- TRUE;
    obj <- CplexSet[["para"]][["obj"]];
    sense <- CplexSet[["para"]][["sense"]];
    sense[sense == "E"] <- "==";
    sense[sense == "L"] <- "<=";
    dir <- sense;
    rhs <- CplexSet[["para"]][["rhs"]];

    # NOTE: "I" is integer, "B" is binary,
    # both of them are running too slowly (>1h)
    # for huge constraints matrix
    types <- "C"
    OptiSolution <- lpsymphony::lpsymphony_solve_LP(obj,
                                        mat,
                                        dir,
                                        rhs,
                                        types = types,
                                        max = max,
                                        verbosity = -1, gap_limit = 1e-3);

    CplexSet <- add_SYM_solution(CplexSet, OptiSolution);

    require("igraph")
    NetworkSet = list();
    NetworkSet = Initiate_networkset(CplexSet,
                                     PeakSet$StructureSet_df,
                                     PeakSet$LibrarySet,
                                     solution = "ilp_solution");

    PeakSet$NetworkSet <- NetworkSet;
    ValidNetwork <-  NetworkSet[["g_all_valid"]];
    v_att_nms <- names(vertex.attributes(ValidNetwork));
    ValidNetwork <- delete_vertex_attr(ValidNetwork, "score_class");
    ValidNetwork <- delete_vertex_attr(ValidNetwork, "score_known_rt");
    ValidNetwork <- delete_vertex_attr(ValidNetwork, "score_RDBE");
    ValidNetwork <- delete_vertex_attr(ValidNetwork, "score_element_ratio");
    ValidNetwork <- delete_vertex_attr(ValidNetwork, "score_missing_isotope_Cl");
    ValidNetwork <- delete_vertex_attr(ValidNetwork, "score_database_origin");
    ValidNetwork <- delete_vertex_attr(ValidNetwork, "score_mz");
    ValidNetwork <- delete_vertex_attr(ValidNetwork, "score_propagation");
    ValidNetwork <- delete_vertex_attr(ValidNetwork, "ilp_solution");
    ValidNetwork <- delete_vertex_attr(ValidNetwork, "log10_inten");
    ValidNetwork <- delete_vertex_attr(ValidNetwork, "struct_set_id");
    ValidNetwork <- delete_vertex_attr(ValidNetwork, "steps");
    PeakSet$NetworkSet$cleanNetwork <- ValidNetwork;

    PeakSet$ILPSet <- CplexSet;
    return(PeakSet)
}

.filterCurrency <- function(nodesVec = NULL) {
  if(is.null(nodesVec)) stop("No nodes procvided!")
  if(all(is.na(nodesVec$KEGGID))) return(nodesVec)

  if(exists(".on.public.web",envir = .GlobalEnv)) {
    .on.public.web <- get(".on.public.web", envir = .GlobalEnv)
  } else {
    .on.public.web <- FALSE;
  }
  if(!.on.public.web){
    file_path <- system.file('db/currency.qs', package = "OmicsNetR")
    curVec <- qs::qread(file_path)
  } else {
    curVec <- qs::qread("../../data/lib/currency.qs");
  }

  res <-
    vapply(nodesVec$KEGGID, FUN = function(x){x %in% curVec$V1},
           FUN.VALUE = vector(mode = "logical",
                              length = 1))
  nodesVec$KEGGID[res] <- "";
  nodesVec;
}

.filterNonSigCMPDs <- function(PeakSet) {
  pval <- PeakSet$ParamSet$p_cutoff;
  if(pval == -1){cat("Initial Analysis -- No filtration !\n")}
  if(pval == -1 | is.null(PeakSet[["Data"]]$pvalue)) {return(PeakSet)} # No or notuse pvalue, use all, not to filter!
  sigPeaks <- PeakSet$Data$id[PeakSet[["Data"]]$pvalue < pval];
  edgeDT <- PeakSet$edges
  res0 <- apply(edgeDT, 1, FUN = function(x) {
    (x[4] %in% sigPeaks & x[5] %in% sigPeaks) #& (x[7] != "Biotransform")
  })
  sigEDT <- PeakSet$edges <- edgeDT[res0,];
  nodes2kepp <- unique(c(sigEDT$from, sigEDT$to));
  nodeDT <- PeakSet$nodes;
  PeakSet$nodes <- nodeDT[nodeDT$name %in% nodes2kepp,]
  PeakSet
}

.cleanPeakSet <- function(PeakSet) {
  PeakSet$ParamSet <-
    PeakSet$raw_data <-
    PeakSet$Data <-
    PeakSet$HMDB_library <-
    PeakSet$empirical_rules <-
    PeakSet$propagation_rule <-
    PeakSet$MS2_library <-
    PeakSet$LibrarySet <-
    PeakSet$NodeSet <-
    PeakSet$EdgeSet <-
    PeakSet$StructureSet <-
    PeakSet$EdgeSet_all <-
    PeakSet$StructureSet_df <-
    PeakSet$ILPSet <-
    PeakSet$NetworkSet <- NULL;
  gc(full = TRUE);
  PeakSet;
}

.enforceNodeAdduction <- function(PeakSet) {
  nodt <- PeakSet$nodes;
  res <- gsub(pattern = "([a-zA-Z]{1})1{1}(\\D|$)",
              replacement = "\\1\\2",
              x = gsub(pattern = "([a-zA-Z]{1})1{1}(\\D|$)",
                       replacement = "\\1\\2",
                       x =nodt$transform));

  PeakSet$nodes <-
    cbind(nodt, adduction =
            vapply(1:length(res), FUN = function(x){
              if(nodt[x,6] == "Heterodimer"){
                return("Heterodimer")
              } else if(nodt[x,6] == "Natural_abundance") {
                if(nodt[x,9] == "[13]C1C-1") return("[M+1(13C)]")
                if(nodt[x,9] == "[18]O1O-1") return("[M+2(18O)]")
                if(nodt[x,9] == "[37]Cl1Cl-1") return("[M+2(35cl)]")
                return("[M+Isotope]")
              } else if (res[x] =="") {
                return("[M]")
              } else {
                return(paste0("[M+", res[x], "]"))
              }},
              FUN.VALUE =
                vector(mode = "character",
                       length = 1)))


  return(PeakSet)
}

# RequiredPackageCheck <- function() {
#     insPkgs <- names(installed.packages()[,1]);
#     checkPkgsList <- c("enviPat", "slam", "readr", "stringi", "pracma", "stringr", "janitor", "lpsymphony")
#     if(!("lc8" %in% insPkgs)) {
#         cat("R package lc8 is required, but not found in this PC/server, is being installed automatically now, don't worry...\n")
#         install.packages("https://www.dropbox.com/s/ef4o314x2caxz7u/Formula_manipulation.tar.gz", repos = NULL, method = "wget")
#     }
#
#     if(any(!(checkPkgsList %in% insPkgs))) {
#         cat("R package ", checkPkgsList[!(checkPkgsList %in% insPkgs)],
#             " are required, but not found in this PC/server, are being installed automatically now,
#             may need some time, don't worry...\n")
#         BiocManager::install(checkPkgsList[!(checkPkgsList %in% insPkgs)])
#     }
# }

### NetID network expansion
## This section is designed to expand the network generated by netID
## by extracting the metabolites (already annotated as seeds compound)
## and then fetching entire potential metabolic network from KEGG
## and then filtering the network by using PCSF to minimize the network
## and then merging the network generated by NetID with the minimal one from the KEGG
## finally outputting the network for visualization

extendMetPeakNetwork <- function(table.nm) {
  if(exists("PeakSet",envir = .GlobalEnv)) {
    PeakSet <- get("PeakSet", envir = .GlobalEnv)
  } else {
    return(0)
  }

  # if(!is.null(edgeu.res.list[["m2m"]])){
  #   cat("m2m already been performed! skipping... \n")
  #   return(1)
  # }
  if(!exists("old_m2mType",envir = .GlobalEnv)){
    old_m2mType <<- table.nm;
  } else if(old_m2mType == table.nm) {
    cat("m2m already been performed! skipping... \n")
    return(1)
  } else {
    old_m2mType <<- table.nm;
  }
  protein.vec <- result.listu[["protein.vec"]];

  kgids1 <-
    PeakSet[["nodes"]]$KEGGID[
      grepl(pattern = "C[0-9][0-9][0-9][0-9][0-9]",
            x = PeakSet[["nodes"]]$KEGGID)];

  kgids2 <-
    protein.vec[
      grepl(pattern = "C[0-9][0-9][0-9][0-9][0-9]",
            x = protein.vec)];

  kgids <- unique(c(kgids1, kgids2));
  netInv <- "direct";
  continueidx <- TRUE; intk <- 0;
  while(continueidx){
    ## this while loop is designed to extract all metabolic components from
    ## KEGG network. This is fast to finish. Usually less than 1 sec.
    cur.num <- length(kgids);
    res <- QueryM2mSQLiteNet(table.nm = table.nm,
                             q.vec = unique(kgids),
                             inv = netInv);
    if(intk < 1) {
      kgids0 <- c(kgids, res$productID);
    }
    intk <- intk + 1;
    kgids <- unique(c(kgids, res$productID));
    if(cur.num == length(kgids)){
      #reach the maximum of the subnetwork -> to stop
      continueidx <- FALSE
    }
  }

  ## Now, construct a graph
  require("igraph")
  require("dplyr")

  df <- data.frame(from = res$sourceID,
                   to = res$productID, weight = 0.2) #TODO: the weight can be optimized
  cmpd <- data.frame(ID = c(res[,1],res[,3]), NM = c(res[,2],res[,4])) %>% distinct()

  ig.obj <- igraph::graph_from_data_frame(df, directed = FALSE, vertices = cmpd)
  ppi <- simplify(ig.obj)

  ## Now, filter with PCSF
  ## Giving prize to the ones originally from peak or directly extended from the metabolites of peaks
  kgids_prize <- unique(c(kgids0, kgids1, kgids2))
  expr.vec <- sapply(kgids_prize, FUN = function(x) {
    rres <- mean(PeakSet[["NetID_output"]][["log10_inten"]][x == PeakSet[["nodes"]][["KEGGID"]]]);
    if(is.nan(rres)){return(1)} else {return(round(rres,2))}
  })
  names(expr.vec) <- unique(c(kgids0, kgids1, kgids2))

  g <- Compute.SteinerForest(ppi, expr.vec, w = 5, b = 100, mu = 0.0005);

  edgeList <- get.data.frame(g, "edges");
  edgeList <- cbind(rownames(edgeList), edgeList);
  colnames(edgeList) <- c("Id", "sourceID", "productID");

  ## Now, merge with original reuslts
  ##
  edge.res <- data.frame(Source=edgeList[,"sourceID"],
                         Target=edgeList[,"productID"],
                         stringsAsFactors=FALSE);

  # Filter out rows with NA values in Source or Target columns
  edge.res <- edge.res[complete.cases(edge.res), ];

  if(nrow(edge.res)!=0){
    row.names(edge.res) <- 1:nrow(edge.res);
  }
  fast.write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
  node.ids <- c(edgeList[,"sourceID"], edgeList[,"productID"])
  node.nms <- sapply(node.ids, FUN = function(x){
    cmpd$NM[cmpd$ID == x]
  });

  ## enhance annotateion by matching formula to KEGG
  node.ids_enhanced <- enhanceKBAnnot(table.nm);

  edgeu.res <<- rbind(edgeu.res, edge.res) %>% distinct();
  nodeu.ids <<- c(nodeu.ids, node.ids, node.ids_enhanced$node.id);
  nodeu.nms <<- c(nodeu.nms, node.nms, node.ids_enhanced$node.nm);

  if(netInv == "direct"){
    if("peak" %in% dataSet$type){
      ## kncmpd is abbreviated from "knowledged compound"
      not.kncmpd.ids <- c(net.info$peak.ids,
                          net.info$met.ids)
      net.info$kncmpd.ids <- c(net.info$kncmpd.ids,
                               node.ids[!node.ids %in% not.kncmpd.ids])
      #net.info$int.ids <- net.info$kncmpd.ids
    } else {
      net.info$met.ids <- unique(protein.vec)
      net.info$kncmpd.ids <- c(net.info$kncmpd.ids,
                               node.ids[!node.ids %in% protein.vec])
    }
  } else {
    net.info$met.ids <- node.ids[!node.ids %in% protein.vec]
    if(!zero){
      net.info$kncmpd.ids <- c(net.info$kncmpd.ids, unique(protein.vec))
    }
  }

  return(list(edge.res = edge.res, net.info = net.info))
}

enhanceKBAnnot <- function(table.nm) {
  ## This function is developed to enhance the annotation by
  ## leveraging the knowledge from Database
  if(exists("edgeu.res",envir = .GlobalEnv)) {
    edgeu.res <- get("edgeu.res", envir = .GlobalEnv)
  } else {
    return(0)
  }
  cmpdDB <- qs::qread("../../data/lib/hmdb_lib.qs")

  enh.idx <- apply(edgeu.res, 1, FUN = function(x) {
    grepl(pattern = "(C|G)[0-9][0-9][0-9][0-9][0-9]", x = x[1]) &
      !grepl(pattern = "(C|G)[0-9][0-9][0-9][0-9][0-9]", x = x[2])
  })
  edge2en <- edgeu.res[enh.idx,] #edge2en means: edge to enhance

  res <- QueryM2mSQLiteNet(table.nm = table.nm,
                           q.vec = unique(unique(edge2en$Source)),
                           inv = "direct")
  edge2en_new <- apply(edge2en, 1, FUN = function(x) {
    sourID <- as.character(x[1]);
    targFM <- as.character(x[2]);
    subres <- res[res$sourceID == sourID,]
    if(nrow(subres) == 0){return(x)}
    kgFM <- cmpdDB[sapply(unique(subres$productID), FUN = function(x) {which(x == cmpdDB$KEGG)[1]}), 6]

    matchRes <- which(kgFM == targFM)
    if(length(matchRes) > 0) {
      x[2] <- unique(subres$productID)[matchRes[1]]
    }
    return(x);
  })
  t(edge2en_new) -> edgeu.res[enh.idx,]
  edgeu.res <<- edgeu.res;
  nd.id <- edge2en_new[2,][grepl(pattern = "(C|G)[0-9][0-9][0-9][0-9][0-9]",
                                 x = edge2en_new[2,])]
  nd.nm <- unlist(sapply(as.character(nd.id), FUN = function(x){
    res$productNM[res$productID == x]
  }));
  return(list(node.id = unname(nd.id),
              node.nm = unname(nd.nm)))
}

redundancyClean <- function(node.res, edge.res) {
  ## This function is designed to clean the expanded "metabolite-metabolite" network
  ## by removing the nodes which only connect with the expanded nodes (degree = 2) [case 1];
  ## or the expanded ones but with degree = 1 [case 2];
  if(exists("net.info",envir = .GlobalEnv)) {
    net.info <- get("net.info", envir = .GlobalEnv)
  } else {
    return(0)
  }

  require("igraph")

  kncmpd.ids <- net.info$kncmpd.ids;
  met.ids <- net.info$met.ids;

  ## case 1
  nonkwdt <- edge.res[apply(edge.res, 1, FUN = function(x) {
    !(x[1] %in% kncmpd.ids) | !(x[2] %in% kncmpd.ids)
  }),];
  nonrm <- c(nonkwdt$Source, nonkwdt$Target);
  kwdt_rmed <- edge.res[!apply(edge.res, 1, FUN = function(x) {
    (x[1] %in% kncmpd.ids) & (x[2] %in% kncmpd.ids) & !(x[1] %in% nonrm) & !(x[2] %in% nonrm)
  }),]

  ## case 2
  temp.graph <- simplify(graph.data.frame(kwdt_rmed, directed=FALSE))
  res_degree <- degree(temp.graph)

  # a kncmpd with degree = 1 should be removed
  node2rm1 <- names(res_degree)[(names(res_degree) %in% kncmpd.ids) & (res_degree == 1)];
  # a kncmpd (degree = 2) linked to a kncmp with degree = 1 should also be removed
  node2rm2 <- names(res_degree)[(names(res_degree) %in% kncmpd.ids) &
                                  sapply(names(res_degree), function(x) {
                                    r1 <- kwdt_rmed[(kwdt_rmed$Source == x | kwdt_rmed$Target == x),]
                                    r2 <- apply(r1, 1, FUN = function(r){
                                      r[1] %in% node2rm1 | r[2] %in% node2rm1
                                    })
                                    return(length(which(r2)) > 0)
                                  }) &
                                  (res_degree >= 2)];
  node2rm <- c(node2rm1, node2rm2)
  cleandt <- kwdt_rmed[!apply(kwdt_rmed, 1, FUN = function(x){
    (x[1] %in% node2rm) | (x[2] %in% node2rm)
  }),]

  ## Organize
  allnodes <- c(cleandt$Source, cleandt$Target)
  node.res <- node.res[node.res$Id %in% allnodes,]
  net.info[["kncmpd.ids"]] <- net.info[["kncmpd.ids"]][net.info[["kncmpd.ids"]] %in% allnodes]

  net.info <<- net.info;
  return(list(node.res = node.res, edge.res = cleandt))

}

# elem_table = read_csv("elem_table.csv")
# save(elem_table, file="elem_table.rda")
# sinew::makeOxygen(elem_table)
#' @title elem_table
#' @description elem_table
#' @format A data frame with 288 rows and 7 variables:
#' \describe{
#'   \item{\code{element}}{character COLUMN_DESCRIPTION}
#'   \item{\code{isotope}}{character COLUMN_DESCRIPTION}
#'   \item{\code{mass}}{double COLUMN_DESCRIPTION}
#'   \item{\code{abundance}}{double COLUMN_DESCRIPTION}
#'   \item{\code{ratioC}}{double COLUMN_DESCRIPTION}
#'   \item{\code{Mass_Dif}}{double COLUMN_DESCRIPTION}
#'   \item{\code{unsaturation}}{double COLUMN_DESCRIPTION}
#'}
#' @details DETAILS
"elem_table"
