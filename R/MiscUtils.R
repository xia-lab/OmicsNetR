##################################################
## R scripts for OmicsNet
## Various utility methods
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# new range [a, b]
rescale2NewRange <- function(qvec, a, b){
    qvec = replace(qvec, qvec == 0, 1)
    q.min <- min(qvec);
    q.max <- max(qvec);
    if(length(qvec) < 50){
        a <- a*2;
    }
    if(q.max == q.min){
        new.vec <- rep(8, length(qvec));
    }else{
        coef.a <- (b-a)/(q.max-q.min);
        const.b <- b - coef.a*q.max;
        new.vec <- coef.a*qvec + const.b;
    }
    return(new.vec);
}

# #FFFFFF to rgb(1, 0, 0)
hex2rgba <- function(cols){
  return(apply(sapply(cols, col2rgb), 2, function(x){paste("rgba(", x[1], ",", x[2], ",", x[3], ",0.8)", sep="")}));
}

# re-arrange one vector elements according to another vector values
# usually src is character vector to be arranged
# target is numberic vector of same length

# given a data with duplicates, dups is the one with duplicates
RemoveDuplicates <- function(data, lvlOpt, quiet=T){

    all.nms <- rownames(data);
    colnms <- colnames(data);
    dup.inx <- duplicated(all.nms);
    dim.orig  <- dim(data);
    data <- suppressWarnings(apply(data, 2, as.numeric)); # force to be all numeric
    dim(data) <- dim.orig; # keep dimension (will lost when only one item)
    rownames(data) <- all.nms;
    colnames(data) <- colnms;
    if(sum(dup.inx) > 0){
        uniq.nms <- all.nms[!dup.inx];
        uniq.data <- data[!dup.inx,,drop=F];

        dup.nms <- all.nms[dup.inx];
        uniq.dupnms <- unique(dup.nms);
        uniq.duplen <- length(uniq.dupnms);

        for(i in 1:uniq.duplen){
            nm <- uniq.dupnms[i];
            hit.inx.all <- which(all.nms == nm);
            hit.inx.uniq <- which(uniq.nms == nm);

            # average the whole sub matrix
            if(lvlOpt == "mean"){
                uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, mean, na.rm=T);
            }else if(lvlOpt == "median"){
                uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, median, na.rm=T);
            }else if(lvlOpt == "max"){
                uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, max, na.rm=T);
            }else{ # sum
                uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, sum, na.rm=T);
            }
        }
        if(!quiet){
            current.msg <<- paste(current.msg, paste("A total of ", sum(dup.inx), " of duplicates were replaced by their ", lvlOpt, ".", sep=""), collapse="\n");
        }
        return(uniq.data);
    }else{
        if(!quiet){
            current.msg <<- paste(current.msg, "All IDs are unique.", collapse="\n");
        }
        return(data);
    }
}

generate_breaks = function(x, n, center = F){
    if(center){
        m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))));
        res = seq(-m, m, length.out = n + 1);
    } else {
        res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1);
    }
    return(res)
}

ComputeColorGradient <- function(nd.vec, background = "black", centered = FALSE) {
  require("RColorBrewer")
  
  # Check if the input is a list of vectors or a vector
  if (is.list(nd.vec) && all(sapply(nd.vec, is.numeric))) {
    # Flatten the list of vectors into a single numerical vector
    flat_vec <- unlist(nd.vec)
  } else if (is.numeric(nd.vec)) {
    flat_vec <- nd.vec
  } else {
    stop("Input must be a vector or a list of numeric vectors")
  }
  
  # Compute min and max values
  minval <- min(flat_vec, na.rm = TRUE)
  maxval <- max(flat_vec, na.rm = TRUE)
  res <- maxval - minval
  
  # Handle the case where all values are the same
  if (res == 0) {
    return(rep("#0b6623", length(flat_vec)))
  }
  
  # Generate the color gradient
  color <- GetColorGradient(background, centered)
  breaks <- generate_breaks(flat_vec, length(color), center = centered)
  color_vec <- scale_vec_colours(flat_vec, col = color, breaks = breaks)
  
  # If input was a list of vectors, map the colors back to the original list structure
  if (is.list(nd.vec) && all(sapply(nd.vec, is.numeric))) {
    color_list <- split(color_vec, ceiling(seq_along(color_vec) / 2))
    return(color_list)
  } else {
    return(color_vec)
  }
}


GetColorGradient <- function(background, center){
    if(background == "black"){
        if(center){
            return(c(colorRampPalette(c("#31A231", "#5BC85B", "#90EE90", "#C1FFC1"))(50), colorRampPalette(c("#FF9DA6", "#FF7783", "#E32636", "#BD0313"))(50)));
        }else{
            return(c(colorRampPalette(c("#edf8fb","#b2e2e2","#66c2a4","#2ca25f","#006d2c"))(100)));
        }
    }else{ # white background
        if(center){
            return(c(colorRampPalette(c("#137B13", "#31A231", "#5BC85B", "#90EE90"))(50), colorRampPalette(c("#FF7783", "#E32636", "#BD0313", "#96000D"))(50)));
        }else{
            return(colorRampPalette(hsv(h = seq(0.72, 1, 0.035), s = 0.72, v = 1))(100));
        }
    }
}

scale_vec_colours = function(x, col = rainbow(10), breaks = NA){
    breaks <- sort(unique(breaks));
    return(col[as.numeric(cut(x, breaks = breaks, include.lowest = T))]);
}

scale_colours = function(mat, col = rainbow(10), breaks = NA){
    mat = as.matrix(mat)
    return(matrix(scale_vec_colours(as.vector(mat), col = col, breaks = breaks), nrow(mat), ncol(mat), dimnames = list(rownames(mat), colnames(mat))))
}

###########
# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
                           capture.output(format(utils::object.size(x), units = "auto")) })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}

# shorthand
ShowMemoryUse <- function(..., n=30) {
    require("pryr");
    sink(); # make sure print to screen
    print(pryr::mem_used());
    print(sessionInfo());
    print(.ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n));
    print(warnings());
}

CleanMemory <- function(){
   gc();
}

# need to obtain the full path to convert (from imagemagik) for cropping images
GetBashFullPath<-function(){
  path <- system("which bash", intern=TRUE);
  if((length(path) == 0) && (typeof(path) == "character")){
    print("Could not find bash in the PATH!");
    return("NA");
  }
  return(path);
}

# private method for all sqlite queries
.query.sqlite <- function(db.con, statement, offline=TRUE){
  rs <- dbSendQuery(db.con, statement);
  res <- fetch(rs, n=-1); # get all records
  dbClearResult(rs);
  if(offline){
    dbDisconnect(db.con);
  }
  cleanMem();
  return(res);
}

fast.write.csv <- function(dat, file, row.names=TRUE){
    tryCatch(
        {
           if(is.data.frame(dat)){
                # there is a rare bug in data.table (R 3.6) which kill the R process in some cases
                data.table::fwrite(dat, file, row.names=row.names);
           }else{
                write.csv(dat, file, row.names=row.names);
           }
        }, error=function(e){
            print(e);
            write.csv(dat, file, row.names=row.names);
        }, warning=function(w){
            print(w);
            write.csv(dat, file, row.names=row.names);
        });
}

cleanMem <- function(n=10) { for (i in 1:n) gc() }

.set.nSet <- function(dataSetObj=NA){
  dataSet <<- dataSetObj;
  if(!exists(".on.public.web")){
    .on.public.web <<- T;
  }
  if(.on.public.web){
    return (1);
  }else{
    return(dataSetObj);
  }
}

.get.nSet <- function(dataSetObj=NA){
  if(.on.public.web){
    return(dataSet)
  }else{
    return(dataSetObj);
  }
}

ReadList <- function(dataSetObj=NA, fullPath, fileNm){
    fullUrl <- url(paste0(fullPath,"/", fileNm))
    all_str <- paste0(readLines(fullUrl),collapse="\n");
    return(all_str);
}


makeReadable <- function(str){
    result <- switch(str,
                 pct = "Percent",
                 abs = "Absolute",
                 log = "Log2",
                 rle = "RLE",
                 array = "Microarray",
                 count= "RNA-Seq",
                 hsa = "H. sapiens (human)",
                 mmu = "M. musculus (mouse)",
                 rno = "R. norvegicus (rat)",
                 cel = "C. elegans (roundworm)",
                 dme = "D. melanogaster (fruitfly)",
                 dre = "D. rerio (zebrafish)",
                 sce = "S. cerevisiae (yeast)",
                 eco = "E. coli",
                 ath = "A. thaliana (Arabidopsis)",
                 bta = "B. taurus (cow)",
                 gga = "G. gallus (chicken)",
                 mun = "M. unguiculatus (Mongolian gerbil)",
                 bsu = "B. subtilis",
                 pae = "P. aeruginosa",
                 mtb = "M. tuberculosis",
                 smm = "S. mansoni (schistosomiasis)",
                 tbr = "T. brucei (trypanosoma)",
                 pfa = "P. falciparum (malaria)",
                 cjo = "C. japonica (japanese quail)",
                 xla = "X. laevis (African clawed frog)",
                 ppr = "P. promelas (fathead minnow; custom)",
                 fhm = "P. promelas (fathead minnow; NCBI)",
                 nlf = "L. pipiens (northern leopard frog)",
                 omk = "O. mykiss (rainbow trout)",
                 ham = "H. americanus (American lobster)",
                 cdi = "C. dilutus",
                 dma = "D. magna",
                 rsu = "R. subcapitata",
                 haz = "H. azteca",
                 fcd = "F. candida",
                 "entrez" = "Entrez ID",
                 "refseq" = "RefSeq ID",
                   "gb" = "Genbank ID",
                   "symbol" = "Official Gene Symbol",
                   "embl_gene" = "Ensembl Gene ID",
                   "embl_transcript" = "Ensemble Transcript ID",
                   "embl_protein" = "Ensembl Protein ID",
                   "uniprot" = "Uniprot Accession ID",
                   "hgu95a" = "Affymetrix Human Genome U95 (chip hgu95a)",
                   "hgu95av2" = "Affymetrix Human Genome U95 (chip hgu95av2)",
                   "hgu95b" = "Affymetrix Human Genome U95 (chip hgu95b)",
                   "hgu95c" = "Affymetrix Human Genome U95 (chip hgu95c)",
                   "hgu95d" = "Affymetrix Human Genome U95 (chip hgu95d)",
                   "hgu95e" = "Affymetrix Human Genome U95 (chip hgu95e)",
                   "hgu133a" = "Affymetrix Human Genome U133 (chip hgu133a)",
                   "hgu133b" = "Affymetrix Human Genome U133 (chip hgu133b)",
                   "hgu133plus2" = "Affymetrix Human Genome U133plus2 (hgu133plus2)",
                   "hgu133plus2pm" = "Affymetrix Human Genome U133plus2_PM (hgu133plus2pm)",
                   "lumiht12v3" = "Illumina HumanHT-12 V3 BeadArray",
                   "lumiht12v4" = "Illumina HumanHT-12 V4 BeadArray",
                   "lumiref8v2" = "Illumina HumanRef-8 V2 BeadArray",
                   "lumiref8v3" = "Illumina HumanRef-8 V3 BeadArray",
                   "lumiwg6v2" = "Illumina HumanWG-6 V2 BeadArray",
                   "lumiwg6v3" = "Illumina HumanWG-6 V3 BeadArray",
                   "agi4100a" = "Agilent Human 1 cDNA Microarray (4100A)",
                   "agi4101a" = "Agilent Human 2 cDNA Microarray (4101A)",
                   "agi4110b" = "Agilent Human 1A cDNA Microarray (4110B)",
                   "agi4111a" = "Agilent Human 1B cDNA Microarray (4111A)",
                   "agi4112a" = "Agilent Human Genome Whole Microarray (4x44k/4112)",
                   "agi4845a" = "Agilent Human AMADID 026652 Microarray (4845A)",
                   "lumiwg6v1" = "Illumina MouseWG-6 v1.0 Bead Array",
                   "lumiwg6v11" = "Illumina MouseWG-6 v1.1 Bead Array",
                   "lumiwg6v2" = "Illumina MouseWG-6 v2.0 Bead Array",
                   "lumiref8v1" = "Illumina MouseRef-8 v1.0 Bead Array",
                   "lumiref8v2" = "Illumina MouseRef-8 v2.0 Bead Array",
                   "mgu74a" = "Affymetrix Murine Genome U74v2 (chip mgu74a)",
                   "mgu74av2" = "Affymetrix Murine Genome U74v2 (chip mgu74av2)",
                   "mgu74b" = "Affymetrix Murine Genome U74v2 (chip mgu74b)",
                   "mgu74bv2" = "Affymetrix Murine Genome U74v2 (chip mgu74bv2)",
                   "mgu74c" = "Affymetrix Murine Genome U74v2 (chip mgu74c)",
                   "mgu74cv2" = "Affymetrix Murine Genome U74v2 (chip mgu74cv2)",
                   "moe430a" = "Affymetrix Mouse Expression Set 430 (chip moe430a)",
                   "moe430b" = "Affymetrix Mouse Expression Set 430 (chip moe430b)",
                   "moe430_2" = "Affymetrix GeneChip Mouse Genome 430 2.0",
                   "mgi_st1" = "Affymetrix Mouse Gene 1.0 ST Array",
                   "mgu4101a" = "Agilent Mouse Array (chip mgug4104a)",
                   "mgu4120a" = "Agilent Mouse Array (chip mgug4120a)",
                   "mgu4121a" = "Agilent Mouse Array (chip mgug4121a)",
                   "mgu4122a" = "Agilent Mouse Array (chip mgug4122a)",
                   "kegg" = "KEGG",
                   "keggp" = "KEGG species specific",
                    "reactome" = "Reactome",
                    "go_bp" = "GO:BP",
                    "go_mf" = "GO:MF",
                    "go_cc" = "GO:CC",
                    "panth" = "PANTHER Slim",
                    "motif_set" = "Motif",
                    "trrust" = "TRRUST",
                    "encode" = "ENCODE",
                    "chea" = "CHEA",
                    "jaspar" = "JASPAR",
                    "innate" = "InnateDB",
                    "string" = "StringDB",
                    "intact" = "IntAct",
                    "huri" = "HuRI",
                    "mirtarbase" = "miRTarBase",
                    "tarbase" = "TarBase",
                    "mirecords" = "miRecords",
                    "mirtarbase" = "miRTarBase",
                    "tarbase" = "TarBase",
                    "recon2" = "Recon2",
                    "recon3" = "Recon3D",
                    "agora" = "AGORA",
                    "embl" ="EMBL",
                    "vep" = "VEP",
                    "admire" = "ADmiRE",
                    "snp2tfbs" = "SNP2TFBS",
                    "hmdb" = "HMDB",
                    "pubchem" = "PubChem",
                 str)
}


GetNodeMat <- function(){
  if(is.null(dataSet$imgSet$node_table)){
    df <- .readDataTable('node_table.csv')
    df[,-c(1:2)] <- lapply(df[,-c(1:2)], function(col) as.numeric(as.character(col)))
    dataSet$imgSet$node_table <<- df;
  }
  return(as.matrix(dataSet$imgSet$node_table[,-c(1:2)]))  # ensure matrix of numerics
}

GetNodeRowNames <- function(){
  if(is.null(dataSet$imgSet$node_table)){
  df <- .readDataTable('node_table.csv')
    dataSet$imgSet$node_table <<- df;

  }
  dataSet$imgSet$node_table$Id;
}

GetNodeGeneSymbols <- function(){
  if(is.null(dataSet$imgSet$node_table)){
  df <- .readDataTable('node_table.csv')
    dataSet$imgSet$node_table <<- df;

  }
  dataSet$imgSet$node_table$Label;
}

GetNodeColNames <- function(){
  if(is.null(dataSet$imgSet$node_table)){
  df <- .readDataTable('node_table.csv')
    dataSet$imgSet$node_table <<- df;

  }
  return(colnames(dataSet$imgSet$node_table[,-c(1:2)]));

}

CleanNumber <-function(bdata){
  if(sum(bdata==Inf)>0){
    inx <- bdata == Inf;
    bdata[inx] <- NA;
    bdata[inx] <- 999999;
  }
  if(sum(bdata==-Inf)>0){
    inx <- bdata == -Inf;
    bdata[inx] <- NA;
    bdata[inx] <- -999999;
  }
  bdata;
}

CheckDetailsTablePerformed <-function(type){
  performed <- T;
  if(type == "node"){
    performed <- file.exists("node_table.csv");
  }else if(type %in% c( "network_enr", "regNetwork_enr", "gba_enr", "module_enr", "defaultEnr")){
    clean_type <- gsub("_enr", "", type);
    performed <- !is.null(dataSet$imgSet$enrTables[[clean_type]]);
  }else if(type %in% c( "peak_anot")){
    performed <- file.exists("peak_annotation.csv")
  }
  print(paste("checkPerformed=", type, "====",performed));

return(performed)
}


GetPeakAnnotMat <- function(){
  if(is.null(dataSet$imgSet$peak_annotation)){
    df <- .readDataTable('peak_annotation.csv')
    dataSet$imgSet$peak_annotation <<- df;
  }
  df <- dataSet$imgSet$peak_annotation
cols_to_exclude <- c("class", "formula", "annotation")
cols_to_convert <- setdiff(names(df), cols_to_exclude)

df[cols_to_convert] <- lapply(df[cols_to_convert], function(col) as.numeric(as.character(col)))
  return(as.matrix(df[cols_to_convert]))  # ensure matrix of numerics
}

GetPeakAnnotRowNames <- function(){
  if(is.null(dataSet$imgSet$peak_annotation)){
  df <- .readDataTable('peak_annotation.csv')
    dataSet$imgSet$peak_annotation <<- df;

  }
  dataSet$imgSet$peak_annotation$annotation;
}

GetPeakAnnotClass <- function(){
  if(is.null(dataSet$imgSet$peak_annotation)){
  df <- .readDataTable('peak_annotation.csv')
    dataSet$imgSet$node_table <<- df;

  }
  dataSet$imgSet$peak_annotation$class;
}


GetPeakAnnotColNames <- function(){
  if(is.null(dataSet$imgSet$peak_annotation)){
  df <- .readDataTable('peak_annotation.csv')
    dataSet$imgSet$peak_annotation <<- df;

  }
  df <- dataSet$imgSet$peak_annotation
cols_to_exclude <- c("class", "formula", "annotation")
cols_to_convert <- setdiff(names(df), cols_to_exclude)

  return(colnames(dataSet$imgSet$peak_annotation[cols_to_convert]));

}