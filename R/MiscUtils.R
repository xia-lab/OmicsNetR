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
        m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
        res = seq(-m, m, length.out = n + 1)
    }
    else{
        res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
    }
    return(res)
}

ComputeColorGradient <- function(nd.vec, background="black", centered){
    require("RColorBrewer");
    minval = min(nd.vec, na.rm=TRUE);
    maxval = max(nd.vec, na.rm=TRUE);
    res = maxval-minval;

    if(res == 0){
        return(rep("#0b6623", length(nd.vec)));
    }

    #if(sum(nd.vec<0, na.rm=TRUE) > 0){
    #    centered <- T;
    #}else{
    #    centered <- F;
    #}
    color <- GetColorGradient(background, centered);
    breaks <- generate_breaks(nd.vec, length(color), center = centered);
    return(scale_vec_colours(nd.vec, col = color, breaks = breaks));
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
    for (i in 1:10){
        gc(reset = T);
    }
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
