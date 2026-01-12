  
my.snp.query <- function(dataSet, q.vec, db.type, netInv, zero=FALSE, snpRegion=FALSE){
    require('RSQLite');
    if(db.type == "PhenoScanner"){

      require("tidyr");
      if(dataSet$snp2gene$phesc.opt=="eqtl"){

        res <- Query.PhenoScanner(snpquery= q.vec,catalogue="eQTL");
        fast.write.csv(res$results, file="phenoscanner.eQtl.csv",row.names=FALSE);
        res <- res$results[,c("rsid","exp_gene")];
        res <- unique(res[which(!(is.na(res$exp_gene)) &res$exp_gene !=""&res$exp_gene !="-" ),]);
        colnames(res)[2] = "symbol";
        res <- unique(data.frame(tidyr::separate_rows(res,symbol,sep =";"),stringsAsFactors = F));

      }else{

        res <- Query.PhenoScanner(snpquery= q.vec);
        res <- res$snps;
        fast.write.csv(res, file="phenoscanner.nearest.csv",row.names=FALSE);
        colnames(res)[which(colnames(res)=="hgnc")] = "symbol";
        res <- unique(data.frame(tidyr::separate_rows(res,symbol,sep =";"),stringsAsFactors = F))
      }

      res$entrez = doGeneIDMapping(res$symbol,type="symbol");
      res$entrez[is.na(res$entrez)] = res$symbol[is.na(res$entrez)];

    }else if(db.type == "vep"){

      if(dataSet$snp2gene$vep.opt=="dis"){
        res <- QueryVEP(q.vec, vepDis=dataSet$snp2gene$vep.dis, snpRegion = snpRegion);
        fast.write.csv(res, file="vep.map.csv",row.names=FALSE);
        res <-  unique(res[which(res$gene_symbol !="NA"),c("rsid","gene_symbol")]);
      }else{
        require("dplyr")
        res<- QueryVEP(q.vec, vepDis=50,snpRegion = snpRegion)
        fast.write.csv(res, file="vep.map.csv",row.names=FALSE);
        res$distance[which(res$distance=="NA")] = 0;
        res <- unique(res[!(is.na(res$gene_symbol)),c("gene_symbol","rsid","distance")]);
        res <- unique(res[!(duplicated(res$gene_symbol,res$rsid)),]);
        res <- data.frame(res %>% arrange(distance) %>%
                           group_by(rsid) %>%
                           mutate(rank = rank(distance)),stringsAsFactors=F);
        res <- res[which(res$rank<(dataSet$snp2gene$vep.num + 1)),];
      }

      res$entrez <- doGeneIDMapping(res$gene_symbol,type="symbol")
      res$entrez[is.na(res$entrez)] <- res$gene_symbol[is.na(res$entrez)]

    }else{
      if(db.type == "admire"){
        file.nm <- "snp2mir";
      }else{
        file.nm <- "snp2tfbs";
      }
      db.path <- paste(sqlite.path, file.nm, sep="");
      table.nm <- "hsa";
      col.nm <- "rsid";
      res <- Query.snpDB(db.path, q.vec, table.nm, col.nm);
    }

    if(db.type == "admire"){
      targetColNm <- "MIRNA_Name";
    }else if(db.type == "PhenoScanner"){
      targetColNm <- "entrez";
    }else if(db.type == "vep"){
      targetColNm <- "entrez";
    }else{
      targetColNm <- "entrez";
    }
    colInx <- which(colnames(res) == targetColNm);

    snp.list <- list();
    snp.type.list <- list();
    snpTable <<- snp.list;

    edge.res <- data.frame(Source=res[,"rsid"],Target=res[,colInx], stringsAsFactors=FALSE);

    # Filter out rows with NA values in Source or Target columns
    edge.res <- edge.res[complete.cases(edge.res), ];

    if(nrow(edge.res)!=0){
      row.names(edge.res) <- 1:nrow(edge.res);
    }
    fast.write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    node.ids <- c(edge.res[,"Source"], edge.res[,"Target"])

    #  print(node.ids)
    if(db.type == "admire"){
      symb <- res[,"MIRNA_Acc"];
    }else if(db.type == "PhenoScanner"){
      symb <- res[,"symbol"];
    }else if(db.type == "vep"){
      symb <- res[,"gene_symbol"];
    }else{
      symb <- doEntrez2SymbolMapping(res[,"entrez"]);
    }
    node.nms <- c(res[,"rsid"], symb);

    if(netInv == "direct"){
      targetIds <- unique(c(net.info$gene.ids, node.ids[!node.ids %in% q.vec]));
      net.info$snp.ids <- unique(q.vec);
    }else{
      net.info$snp.ids <- unique(res[,"rsid"]);
      if(!zero){
        targetIds <- c(unique(q.vec));
      }
    }

    if(db.type == "admire"){
      net.info$mir.ids <- targetIds;
    }else if(db.type == "snp2tfbs"){
      net.info$tf.ids <- targetIds;
    }else if(db.type == "PhenoScanner"){
      net.info$gene.ids <- targetIds;
    }else{
      net.info$gene.ids <- targetIds;
    }

    net.info$snpi.ids <- unique(node.ids);
    
    res <- list(
        edge.res = edge.res,
        node.ids = node.ids,
        node.nms = node.nms,
        net.info= net.info
    );
    return(res);
}