##################################################
## Compound Name Approximate Matching Utilities
## Adapted from MetaboAnalyst for OmicsNet
## Performs fuzzy matching using synonym database
##################################################

#' Perform approximate compound name matching using synonym database
#' @description Given a query compound name, performs fuzzy matching against
#' a synonym database to find the best matches
#' @param q Input query compound name
#' @param cmpd.db Compound database (from compound_db.rds)
#' @param syn.db Synonym database (from syn_nms.qs)
#' @return A data frame with matched candidates (index, value, score) or NULL if no matches
#' @export
PerformCompoundApproxMatch <- function(q, cmpd.db, syn.db) {

  # Check if databases match in size
  if(length(syn.db$syns.vec) != nrow(cmpd.db)) {
    # Synonym database size doesn't match - use subset that matches
    min.size <- min(length(syn.db$syns.vec), nrow(cmpd.db));
    com.nms <- cmpd.db$name[1:min.size];
    syns.vec <- syn.db$syns.vec[1:min.size];
    syns.list <- syn.db$syns.list[1:min.size];
  } else {
    # Filter out lipids if lipid column exists
    if("lipid" %in% colnames(cmpd.db)) {
      nonLipidInx <- cmpd.db$lipid == 0;
      com.nms <- cmpd.db$name[nonLipidInx];
      syns.vec <- syn.db$syns.vec[nonLipidInx];
      syns.list <- syn.db$syns.list[nonLipidInx];
    } else {
      # No lipid column - use all entries
      com.nms <- cmpd.db$name;
      syns.vec <- syn.db$syns.vec;
      syns.list <- syn.db$syns.list;
    }
  }

  q.length <- nchar(q);
  s <- c(0, 0.1, 0.2);  # tolerance levels for fuzzy matching

  candidates <- NULL;

  # Try increasingly fuzzy matching
  for (j in s) {
    new.q <- q;
    if(q.length > 32){ # agrep fails for exact match when length over 32 characters
      new.q <- substr(q, 1, 32);
    }

    # Fuzzy match against concatenated synonym strings
    matched.inx <- agrep(new.q, syns.vec, ignore.case=TRUE, max.distance=j, useBytes=TRUE);

    if(length(matched.inx) > 0) {
      # Found matches, create candidates data frame
      candidates <- data.frame(
        index = vector(mode = "numeric", length=length(matched.inx)),
        value = vector(mode = "character", length=length(matched.inx)),
        score = vector(mode = "numeric", length=length(matched.inx)),
        stringsAsFactors = FALSE
      );

      for(n in 1:length(matched.inx)){
        nm.vec <- syns.list[[matched.inx[n]]];

        # Try approximate match against individual synonyms
        hit3.inx <- agrep(q, nm.vec, ignore.case=TRUE, max.distance=j, useBytes=TRUE);

        if(length(hit3.inx) > 0){
          hit3.nm <- vector(mode = "character", length=length(hit3.inx));
          hit3.score <- vector(mode = "numeric", length=length(hit3.inx));

          for(k in 1:length(hit3.inx)){
            idx <- hit3.inx[k];
            hit3.nm[k] <- nm.vec[idx];
            # Score based on edit distance and length difference
            hit3.score[k] <- j + abs(nchar(nm.vec[idx]) - nchar(q)) / (10 * nchar(q));
          }

          # Bonus for matching first 1-2 characters
          matches2 <- c();
          if(length(grep("^[1-9a-z]{2}", q, ignore.case=TRUE)) > 0){
            matches2 <- grep(paste("^", substr(q, 1, 2), sep=""), hit3.nm, ignore.case=TRUE);
          } else if (length(grep("^[1-9a-z]", q, ignore.case=TRUE)) > 0){
            matches2 <- grep(paste("^", substr(q, 1, 1), sep=""), hit3.nm, ignore.case=TRUE);
          }

          if(length(matches2) > 0){
            hit3.score[matches2] <- hit3.score[matches2] - 0.05;
          }

          # Select best match for this compound
          best.inx <- which(hit3.score == min(hit3.score))[1];
          candidates[n, 1] <- matched.inx[n];
          candidates[n, 2] <- com.nms[matched.inx[n]]; # show common name
          candidates[n, 3] <- hit3.score[best.inx];
        }
      }

      # Remove NA entries
      rm.inx <- is.na(candidates[,2]) | candidates[,2] == "NA" | candidates[,2] == "";
      candidates <- candidates[!rm.inx, , drop=FALSE];

      # Sort by score (best matches first)
      candidates <- candidates[order(candidates[,3], decreasing=FALSE), , drop=FALSE];

      # Limit to top 10 matches
      if(nrow(candidates) > 10){
        candidates <- candidates[1:10, , drop=FALSE];
      }

      # If we found matches, return them
      if(nrow(candidates) > 0) {
        return(candidates);
      }
    }
  }

  # No matches found
  return(NULL);
}

#' Match compound name with fuzzy matching
#' @description Wrapper function that performs compound name matching
#' @param query.vec Vector of query compound names
#' @param cmpd.db Compound database
#' @param syn.db Synonym database
#' @return List with hit.inx (indices), hit.values (matched names), and match.state (1=match, 0=no match)
#' @export
MatchCompoundNames <- function(query.vec, cmpd.db, syn.db) {

  n <- length(query.vec);
  hit.inx <- rep(0, n);
  hit.values <- rep("", n);
  match.state <- rep(0, n);

  # Determine if we need to adjust indices for database size mismatch
  use.subset <- FALSE;
  subset.size <- nrow(cmpd.db);

  if(length(syn.db$syns.vec) != nrow(cmpd.db)) {
    use.subset <- TRUE;
    subset.size <- min(length(syn.db$syns.vec), nrow(cmpd.db));
  }

  for(i in 1:n) {
    q <- query.vec[i];

    # Try approximate matching with synonyms
    candidates <- PerformCompoundApproxMatch(q, cmpd.db, syn.db);

    if(!is.null(candidates) && nrow(candidates) > 0) {
      # Take the best match (first row after sorting)
      best.idx <- candidates[1, 1];  # index in working database

      # Validate index is within bounds
      if(best.idx > 0 && best.idx <= nrow(cmpd.db)) {
        hit.inx[i] <- best.idx;
        hit.values[i] <- cmpd.db$name[best.idx];
        match.state[i] <- 1;
      }
    }
  }

  return(list(
    hit.inx = hit.inx,
    hit.values = hit.values,
    match.state = match.state
  ));
}
