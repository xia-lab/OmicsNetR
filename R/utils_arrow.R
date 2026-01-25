##################################################
## R script for OmicsNet
## Description: Arrow utilities for zero-copy data exchange with Java
## Author: OmicsNet Team
## Part of the Rserve/qs to Apache Arrow migration
###################################################

#' Safe column extraction with name-first, index-fallback strategy
#'
#' This function provides a safe way to extract columns from data frames
#' during Arrow migration. It ensures data integrity by:
#' 1. Trying named column access first (preferred)
#' 2. Falling back to index if name doesn't exist (with warning)
#' 3. Returning NA vector with error if both fail
#'
#' @param tab Data frame or matrix to extract from
#' @param name Column name to try first (character)
#' @param idx Fallback column index (1-based integer)
#' @param context Optional context string for logging (e.g., "GetNodeTable:gene")
#' @return Column values as vector, or NA vector if extraction fails
#' @export
safeGetCol <- function(tab, name, idx = NULL, context = "") {
  nrows <- nrow(tab)
  if (is.null(nrows) || nrows == 0) {
    return(character(0))
  }

  # Strategy 1: Try by name first (preferred - most reliable)
  if (!is.null(name) && name %in% colnames(tab)) {
    return(tab[[name]])
  }

  # Strategy 2: Fall back to index (with warning)
  if (!is.null(idx) && is.numeric(idx) && idx > 0 && ncol(tab) >= idx) {
    actualName <- colnames(tab)[idx]
    warning(sprintf("[%s] Column '%s' not found, using index %d (actual column: '%s'). Consider updating column name mapping.",
                    context, name, idx, actualName))
    return(tab[, idx])
  }

  # Strategy 3: Both failed - return NA with error
  warning(sprintf("[%s] FAILED to extract column: name='%s', idx=%s. Available columns: %s",
                  context, name, ifelse(is.null(idx), "NULL", idx),
                  paste(colnames(tab), collapse=", ")))
  return(rep(NA, nrows))
}

#' Safe column extraction for multiple columns at once
#'
#' Extracts multiple columns using safeGetCol, returning a data frame.
#'
#' @param tab Source data frame
#' @param mapping Named list where names are output column names and values are
#'                lists with 'name' (preferred column name) and 'idx' (fallback index)
#' @param context Context string for logging
#' @return Data frame with extracted columns
#' @export
safeGetCols <- function(tab, mapping, context = "") {
  df <- data.frame(row.names = rownames(tab))
  for (outCol in names(mapping)) {
    spec <- mapping[[outCol]]
    df[[outCol]] <- safeGetCol(tab, spec$name, spec$idx, context)
  }
  return(df)
}

#' Validate that required columns exist in a data frame
#'
#' @param tab Data frame to validate
#' @param required Character vector of required column names
#' @param context Context string for error messages
#' @return TRUE if all columns exist, FALSE otherwise (with warnings)
#' @export
validateColumns <- function(tab, required, context = "") {
  missing <- setdiff(required, colnames(tab))
  if (length(missing) > 0) {
    warning(sprintf("[%s] Missing required columns: %s. Available: %s",
                    context, paste(missing, collapse=", "),
                    paste(colnames(tab), collapse=", ")))
    return(FALSE)
  }
  return(TRUE)
}

#' Shadow save function for mixed-type data frames
#'
#' Saves data as both .qs (for R) and .arrow (for Java zero-copy access).
#' ALWAYS preserves rownames as the first column (row_names_id) in Arrow files.
#'
#' @param obj The object to save (matrix or data.frame)
#' @param file The .qs file path (will auto-generate .arrow path)
#' @param compress Compression for Arrow (default: "uncompressed" for memory-mapping)
#' @return Invisible NULL
#' @export
shadow_save_mixed <- function(obj, file, compress = "uncompressed") {
    # Always save to qs for R compatibility
    qs::qsave(obj, file)

    # Generate Arrow path
    arrow_path <- sub("\\.qs$", ".arrow", file)

    tryCatch({
        if (is.matrix(obj) || is.data.frame(obj)) {
            df <- as.data.frame(obj)

            # Preserve rownames as first column
            rn <- rownames(obj)
            if (!is.null(rn) && length(rn) > 0 && !all(rn == as.character(1:nrow(df)))) {
                df <- cbind(row_names_id = as.character(rn), df)
            }

            # Ensure all character columns are properly typed
            for (col in names(df)) {
                if (is.factor(df[[col]])) {
                    df[[col]] <- as.character(df[[col]])
                }
            }

            arrow::write_feather(df, arrow_path, compression = compress)
        }
    }, error = function(e) {
        warning(paste("Arrow shadow save failed:", e$message))
    })

    invisible(NULL)
}

#' Shadow save for simple data frames (no special rownames handling)
#'
#' @param obj The object to save
#' @param file The .qs file path
#' @param compress Compression for Arrow
#' @return Invisible NULL
#' @export
shadow_save <- function(obj, file, compress = "uncompressed") {
    qs::qsave(obj, file)

    arrow_path <- sub("\\.qs$", ".arrow", file)

    tryCatch({
        if (is.data.frame(obj)) {
            # Ensure all character columns are properly typed
            df <- obj
            for (col in names(df)) {
                if (is.factor(df[[col]])) {
                    df[[col]] <- as.character(df[[col]])
                }
            }
            arrow::write_feather(df, arrow_path, compression = compress)
        }
    }, error = function(e) {
        warning(paste("Arrow shadow save failed:", e$message))
    })

    invisible(NULL)
}

#' Save network node table to Arrow format
#'
#' Specialized function for saving network node data for Java visualization.
#'
#' @param nodes Data frame with node data (id, label, type, degree, etc.)
#' @param file_prefix File prefix (will create .qs and .arrow files)
#' @return The Arrow file path
#' @export
save_node_table_arrow <- function(nodes, file_prefix) {
    qs_path <- paste0(file_prefix, ".qs")
    arrow_path <- paste0(file_prefix, ".arrow")

    qs::qsave(nodes, qs_path)

    tryCatch({
        df <- nodes
        # Ensure all columns are properly typed
        for (col in names(df)) {
            if (is.factor(df[[col]])) {
                df[[col]] <- as.character(df[[col]])
            }
        }
        arrow::write_feather(df, arrow_path, compression = "uncompressed")
        return(arrow_path)
    }, error = function(e) {
        warning(paste("Node table Arrow save failed:", e$message))
        return(NULL)
    })
}

#' Save enrichment results to Arrow format
#'
#' Specialized function for saving enrichment analysis results.
#'
#' @param results Data frame with enrichment results
#' @param file_prefix File prefix (will create .qs and .arrow files)
#' @return The Arrow file path
#' @export
save_enrichment_arrow <- function(results, file_prefix) {
    qs_path <- paste0(file_prefix, ".qs")
    arrow_path <- paste0(file_prefix, ".arrow")

    qs::qsave(results, qs_path)

    tryCatch({
        df <- results
        # Preserve rownames if present
        rn <- rownames(results)
        if (!is.null(rn) && length(rn) > 0 && !all(rn == as.character(1:nrow(df)))) {
            df <- cbind(row_names_id = as.character(rn), df)
        }
        # Ensure all columns are properly typed
        for (col in names(df)) {
            if (is.factor(df[[col]])) {
                df[[col]] <- as.character(df[[col]])
            }
        }
        arrow::write_feather(df, arrow_path, compression = "uncompressed")
        return(arrow_path)
    }, error = function(e) {
        warning(paste("Enrichment Arrow save failed:", e$message))
        return(NULL)
    })
}

#' Write data frame directly to Arrow file
#'
#' Simple function to write a data frame to Arrow format.
#'
#' @param df Data frame to write
#' @param arrow_path Path for the Arrow file
#' @param preserve_rownames Whether to preserve rownames as first column
#' @return The Arrow file path if successful, NULL otherwise
#' @export
write_arrow <- function(df, arrow_path, preserve_rownames = TRUE) {
    tryCatch({
        # Preserve rownames if requested
        if (preserve_rownames) {
            rn <- rownames(df)
            if (!is.null(rn) && length(rn) > 0 && !all(rn == as.character(1:nrow(df)))) {
                df <- cbind(row_names_id = as.character(rn), df)
            }
        }

        # Convert factors to characters
        for (col in names(df)) {
            if (is.factor(df[[col]])) {
                df[[col]] <- as.character(df[[col]])
            }
        }

        arrow::write_feather(df, arrow_path, compression = "uncompressed")
        return(arrow_path)
    }, error = function(e) {
        warning(paste("Arrow write failed:", e$message))
        return(NULL)
    })
}
