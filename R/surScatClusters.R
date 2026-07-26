## extractClusters ----
# Extract cluster assignments from various clustering objects
# Supports: kmeans, poLCA, tidyLPA, hclust (via cutree), or direct vector/factor
extractClusters <- function(x, sourceData=NULL) {
  # If already a vector or factor, return as-is
  if(is.vector(x) || is.factor(x)) {
    return(as.vector(as.integer(as.factor(as.character(x)))))
  }

  # kmeans object
  if(inherits(x, "kmeans")) {
    return(as.vector(x$cluster))
  }

  # poLCA object
  if(inherits(x, "poLCA")) {
    if(!is.null(x$predclass)) return(as.vector(x$predclass))
    stop("poLCA object must have predclass attribute")
  }

  # tidyLPA result
  if(inherits(x, "tidyLPA")) {
    data <- tidyLPA::get_data(x)
    if(!is.null(data$Class)) return(as.vector(as.integer(data$Class)))
    stop("tidyLPA object must have Class column in returned data")
  }

  # hclust result (should be cut first, but allow passing pre-cut vector)
  if(inherits(x, "hclust")) {
    stop("hclust object should be cut first with cutree(). Pass the result of cutree() instead.")
  }

  # data.frame with single column
  if(is.data.frame(x)) {
    if(ncol(x) == 1) return(as.vector(as.integer(as.factor(as.character(x[[1]])))))
    stop("clusters data.frame must have exactly one column")
  }

  # Unknown type
  stop("Don't know how to extract clusters from object of class ", paste(class(x), collapse=", "))
}


## validateClusterOrder ----
# Validate that clusters are in the same order/subset as sourceData
# Returns invisibly TRUE if valid, stops with error otherwise
validateClusterOrder <- function(clusters, sourceData, nCases) {
  if(is.null(sourceData)) return(invisible(TRUE))

  if(!is.data.frame(sourceData)) {
    stop("sourceData must be a data.frame")
  }

  if(nrow(sourceData) != nCases) {
    stop("sourceData has ", nrow(sourceData), " rows but clusters expect ", nCases,
         " cases. Make sure sourceData is the same data used to compute the clusters.")
  }

  invisible(TRUE)
}


## addClusters ----
# Add a cluster column to a netCoin object from surScat
# clusters: vector, factor, or clustering object (kmeans, poLCA, tidyLPA)
# sourceData: optional data.frame for validation that clusters are in correct order
addClusters <- function(scatObj, clusters, name=NULL, sort=TRUE, weight=NULL, sourceData=NULL) {
  if(!inherits(scatObj, "netCoin"))
    stop("scatObj must be a netCoin object returned by surScat")

  if(is.null(clusters))
    stop("clusters must be provided")

  # Extract clusters from clustering object if needed
  clusters <- extractClusters(clusters, sourceData)

  # Validate that clusters are in correct order if sourceData provided
  validateClusterOrder(clusters, sourceData, length(clusters))

  # Auto-collapse clusters if they are case-level but object has pattern-level nodes
  collapseReport <- NULL
  clusterDisplay <- NULL  # Will hold display values (e.g., "1|2") for visualization
  if(length(clusters) != nrow(scatObj$nodes)) {
    idx <- attr(scatObj, "caseToPattern")
    if(!is.null(idx) && length(clusters) == length(idx)) {
      # Collapse clusters to pattern level using mode (most frequent value per pattern)
      groups <- split(seq_along(idx), idx)
      collapsedClust <- character(length(groups))  # Modal cluster for sorting
      clusterDisplay <- character(length(groups))  # Display value with all clusters
      conflicts <- data.frame(pattern=integer(), nCases=integer(), clusters=character(), stringsAsFactors=FALSE)

      for(i in seq_along(groups)) {
        casesInPattern <- groups[[i]]
        vals <- as.character(clusters[casesInPattern])
        freq <- table(vals)
        # Order clusters by frequency (descending), ties broken alphabetically
        orderedClust <- names(freq)[order(-freq, names(freq))]

        # Modal (most frequent) cluster for sorting
        collapsedClust[i] <- orderedClust[1]
        # Display all clusters ordered by frequency, joined by "|"
        clusterDisplay[i] <- paste(as.character(orderedClust), collapse="|")

        # Record if there were conflicts (multiple distinct clusters in this pattern)
        if(length(orderedClust) > 1) {
          conflicts <- rbind(conflicts, data.frame(
            pattern=i,
            nCases=length(casesInPattern),
            clusters=paste(as.character(orderedClust), collapse="|"),
            stringsAsFactors=FALSE
          ))
        }
      }

      collapseReport <- list(
        nPatterns = length(groups),
        nCases = length(clusters),
        nConflicts = nrow(conflicts),
        conflictRate = nrow(conflicts) / length(groups),
        conflicts = conflicts
      )

      clusters <- collapsedClust  # Use modal for sorting
      clusterDisplay <- clusterDisplay  # Use display values later
    } else {
      stop(paste0("clusters has length ", length(clusters), " but scatObj has ",
                  nrow(scatObj$nodes), " nodes. If clusters are case-level and scatObj has ",
                  "pattern-level nodes, surScat must have been called with vPatterns."))
    }
  }

  # Default name
  if(is.null(name)) {
    name <- make.unique(c(names(scatObj$nodes), "clusters"))[length(names(scatObj$nodes))+1]
  }

  # Get unique clusters (use modal values for sorting if collapse happened)
  cl <- as.character(clusters)
  unique_cl <- unique(cl)
  n_clusters <- length(unique_cl)

  # If there's clusterDisplay (conflicts exist), get ALL unique values from it
  all_cluster_values <- unique_cl
  if(!is.null(clusterDisplay)) {
    # Extract all unique cluster values that appear in display (including non-modal)
    all_vals_in_display <- unique(unlist(strsplit(clusterDisplay, "\\|")))
    all_cluster_values <- union(all_cluster_values, all_vals_in_display)
  }

  # Sort clusters if requested (by mean coordinate on first layout axis)
  if(sort && !is.null(scatObj$layout) && ncol(scatObj$layout) > 0) {
    cf <- factor(cl, levels=unique_cl)
    layout_first <- scatObj$layout[,1]
    if(is.null(weight)) {
      mu <- tapply(layout_first, cf, mean)
    } else {
      if(length(weight) != length(cl))
        stop("weight must have same length as clusters")
      mu <- tapply(seq_along(cl), cf, function(k) weighted.mean(layout_first[k], weight[k]))
    }
    # Renumber clusters based on sorted order (use all_cluster_values so all values have a mapping)
    order_idx <- order(mu[as.character(all_cluster_values)])
    all_cluster_values_sorted <- all_cluster_values[order_idx]

    # Map cl values: replace old values with new numbering
    cl_new <- rep(NA_character_, length(cl))
    for(i in seq_along(all_cluster_values_sorted)) {
      cl_new[cl == all_cluster_values_sorted[i]] <- as.character(i)
    }
    cl <- cl_new

    if(!is.null(clusterDisplay)) {
      # Map display values by replacing each cluster value with its new number
      clusterDisplay_new <- character(length(clusterDisplay))
      for(j in seq_along(clusterDisplay)) {
        vals <- strsplit(clusterDisplay[j], "\\|")[[1]]
        new_vals <- character(length(vals))
        for(k in seq_along(vals)) {
          idx_in_sorted <- which(all_cluster_values_sorted == vals[k])
          if(length(idx_in_sorted) > 0) {
            new_vals[k] <- as.character(idx_in_sorted)
          }
        }
        clusterDisplay_new[j] <- paste(new_vals, collapse="|")
      }
      clusterDisplay <- clusterDisplay_new
    }
    n_clusters <- length(all_cluster_values_sorted)
    unique_cl <- as.character(seq_len(n_clusters))
  }

  # Use display values (with multiple clusters) if available, otherwise use cl values
  cl_values <- if(!is.null(clusterDisplay)) clusterDisplay else cl
  unique_cl_values <- unique(cl_values)
  n_final <- length(unique_cl_values)

  # Create ordered factor with "Group: N" labels based on actual unique values
  labels <- paste0("Group", ": ", sprintf(paste0("%0", nchar(n_final), "d"), seq_len(n_final)))
  cl_factor <- factor(cl_values, levels=unique_cl_values, labels=labels, ordered=TRUE)

  # Add the new cluster column
  scatObj$nodes[[name]] <- cl_factor

  # Store collapse report if applicable
  if(!is.null(collapseReport)) {
    attr(scatObj, paste0("collapse_", name)) <- collapseReport
  }

  return(scatObj)
}


## replaceClusters ----
# Replace cluster columns in a netCoin object from surScat
# Removes all existing cluster columns and adds the new one
replaceClusters <- function(scatObj, clusters, name=NULL, sort=TRUE, weight=NULL, sourceData=NULL) {
  if(!inherits(scatObj, "netCoin"))
    stop("scatObj must be a netCoin object returned by surScat")

  # Find and remove existing cluster columns
  # Look for column names that look like clustering results (e.g., "Groups(2)", "LCA(3)", etc.)
  # Pattern: word, opening paren, number, closing paren - or just start with "Group"
  clusterPatterns <- c("^Groups?\\([0-9]", "^[A-Za-z]+\\([0-9]", "^Groups?:")
  clusterCols <- character(0)
  for(pat in clusterPatterns) {
    clusterCols <- c(clusterCols, grep(pat, names(scatObj$nodes), value=TRUE))
  }
  clusterCols <- unique(clusterCols)

  if(length(clusterCols) > 0) {
    scatObj$nodes[clusterCols] <- NULL
  }

  # If name not specified, try to reuse the first removed column name
  if(is.null(name) && length(clusterCols) > 0) {
    name <- clusterCols[1]
  }

  # Use addClusters for the rest
  addClusters(scatObj, clusters, name=name, sort=sort, weight=weight, sourceData=sourceData)
}
