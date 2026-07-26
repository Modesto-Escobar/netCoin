## addClusters ----
# Add a cluster column to a netCoin object from surScat
addClusters <- function(scatObj, clusters, name=NULL, sort=TRUE, weight=NULL) {
  if(!inherits(scatObj, "netCoin"))
    stop("scatObj must be a netCoin object returned by surScat")

  if(is.null(clusters))
    stop("clusters must be provided")

  if(is.data.frame(clusters)) {
    if(ncol(clusters) == 1) clusters <- clusters[[1]]
    else stop("clusters data.frame must have exactly one column")
  }

  clusters <- as.vector(clusters)

  if(length(clusters) != nrow(scatObj$nodes))
    stop(paste0("clusters has length ", length(clusters), " but scatObj has ",
                nrow(scatObj$nodes), " nodes"))

  # Default name
  if(is.null(name)) {
    name <- make.unique(c(names(scatObj$nodes), "clusters"))[length(names(scatObj$nodes))+1]
  }

  # Get unique clusters
  cl <- as.character(clusters)
  unique_cl <- unique(cl)
  n_clusters <- length(unique_cl)

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
    unique_cl <- unique_cl[order(mu)]
  }

  # Create ordered factor with "Group: N" labels
  labels <- paste0("Group", ": ", sprintf(paste0("%0", nchar(n_clusters), "d"), seq_len(n_clusters)))
  cl_factor <- factor(cl, levels=unique_cl, labels=labels, ordered=TRUE)

  # Add the new cluster column
  scatObj$nodes[[name]] <- cl_factor

  return(scatObj)
}


## replaceClusters ----
# Replace cluster columns in a netCoin object from surScat
replaceClusters <- function(scatObj, clusters, name=NULL, sort=TRUE, weight=NULL) {
  if(!inherits(scatObj, "netCoin"))
    stop("scatObj must be a netCoin object returned by surScat")

  # Remove existing cluster columns
  groupCols <- grep("^Group", names(scatObj$nodes), value=TRUE)
  scatObj$nodes[groupCols] <- NULL

  # If name not specified, try to reuse the first removed column name
  if(is.null(name) && length(groupCols) > 0) {
    name <- groupCols[1]
  }

  # Use addClusters for the rest
  addClusters(scatObj, clusters, name=name, sort=sort, weight=weight)
}
