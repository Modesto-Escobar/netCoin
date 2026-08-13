## selectL4CReductions ----
# Pick the dimensionality reductions out of a looking4clusters object.
# As with the clusterizations, the object is read through its plain list structure, so
# looking4clusters (a GitHub package, not on CRAN) needs not be declared nor installed.
# Its $reductions slot is a named list holding one two-column matrix per reduction
# ("pca", "tsne", "cmds", "nmf", "umap", plus any added through addreduction), each with
# one row per case.
# which: NULL (every reduction), indices, or names.
selectL4CReductions <- function(x, which=NULL) {
  entries <- x$reductions
  if(is.null(entries) || !length(entries))
    stop("looking4clusters object holds no dimensionality reduction")

  if(is.null(which)) {
    keys <- names(entries)
  } else if(is.numeric(which)) {
    if(any(is.na(which) | which < 1 | which > length(entries)))
      stop("which must be indices between 1 and ", length(entries))
    keys <- names(entries)[which]
  } else if(is.character(which)) {
    unknown <- setdiff(which, names(entries))
    if(length(unknown))
      stop("which=\"", paste(unknown, collapse="\", \""),
           "\" is not among the reductions of the looking4clusters object: ",
           paste(names(entries), collapse=", "))
    keys <- which
  } else {
    stop("which must be NULL, indices or names")
  }

  lapply(entries[keys], function(entry) {
    m <- as.matrix(entry)
    if(!is.numeric(m)) stop("the reductions of a looking4clusters object must be numeric")
    m
  })
}


## splitAxisPairs ----
# A matrix or data.frame of an even number of columns holds one plane per consecutive
# pair of columns: the first is drawn horizontally and the second vertically.
# Each plane is named after its two columns, as in "PC1-PC2".
splitAxisPairs <- function(axes) {
  if(is.data.frame(axes)) axes <- as.matrix(axes)
  if(!is.matrix(axes)) axes <- as.matrix(axes)
  if(!is.numeric(axes))
    stop("axes must hold numeric coordinates")
  if(ncol(axes) < 2 || ncol(axes) %% 2 != 0)
    stop("axes must have an even number of columns, a pair of coordinates per plane, but it has ",
         ncol(axes))

  cn <- colnames(axes)
  if(is.null(cn)) cn <- paste0("axis", seq_len(ncol(axes)))
  out <- list()
  for(j in seq(1, ncol(axes), by=2)) {
    pair <- axes[, c(j, j+1), drop=FALSE]
    colnames(pair) <- cn[c(j, j+1)]
    out[[paste(cn[c(j, j+1)], collapse="-")]] <- pair
  }
  out
}


## extractAxesList ----
# Reduce any accepted input to a named list of two-column numeric matrices, so that
# holders of several planes (a matrix of many columns, a looking4clusters object, a list)
# and single ones all travel the same path. An empty name means the plane has none of its
# own, and the caller supplies it.
extractAxesList <- function(axes, which=NULL) {
  if(inherits(axes, "looking4clusters"))
    return(selectL4CReductions(axes, which))

  # a list of matrices, one or several planes each
  if(is.list(axes) && !is.data.frame(axes)) {
    if(!length(axes)) stop("axes holds no coordinates")
    out <- list()
    for(i in seq_along(axes)) {
      own   <- names(axes)[i]
      if(is.null(own) || is.na(own)) own <- ""
      pairs <- splitAxisPairs(axes[[i]])
      # the name of the entry prevails when it holds a single plane, and prefixes the
      # names of its columns when it holds several
      if(nzchar(own))
        names(pairs) <- if(length(pairs) == 1) own else paste0(own, ".", names(pairs))
      out <- c(out, pairs)
    }
    return(out)
  }

  splitAxisPairs(axes)
}


## alignAxes ----
# Bring a pair of coordinates to the rows of the node table. Coordinates given case by
# case are collapsed into patterns, by their (weighted) mean, when surScat collapsed them.
alignAxes <- function(scatObj, pair, weight=NULL, label="") {
  what <- if(nzchar(label)) paste0("axes \"", label, "\"") else "axes"
  n <- nrow(scatObj$nodes)
  if(nrow(pair) == n) return(pair)

  idx <- attr(scatObj, "caseToPattern")
  if(is.null(idx) || nrow(pair) != length(idx))
    stop(what, " has ", nrow(pair), " rows but scatObj has ", n,
         " nodes. If the coordinates are case-level and scatObj has pattern-level nodes,",
         " surScat must have been called with vPatterns.")

  # surScat draws a sample of the patterns when there are more than maxN of them, and does
  # not record which ones, so case-level coordinates can no longer be matched to the nodes
  if(max(idx) != n)
    stop(what, " is case-level, but the ", max(idx), " patterns of scatObj were subsampled",
         " down to ", n, " nodes by maxN. Raise maxN in surScat to add axes case by case.")

  if(!is.null(weight) && length(weight) != nrow(pair))
    stop("weight must have one value per row of ", what)

  groups <- split(seq_along(idx), idx)
  out <- vapply(groups, function(i)
    if(is.null(weight)) colMeans(pair[i, , drop=FALSE])
    else colSums(pair[i, , drop=FALSE]*weight[i])/sum(weight[i]),
    numeric(2))
  out <- t(out)
  colnames(out) <- colnames(pair)
  rownames(out) <- NULL
  out
}


## currentLayouts ----
# The planes the object already holds, as the named list of matrices rD3plot draws a
# selector from. An object that has not been through addAxes keeps a single plane in the
# fx/fy columns of its node table, which becomes the first entry, named after the axes it
# was drawn with (their percentage of variance dropped, as it belongs to the label and not
# to the name of the plane).
currentLayouts <- function(scatObj) {
  if(length(scatObj$layouts)) return(scatObj$layouts)

  nodes <- scatObj$nodes
  if(!all(c("fx","fy") %in% names(nodes)))
    stop("scatObj holds no coordinates. Only a scattergram built by surScat, or a netCoin",
         " object given a layout, can be given further axes.")

  first <- cbind(nodes$fx, nodes$fy)
  labs  <- scatObj$options$axesLabels
  if(length(labs) >= 2) {
    colnames(first) <- sub(" *\\([^()]*\\) *$", "", as.character(labs[1:2]))
    key <- paste(colnames(first), collapse="-")
  } else {
    colnames(first) <- c("fx","fy")
    key <- "layout1"
  }
  out <- list(first)
  names(out) <- key
  out
}


## addAxes ----
# Add one or several planes to a netCoin object from surScat, which the viewer then offers
# in a selector next to the graph. The coordinates already drawn become the first plane.
# axes: a matrix or data.frame with an even number of columns, a pair per plane, a list of
#       such matrices, or a looking4clusters object, whose reductions (pca, tsne, cmds,
#       nmf, umap) are taken as planes.
# name: name of the new plane. Each plane keeps its own name when several are added at
#       once, and name becomes a prefix.
# which: reduction to take from a looking4clusters object, as a name or an index. When it
#        is left empty, every reduction held by the object is added.
# weight: weights of the cases, used when collapsing them into patterns.
addAxes <- function(scatObj, axes, name=NULL, which=NULL, weight=NULL) {
  if(!inherits(scatObj, "netCoin"))
    stop("scatObj must be a netCoin object returned by surScat")

  if(is.null(axes))
    stop("axes must be provided")

  sets    <- extractAxesList(axes, which)
  several <- length(sets) > 1
  layouts <- currentLayouts(scatObj)

  for(i in seq_along(sets)) {
    own <- names(sets)[i]
    if(is.null(own) || is.na(own)) own <- ""

    # its own name, prefixed by name when it was given; name alone when there is only one
    key <- if(is.null(name)) {
             if(nzchar(own)) own else "axes"
           } else if(several) {
             paste0(name, ".", own)
           } else name
    key <- make.unique(c(names(layouts), key))[length(layouts)+1]

    layouts[[key]] <- alignAxes(scatObj, sets[[i]], weight=weight, label=key)
  }

  scatObj$layouts <- layouts
  return(scatObj)
}
