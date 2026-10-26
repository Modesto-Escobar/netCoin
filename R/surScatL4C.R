## surScat and looking4clusters ----
# surScat hands looking4clusters the very cases its own k-means groups, but as they come: the
# quantitative active variables unstandardized, or the 0/1 indicators of the qualitative ones.
# Its planes and its groups are thus a second look at the data, not a copy of the first, and
# they are added to the plane and to the groups surScat draws, never put in their place.

l4cReductions <- c("pca", "tsne", "cmds", "nmf", "umap")
l4cMethods    <- c("kmeans", "pam", "hclust")

# Beyond this number of cases looking4clusters leaves out the methods that need the whole
# matrix of distances between them (pam, hclust and cmds), unless force_execution is stated.
l4cMaxCases <- 5000


## l4cLayout ----
# The reductions surScat's layouts asks for, in the order asked. TRUE asks for all of them, and
# "mds" is read as "cmds". A scattergram has no links, so the network layouts netCoin takes
# ("fo", "la"...) mean nothing to it: whatever is not a reduction is dropped with a warning.
l4cLayout <- function(layouts) {
  if(is.null(layouts) || isFALSE(layouts)) return(character(0))
  if(isTRUE(layouts)) return(l4cReductions)
  if(!is.character(layouts)) {
    warning("surScat has no links, so a network layout does not apply to it and was ignored. ",
            "layouts takes the names of dimensionality reductions: ",
            paste(l4cReductions, collapse=", "), call.=FALSE)
    return(character(0))
  }
  asked <- rotationName(layouts)
  asked[asked == "mds"] <- "cmds"
  if("all" %in% asked) { # every reduction too, next to whatever else was stated
    layouts <- c(layouts[asked != "all"], l4cReductions)
    asked   <- c(asked[asked != "all"],   l4cReductions)
  }
  own <- asked %in% rotationLayouts # rotations of surScat's own plane, which rotationAxes takes
  layouts <- layouts[!own]
  asked   <- asked[!own]
  bad <- layouts[!(asked %in% l4cReductions)]
  if(length(bad))
    warning("surScat has no links, so network layouts do not apply to it: ",
            paste0("\"", bad, "\"", collapse=", "), " ignored. layouts takes the names of ",
            "dimensionality reductions: ", paste(l4cReductions, collapse=", "),
            ", and the rotations ", paste(rotationLayouts, collapse=", "),
            ", or \"all\" for all of them", call.=FALSE)
  unique(asked[asked %in% l4cReductions])
}


## l4cClusters ----
# The clustering methods surScat's clusters asks for. TRUE asks for all of them, as layouts
# does: the k-means of looking4clusters works on the active variables as they come, so it is
# not a copy of surScat's own, which works on the standardized factorial coordinates.
# The latent models "lpa" and "lca" are taken here too, TRUE standing for whichever of them
# suits the type of analysis ("latent"), which latentClusters settles once type is known.
l4cClusters <- function(clusters) {
  if(is.null(clusters) || isFALSE(clusters)) return(character(0))
  if(isTRUE(clusters)) return(c(l4cMethods, "latent"))
  known <- c(l4cMethods, "lpa", "lca")
  if(!is.character(clusters))
    stop("clusters must be TRUE, FALSE or the names of clustering methods: ",
         paste(known, collapse=", "))
  asked <- tolower(clusters)
  bad <- clusters[!(asked %in% known)]
  if(length(bad))
    stop("clusters takes the names of clustering methods (", paste(known, collapse=", "),
         "), not ", paste0("\"", bad, "\"", collapse=", "))
  unique(asked)
}


## l4cRun ----
# Run looking4clusters on the cases, one row each, and hand back the planes and the
# clusterizations asked for, case by case: list(axes, clusters).
# Only the numbers of groups in nclusters are kept, named as surScat names its own columns
# ("pam(3)"), method by method and, within each, in the order nclusters states them.
# looking4clusters computes every number of groups between selectedk-5 and selectedk+5, so it
# is called as many times as needed for that span to cover nclusters.
# seed is set again before every method, so that each result is the same whatever else has
# been asked for alongside it.
l4cRun <- function(data, reductions, methods, nclusters, seed=NULL, force_execution=FALSE) {
  n <- nrow(data)
  object  <- create_l4c(data)
  big     <- n > l4cMaxCases && !force_execution
  skipped <- character(0)
  reseed  <- function() if(!is.null(seed)) set.seed(seed)
  # a method that fails is reported and skipped, so that it does not cost the scattergram
  attempt <- function(what, expr) tryCatch(expr, error=function(e) {
    warning("looking4clusters could not compute ", what, ": ", conditionMessage(e), call.=FALSE)
    object
  })

  ks <- unique(nclusters)
  ks <- ks[ks >= 2 & ks < n]
  if(length(methods) && length(ks)) {
    pamHclust <- intersect(methods, c("pam", "hclust"))
    if(big && length(pamHclust)) {
      skipped   <- c(skipped, pamHclust)
      pamHclust <- character(0)
    }
    rest <- sort(ks)
    while(length(rest)) {
      sk <- rest[1] + 5 # whose span starts at the smallest number of groups still missing
      if("kmeans" %in% methods) {
        reseed()
        object <- attempt("kmeans", run_kmeans(object, selectedk=sk))
      }
      if(length(pamHclust)) {
        reseed()
        object <- attempt(paste(pamHclust, collapse=" and "),
                          run_pam_hclust(object, selectedk=sk))
      }
      rest <- rest[rest > sk + 5]
    }
  }

  clusters <- list()
  if(length(object$clusters)) {
    found <- listL4CClusters(object)
    width <- nchar(max(ks))
    for(meth in intersect(methods, setdiff(l4cMethods, skipped)))
      for(k in ks) {
        key <- paste0(meth, "(", k, ")")
        if(key %in% names(found))
          clusters[[paste0(meth, "(", sprintf(paste0("%0", width, "d"), k), ")")]] <- found[[key]]
      }
  }

  for(r in reductions) {
    if(big && r == "cmds") {
      skipped <- c(skipped, r)
      next
    }
    reseed()
    object <- attempt(r, switch(r,
      pca  = run_pca(object),
      tsne = run_tsne(object),
      cmds = run_mds(object),
      nmf  = run_nmf(object),
      umap = run_umap(object)))
  }
  axes <- object$reductions[intersect(reductions, names(object$reductions))]
  # the pca of looking4clusters works on the covariances of the variables as they come, not on
  # the standardized ones surScat's own PCA plane comes from, so its plane and axes are told apart
  if(!is.null(axes$pca)) {
    colnames(axes$pca) <- c("cPC1", "cPC2")
    names(axes)[names(axes) == "pca"] <- "cPC"
  }
  # the factorial axes (pca, cmds) are labelled with their percentage of variance or inertia,
  # as surScat labels its own; the axes of tsne, nmf and umap have none
  for(r in names(axes)) {
    pct <- attr(axes[[r]], "pctvar")
    if(length(pct) == 2 && all(is.finite(pct)))
      colnames(axes[[r]]) <- sprintf("%s (%.1f%%)", colnames(axes[[r]]), pct)
  }

  if(length(skipped))
    warning("with ", n, " cases, more than ", l4cMaxCases, ", ",
            paste(skipped, collapse=", "), " would take too long and ",
            if(length(skipped) == 1) "was" else "were",
            " left out. State force_execution=TRUE to compute ",
            if(length(skipped) == 1) "it" else "them", " anyway.", call.=FALSE)

  list(axes=axes, clusters=clusters)
}


## addL4C ----
# Add what l4cRun found, case by case, to the scattergram surScat has just built, whose nodes
# may be patterns of cases, a sample of the cases (maxN), or a sample of the patterns.
# addClusters and addAxes collapse cases into patterns through the caseToPattern map, which
# surScat drops when maxN has drawn a sample. It is laid here for the cases of the nodes
# drawn, and the object then gets back the attributes it had.
# idx: the case-to-pattern map (NULL when the nodes are cases). kept: the nodes maxN drew.
# suffix: NULL, or the suffix under which every added plane is also shown aligned by
# Procrustes to the plane surScat drew (see layoutsAlign).
# alignTo: NULL, or the name or index of the plane every plane is brought to, in place of the
# plane surScat drew. A name is read as the plane's own, or as a rotation ("varimax" stands for
# "varimax 1-2", the first pair of its components), whatever the case.
# axesLabels: passed to addAxes; NA keeps the labels the caller stated.
addL4C <- function(xnc, found, idx=NULL, kept=NULL, caseWeight=NULL, sort=TRUE,
                   suffix=NULL, alignTo=NULL, axesLabels=NULL) {
  if(!length(found$axes) && !length(found$clusters)) return(xnc)

  sel <- seq_len(if(length(found$clusters)) length(found$clusters[[1]])
                 else nrow(found$axes[[1]]))
  map <- idx
  if(!is.null(kept)) {
    if(is.null(idx)) sel <- kept
    else {
      sel <- which(idx %in% kept)
      map <- match(idx[sel], kept)
    }
  }

  saved <- list(map=attr(xnc, "caseToPattern"), weight=attr(xnc, "caseWeight"))
  attr(xnc, "caseToPattern") <- map
  attr(xnc, "caseWeight")    <- if(!is.null(map) && !is.null(caseWeight)) caseWeight[sel]

  if(length(found$clusters))
    xnc <- addClusters(xnc, as.data.frame(lapply(found$clusters, function(cl) cl[sel]),
                                          check.names=FALSE),
                       sort=sort, maxGroups=NULL)
  if(length(found$axes))
    xnc <- addAxes(xnc, lapply(found$axes, function(a) a[sel, , drop=FALSE]),
                   axesLabels=axesLabels)

  attr(xnc, "caseToPattern") <- saved$map
  attr(xnc, "caseWeight")    <- saved$weight
  # the centroids are worked out again with the attributes the object keeps, as any later
  # addClusters or addAxes will
  if(length(attr(xnc, "clusterColumns"))) {
    xnc$clusters <- currentClusters(xnc)
    xnc$clusterMeans <- currentClusterMeans(xnc)
  }

  if(!is.null(suffix) && length(found$axes)) { # every plane, so that the reference is not repeated
    planes <- names(currentLayouts(xnc))
    to <- alignPlane(alignTo, planes)
    xnc <- suppressMessages(layoutsAlign(xnc, "procrustes", to=to, suffix=suffix))
  }
  xnc
}


## alignPlane ----
# The plane of planes that alignTo points to, as the index layoutsAlign takes: 1, the plane
# surScat drew, when alignTo is NULL or points to none (with a warning in this last case).
alignPlane <- function(alignTo, planes) {
  if(is.null(alignTo)) return(1)
  hit <- integer(0)
  if(is.numeric(alignTo)) hit <- if(alignTo %in% seq_along(planes)) as.integer(alignTo) else integer(0)
  else {
    key  <- tolower(c(alignTo, rotationName(alignTo)))
    low  <- tolower(planes)
    hit  <- which(low %in% key)
    if(!length(hit)) hit <- which(Reduce(`|`, lapply(key, function(k) startsWith(low, paste0(k, " ")))))
  }
  if(!length(hit)) {
    warning("alignTo=", if(is.character(alignTo)) paste0("\"", alignTo, "\"") else alignTo,
            " is none of the planes (", paste(planes, collapse=", "),
            "), so they are aligned to the main plane", call.=FALSE)
    return(1)
  }
  hit[1]
}
