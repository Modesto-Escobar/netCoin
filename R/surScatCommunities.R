## communityAlgorithms ----
# The community detection methods addCommunities offers, keyed by the same two-letter code
# the community argument of netCoin takes, so that both name them alike. The full name is
# what the cluster column is called, in the "method(k)" style the other clusterizations of a
# scattergram follow, so that a louvain of five groups reads "louvain(5)".
# maxNodes is how large a network each one can still be asked about. A case-by-case
# correlation network holds one node per case, that is thousands of them, and the cost of
# these methods does not grow alike: cluster_louvain answers on thousands in seconds, while
# cluster_optimal is exponential and gives up inside igraph well below a thousand. The cap
# is stated per method so that the refusal names a number the method can actually meet,
# rather than letting igraph fail deep inside with a message about none of this.
# random marks the three which draw on the random number generator, and so answer differently
# on the same network from one run to the next: louvain and label propagation sweep the nodes
# in a random order, and spinglass anneals from a random state. The other five are determinate.
# Only these three need the seed, but it is set for all of them alike, since a method may be
# given the generator by igraph for a reason this list does not know about.
communityAlgorithms <- list(
  ed = list(name="edgeBetweenness", maxNodes=2000, random=FALSE,
            fun=function(g, w) igraph::cluster_edge_betweenness(g, weights=w)),
  fa = list(name="fastGreedy",      maxNodes=Inf,  random=FALSE,
            fun=function(g, w) igraph::cluster_fast_greedy(g, weights=w)),
  la = list(name="labelProp",       maxNodes=Inf,  random=TRUE,
            fun=function(g, w) igraph::cluster_label_prop(g, weights=w)),
  le = list(name="leadingEigen",    maxNodes=Inf,  random=FALSE,
            fun=function(g, w) igraph::cluster_leading_eigen(g, weights=w)),
  lo = list(name="louvain",         maxNodes=Inf,  random=TRUE,
            fun=function(g, w) igraph::cluster_louvain(g, weights=w)),
  op = list(name="optimal",         maxNodes=100,  random=FALSE,
            fun=function(g, w) igraph::cluster_optimal(g, weights=w)),
  sp = list(name="spinglass",       maxNodes=2000, random=TRUE,
            fun=function(g, w) igraph::cluster_spinglass(g, weights=w)),
  wa = list(name="walktrap",        maxNodes=Inf,  random=FALSE,
            fun=function(g, w) igraph::cluster_walktrap(g, weights=w))
)


## communityCode ----
# Reduce one community method to its two-letter code, accepting the code itself ("lo"), the
# igraph name ("louvain", "cluster_louvain") or the name this file gives it. Unlike the
# congloControl of rD3plot, which warns and falls back to no communities at all, an
# unreadable method stops: it was asked for by name and dropping it silently would leave
# the caller with fewer columns than requested and no clue why.
communityCode <- function(method) {
  if(!is.character(method) || length(method) != 1 || is.na(method))
    stop("each community method must be a single name")
  key <- tolower(sub("^cluster_", "", method))
  full <- vapply(communityAlgorithms, function(a) tolower(a$name), character(1))
  if(key %in% full) return(names(full)[match(key, full)])
  key <- substr(key, 1, 2)
  if(key %in% names(communityAlgorithms)) return(key)
  stop("\"", method, "\" is not a community detection method. Available: ",
       paste0(names(communityAlgorithms), " (", vapply(communityAlgorithms,
              function(a) a$name, character(1)), ")", collapse=", "))
}


## lumpCommunities ----
# Gather into one residual group the communities holding fewer cases than minSize, read as a
# proportion of the cases when below 1 and as a count of them otherwise. A correlation
# network of cases leaves isolated nodes behind, and each of them comes out of igraph as a
# community of its own, so without this a scattergram of three readable groups reports
# dozens and runs into maxGroups.
# maxCommunities caps how many groups come out of all this, by size rather than by rule: what
# minSize cannot foresee is how many communities will clear it, and a legend has room for so
# many. The two criteria gather into the same residual group, and a community falling to
# either of them ends up there.
# Answers the renumbered assignments, the number of the residual group (NA when every
# community was large enough) and what was gathered, for the message to the caller.
lumpCommunities <- function(membership, minSize, maxCommunities=NULL) {
  n <- length(membership)
  # Read as a proportion, the threshold falls to one on few cases, and a threshold of one
  # gathers nothing, since every community holds at least one case. That is where the
  # isolated cases are most of them, so the proportional reading keeps a floor of two: a
  # community of a single case is not a group. An outright count is taken as it is, so that
  # minSize=1 still means that every community is to be kept apart.
  # Both readings round up, so that a community of fewer cases than minSize is gathered
  # whatever minSize is: truncating a count instead kept communities of two on minSize=2.9.
  threshold <- if(minSize < 1) max(2L, as.integer(ceiling(minSize*n))) else as.integer(ceiling(minSize))
  sizes <- table(membership)
  big   <- names(sizes)[sizes >= threshold]
  if(!length(big))
    stop("every community holds fewer than ", threshold, " cases, so all of them would be ",
         "gathered into the residual group. Lower minSize, or raise minL so that the ",
         "correlation network is sparser and its communities larger.")

  # Big ones first, by decreasing size; sort=TRUE renumbers them by coordinate afterwards
  ordered <- names(sort(sizes[big], decreasing=TRUE))
  nSmall  <- length(sizes) - length(big)

  # The cap counts the groups of the legend, the residual one among them: asked for five with
  # something left over, four communities are kept apart and the fifth group is the residual.
  # Nothing is cut when the groups already fit, so a cap of five over four communities and no
  # residual leaves the four alone rather than making a fifth group out of nothing. The cut
  # runs over the communities minSize already found large enough, which reach this point
  # sorted by decreasing size, so what it drops is always the smallest of them.
  nRanked <- 0L
  if(!is.null(maxCommunities) && length(ordered) + (nSmall > 0L) > maxCommunities) {
    keep    <- as.integer(maxCommunities) - 1L
    nRanked <- length(ordered) - keep
    ordered <- ordered[seq_len(keep)]
  }

  map <- stats::setNames(seq_along(ordered), ordered)
  out <- unname(map[as.character(membership)])
  lumped <- is.na(out)
  restNum <- NA_integer_
  if(any(lumped)) {
    restNum <- length(ordered) + 1L
    out[lumped] <- restNum
  }
  list(clusters = out,
       rest     = restNum,
       nGroups  = max(out),
       nLumped  = sum(lumped),
       nSmall   = nSmall,
       nRanked  = nRanked,
       threshold= threshold)
}


## addCommunities ----
# Add to a scattergram one cluster column per community detection method, taken from the
# network of correlations between its cases.
#
# The cases are correlated with one another over the variables, which is the transposition of
# the usual reading of a data frame: each case becomes a column and each variable a row, so
# that a link joins two cases whose profiles run parallel. The communities of that network
# are then a grouping of cases, of the same shape as any other clusterization, and reach the
# node table through addClusters.
#
# scatObj: a netCoin object returned by surScat
# data: the data frame the scattergram was drawn from
# variables: the variables the cases are compared over
# community: one or several methods, by code ("lo") or by name ("louvain")
# scale: standardize each variable before transposing. On by default, since a correlation
#        between two cases over variables of unlike ranges measures the ranges, not the profiles
# minSize: communities below this are gathered into a residual group. See lumpCommunities
# maxCommunities: how many groups may come out at most, the residual one among them. NULL
#                 leaves as many as clear minSize. See lumpCommunities
# maxNodes: how large a network a method may be asked about, overriding the ceiling each one
#           carries. NULL keeps those; a number waits for a method past what it can bear
# weight: passed to addClusters, where it weighs the mean coordinate each group is sorted by.
#         It does not weigh the correlations, which cor() computes unweighted
# seed: the random seed the methods are run from, so that a scattergram comes out the same
#       twice. 2020 as in surScat; NULL leaves the generator alone. See below
addCommunities <- function(scatObj, data, variables=names(data), community="lo",
                           scale=TRUE, method=c("pearson", "kendall", "spearman"),
                           criteria="value", minL=0.95, maxL=Inf, pairwise=FALSE,
                           minSize=0.01, maxCommunities=NULL, name=NULL, sort=TRUE, weight=NULL,
                           maxGroups=25, maxNodes=NULL, seed=2020) {
  if(!inherits(scatObj, "netCoin"))
    stop("scatObj must be a netCoin object returned by surScat")
  if(missing(data) || is.null(data))
    stop("data must be the data frame the scattergram was drawn from")
  if(length(minSize) != 1 || is.na(minSize) || minSize <= 0)
    stop("minSize must be a proportion of the cases (below 1) or a count of them (1 or more)")
  if(!is.null(maxCommunities)) {
    if(length(maxCommunities) != 1 || !is.numeric(maxCommunities) || is.na(maxCommunities) ||
       maxCommunities < 2)
      stop("maxCommunities must be a count of two groups or more, the residual one among ",
           "them, or NULL to leave as many groups as clear minSize")
    # A cap of five and a half groups is a cap of five: unlike minSize, which rounds up so
    # that a community smaller than it is always gathered, a bound on a count of groups is
    # only honoured by rounding down
    maxCommunities <- as.integer(maxCommunities)
  }

  # Three of the methods draw on the random number generator, so without a seed the same
  # call answers differently from one run to the next: the communities move, their sizes
  # move with them, and which of them minSize gathers moves too, which reads as if minSize
  # itself were behaving erratically. The generator is left as it was found, as surScat
  # leaves it, so that seeding a scattergram does not seed whatever the caller does next.
  if(!is.null(seed)) {
    if(length(seed) != 1 || !is.numeric(seed) || is.na(seed))
      stop("seed must be a single number, or NULL to leave the random state as it is")
    if(exists(".Random.seed", envir=globalenv())) {
      oldseed <- get(".Random.seed", envir=globalenv())
      on.exit(assign(".Random.seed", oldseed, envir=globalenv()), add=TRUE)
    } else {
      on.exit(rm(".Random.seed", envir=globalenv()), add=TRUE)
    }
  }

  codes <- vapply(community, communityCode, character(1), USE.NAMES=FALSE)
  if(anyDuplicated(codes))
    stop("community holds the same method more than once: ",
         paste(unique(codes[duplicated(codes)]), collapse=", "))

  data <- as.data.frame(data)
  missingVars <- setdiff(variables, names(data))
  if(length(missingVars))
    stop("variables not found in data: ", paste(missingVars, collapse=", "))
  if(length(variables) < 3)
    warning("the cases are correlated over ", length(variables), " variables, so each ",
            "correlation rests on ", length(variables), " points and says very little. ",
            "Compare the cases over more variables.", call.=FALSE)

  # Each variable as a number, in the same spirit as the axes of surScat: an ordered factor
  # is read by its level order, a nominal one has none to read and stops
  D <- data[, variables, drop=FALSE]
  D <- as.data.frame(lapply(stats::setNames(variables, variables), function(v) {
    x <- D[[v]]
    if(is.ordered(x)) return(as.numeric(x))
    if(is.factor(x))
      stop(v, " is a nominal variable, so correlating the cases over it would impose an ",
           "arbitrary order on its categories. Turn it into an ordered factor or a numeric one.")
    x <- unclass(x) # labelled vectors keep their underlying codes
    if(!is.numeric(x))
      stop(v, " must be numeric or an ordered factor to compare the cases over it")
    as.numeric(x)
  }))
  D <- D[complete.cases(D), , drop=FALSE] # the same filter surScat applies to its own cases

  # Line the cases up with the nodes that were drawn. Three shapes reach this point: nodes
  # that are the cases themselves, nodes that are patterns of them (addClusters collapses the
  # case-level communities afterwards), and nodes that are a sample drawn by maxN.
  nNodes  <- nrow(scatObj$nodes)
  idx     <- attr(scatObj, "caseToPattern")
  sampled <- attr(scatObj, "sampledNodes")
  from    <- attr(scatObj, "sampledFrom")
  if(!is.null(sampled)) {
    # maxN capped the nodes. Which cases they are is recorded, so the network can be built on
    # them rather than on the whole file, and nothing has to be guessed. That only holds when
    # the nodes were cases: a sample of patterns cannot be traced back to the cases behind it,
    # since the map that would do it is not kept once a sample is drawn.
    if(!identical(as.integer(from), nrow(D)))
      stop("the ", nNodes, " nodes are a sample of ", from, " patterns, not of the ", nrow(D),
           " complete cases of data, so which cases each node stands for is no longer ",
           "recorded and the communities cannot be lined up with them. Raise maxN in surScat ",
           "so that every pattern is drawn, or drop patterns.")
    D <- D[sampled, , drop=FALSE]
  } else if(!is.null(idx)) {
    if(length(idx) != nrow(D))
      stop("data holds ", nrow(D), " complete cases but the scattergram was drawn from ",
           length(idx), ". Pass the same data frame and the same variables surScat was given.")
  } else if(nrow(D) != nNodes) {
    stop("data holds ", nrow(D), " complete cases but the scattergram has ", nNodes,
         " nodes. Pass the same data frame and the same variables surScat was given.",
         sampledNote(scatObj))
  }
  if(nrow(D) < 3) stop("at least three cases are needed to build a network of them")

  # A correlation between two cases centres each of them on its own mean across the
  # variables, so it compares the shape of their profiles. Standardizing first puts every
  # variable on one footing; without it the variable of widest range writes the profile.
  if(scale) {
    constant <- vapply(D, function(x) stats::sd(x) == 0, logical(1))
    if(any(constant))
      stop("these variables take a single value across the cases, so they cannot be ",
           "standardized: ", paste(variables[constant], collapse=", "),
           ". Drop them, or set scale=FALSE.")
    D <- as.data.frame(base::scale(as.matrix(D)))
  }

  # Transposed, each case is a column and netCorr correlates the columns, so the nodes of the
  # network are the cases. Names are the position of each case among the ones drawn, padded
  # so that they read in order, and unique, which netCorr requires of its node names.
  tD <- as.data.frame(t(as.matrix(D)))
  names(tD) <- sprintf(paste0("%0", nchar(nrow(D)), "d"), seq_len(nrow(D)))

  N <- netCorr(tD, method=method, criteria=criteria, minL=minL, maxL=maxL, pairwise=pairwise)
  g <- toIgraph(N)
  if(igraph::vcount(g) != nrow(D))
    stop("the correlation network came back with ", igraph::vcount(g), " nodes for ",
         nrow(D), " cases") # not reachable through the arguments; guards a silent misalignment

  # netCorr hands the correlation to the links as their weight, which is what the igraph
  # methods read. A negative one has no meaning for them, so those are dropped rather than
  # letting igraph fail deep inside with a message about the graph.
  w <- if("weight" %in% igraph::edge_attr_names(g)) igraph::E(g)$weight else NULL
  if(!is.null(w) && any(w < 0)) {
    warning("the network holds links of negative correlation, which the community methods ",
            "cannot weigh, so the communities were found without weighing them. Raise minL ",
            "above 0 to keep only links of positive correlation.", call.=FALSE)
    w <- NULL
  }

  ids <- igraph::V(g)$name # membership comes in this order, which need not be the node order
  reports <- character(0)
  ranked  <- FALSE
  for(code in codes) {
    algo <- communityAlgorithms[[code]]
    limit <- if(is.null(maxNodes)) algo$maxNodes else maxNodes
    if(igraph::vcount(g) > limit)
      stop(algo$name, " does not answer on a network of ", igraph::vcount(g),
           " nodes; it bears ", limit, ". Choose one of ",
           paste(vapply(communityAlgorithms[vapply(communityAlgorithms,
                 function(a) is.infinite(a$maxNodes), logical(1))],
                 function(a) a$name, character(1)), collapse=", "),
           ", cap the nodes with maxN in surScat, or raise maxNodes to wait for it anyway.")

    # Seeded once per method rather than once for the call, so that each method answers the
    # same whichever others were asked for alongside it: otherwise a method would inherit
    # whatever state the ones before it had left, and adding a method to community would
    # move the communities of the ones after it.
    if(!is.null(seed)) set.seed(seed)
    comm <- tryCatch(algo$fun(g, w), error=function(e)
      stop(algo$name, " could not find communities in this network: ", conditionMessage(e),
           if(sum(igraph::degree(g) == 0)) paste0(" It holds ", sum(igraph::degree(g) == 0),
           " isolated cases, which some methods refuse. Lower minL so that fewer cases are ",
           "left without links, or choose another method."), call.=FALSE))

    memb <- igraph::membership(comm)
    memb <- as.integer(memb[match(ids, names(memb))]) # back to the order of the cases
    lump <- lumpCommunities(memb, minSize, maxCommunities)

    # The column is named in the "method(k)" style the rest of the scattergram follows, so
    # that addOneCluster reads its prefix off the name and the legend says "Louvain: 1"
    own <- paste0(algo$name, "(", lump$nGroups, ")")
    column <- if(is.null(name)) own else if(length(codes) > 1) paste0(name, ".", own) else name

    # Padded to a fixed width, so that the character ordering addClusters puts the values
    # through keeps them in numeric order and the residual group stays the last of them
    codesOut <- sprintf(paste0("%0", nchar(lump$nGroups), "d"), lump$clusters)
    frame <- data.frame(codesOut, stringsAsFactors=FALSE, check.names=FALSE)
    names(frame) <- column

    scatObj <- addClusters(scatObj, frame, sort=sort, weight=weight,
                           maxGroups=maxGroups, rest=lump$rest)

    # Which of the two criteria gathered what, since the way out differs: minSize is lowered,
    # maxCommunities is raised, and the caller cannot tell which is in his way from a count
    if(lump$nLumped) {
      why <- c(if(lump$nSmall)  paste0(lump$nSmall, " of fewer than ", lump$threshold, " cases"),
               if(lump$nRanked) paste0(lump$nRanked, " past the ", maxCommunities-1L, " largest"))
      reports <- c(reports, paste0(column, ": ", lump$nSmall + lump$nRanked, " communities (",
                   paste(why, collapse=", "), "), ", lump$nLumped, " cases in all"))
      ranked  <- ranked || lump$nRanked > 0L
    }
  }

  if(length(reports))
    message("Gathered into the residual group -- ", paste(reports, collapse="; "), ". Set ",
            if(ranked) "minSize=1 and maxCommunities=NULL " else "minSize=1 ",
            "to keep every community apart.")

  return(scatObj)
}
