## selectL4CClusters ----
# Pick one clusterization out of a looking4clusters object.
# The object is read through its plain list structure, so looking4clusters
# (a GitHub package, not on CRAN) needs not be declared nor installed.
# Its $clusters slot mixes single factors (added through addcluster) with
# data.frames holding one column per number of groups ("levels_2" ... "levels_10"),
# tagged with an "optim_cluster" attribute that stores the optimal k itself.
# which: NULL (options$myGroups, else the first entry), an index, a name, or
# "name:k" to ask for a given number of groups.
selectL4CClusters <- function(x, which=NULL) {
  entries <- x$clusters
  if(is.null(entries) || !length(entries))
    stop("looking4clusters object has no clusters")

  k <- NULL
  if(is.null(which)) {
    key <- x$options$myGroups
    if(is.null(key) || !(key %in% names(entries))) key <- names(entries)[1]
  } else if(is.numeric(which)) {
    if(length(which) != 1 || which < 1 || which > length(entries))
      stop("which must be a single index between 1 and ", length(entries))
    key <- names(entries)[which]
  } else if(is.character(which) && length(which) == 1) {
    parts <- strsplit(which, ":", fixed=TRUE)[[1]]
    key <- parts[1]
    if(length(parts) > 1) {
      k <- suppressWarnings(as.integer(parts[2]))
      if(is.na(k)) stop("the number of groups in which=\"", which, "\" must be an integer")
    }
    if(!(key %in% names(entries)))
      stop("which=\"", key, "\" is not among the clusterizations of the looking4clusters object: ",
           paste(names(entries), collapse=", "))
  } else {
    stop("which must be NULL, a single index or a single name")
  }

  entry <- entries[[key]]

  # Single clusterization: added through addcluster()
  if(!is.data.frame(entry)) {
    if(!is.null(k))
      stop("clusterization \"", key, "\" holds a single solution, so a number of groups cannot be chosen")
    return(as.vector(as.integer(as.factor(entry))))
  }

  # Several solutions: one column per number of groups
  if(is.null(k)) {
    k <- attr(entry, "optim_cluster")
    if(is.null(k)) return(as.vector(as.integer(as.factor(entry[[1]]))))
  }
  column <- paste0("levels_", k)
  if(!(column %in% names(entry)))
    stop("clusterization \"", key, "\" has no solution with ", k, " groups. Available: ",
         paste(sub("^levels_", "", names(entry)), collapse=", "))

  as.vector(as.integer(as.factor(entry[[column]])))
}


## listL4CClusters ----
# Enumerate every clusterization held by a looking4clusters object, as a named
# list of cluster vectors. Methods run for several numbers of groups yield one
# entry per solution. Names follow the "method(k)" style of surScat columns.
listL4CClusters <- function(x) {
  entries <- x$clusters
  if(is.null(entries) || !length(entries))
    stop("looking4clusters object has no clusters")

  out <- list()
  for(key in names(entries)) {
    entry <- entries[[key]]
    if(!is.data.frame(entry)) {
      cl <- as.vector(as.integer(as.factor(entry)))
      out[[paste0(key, "(", length(unique(cl)), ")")]] <- cl
    } else {
      for(column in names(entry)) {
        cl <- as.vector(as.integer(as.factor(entry[[column]])))
        out[[paste0(key, "(", sub("^levels_", "", column), ")")]] <- cl
      }
    }
  }
  out
}


## tidyLPAName ----
# Name of a tidyLPA solution, in the "method(k)" style the other clusterizations follow, so
# that a latent profile analysis of four profiles reads "LPA(4)" and is labelled "LPA: 1" ...
# "LPA: 4". A tidyLPA object estimated for several model specifications at once tells them
# apart by the number tidyLPA gives each, as in "LPA1(4)" and "LPA6(4)", since the analysis
# is then no longer a single one; oneModel says whether that is the case.
tidyLPAName <- function(model, classes, oneModel=TRUE) {
  if(is.null(model) || is.na(model)) model <- 1
  stem <- if(oneModel) "LPA" else paste0("LPA", model)
  if(is.null(classes) || is.na(classes)) return(stem)
  paste0(stem, "(", classes, ")")
}


## listTidyLPAClusters ----
# Enumerate every solution held by a tidyLPA object, as a named list of cluster
# vectors, one entry per (model_number, classes_number) pair.
# get_data() answers in two shapes: a single solution comes wide, one row per case,
# carrying no id; several solutions come long, one row per case and class, stacked
# in one block per pair, so each block is reduced back to one row per case through
# its id column. Reading the long shape as if it were the wide one is what made
# addClusters see n*sum(k) assignments instead of n.
listTidyLPAClusters <- function(x) {
  data <- as.data.frame(tidyLPA::get_data(x)) # a tibble warns on absent columns
  has <- function(column) column %in% names(data)
  if(!has("Class"))
    stop("tidyLPA object must have Class column in returned data")
  model   <- if(has("model_number"))   data$model_number   else rep(1L, nrow(data))
  classes <- if(has("classes_number")) data$classes_number else rep(NA_integer_, nrow(data))

  if(!has("id")) # a single solution already holds one row per case
    return(stats::setNames(list(as.vector(as.integer(as.factor(data$Class)))),
                           tidyLPAName(model[1], classes[1])))

  combos <- unique(data.frame(model=model, classes=classes))
  oneModel <- length(unique(combos$model)) == 1 # one specification: its number adds nothing
  out <- list()
  for(i in seq_len(nrow(combos))) {
    block <- data[model == combos$model[i] & classes == combos$classes[i], , drop=FALSE]
    block <- block[!duplicated(block$id), , drop=FALSE] # k rows per case, Class constant
    block <- block[order(block$id), , drop=FALSE]       # back to the original case order
    out[[tidyLPAName(combos$model[i], combos$classes[i], oneModel)]] <-
      as.vector(as.integer(as.factor(block$Class)))
  }
  out
}


## selectTidyLPAClusters ----
# Pick one solution out of a tidyLPA object.
# which: NULL (the first one), an index, or a name such as "model_1(3)".
selectTidyLPAClusters <- function(x, which=NULL) {
  sets <- listTidyLPAClusters(x)
  if(is.null(which)) {
    if(length(sets) > 1)
      warning("the tidyLPA object holds ", length(sets),
              " solutions and the first one was taken; state which to choose another: ",
              paste(names(sets), collapse=", "), call.=FALSE)
    return(sets[[1]])
  }
  if(is.numeric(which)) {
    if(length(which) != 1 || which < 1 || which > length(sets))
      stop("which must be a single index between 1 and ", length(sets))
    return(sets[[which]])
  }
  if(is.character(which) && length(which) == 1) {
    if(!(which %in% names(sets)))
      stop("which=\"", which, "\" is not among the solutions of the tidyLPA object: ",
           paste(names(sets), collapse=", "))
    return(sets[[which]])
  }
  stop("which must be NULL, a single index or a single name")
}


## extractClusters ----
# Extract cluster assignments from various clustering objects
# Supports: kmeans, poLCA, tidyLPA, looking4clusters, hclust (via cutree), or direct vector/factor
extractClusters <- function(x, sourceData=NULL, which=NULL) {
  # looking4clusters object (checked first: it is a classed list)
  if(inherits(x, "looking4clusters")) {
    return(selectL4CClusters(x, which))
  }

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
    return(selectTidyLPAClusters(x, which))
  }

  # hclust result (should be cut first, but allow passing pre-cut vector)
  if(inherits(x, "hclust")) {
    stop("hclust object should be cut first with cutree(). Pass the result of cutree() instead.")
  }

  # data.frame with single column (several columns are handled by extractClusterList)
  if(is.data.frame(x)) {
    if(ncol(x) == 1) return(as.vector(as.integer(as.factor(as.character(x[[1]])))))
    stop("clusters data.frame must have exactly one column")
  }

  # Unknown type
  stop("Don't know how to extract clusters from object of class ", paste(class(x), collapse=", "))
}


## extractClusterList ----
# Reduce any accepted input to a named list of cluster vectors, so that holders of
# several clusterizations (a data.frame, a looking4clusters object, a tidyLPA object
# estimated for several numbers of profiles) and single ones (a vector, kmeans,
# poLCA, a tidyLPA object of one solution) all travel the same path.
# Names are the ones each clusterization carries; an empty name means it has none of
# its own, and the caller supplies it.
extractClusterList <- function(x, sourceData=NULL, which=NULL) {
  # looking4clusters with no choice made: every clusterization it holds
  if(inherits(x, "looking4clusters") && is.null(which)) {
    return(listL4CClusters(x))
  }

  # tidyLPA with no choice made: every solution it holds
  if(inherits(x, "tidyLPA") && is.null(which)) {
    return(listTidyLPAClusters(x))
  }

  # data.frame of several columns, one clusterization each
  if(is.data.frame(x) && ncol(x) > 1) {
    out <- lapply(x, function(column) as.vector(as.integer(as.factor(as.character(column)))))
    names(out) <- names(x)
    return(out)
  }

  # data.frame of one column: keeps its name, which the caller may override
  if(is.data.frame(x) && ncol(x) == 1) {
    out <- list(as.vector(as.integer(as.factor(as.character(x[[1]])))))
    names(out) <- names(x)
    return(out)
  }

  # anything else resolves to a single unnamed clusterization
  out <- list(extractClusters(x, sourceData, which))
  names(out) <- ""
  out
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
# Add one or several cluster columns to a netCoin object from surScat
# clusters: vector, factor, data.frame of one clusterization per column, or clustering
#           object (kmeans, poLCA, tidyLPA, looking4clusters)
# sourceData: optional data.frame for validation that clusters are in correct order
# which: clusterization to take from a looking4clusters object. When it is left
#        empty, every clusterization held by the object is added as its own column.
# maxGroups: clusterizations with more groups than this are skipped with a warning,
#            which catches continuous variables passed in by mistake. NULL disables it.
addClusters <- function(scatObj, clusters, name=NULL, sort=TRUE, weight=NULL, sourceData=NULL,
                        which=NULL, maxGroups=25) {
  if(!inherits(scatObj, "netCoin"))
    stop("scatObj must be a netCoin object returned by surScat")

  if(is.null(clusters))
    stop("clusters must be provided")

  sets <- extractClusterList(clusters, sourceData, which)
  several <- length(sets) > 1
  renamed <- renamedTo <- character(0) # clusterizations the node table already held

  for(i in seq_along(sets)) {
    own <- names(sets)[i]
    if(is.na(own)) own <- ""

    # Skip what has too many groups to be a clusterization
    nGroups <- length(unique(sets[[i]]))
    if(!is.null(maxGroups) && nGroups > maxGroups) {
      warning("clusterization ", if(nzchar(own)) paste0("\"", own, "\" ") else "",
              "has ", nGroups, " groups, more than maxGroups=", maxGroups,
              ", and was not added. Raise maxGroups to add it anyway.", call.=FALSE)
      next
    }

    # Its own name, prefixed by name when it was given; name alone when there is only one
    column <- if(is.null(name)) {
                if(nzchar(own)) own else NULL
              } else if(several) {
                paste0(name, ".", own)
              } else name
    if(!is.null(column)) {
      asked  <- column
      column <- make.unique(c(names(scatObj$nodes), column))[length(names(scatObj$nodes))+1]
      if(!identical(column, asked)) { # a column of that name is already there
        renamed   <- c(renamed, asked)
        renamedTo <- c(renamedTo, column)
      }
    }

    scatObj <- addOneCluster(scatObj, sets[[i]], name=column, sort=sort,
                             weight=weight, sourceData=sourceData)
  }

  if(length(renamed))
    warning("the node table already held ", nameClash(renamed, renamedTo),
            ". State name to tell them apart, or drop the previous columns to avoid ",
            "repeating a clusterization.", call.=FALSE)

  return(scatObj)
}


## addOneCluster ----
# Add a single cluster column to a netCoin object. Expects clusters already
# reduced to a plain vector by extractClusters.
addOneCluster <- function(scatObj, clusters, name=NULL, sort=TRUE, weight=NULL, sourceData=NULL) {
  # Extract language setting for labels
  language <- scatObj$options$language
  if(is.null(language)) language <- "en"

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
        # Store cluster values without prefix (will add prefix as labels later)
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
                  nrow(scatObj$nodes), " nodes. Clusters must hold one assignment per node, ",
                  "or one per case when the nodes are patterns collapsed by vPatterns. ",
                  "A clustering object holding several solutions at once, such as a tidyLPA ",
                  "estimated for several numbers of profiles, is taken apart into one ",
                  "clusterization per solution, so it need not be reduced beforehand.",
                  sampledNote(scatObj)))
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

  # Sort clusters if requested, by their mean coordinate on the first axis as drawn.
  # The layout is read exactly, as $layout would partially match the $layouts that
  # addAxes adds; and a surScat object keeps no layout of its own, since its
  # coordinates reach the nodes as the fx/fy columns instead.
  lay <- scatObj[["layout"]]
  layout_first <- if(!is.null(lay) && !is.null(ncol(lay)) && ncol(lay) > 0) lay[,1]
                  else scatObj$nodes$fx
  if(sort && !is.null(layout_first) && length(layout_first) == length(cl)) {
    # For sorting, use ALL cluster values (including non-modal), but compute mean only for modals
    # Non-modal clusters inherit the mean of their pattern's modal cluster
    cf <- factor(cl, levels=unique_cl)
    if(is.null(weight)) {
      mu_modal <- tapply(layout_first, cf, mean)
    } else {
      if(length(weight) != length(cl))
        stop("weight must have same length as clusters")
      mu_modal <- tapply(seq_along(cl), cf, function(k) weighted.mean(layout_first[k], weight[k]))
    }

    # Assign modal values' means to all cluster values (for sorting purposes)
    mu_all <- numeric(length(all_cluster_values))
    names(mu_all) <- all_cluster_values
    for(i in seq_along(unique_cl)) {
      # Find all patterns with this modal value
      patterns_with_modal <- which(cl == unique_cl[i])
      if(length(patterns_with_modal) > 0) {
        mu_all[unique_cl[i]] <- mu_modal[unique_cl[i]]
      }
    }

    # Renumber clusters based on sorted order (use all_cluster_values so all values have a mapping)
    order_idx <- order(mu_all[as.character(all_cluster_values)])
    all_cluster_values_sorted <- all_cluster_values[order_idx]

    # Map cl values: replace old values with new numbering
    cl_new <- rep(NA_character_, length(cl))
    for(i in seq_along(all_cluster_values_sorted)) {
      cl_new[cl == all_cluster_values_sorted[i]] <- as.character(i)
    }
    cl <- cl_new

    if(!is.null(clusterDisplay)) {
      # Map display values: extract cluster numbers and remap them (without prefix)
      clusterDisplay_new <- character(length(clusterDisplay))

      for(j in seq_along(clusterDisplay)) {
        # Split the cluster values (they are just numbers without prefix: "1|2|3")
        cluster_nums <- strsplit(clusterDisplay[j], "\\|")[[1]]

        # Remap the numbers based on sorting
        new_nums <- character(length(cluster_nums))
        for(k in seq_along(cluster_nums)) {
          idx_in_sorted <- which(all_cluster_values_sorted == cluster_nums[k])
          if(length(idx_in_sorted) > 0) {
            new_nums[k] <- as.character(idx_in_sorted)
          }
        }

        # Reconstruct without prefix (will add prefix as labels later)
        clusterDisplay_new[j] <- paste(new_nums, collapse="|")
      }
      clusterDisplay <- clusterDisplay_new
    }
    n_clusters <- length(all_cluster_values_sorted)
    unique_cl <- as.character(seq_len(n_clusters))
  }

  # Use display values (with multiple clusters) if available, otherwise use cl values
  cl_values <- if(!is.null(clusterDisplay)) clusterDisplay else cl

  # Create ordered factor with localized labels
  # Values are cluster numbers (e.g., "1|2|3"), labels have localized prefix (e.g., "Group: 1|Group: 2|Group: 3")
  # Group numbers are padded with zeros to the width of the largest one, so that
  # sorting them as text, here or further downstream, keeps them in numeric order.
  unique_cl_values <- unique(cl_values)

  # Labels are prefixed with the name of the clusterization they come from, dropping the
  # trailing number of groups and capitalized, so "hclust(10)" reads "Hclust: 01".
  # A clusterization added twice reaches this point renamed by make.unique, so the
  # suffix it appends is dropped as well and "hclust(10).1" reads "Hclust: 01" too.
  # The localized word is kept for the case of a column with no name of its own.
  groupWord <- sub("\\([0-9]+\\)(\\.[0-9]+)?$", "", name)
  if(nzchar(groupWord))
    groupWord <- paste0(toupper(substring(groupWord, 1, 1)), substring(groupWord, 2))
  else
    groupWord <- getByLanguage(groupList, language)

  splitGroups <- function(val) as.integer(strsplit(val, "\\|")[[1]])
  width <- max(nchar(as.character(unlist(lapply(unique_cl_values, splitGroups)))), 1L)
  padded <- function(val) sprintf(paste0("%0", width, "d"), splitGroups(val))

  # Levels must follow the group number, not the order the values happen to appear in
  unique_cl_values <- unique_cl_values[order(vapply(unique_cl_values,
                                                    function(val) paste(padded(val), collapse="|"),
                                                    character(1)))]

  # Create labels by adding prefix to each cluster number
  labels_with_prefix <- vapply(unique_cl_values, function(val) {
    # Pad each number, add prefix to each, rejoin
    paste(paste0(groupWord, ": ", padded(val)), collapse="|")
  }, character(1), USE.NAMES=FALSE)

  cl_factor <- factor(cl_values, levels=unique_cl_values, labels=labels_with_prefix, ordered=TRUE)

  # Add the new cluster column
  scatObj$nodes[[name]] <- cl_factor

  # Record it as a clusterization, so that replaceClusters knows what to remove
  # without having to guess it from the column names
  attr(scatObj, "clusterColumns") <- union(attr(scatObj, "clusterColumns"), name)

  # Store collapse report if applicable
  if(!is.null(collapseReport)) {
    attr(scatObj, paste0("collapse_", name)) <- collapseReport
  }

  return(scatObj)
}


## replaceClusters ----
# Replace cluster columns in a netCoin object from surScat
# Removes all existing cluster columns and adds the new one
replaceClusters <- function(scatObj, clusters, name=NULL, sort=TRUE, weight=NULL, sourceData=NULL,
                            which=NULL, maxGroups=25) {
  if(!inherits(scatObj, "netCoin"))
    stop("scatObj must be a netCoin object returned by surScat")

  # Find and remove existing cluster columns. Both surScat and addClusters record the
  # columns they create, so there is no need to tell them apart by their name: guessing
  # it from a pattern used to delete variables of the study that happened to be named
  # like "Income(2020)".
  clusterCols <- intersect(attr(scatObj, "clusterColumns"), names(scatObj$nodes))

  # Objects built before clusterColumns was recorded fall back to the name of the
  # column surScat used to create
  if(!length(clusterCols))
    clusterCols <- grep("^Groups?\\([0-9]", names(scatObj$nodes), value=TRUE)

  if(length(clusterCols) > 0) {
    scatObj$nodes[clusterCols] <- NULL
    attr(scatObj, "clusterColumns") <- setdiff(attr(scatObj, "clusterColumns"), clusterCols)
  }

  # If name not specified, try to reuse the first removed column name. Not when several
  # clusterizations are added at once, since each one names its own column.
  addingSeveral <- (inherits(clusters, "looking4clusters") && is.null(which)) ||
                   (is.data.frame(clusters) && ncol(clusters) > 1)
  if(is.null(name) && length(clusterCols) > 0 && !addingSeveral) {
    name <- clusterCols[1]
  }

  # Use addClusters for the rest
  addClusters(scatObj, clusters, name=name, sort=sort, weight=weight, sourceData=sourceData,
              which=which, maxGroups=maxGroups)
}
