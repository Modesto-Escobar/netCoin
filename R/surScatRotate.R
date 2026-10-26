## surScat and the rotations of its own plane ----
# layouts also takes the rotations of the components surScat extracts. They are not
# reductions of looking4clusters: they come from the very factors of the main plane, which
# they leave as they are and complement with another plane, added as any other.

rotationLayouts <- c("varimax", "quartimax", "equamax")

# "equimax" is read as "equamax", as some programs spell it
rotationName <- function(x) {
  x <- tolower(x)
  x[x == "equimax"] <- "equamax"
  x
}


## rotationLayout ----
# The rotations surScat's layouts asks for, in the order asked. They are not part of
# layouts=TRUE, which stands for the reductions alone, but of layouts="all", which stands
# for the reductions and the rotations.
rotationLayout <- function(layouts) {
  if(!is.character(layouts)) return(character(0))
  asked <- rotationName(layouts)
  if("all" %in% asked) asked <- c(asked[asked != "all"], rotationLayouts)
  unique(asked[asked %in% rotationLayouts])
}


## rotationAxes ----
# The plane of the first nfactors components of a principal components analysis once rotated,
# one row per case, as list(varimax=matrix) ready for addAxes, or NULL when it cannot be made.
# ff: the components surScat extracted; vm: the standardized active variables they come from;
# eig: the full spectrum of eigenvalues; weight: that of the cases (NULL when there is none).
# The components are brought to unit variance before rotating (the loadings are then the
# correlations of the variables with them), as the rotated components of the usual packages
# are. The rotated ones are ordered by the variance they account for, and signed so that
# their loadings add up to a positive number, as the rotation fixes neither.
# With nfactors=2 the rotation turns the plane of the first two components within itself,
# so the cloud is the same one at another angle: it is from three components on that the
# plane differs from the main one. Then one pair of components is only a projection of the
# solution, so there is a plane for every pair ("varimax 1-2", "varimax 1-3"...).
# quiet: layouts="all" asks for every plane that can be made, so those that cannot are
# left out without a word.
# With two components every rotation turns the same plane, so layouts="all" keeps varimax, the
# most common, and the others only if layouts stated them (asked).
rotationAxes <- function(ff, vm, eig, weight, nfactors, rotations, type, quiet=FALSE,
                         asked=character(0)) {
  if(!length(rotations)) return(NULL)
  if(type != "pca") {
    if(!quiet)
      warning("rotations are only made of principal components: ",
              paste(rotations, collapse=", "), " ignored with type=\"", type, "\"", call.=FALSE)
    return(NULL)
  }
  k <- min(nfactors, ncol(ff))
  if(k < 2) {
    if(!quiet)
      warning("a rotation needs at least two components: ",
              paste(rotations, collapse=", "), " ignored with nfactors=", nfactors, call.=FALSE)
    return(NULL)
  }

  if(quiet && k == 2) rotations <- intersect(rotations, c("varimax", asked))

  Z <- sweep(ff[, seq_len(k), drop=FALSE], 2, sqrt(eig[seq_len(k)]), "/")
  L <- corr(vm, Z, weight=weight)
  out <- list()
  for(rot in rotations) {
    # all three are orthogonal rotations of the loadings with Kaiser's normalization, and L %*% Th
    # gives the rotated loadings, so the scores are rotated by the same matrix
    R <- tryCatch(switch(rot,
           varimax   = stats::varimax(L, normalize=TRUE),
           quartimax = GPArotation::quartimax(L, normalize=TRUE),
           # Crawford-Ferguson family with kappa = components/(2*variables)
           equamax   = GPArotation::cfT(L, kappa=k/(2*nrow(L)), normalize=TRUE)),
         error=function(e) {
           warning("the ", rot, " rotation could not be made: ", conditionMessage(e), call.=FALSE)
           NULL
         })
    if(is.null(R)) next
    Th <- if(rot == "varimax") R$rotmat else R$Th
    scores <- Z %*% Th
    load   <- unclass(R$loadings)
    ord    <- order(colSums(load^2), decreasing=TRUE)
    scores <- scores[, ord, drop=FALSE]
    load   <- load[, ord, drop=FALSE]
    flip   <- ifelse(colSums(load) < 0, -1, 1)
    scores <- sweep(scores, 2, flip, "*")
    pct    <- colSums(load^2)/ncol(vm)*100 # of the variance of the variables, as for the main plane
    colnames(scores) <- sprintf("RC%d (%.1f%%)", seq_len(k), pct)
    rownames(scores) <- rownames(ff)
    if(k == 2) out[[rot]] <- scores
    else { # a plane for every pair of rotated components, as one pair alone is a mere projection
      pairs <- utils::combn(k, 2)
      for(j in seq_len(ncol(pairs)))
        out[[paste0(rot, " ", pairs[1, j], "-", pairs[2, j])]] <- scores[, pairs[, j], drop=FALSE]
    }
  }
  out
}
