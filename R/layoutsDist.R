## asPlane ----
# Any entry of the layouts of a scattergram (a matrix or a data.frame) left as a numeric
# matrix. The planes of surScat always hold two columns, but nothing that follows depends on
# it, so no more is asked than that the two of them have the same shape.
asPlane <- function(z, nm, who) {
  if(is.data.frame(z)) z <- as.matrix(z)
  if(!is.matrix(z)) z <- as.matrix(z)
  if(!is.numeric(z)) stop(who, ": ", nm, " does not hold numeric coordinates")
  if(!ncol(z)) stop(who, ": ", nm, " holds no column of coordinates")
  z
}


## normPlane ----
# Bring the centroid to the origin and divide by the root mean square radius, which is the
# typical distance of the points to the centre. The divisor is one for the whole plane rather
# than one per column: dividing every axis by its own spread would stretch the second one
# until it matched the first, and in a reduction such as pca that undoes precisely what the
# reduction states, that PC1 carries more variation than PC2. The plane keeps its shape and
# only changes size. With weights, centroid and radius are computed weighted.
normPlane <- function(z, weights=NULL, who="normPlane") {
  if(is.null(weights)) {
    ctr <- colMeans(z)
    z   <- sweep(z, 2, ctr)
    s   <- sqrt(sum(z^2)/nrow(z))
  } else {
    sw  <- sum(weights)
    if(sw <= 0) stop(who, ": the weights add up to ", sw)
    ctr <- colSums(z*weights)/sw
    z   <- sweep(z, 2, ctr)
    s   <- sqrt(sum(rowSums(z^2)*weights)/sw)
  }
  if(!is.finite(s) || s <= 0)
    stop(who, ": the plane has no size (every node sits on the same point), so it cannot be",
         " brought to a common scale")
  z <- z/s
  attr(z, "center") <- ctr
  attr(z, "scale")  <- s
  z
}


## canvasPlane ----
# Every column brought to [0, 1] by its minimum and its maximum, which is what rD3plot does
# when it sets the domain of each axis to the range of the coordinates and its own range to
# the width or the height of the canvas. The 4% margin the viewer adds is an equal dilation at
# both ends of each axis, so it alters none of this.
# The canvas is taken square, [0,1] x [0,1], because the shape of the window is not known from
# R; worth remembering that on screen the cloud is stretched further by the proportions of the
# window, and changes as it is resized.
# Weights play no part: a node of little weight takes up the canvas like any other, and the
# viewer counts it in for the range of the axis. A single extreme point therefore sets the
# scale of the whole plane, here as there.
canvasPlane <- function(z, who="canvasPlane") {
  rng <- apply(z, 2, range)
  d   <- rng[2,] - rng[1,]
  if(any(!is.finite(d)) || any(d <= 0))
    stop(who, ": one axis has no range (every node holds the same value), so it cannot be",
         " brought to the canvas")
  z <- sweep(sweep(z, 2, rng[1,]), 2, d, "/")
  attr(z, "scale") <- d
  z
}


## procrustesRotate ----
# Orthogonal Procrustes: with the two planes already centred and at the same scale, find the
# orthogonal matrix R minimizing the distance from a*R to b, and hand back a rotated.
# Reflection is allowed alongside rotation, because in a reduction the mirror says nothing: a
# reflected map shows the same structure.
# No scaling factor of its own is added, though the classical Procrustes fit carries one:
# normPlane has already levelled the sizes, and that factor would always be a shrinking
# towards the centre which would break the symmetry of the measure, and a matrix of distances
# between planes needs it.
procrustesRotate <- function(a, b, weights=NULL) {
  m  <- if(is.null(weights)) crossprod(a, b) else crossprod(a*weights, b)
  sv <- svd(m)
  r  <- sv$u %*% t(sv$v)                 # a %*% r is as close to b as a can get
  out <- a %*% r
  attr(out, "rotation") <- r
  out
}


## layoutDist ----
# The distance between two planes of a scattergram: how far one and the same node travels
# from one to the other, added up over the nodes and divided by how many they are.
layoutDist <- function(a, b, align=c("cloud", "procrustes", "canvas", "none"), weights=NULL,
                       normalize=TRUE, who="layoutDist") {
  align <- match.arg(align)
  a <- asPlane(a, "a", who)
  b <- asPlane(b, "b", who)
  if(nrow(a) != nrow(b))
    stop(who, ": the planes hold ", nrow(a), " and ", nrow(b),
         " rows; the same ones are needed, one per node")
  if(ncol(a) != ncol(b))
    stop(who, ": the planes hold ", ncol(a), " and ", ncol(b), " columns; the same ones are needed")

  ok <- complete.cases(a) & complete.cases(b)
  if(!is.null(weights)) {
    if(!is.numeric(weights)) stop(who, ": the weights must be numeric")
    if(length(weights) == 1L) weights <- rep(weights, nrow(a))
    if(length(weights) != nrow(a))
      stop(who, ": weights holds ", length(weights), " values and the planes ", nrow(a), " rows")
    if(any(weights < 0, na.rm=TRUE)) stop(who, ": there are negative weights")
    ok <- ok & !is.na(weights)
  }
  if(!any(ok)) stop(who, ": no comparable row is left")

  # the alignment is computed on the rows about to be compared, not on the ones dropped
  a <- a[ok, , drop=FALSE]
  b <- b[ok, , drop=FALSE]
  w <- if(is.null(weights)) NULL else weights[ok]

  scales <- NULL
  rot    <- NULL
  if(align == "canvas") {
    a <- canvasPlane(a, who)
    b <- canvasPlane(b, who)
    scales <- list(a=attr(a, "scale"), b=attr(b, "scale"))
  } else if(align != "none") {
    a <- normPlane(a, w, who)
    b <- normPlane(b, w, who)
    scales <- list(a=attr(a, "scale"), b=attr(b, "scale"))
    if(align == "procrustes") {
      a   <- procrustesRotate(a, b, w)
      rot <- attr(a, "rotation")
    }
  }

  d   <- sqrt(rowSums((a-b)^2))
  out <- if(is.null(w)) sum(d) else sum(d*w)
  if(normalize) {
    div <- if(is.null(w)) sum(ok) else sum(w)
    if(div <= 0) stop(who, ": there is nothing to divide by (the weights add up to ", div, ")")
    out <- out/div
  }

  attr(out, "rows")    <- sum(ok)
  attr(out, "dropped") <- sum(!ok)
  attr(out, "align")   <- align
  if(!is.null(scales)) attr(out, "scales") <- scales
  if(!is.null(rot) && ncol(a) == 2L) {
    attr(out, "angle")      <- atan2(rot[2,1], rot[1,1])*180/pi
    attr(out, "reflection") <- det(rot) < 0
  }
  out
}


## layoutList ----
# The planes to be compared, as a named list. A netCoin object is read through
# currentLayouts, which also covers the object that has not been through addAxes and holds a
# single plane in the fx/fy columns of its node table; a plain list of matrices is taken as
# it comes, so that planes held in no object can be compared as well.
# which: NULL (every one of them), indices or names.
layoutList <- function(x, which=NULL, who="layoutList") {
  if(inherits(x, "netCoin")) L <- currentLayouts(x)
  else if(is.list(x) && !is.data.frame(x)) L <- x
  else stop(who, ": x must be a netCoin object from surScat or a list of planes")

  if(!length(L)) stop(who, ": there is no plane at all")
  if(is.null(names(L)) || any(!nzchar(names(L)))) names(L) <- paste0("layout", seq_along(L))

  if(!is.null(which)) {
    if(is.numeric(which)) {
      if(any(is.na(which) | which < 1 | which > length(L)))
        stop(who, ": which must hold indices between 1 and ", length(L))
      L <- L[which]
    } else if(is.character(which)) {
      bad <- setdiff(which, names(L))
      if(length(bad))
        stop(who, ": there is no plane called \"", paste(bad, collapse="\", \""),
             "\". Available: ", paste(names(L), collapse=", "))
      L <- L[which]
    } else stop(who, ": which must be NULL, indices or names")
  }
  if(length(L) < 2L)
    stop(who, ": at least two planes are needed to compare them, and there ",
         if(length(L) == 1L) "is 1" else paste("are", length(L)))
  L
}


## layoutWeights ----
# The weight of every node, settled just as matchSetup settles it in clusterMatch, so that
# both faces of the object -clusters and layouts- count alike: where the nodes are patterns,
# each one weighs by default the cases it stands for.
# weights: NULL (the criterion above), the name of a numeric node column, or a vector holding
# one value per node or per case. weights=1 forces the plain count of nodes.
layoutWeights <- function(x, weights, n, who) {
  nodes  <- if(is.list(x) && is.data.frame(x$nodes)) x$nodes else NULL
  idx    <- attr(x, "caseToPattern")
  nCases <- if(is.null(idx)) NULL else length(idx)
  caseW  <- attr(x, "caseWeight")
  if(!is.null(caseW) && !is.null(nCases) && length(caseW) != nCases) caseW <- NULL

  if(is.null(weights)) {
    if(is.null(idx)) return(list(w=NULL, source="node"))
    cw <- if(is.null(caseW)) rep(1, nCases) else caseW
    w  <- as.vector(tapply(cw, factor(idx, levels=seq_len(n)), sum))
    w[is.na(w)] <- 0
    return(list(w=w, source=if(is.null(caseW)) "cases per pattern"
                            else "cases per pattern (caseWeight)"))
  }

  if(is.character(weights) && length(weights) == 1L) {
    if(is.null(nodes) || !(weights %in% names(nodes)))
      stop(who, ": there is no node column called \"", weights, "\"")
    w   <- nodes[[weights]]
    src <- paste0("column \"", weights, "\"")
  } else {
    w   <- weights
    src <- "weights"
  }
  if(!is.numeric(w)) stop(who, ": the weights must be numeric")
  if(length(w) == 1L) w <- rep(w, n)
  if(length(w) != n) {
    # weights given case by case but a comparison by nodes: they are added up
    if(!is.null(nCases) && length(w) == nCases) {
      w   <- as.vector(tapply(w, factor(idx, levels=seq_len(n)), sum))
      w[is.na(w)] <- 0
      src <- paste0(src, ", added up by pattern")
    } else
      stop(who, ": weights holds ", length(w), " values and there are ", n, " nodes",
           if(!is.null(nCases)) paste0(" (or ", nCases, " cases)"), "", sampledNote(x))
  }
  # ones amount to no weighting at all: it is said so and the unweighted path is taken, which
  # yields the same and announces no weighting where there is none. This is what weights=1
  # achieves on pattern nodes, counting nodes rather than cases.
  if(all(w == 1, na.rm=TRUE)) return(list(w=NULL, source="node"))
  list(w=w, source=src)
}


## layoutsDist ----
# The distance between every pair of planes of a scattergram, as the dist object hclust,
# cmdscale and the like expect.
layoutsDist <- function(x, align=c("cloud", "procrustes", "canvas", "none"), weights=NULL,
                        which=NULL, normalize=TRUE, who="layoutsDist") {
  align <- match.arg(align)
  L  <- layoutList(x, which, who)
  nm <- names(L)
  k  <- length(L)

  n <- unique(vapply(L, nrow, 0L))
  if(length(n) > 1L)
    stop(who, ": the planes hold a different number of rows (", paste(sort(n), collapse=", "),
         "); every one of them must hold one per node")
  wl <- layoutWeights(x, weights, n, who)

  v    <- numeric(k*(k-1)/2)
  rows <- numeric(length(v))
  ang  <- matrix(NA_real_, k, k, dimnames=list(nm, nm))
  ref  <- matrix(NA, k, k, dimnames=list(nm, nm))
  diag(ang) <- 0
  diag(ref) <- FALSE

  # the order of a dist is that of the lower triangle swept column by column:
  # (2,1), (3,1), ..., (k,1), (3,2), ...
  p <- 0L
  for(j in seq_len(k-1L)) for(i in (j+1L):k) {
    p <- p+1L
    z <- layoutDist(L[[i]], L[[j]], align=align, weights=wl$w, normalize=normalize, who=who)
    v[p]    <- z
    rows[p] <- attr(z, "rows")
    if(!is.null(attr(z, "angle"))) {
      ang[i,j] <-  attr(z, "angle")
      ang[j,i] <- -attr(z, "angle")
      ref[i,j] <- ref[j,i] <- attr(z, "reflection")
    }
  }

  out <- structure(v, class="dist", Size=k, Labels=nm, Diag=FALSE, Upper=FALSE, method=align)
  attributes(rows) <- attributes(out)
  attr(out, "rows")      <- rows
  attr(out, "align")     <- align
  attr(out, "normalize") <- normalize
  attr(out, "weights")   <- wl$w
  attr(out, "wSource")   <- wl$source
  if(align == "procrustes") {
    attr(out, "angle")      <- ang
    attr(out, "reflection") <- ref
  }
  # the weighting changes what the figures stand for, so it is stated
  if(!is.null(wl$w))
    message(who, ": weighted by ", wl$source, " (N = ", format(sum(wl$w), big.mark=""), ")")
  out
}
