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
# nperm: how many times the rows of one plane of each pair are shuffled, which breaks every
#     link between the nodes of the two planes and so shows what distance chance alone yields.
#     From it come the mean distance by chance, the distance relative to it (0 the same
#     plane, 1 no more alike than by chance) and the probability of a distance as small as the
#     one observed if the planes had nothing to do with each other.
# value: which of those the dist holds; the others go along as attributes of the same shape.
# digits: decimals shown when printed. The figures themselves are kept whole.
layoutsDist <- function(x, align=c("cloud", "procrustes", "canvas", "none"), weights=NULL,
                        which=NULL, normalize=TRUE, nperm=0,
                        value=c("distance", "relative", "p", "chance"), digits=3,
                        who="layoutsDist") {
  align <- match.arg(align)
  value <- match.arg(value)
  if(!is.numeric(digits) || length(digits) != 1L || is.na(digits) || digits < 0)
    stop(who, ": digits must be a whole number, 0 or more")
  if(!is.numeric(nperm) || length(nperm) != 1L || is.na(nperm) || nperm < 0 ||
     nperm != round(nperm))
    stop(who, ": nperm must be a whole number, 0 or more")
  if(value != "distance" && nperm == 0)
    stop(who, ": value=\"", value, "\" is worked out by shuffling the planes, so it needs",
         " nperm; 999 is a usual choice")
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
  pv   <- chance <- rep(NA_real_, length(v))
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
      # a rotation is undone by turning back, but a reflection is its own inverse: the line it
      # mirrors across is the same whichever plane is laid over the other, and so is its angle
      ang[i,j] <- attr(z, "angle")
      ang[j,i] <- if(attr(z, "reflection")) attr(z, "angle") else -attr(z, "angle")
      ref[i,j] <- ref[j,i] <- attr(z, "reflection")
    }

    if(nperm > 0) {
      # shuffled among the rows the observed distance was measured on, so that every
      # permutation compares the same nodes; each weight stays with the node of the plane
      # left in place
      a  <- asPlane(L[[i]], nm[i], who)
      b  <- asPlane(L[[j]], nm[j], who)
      ok <- complete.cases(a) & complete.cases(b)
      if(!is.null(wl$w)) ok <- ok & !is.na(wl$w)
      a  <- a[ok, , drop=FALSE]
      b  <- b[ok, , drop=FALSE]
      w  <- if(is.null(wl$w)) NULL else wl$w[ok]
      sim <- vapply(seq_len(nperm), function(r)
               c(layoutDist(a[sample.int(nrow(a)), , drop=FALSE], b, align=align, weights=w,
                            normalize=normalize, who=who)), 0)
      chance[p] <- mean(sim)
      pv[p]     <- (1 + sum(sim <= z))/(nperm + 1)
    }
  }

  # every figure a dist, which hclust and cmdscale take as it is, and printed with a fixed
  # number of decimals, as print.dist falls back on scientific notation once a p-value is small
  asDist <- function(y, class=c("layoutsDist", "dist"))
    structure(y, class=class, Size=k, Labels=nm, Diag=FALSE, Upper=FALSE, method=align,
              digits=digits)
  dists  <- list(distance=asDist(v))
  if(nperm > 0) {
    dists$relative <- asDist(v/chance)
    dists$p        <- asDist(pv)
    dists$chance   <- asDist(chance)
  }
  out <- dists[[value]]
  attr(out, "value") <- value
  for(s in setdiff(names(dists), value)) attr(out, s) <- dists[[s]]
  if(nperm > 0) attr(out, "nperm") <- nperm
  attributes(rows) <- attributes(asDist(0, "dist"))
  attr(rows, "digits") <- NULL
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


## print.layoutsDist ----
# The lower triangle, as print.dist lays it out, with every figure given the same number of
# decimals and never in scientific notation. digits defaults to what layoutsDist was told.
print.layoutsDist <- function(x, digits=attr(x, "digits"), ...) {
  if(is.null(digits)) digits <- 3
  m   <- as.matrix(x)
  txt <- formatC(m, format="f", digits=digits)
  txt[is.na(m)] <- "NA"
  txt[upper.tri(txt, diag=TRUE)] <- ""
  txt <- txt[-1, -ncol(txt), drop=FALSE]
  print(noquote(txt), right=TRUE)
  invisible(x)
}


## planeName ----
# The name of one plane among those available, given as its name or as its index.
planeName <- function(to, available, who) {
  if(length(to) != 1L || is.na(to)) stop(who, ": to must be a single name or index")
  if(is.numeric(to)) {
    if(to < 1 || to > length(available))
      stop(who, ": to must be an index between 1 and ", length(available))
    return(available[to])
  }
  if(!is.character(to)) stop(who, ": to must be a name or an index")
  if(!(to %in% available))
    stop(who, ": there is no plane called \"", to, "\". Available: ",
         paste(available, collapse=", "))
  to
}


## layoutsAlign ----
# The planes of a scattergram laid over one another, as layoutsDist lays them to measure how
# far apart they are, but handing back the coordinates themselves.
# to: the plane the rest are brought to, whose units and orientation the result is read in.
#     It is left as it was, and need not be among which.
# gpa: with align="procrustes", rotate every plane towards the mean of them all (generalized
#     Procrustes analysis) rather than towards to, which then only sets how the whole set is
#     finally turned and measured. No plane is favoured over the others this way.
# suffix: NULL replaces the planes of a netCoin object by the aligned ones; a string keeps the
#     originals and adds the aligned ones, their names followed by it (as in "sepal'").
# The alignment is estimated on the rows complete in every plane involved, and then applied
# to every row of each plane, so a node missing from one plane keeps its place in the others.
layoutsAlign <- function(x, align=c("procrustes", "cloud", "canvas", "none"), to=1, gpa=FALSE,
                         weights=NULL, which=NULL, suffix=NULL, tol=1e-10, maxit=100,
                         who="layoutsAlign") {
  align <- match.arg(align)
  if(!is.null(suffix) && (!is.character(suffix) || length(suffix) != 1L || !nzchar(suffix)))
    stop(who, ": suffix must be a single non-empty string")
  if(gpa && align != "procrustes")
    stop(who, ": gpa rotates the planes, so it needs align=\"procrustes\"")

  Lall <- layoutList(x, NULL, who)
  sel  <- names(layoutList(x, which, who))
  ref  <- planeName(to, names(Lall), who)
  use  <- union(sel, ref)
  L    <- setNames(lapply(use, function(nm) asPlane(Lall[[nm]], nm, who)), use)

  n <- unique(vapply(L, nrow, 0L))
  if(length(n) > 1L)
    stop(who, ": the planes hold a different number of rows (", paste(sort(n), collapse=", "),
         "); every one of them must hold one per node")
  p <- unique(vapply(L, ncol, 0L))
  if(length(p) > 1L)
    stop(who, ": the planes hold a different number of columns (", paste(sort(p), collapse=", "),
         "); they can only be laid over one another with the same ones")
  wl <- layoutWeights(x, weights, n, who)

  # the rows every fit is estimated on: the same for all planes, so that they are comparable
  ok <- Reduce(`&`, lapply(L, complete.cases))
  if(!is.null(wl$w)) ok <- ok & !is.na(wl$w)
  if(!any(ok)) stop(who, ": no row is complete in every plane")
  w <- if(is.null(wl$w)) NULL else wl$w[ok]

  # what takes each plane to the common footing: u = ((z - center)/scale) %*% rotation, with
  # the scale per column under canvas and one for the whole plane otherwise
  center <- scale <- rot <- setNames(vector("list", length(use)), use)
  for(nm in use) {
    z <- L[[nm]][ok, , drop=FALSE]
    if(align == "canvas") {
      rng <- apply(z, 2, range)
      d   <- rng[2,] - rng[1,]
      if(any(!is.finite(d)) || any(d <= 0))
        stop(who, ": one axis of \"", nm, "\" has no range, so it cannot be brought to the canvas")
      center[[nm]] <- rng[1,]
      scale[[nm]]  <- d
    } else if(align != "none") {
      z <- normPlane(z, w, who)
      center[[nm]] <- attr(z, "center")
      scale[[nm]]  <- attr(z, "scale")
    } else {
      center[[nm]] <- rep(0, p)
      scale[[nm]]  <- 1
    }
    rot[[nm]] <- diag(p)
  }
  common <- function(nm, z=L[[nm]]) sweep(sweep(z, 2, center[[nm]]), 2, scale[[nm]], "/")

  iter <- converged <- consensus <- NULL
  if(align == "procrustes") {
    Z <- lapply(setNames(use, use), function(nm) common(nm)[ok, , drop=FALSE])
    if(!gpa) {
      for(nm in setdiff(use, ref))
        rot[[nm]] <- attr(procrustesRotate(Z[[nm]], Z[[ref]], w), "rotation")
    } else {
      # every plane turned towards the mean of them all, and the mean worked out again, until
      # the scatter around it stops shrinking. Only the planes asked for make up the mean.
      M    <- Z[[ref]]
      prev <- Inf
      converged <- FALSE
      for(iter in seq_len(maxit)) {
        for(nm in sel) rot[[nm]] <- attr(procrustesRotate(Z[[nm]], M, w), "rotation")
        Y  <- lapply(sel, function(nm) Z[[nm]] %*% rot[[nm]])
        M  <- Reduce(`+`, Y)/length(Y)
        ss <- sum(vapply(Y, function(y) {
                d2 <- rowSums((y - M)^2)
                if(is.null(w)) sum(d2) else sum(d2*w)
              }, 0))
        if(prev - ss <= tol*max(1, ss)) { converged <- TRUE; break }
        prev <- ss
      }
      if(!converged)
        warning(who, ": the generalized Procrustes fit did not settle in ", maxit,
                " iterations; raise maxit", call.=FALSE)
      # the mean has no orientation of its own: the whole set is turned as one, which alters
      # none of the fits, so that it reads like the reference plane
      q <- attr(procrustesRotate(M, Z[[ref]], w), "rotation")
      for(nm in sel) rot[[nm]] <- rot[[nm]] %*% q
      consensus <- matrix(NA_real_, n, p, dimnames=list(rownames(L[[ref]]), colnames(L[[ref]])))
      consensus[ok, ] <- M %*% q
    }
  }

  # back to the units of the reference plane, which is thereby left just as it was
  toRef <- function(u) sweep(sweep(u, 2, scale[[ref]], "*"), 2, center[[ref]], "+")
  out <- lapply(setNames(sel, sel), function(nm) {
    a <- toRef(common(nm) %*% rot[[nm]])
    dimnames(a) <- dimnames(L[[nm]])
    a
  })
  if(!is.null(consensus)) consensus <- toRef(consensus)

  angle <- reflection <- NULL
  if(align == "procrustes" && p == 2L) {
    angle      <- vapply(rot[sel], function(r) atan2(r[2,1], r[1,1])*180/pi, 0)
    reflection <- vapply(rot[sel], function(r) det(r) < 0, NA)
  }
  info <- list(align=align, to=ref, gpa=gpa, rows=sum(ok), dropped=sum(!ok),
               center=center[sel], scale=scale[sel], rotation=rot[sel],
               angle=angle, reflection=reflection, weights=wl$w, wSource=wl$source)
  if(gpa) info <- c(info, list(consensus=consensus, iterations=iter, converged=converged))
  info <- info[!vapply(info, is.null, NA)]

  if(!is.null(wl$w))
    message(who, ": weighted by ", wl$source, " (N = ", format(sum(wl$w), big.mark=""), ")")

  if(!inherits(x, "netCoin")) {
    attributes(out) <- c(attributes(out), info)
    return(out)
  }

  # a netCoin object comes back with the aligned planes in place of the originals, the first
  # of them also as the fx/fy it is drawn with, and the cluster centroids worked out again.
  # With a suffix, the originals stay and each aligned plane is placed right after its own,
  # under its name and the suffix, so that the selector offers them side by side; the
  # reference, which only gpa moves, is not repeated.
  lay <- currentLayouts(x)
  if(is.null(suffix)) {
    lay[sel] <- out
    x$nodes$fx <- lay[[1]][,1]
    x$nodes$fy <- lay[[1]][,2]
  } else {
    add   <- if(gpa) out else out[setdiff(sel, ref)]
    # the titles of their axes take the suffix too, as a rotated axis is no longer the one
    # of the original plane, nor that of the plane it was brought to
    add   <- lapply(add, function(a) {
      if(!is.null(colnames(a))) colnames(a) <- ifelse(nzchar(colnames(a)),
                                                      paste0(colnames(a), suffix), "")
      a
    })
    asked <- paste0(names(add), suffix)
    given <- make.unique(c(names(lay), asked))[length(lay) + seq_along(asked)]
    if(any(given != asked))
      warning(who, ": the object already held ", nameClash(asked[given != asked],
              given[given != asked]), ". Choose another suffix to tell them apart.", call.=FALSE)
    names(given) <- names(add)
    merged <- list()
    for(nm in names(lay)) {
      merged[[nm]] <- lay[[nm]]
      if(nm %in% names(add)) merged[[given[[nm]]]] <- add[[nm]]
    }
    lay <- merged
  }
  x$layouts <- lay
  # every plane, the aligned ones included, with the titles of its axes, unless the object
  # states one pair for all of them
  if(length(lay) > 1 && is.null(x$options$axesLabels)) x <- planesAxesLabels(x)
  if(length(attr(x, "clusterColumns"))) x$clusters <- currentClusters(x)
  attr(x, "alignment") <- info
  x
}
