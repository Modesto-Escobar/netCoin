## optimalAssign ----
# The optimal one-to-one assignment (Hungarian algorithm): for every row of x, the column it
# is given, maximizing the total. clue::solve_LSAP does it when the package is at hand;
# otherwise an exact search with pruning, which is enough for the handful of categories a
# clusterization holds.
optimalAssign <- function(x) {
  x <- as.matrix(x)
  if(nrow(x) != ncol(x)) stop("optimalAssign: the matrix must be square")
  if(requireNamespace("clue", quietly=TRUE))
    return(as.integer(clue::solve_LSAP(x, maximum=TRUE)))

  k <- nrow(x)
  best <- list(value=-Inf, perm=seq_len(k))
  # suffix[i]: upper bound of what rows i..k can still contribute
  suffix <- c(rev(cumsum(rev(apply(x, 1, max)))), 0)
  recurse <- function(i, used, value, perm) {
    if(i > k) {
      if(value > best$value) best <<- list(value=value, perm=perm)
      return(invisible(NULL))
    }
    if(value + suffix[i] <= best$value) return(invisible(NULL))
    for(j in which(!used)) {
      used[j] <- TRUE
      perm[i] <- j
      recurse(i+1L, used, value+x[i,j], perm)
      used[j] <- FALSE
    }
    invisible(NULL)
  }
  recurse(1L, logical(k), 0, integer(k))
  best$perm
}


## pairTable ----
# The contingency table of one pair of classifications.
# weights: what every unit counts for (the cases it stands for, a sampling weight...); NULL
# amounts to weights of one, that is, to counting units.
pairTable <- function(pred, ref, weights=NULL, who="pairTable") {
  ok <- !is.na(pred) & !is.na(ref)
  if(!is.null(weights)) {
    if(length(weights) != length(pred))
      stop(who, ": weights holds ", length(weights), " values and the classifications ",
           length(pred))
    if(any(weights < 0, na.rm=TRUE)) stop(who, ": there are negative weights")
    ok <- ok & !is.na(weights)
  }
  pred <- droplevels(as.factor(pred[ok]))
  ref  <- droplevels(as.factor(ref[ok]))
  if(is.null(weights)) {
    m <- table(pred, ref, dnn=c("Predicted", "Model"))
  } else {
    # a weighted contingency table: every cell adds up the weights rather than counting
    m <- tapply(weights[ok], list(Predicted=pred, Model=ref), sum)
    m[is.na(m)] <- 0
    m <- as.table(m)
  }
  attr(m, "rows")    <- sum(ok)          # units actually compared
  attr(m, "dropped") <- sum(!ok)
  m
}


## classMetrics ----
# The one-against-the-rest measures of a confusion matrix, following the definitions of
# caret::confusionMatrix (prediction down the rows, reference along the columns) but computed
# here so as not to depend on caret.
# m: a square table, prediction in rows and reference in columns.
classMetrics <- function(m) {
  n  <- sum(m)
  tp <- diag(m)
  fp <- rowSums(m) - tp
  fn <- colSums(m) - tp
  tn <- n - tp - fp - fn
  sens <- tp/(tp+fn)                     # sensitivity = recall
  spec <- tn/(tn+fp)                     # specificity
  prec <- tp/(tp+fp)                     # precision
  # a form equivalent to 2*P*S/(P+S) which returns 0 -rather than NaN- when a category gets
  # no hit at all; it is NaN only where the category does not exist
  f1   <- 2*tp/(2*tp+fp+fn)
  out <- cbind(Sensitivity=sens, Specificity=spec, Precision=prec, F1=f1,
               Prevalence=colSums(m)/n, N=colSums(m))
  rownames(out) <- colnames(m)
  out
}


## confusionAlign ----
# Compare one pair of classifications: pred is the one under scrutiny and ref the model.
# Categories carry no name of their own, so before comparing they are realigned with those of
# the model through the optimal assignment, which maximizes the diagonal of the table.
# Returns the realigned table, with the measures and the mapping as attributes.
confusionAlign <- function(pred, ref, weights=NULL, align=TRUE) {
  m <- pairTable(pred, ref, weights, who="confusionAlign")
  if(nrow(m) != ncol(m))
    stop("confusionAlign: different number of categories (", nrow(m), " against ", ncol(m), ")")

  rows <- attr(m, "rows"); dropped <- attr(m, "dropped")
  if(align) {
    asig <- optimalAssign(m)             # row i -> column asig[i]
    m    <- m[order(asig), , drop=FALSE]
    rownames(m) <- colnames(m)
  }
  mapa <- data.frame(Model=colnames(m), Predicted=rownames(m), stringsAsFactors=FALSE)

  attr(m, "byClass") <- classMetrics(m)
  attr(m, "mapping") <- mapa
  attr(m, "n")       <- sum(m)           # sum of the weights (units where there are none)
  attr(m, "rows")    <- rows
  attr(m, "dropped") <- dropped
  m
}


## agreeIndices ----
# Indices which need no matching of categories, and hold therefore even where the two
# classifications have a different number of them.
#   ARI: adjusted Rand index. It walks through the pairs of units, counting whether the two
#        classifications put them together or keep them apart, and corrects for chance: 0 is
#        what independent partitions yield, 1 an exact match. It may go slightly negative.
#   NMI: mutual information normalized by the mean of the two entropies. From 0 (independent)
#        to 1 (exact match). Uncorrected for chance, so it tends to rise with the number of
#        categories.
#   VI : variation of information, H(U)+H(V)-2*MI, a distance in nats: 0 is an exact match and
#        it grows with the discrepancy. VI.n divides it by log(N) to keep it between 0 and 1.
# With weights the count of pairs is extended continuously, which takes for granted that the
# weights are on a scale of cases (which is what the default weighting yields); weights given
# as proportions would make ARI meaningless.
agreeIndices <- function(m) {
  n  <- sum(m)
  a  <- rowSums(m); b <- colSums(m)
  ch2 <- function(z) z*(z-1)/2           # combinations two at a time
  idxObs <- sum(ch2(m))
  idxExp <- sum(ch2(a))*sum(ch2(b))/ch2(n)
  idxMax <- (sum(ch2(a))+sum(ch2(b)))/2
  ari <- if(isTRUE(all.equal(idxMax, idxExp))) NA_real_ else
         (idxObs-idxExp)/(idxMax-idxExp)

  ent <- function(p) { p <- p[p > 0]; -sum(p*log(p)) }
  pu <- a/n; pv <- b/n; pm <- m/n
  hu <- ent(pu); hv <- ent(pv)
  outer_p <- outer(pu, pv)
  keep <- pm > 0
  mi  <- sum(pm[keep]*log(pm[keep]/outer_p[keep]))
  nmi <- if((hu+hv) <= 0) NA_real_ else 2*mi/(hu+hv)
  vi  <- hu+hv-2*mi
  c(ARI=ari, NMI=nmi, VI=vi, VI.n=if(n > 1) vi/log(n) else NA_real_)
}


## matchSetup ----
# What clusterMatch and clusterAgree both settle before comparing anything: the model, the
# unit of analysis (nodes or cases), the weighting and the list of candidate classifications.
# Gathering it here keeps the two functions from drifting apart.
matchSetup <- function(x, model, modelName, columns, weights, modal, sameK, who) {
  nodes <- if(is.data.frame(x)) x else x$nodes
  if(!is.data.frame(nodes))
    stop(who, ": x must be a netCoin object or a data.frame")
  nNodes <- nrow(nodes)
  idx    <- attr(x, "caseToPattern")     # the pattern (row of nodes) every case fell into
  nCases <- if(is.null(idx)) NULL else length(idx)

  # Where the nodes are patterns, surScat keeps in every cluster column all the groups the
  # pattern spans ("Group: 2|Group: 1"), ordered from the commonest down, and does not keep
  # which group each case went to. The only comparable thing is then the modal group, which
  # is the first of the list.
  spanning <- integer(0)
  modalize <- function(z, nm) {
    if(is.null(idx) || !is.character(z)) return(z)
    n <- sum(grepl("|", z, fixed=TRUE))
    if(n) {
      spanning[nm] <<- n
      if(modal) z <- sub("\\|.*$", "", z)
    }
    z
  }

  # the model, and with it the unit of analysis
  if(is.character(model) && length(model) == 1L && model %in% names(nodes)) {
    modelName <- model
    modelVec  <- modalize(nodes[[model]], model)
    unit      <- "nodes"
  } else {
    # a single name that is no column of the node table would fall through as a vector of one
    # element and be reported as a model of the wrong length, which says nothing of the
    # misspelling behind it
    if(is.character(model) && length(model) == 1L) {
      cand <- attr(x, "clusterColumns")
      if(is.null(cand))
        cand <- names(nodes)[vapply(nodes, function(z) is.factor(z) || is.character(z), TRUE)]
      stop(who, ": there is no node column called \"", model, "\"",
           if(length(cand)) paste0(". Classifications available: ",
                                   paste(utils::head(cand, 10), collapse=", "),
                                   if(length(cand) > 10)
                                     paste0(" and ", length(cand)-10, " more") else ""))
    }
    modelVec <- model
    if(length(modelVec) == nNodes) unit <- "nodes"
    else if(!is.null(nCases) && length(modelVec) == nCases) unit <- "cases"
    else stop(who, ": the model holds ", length(modelVec), " elements; one per node (",
              nNodes, ")", if(!is.null(nCases)) paste0(" or one per case (", nCases, ")"),
              " was expected.", sampledNote(x))
  }
  modelVec <- droplevels(as.factor(modelVec))
  k <- nlevels(modelVec)
  if(k < 2L) stop(who, ": the model must hold at least 2 categories")

  # the weighting
  caseW <- attr(x, "caseWeight")
  if(!is.null(caseW) && !is.null(nCases) && length(caseW) != nCases) caseW <- NULL
  wSource <- "unit"
  if(is.null(weights)) {
    if(unit == "cases") {
      if(!is.null(caseW)) { w <- caseW; wSource <- "caseWeight" }
      else w <- NULL
    } else if(!is.null(idx)) {
      # pattern nodes: each one weighs the cases it stands for
      cw <- if(is.null(caseW)) rep(1, nCases) else caseW
      w  <- as.vector(tapply(cw, factor(idx, levels=seq_len(nNodes)), sum))
      w[is.na(w)] <- 0
      wSource <- if(is.null(caseW)) "cases per pattern" else "cases per pattern (caseWeight)"
    } else w <- NULL
  } else {
    if(is.character(weights) && length(weights) == 1L) {
      if(!(weights %in% names(nodes)))
        stop(who, ": there is no node column called \"", weights, "\"")
      w <- nodes[[weights]]
      wSource <- paste0("column \"", weights, "\"")
    } else {
      w <- weights
      wSource <- "weights"
    }
    if(!is.numeric(w)) stop(who, ": the weights must be numeric")
    n <- if(unit == "cases") nCases else nNodes
    if(length(w) == 1L) w <- rep(w, n)
    if(length(w) != n) {
      if(unit == "nodes" && !is.null(nCases) && length(w) == nCases) {
        # weights given case by case but a comparison by nodes: they are added up
        w <- as.vector(tapply(w, factor(idx, levels=seq_len(nNodes)), sum))
        w[is.na(w)] <- 0
        wSource <- paste0(wSource, ", added up by pattern")
      } else
        stop(who, ": weights holds ", length(w), " values and the comparison goes by ", unit,
             " (", n, ")")
    }
  }
  if(!is.null(w) && any(w < 0, na.rm=TRUE)) stop(who, ": there are negative weights")

  # the candidates
  if(is.null(columns)) {
    columns <- attr(x, "clusterColumns")
    if(is.null(columns))
      columns <- names(nodes)[vapply(nodes, function(z) is.factor(z) || is.character(z), TRUE)]
  }
  columns <- intersect(columns, names(nodes))
  columns <- setdiff(columns, modelName)
  # the classifications live in the nodes: where the comparison goes case by case, every node
  # is spread over the cases it stands for
  getCol <- function(cl) {
    z <- modalize(nodes[[cl]], cl)
    if(unit == "cases") z[idx] else z
  }
  # this sweep has already seen every candidate column, so spanning is complete by the time
  # it is handed back
  nCat <- vapply(columns, function(cl) {
    z <- getCol(cl); nlevels(droplevels(as.factor(z[!is.na(z)]))) }, 0L)
  if(sameK) {
    columns <- columns[nCat == k]
    if(!length(columns))
      stop(who, ": no classification holds ", k, " categories.",
           if(length(spanning) && !modal)
             paste0(" There are ", length(spanning), " columns whose patterns span several",
                    " groups; with modal=TRUE the modal group of each pattern is compared."))
  } else if(!length(columns))
    stop(who, ": there is no classification to compare with")

  list(nodes=nodes, modelVec=modelVec, modelName=modelName, k=k, unit=unit, w=w,
       wSource=wSource, columns=columns, nCat=nCat[columns], getCol=getCol, spanning=spanning)
}


## sortResult ----
# sort: NULL keeps the order the classifications came in; the name of a column sorts by it
# from the largest down; a vector of names sorts by all of them at once, the first as the main
# criterion and the rest to break its ties. The prefix "-" turns one column around, from the
# smallest up: sort="-VI" puts the most concordant on top, and sort=c("-k","F1") the simplest
# classifications first.
sortResult <- function(res, sort, who) {
  if(is.null(sort)) return(res)
  desc <- !startsWith(sort, "-")
  nm   <- sub("^-", "", sort)
  bad  <- setdiff(nm, names(res))
  if(length(bad))
    stop(who, ": there is no column ", paste0("\"", bad, "\"", collapse=", "),
         " to sort by. Available: ", paste(names(res), collapse=", "),
         ". The prefix \"-\" sorts from the smallest up.")
  # xtfrm carries any column over to numbers ordering just as it does, text included, so
  # negating them turns that column around and only that one, without depending on order()
  # taking one decreasing per criterion
  keys <- Map(function(col, d) { z <- xtfrm(col); if(d) z else -z }, res[nm], desc)
  res <- res[do.call(order, c(unname(keys), list(decreasing=TRUE))), , drop=FALSE]
  rownames(res) <- NULL
  res
}


## reportSetup ----
# What both functions tell the caller once they are through.
reportSetup <- function(s, modal, who) {
  # the weighting changes what the figures stand for, so it is stated
  if(!is.null(s$w))
    message(who, ": comparison by ", s$unit, ", weighted by ", s$wSource, " (N = ",
            format(sum(s$w), big.mark=""), ")")
  # so are the patterns spanning several groups, because the object no longer holds which
  # group each case went to and the comparison is made with the modal one
  if(length(s$spanning))
    warning(who, ": ", paste0(names(s$spanning), " (", s$spanning, " patterns)", collapse=", "),
            if(length(s$spanning) > 1) " span" else " spans", " more than one group. ",
            if(modal)
              paste0("The modal group of each pattern has been compared, which is all the",
                     " object keeps: the cases of the pattern going to another group count",
                     " as hits.")
            else "With modal=FALSE they are left out of the comparison.", call.=FALSE)
}


## clusterMatch ----
# Compare a model classification with the other classifications of a netCoin object from
# surScat, category by category, once they have been matched one to one.
clusterMatch <- function(x, model, columns=NULL, weights=NULL, align=TRUE, weighted=TRUE,
                         modal=TRUE, sort=NULL, digits=3) {

  s <- matchSetup(x, model, deparse(substitute(model)), columns, weights, modal, sameK=TRUE,
                  who="clusterMatch")
  w <- s$w

  tables <- byClass <- mapping <- vector("list", length(s$columns))
  names(tables) <- names(byClass) <- names(mapping) <- s$columns
  res <- data.frame()

  for(cl in s$columns) {
    m    <- confusionAlign(s$getCol(cl), s$modelVec, weights=w, align=align)
    bc   <- attr(m, "byClass")
    prev <- bc[,"Prevalence"]            # relative size of every category
    n    <- attr(m, "n")
    po <- sum(diag(m))/n
    pe <- sum(rowSums(m)*colSums(m))/n^2

    vals <- list(
      Classification = cl,
      k              = s$k,
      N              = n,
      Accuracy       = po,
      Kappa          = (po-pe)/(1-pe),
      # sensitivity and precision, the two F1 is made of, and behind them specificity, which
      # plays no part in it
      Sensitivity    = mean(bc[,"Sensitivity"], na.rm=TRUE),
      Precision      = mean(bc[,"Precision"],   na.rm=TRUE),
      F1             = mean(bc[,"F1"],          na.rm=TRUE),
      Specificity    = mean(bc[,"Specificity"], na.rm=TRUE)
    )
    # every weighted mean next to its plain one. There is no Sensitivity.w because weighting
    # recall by prevalence gives the proportion of hits: it would be Accuracy over again.
    if(weighted) {
      after <- function(l, nm, new) append(l, new, after=match(nm, names(l)))
      vals <- after(vals, "Specificity",
                    list(Specificity.w=sum(prev*bc[,"Specificity"], na.rm=TRUE)))
      vals <- after(vals, "F1", list(F1.w=sum(prev*bc[,"F1"], na.rm=TRUE)))
    }
    if(!is.null(w)) vals <- append(vals, list(Rows=attr(m, "rows")), after=2)
    res <- rbind(res, as.data.frame(vals, stringsAsFactors=FALSE))

    tables[[cl]]  <- m
    byClass[[cl]] <- bc
    mapping[[cl]] <- attr(m, "mapping")
  }

  res <- sortResult(res, sort, "clusterMatch")
  num <- vapply(res, is.numeric, TRUE)
  num["k"] <- FALSE
  if(!is.null(w)) num["Rows"] <- FALSE
  res[num] <- lapply(res[num], round, digits)

  attr(res, "model")    <- s$modelName
  attr(res, "unit")     <- s$unit
  attr(res, "weights")  <- w
  attr(res, "spanning") <- if(length(s$spanning)) s$spanning else NULL
  attr(res, "tables")   <- tables
  attr(res, "byClass")  <- byClass
  attr(res, "mapping")  <- mapping
  reportSetup(s, modal, "clusterMatch")
  res
}


## clusterAgree ----
# The same comparison read as agreement between partitions, which matches no categories and
# so takes classifications with any number of them. The k column is that of each
# classification, not that of the model, which is left in attr(res, "kModel").
clusterAgree <- function(x, model, columns=NULL, weights=NULL, modal=TRUE, sort=NULL,
                         digits=3) {

  s <- matchSetup(x, model, deparse(substitute(model)), columns, weights, modal, sameK=FALSE,
                  who="clusterAgree")
  w <- s$w

  tables <- vector("list", length(s$columns))
  names(tables) <- s$columns
  res <- data.frame()

  for(cl in s$columns) {
    m <- pairTable(s$getCol(cl), s$modelVec, w, who="clusterAgree")
    vals <- c(list(Classification=cl, k=s$nCat[[cl]], N=sum(m)), as.list(agreeIndices(m)))
    if(!is.null(w)) vals <- append(vals, list(Rows=attr(m, "rows")), after=2)
    res <- rbind(res, as.data.frame(vals, stringsAsFactors=FALSE))
    tables[[cl]] <- m
  }

  res <- sortResult(res, sort, "clusterAgree")
  num <- vapply(res, is.numeric, TRUE)
  num["k"] <- FALSE
  if(!is.null(w)) num["Rows"] <- FALSE
  res[num] <- lapply(res[num], round, digits)

  attr(res, "model")    <- s$modelName
  attr(res, "kModel")   <- s$k
  attr(res, "unit")     <- s$unit
  attr(res, "weights")  <- w
  attr(res, "spanning") <- if(length(s$spanning)) s$spanning else NULL
  attr(res, "tables")   <- tables
  reportSetup(s, modal, "clusterAgree")
  res
}
