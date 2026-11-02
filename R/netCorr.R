# Program to apply nets to correlations
## netCorr ----
netCorr<-function(variables, weight=NULL, pairwise=FALSE,
                  sort=FALSE, decreasing=TRUE,
                  frequency=pairwise, means=TRUE,
                  method="pearson", criteria="p", Bonferroni=FALSE,
                  minL=0, maxL=Inf,
                  sortL=NULL, decreasingL=TRUE,
                  igraph=FALSE, ...)
{
  arguments <- list(...)
  arguments$name <- nameByLanguage(arguments$name,arguments$language,arguments$nodes)
  methods <- corrMethods(method)
  # Correlations have a sign: links are bipolar and coloured by their value unless told otherwise
  if(!("linkBipolar" %in% names(arguments))) arguments$linkBipolar <- TRUE
  if(!("lcolor" %in% names(arguments))) arguments$lcolor <- "value"
  bipolar <- isTRUE(arguments[["linkBipolar"]])
  if(!("lwidth" %in% names(arguments))) arguments$lwidth <- "value"
  if(!("lweight" %in% names(arguments))) arguments$lweight <- "value"

  # Weights: the name of a column of variables, or a numeric vector of one weight per case
  variables <- as.data.frame(variables)
  w <- NULL
  if(!is.null(weight)) {
    if(is.character(weight)) {
      if(length(weight)!=1 || !(weight %in% names(variables)))
        stop("weight must be the name of a column of variables or a numeric vector with a weight per case")
      w <- variables[[weight]]
      variables <- variables[, setdiff(names(variables), weight), drop=FALSE]
    } else w <- weight
    if(!is.numeric(w) || length(w)!=nrow(variables))
      stop("weight must be the name of a column of variables or a numeric vector with a weight per case")
    if(any(w<0, na.rm=TRUE)) stop("weights cannot be negative")
    if(any(methods %in% c("kendall","dcor")))
      stop("weights are not available for the methods kendall and dcor")
  }
  noNumeric <- !vapply(variables, is.numeric, logical(1))
  if(any(noNumeric)) stop("variables must be quantitative: ", paste(names(variables)[noNumeric], collapse=", "))

  # Missing values: listwise by default, pairwise on request
  keep <- if(pairwise) rep(TRUE, nrow(variables)) else stats::complete.cases(variables)
  if(!is.null(w)) keep <- keep & !is.na(w)
  variables <- variables[keep, , drop=FALSE]
  w <- w[keep]
  X <- as.matrix(variables)
  storage.mode(X) <- "double"
  m <- ncol(X)
  if(m<2) stop("at least two variables are needed")
  if(nrow(X)<3) stop("at least three complete cases are needed")
  if(is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(m))
  cn <- colnames(X)

  # Nodes: weighted when weights are given
  nodeStat <- function(j) {
    ok <- !is.na(X[,j])
    x <- X[ok,j]
    if(is.null(w)) {
      mu <- mean(x)
      s <- stats::sd(x)
    } else {
      ww <- w[ok]
      mu <- sum(ww*x)/sum(ww)
      s <- sqrt(sum(ww*(x-mu)^2)/(sum(ww)-sum(ww^2)/sum(ww)))
    }
    c(n=length(x), mean=mu, std=s, min=suppressWarnings(min(x)), max=suppressWarnings(max(x)))
  }
  S <- t(vapply(seq_len(m), nodeStat, numeric(5)))
  # sort: the variables go by frequency (valid cases) and, within equal frequency, by mean
  if(isTRUE(sort)) {
    o <- order(S[,"n"], S[,"mean"], decreasing=decreasing)
    X <- X[, o, drop=FALSE]
    S <- S[o, , drop=FALSE]
    cn <- colnames(X)
  }
  statistics <- data.frame(name=cn)
  if(frequency) statistics[[getByLanguage(frequencyList, arguments$language)]] <- S[,"n"]
  if(means) {
    statistics$mean <- round(S[,"mean"],2)
    statistics$std <- round(S[,"std"],2)
    statistics$min <- S[,"min"]
    statistics$max <- S[,"max"]
  }
  colnames(statistics)[1] <- arguments$name
  if(!is.null(arguments$nodes)) arguments$nodes <- merge(statistics, arguments$nodes, by=arguments$name, all.x=TRUE, sort=FALSE)
  else arguments$nodes <- statistics

  # Statistics of every method over the pairs of variables, in the order of the lower triangle
  pw <- pairwise && anyNA(X)
  stats.l <- lapply(methods, function(mt) corrStat(X, w, mt, pw, bipolar))
  ij <- which(lower.tri(diag(m)), arr.ind=TRUE)
  E <- data.frame(Source=cn[ij[,2]], Target=cn[ij[,1]], value=stats.l[[1]]$R[ij],
                  z=stats.l[[1]]$z, p=stats.l[[1]]$p, stringsAsFactors=FALSE)
  if(pw) E$n <- round(stats.l[[1]]$n, 2)
  for(k in seq_along(methods)[-1]) {
    E[[methods[k]]] <- stats.l[[k]]$R[ij]
    E[[paste0("z(", methods[k], ")")]] <- stats.l[[k]]$z
    E[[paste0("p(", methods[k], ")")]] <- stats.l[[k]]$p
  }
  R <- stats.l[[1]]$R

  # Selection criterion. The statistics of the first method are value, z and p
  criteria <- corrCriteria(criteria, methods)
  isP <- criteria=="p" || grepl("^p\\(", criteria)
  if(isP && maxL==Inf) maxL <- .5*(1+bipolar)
  if(isP && Bonferroni) maxL <- maxL/choose(m,2) # one contrast for each pair of variables
  sel <- E[[criteria]]
  if(bipolar && !isP) sel <- abs(sel) # the strength of a bipolar link is its absolute value
  E <- E[!is.na(sel) & sel>=minL & sel<=maxL, , drop=FALSE]
  if (!is.null(sortL)) E<-E[order((-1*decreasingL+!decreasingL)*E[[sortL]]),]
  rownames(E) <- NULL
  arguments$links <- E
  if(!is.null(arguments[["layout"]])) arguments$layout <- layoutSpecial(arguments$layout, list(pc=function() layoutPCA(R)))
  xNx <- do.call(netCoin,arguments)
  if (igraph) return(toIgraph(xNx))
  else return(xNx)
}

# Methods of association netCorr knows
corrMethods <- function(method) {
  known <- c("pearson", "spearman", "kendall", "partial", "dcor")
  method <- tolower(as.character(method))
  if(!length(method) || anyNA(method)) stop("method must name at least one of: ", paste(known, collapse=", "))
  out <- vapply(method, function(mt) {
    if(mt=="p") return("pearson") # "p" has always meant Pearson
    k <- pmatch(mt, known)
    if(is.na(k)) stop("method '", mt, "' is not available. Choose among: ", paste(known, collapse=", "))
    known[k]
  }, character(1))
  unique(unname(out))
}

# The column of the links a criterion refers to: value, z, p, the name of a method (its
# coefficient) or z(method) / p(method)
corrCriteria <- function(criteria, methods) {
  cr <- tolower(criteria[1])
  cols <- c("value", "z", "p")
  if(length(methods)>1) cols <- c(cols, methods[-1], paste0("z(", methods[-1], ")"), paste0("p(", methods[-1], ")"))
  if(cr==methods[1]) return("value")
  inner <- sub("^[zp]\\((.*)\\)$", "\\1", cr)
  if(inner!=cr && inner==methods[1]) return(substr(cr, 1, 1))
  if(!(cr %in% cols))
    stop("criteria '", criteria[1], "' is not available. Choose among: ", paste(unique(c(methods[1], cols)), collapse=", "))
  cr
}

## Statistics of association ----
# Each of them gives the matrix of coefficients R and, for the pairs of the lower triangle
# (column by column), the statistic z of the contrast, its p value and the number of cases n
# on which it rests (the effective number of Kish when there are weights). The contrast is
# one-sided (positive association) unless the links are bipolar.
corrStat <- function(X, w, method, pairwise, bipolar) {
  m <- ncol(X)
  M <- !is.na(X)
  lt <- lower.tri(diag(m))
  ij <- which(lt, arr.ind=TRUE)
  use <- if(pairwise) "pairwise.complete.obs" else "everything"
  if(is.null(w)) Ne <- crossprod(M)
  else Ne <- crossprod(M*w, M)^2 / crossprod(M*w^2, M)
  df <- NULL
  R <- switch(method,
    pearson = if(is.null(w)) stats::cor(X, use=use) else wpearson(X, w),
    spearman = if(is.null(w)) stats::cor(X, method="spearman", use=use)
               else if(!pairwise) wpearson(apply(X, 2, wrank, w=w), w)
               else pairMatrix(X, w, function(x, y, ww) wcor2(wrank(x, ww), wrank(y, ww), ww)),
    kendall = stats::cor(X, method="kendall", use=use),
    partial = {
      cc <- stats::complete.cases(X)
      if(!all(cc)) warning("partial correlations rest on the cases complete in all the variables", call.=FALSE)
      Xc <- X[cc, , drop=FALSE]
      wc <- w[cc]
      Rc <- if(is.null(wc)) stats::cor(Xc) else wpearson(Xc, wc)
      nn <- if(is.null(wc)) nrow(Xc) else sum(wc)^2/sum(wc^2)
      if(nn-m<=0) stop("partial correlations need more complete cases than variables")
      Th <- tryCatch(solve(Rc), error=function(e)
        stop("partial correlations cannot be computed: the correlation matrix is singular", call.=FALSE))
      P <- -Th/sqrt(outer(diag(Th), diag(Th)))
      diag(P) <- 1
      Ne <- matrix(nn, m, m)
      df <- matrix(nn-m, m, m)
      P
    },
    dcor = if(!pairwise) dcorMatrix(X) else pairMatrix(X, NULL, dcor2),
    stop("method '", method, "' is not available"))
  R <- pmin(pmax(R, -1), 1)
  dimnames(R) <- list(colnames(X), colnames(X))
  r <- R[lt]
  n <- Ne[lt]
  if(is.null(df)) df <- n-2 else df <- df[lt]
  pT <- function(t, df) if(bipolar) 2*stats::pt(-abs(t), df) else stats::pt(t, df, lower.tail=FALSE)
  pZ <- function(z) if(bipolar) 2*stats::pnorm(-abs(z)) else stats::pnorm(z, lower.tail=FALSE)

  if(method %in% c("pearson", "spearman", "partial")) {
    z <- r*sqrt(df/(1-r^2))
    p <- pT(z, df)
    if(method=="spearman" && is.null(w)) {
      # exact contrast for the small samples without ties
      small <- which(n>=3 & n<=9)
      if(length(small)) {
        D <- spearmanD(X, M, ij[small, , drop=FALSE], n[small], pairwise)
        ok <- !is.na(D)
        if(any(ok)) {
          k <- small[ok]
          for(nk in unique(n[k])) {
            kk <- k[n[k]==nk]
            p[kk] <- spearmanExactP(D[ok][n[k]==nk], nk, bipolar)
          }
        }
      }
    }
  }
  else if(method=="kendall") {
    ties <- kendallTies(X, M, ij, pairwise)
    n0 <- n*(n-1)/2
    Sc <- r*sqrt((n0-ties$x[1,]/2)*(n0-ties$y[1,]/2))
    V <- (n*(n-1)*(2*n+5) - ties$x[2,] - ties$y[2,])/18 +
         ties$x[1,]*ties$y[1,]/(2*n*(n-1)) +
         ties$x[3,]*ties$y[3,]/(9*n*(n-1)*(n-2))
    z <- Sc/sqrt(V)
    z[n<3] <- NA
    p <- pZ(z)
    # exact contrast for the small samples without ties
    exact <- which(n>=3 & n<50 & ties$x[1,]==0 & ties$y[1,]==0 & !is.na(Sc))
    for(nk in unique(n[exact])) {
      kk <- exact[n[exact]==nk]
      cdf <- kendallCdf(nk)
      S <- round(if(bipolar) abs(Sc[kk]) else Sc[kk])
      pr <- cdf[pmin(pmax((nk*(nk-1)/2-S)/2, 0), length(cdf)-1) + 1] # P(S0 >= S)
      p[kk] <- if(bipolar) pmin(1, 2*pr) else pr
    }
  }
  else if(method=="dcor") {
    v <- n*(n-3)/2
    z <- sqrt(v-1)*r/sqrt(1-r^2)
    z[n<4] <- NA
    p <- pT(z, v-1)
  }
  list(R=R, z=z, p=p, n=n)
}

# Weighted Pearson correlations, each pair resting on the cases observed in both variables
wpearson <- function(X, w) {
  n <- nrow(X)
  M <- !is.na(X)
  X[!M] <- 0
  mu <- colSums(X*w)/colSums(M*w)
  X <- (X - rep(mu, each=n)) * M # centered, and zero where missing
  W <- crossprod(M*w, M)
  A <- crossprod(X*w, M) # A[i,j]: sum of w*x_i over the cases observed in j
  B <- crossprod(X^2*w, M)
  C <- crossprod(X*w, X)
  V <- B/W - (A/W)^2 # variance of i over the cases shared with j
  R <- (C/W - A*t(A)/W^2)/sqrt(V*t(V))
  diag(R) <- 1
  dimnames(R) <- list(colnames(X), colnames(X))
  R
}

wcor2 <- function(x, y, w) {
  W <- sum(w)
  x <- x - sum(w*x)/W
  y <- y - sum(w*y)/W
  sum(w*x*y)/sqrt(sum(w*x^2)*sum(w*y^2))
}

# Ranks of the values weighted by w; ties share the middle of their block. With equal weights
# they are the usual mean ranks, up to a constant
wrank <- function(x, w) {
  o <- order(x)
  xs <- x[o]
  g <- cumsum(c(TRUE, xs[-1]!=xs[-length(xs)]))
  gs <- as.vector(rowsum(w[o], g))
  out <- numeric(length(x))
  out[o] <- (cumsum(gs) - gs/2)[g]
  out
}

# Matrix of a function of two variables over the cases observed in both
pairMatrix <- function(X, w, f) {
  m <- ncol(X)
  M <- !is.na(X)
  R <- diag(m)
  for(j in seq_len(m-1)) for(i in (j+1):m) {
    idx <- M[,i] & M[,j]
    v <- if(sum(idx)<2) NA else f(X[idx,i], X[idx,j], if(is.null(w)) NULL else w[idx])
    R[i,j] <- R[j,i] <- v
  }
  dimnames(R) <- list(colnames(X), colnames(X))
  R
}

## Distance correlation ----
# Distances centered "U" (unbiased) way of Szekely & Rizzo (2014), as a vector
ucenter <- function(x) {
  n <- length(x)
  a <- abs(outer(x, x, "-"))
  rs <- rowSums(a)
  a <- a - outer(rs, rs, "+")/(n-2) + sum(a)/((n-1)*(n-2))
  diag(a) <- 0
  as.vector(a)
}

# Bias-corrected distance correlation between all the columns (no missing values)
dcorMatrix <- function(X) {
  n <- nrow(X)
  m <- ncol(X)
  if(n<4) stop("distance correlation needs at least four cases")
  if(as.numeric(n)^2*m>6e7) stop("too many cases for distance correlation (it keeps n^2 distances per variable)")
  U <- vapply(seq_len(m), function(j) ucenter(X[,j]), numeric(n*n))
  G <- crossprod(U)/(n*(n-3))
  d <- diag(G)
  R <- G/sqrt(outer(d, d))
  diag(R) <- 1
  R
}

dcor2 <- function(x, y, w=NULL) {
  n <- length(x)
  if(n<4) return(NA_real_)
  ux <- ucenter(x)
  uy <- ucenter(y)
  sum(ux*uy)/sqrt(sum(ux*ux)*sum(uy*uy))
}

## Exact contrasts ----
corrCache <- new.env()

# Distribution of the number of inversions of a permutation of n (Mahonian numbers); the
# probability that Kendall's S of n cases without ties reaches a value is P(inversions <= (n0 - S)/2)
kendallCdf <- function(n) {
  key <- paste0("k", n)
  if(is.null(corrCache[[key]])) {
    f <- 1
    for(k in seq_len(n)[-1]) {
      g <- numeric(k*(k-1)/2 + 1)
      for(j in 0:(k-1)) {
        idx <- (j+1):(j+length(f))
        g[idx] <- g[idx] + f
      }
      f <- g
    }
    corrCache[[key]] <- cumsum(f)/sum(f)
  }
  corrCache[[key]]
}

# Distribution of the sum of squared differences of ranks D of n cases without ties
spearmanNull <- function(n) {
  key <- paste0("s", n)
  if(is.null(corrCache[[key]])) {
    perms <- function(n) {
      if(n==1) return(matrix(1L, 1, 1))
      p <- perms(n-1)
      do.call(rbind, lapply(seq_len(n), function(i) cbind(i, matrix(ifelse(p>=i, p+1L, p), nrow(p)))))
    }
    P <- perms(n)
    D <- rowSums((P - rep(seq_len(n), each=nrow(P)))^2)
    tab <- table(D)
    corrCache[[key]] <- list(vals=as.numeric(names(tab)), cum=cumsum(as.vector(tab))/sum(tab))
  }
  corrCache[[key]]
}

spearmanExactP <- function(D, n, bipolar) {
  nl <- spearmanNull(n)
  cum <- c(0, nl$cum)
  le <- cum[findInterval(D, nl$vals) + 1] # P(D0 <= D): as strong a positive association
  if(!bipolar) return(le)
  lt <- cum[findInterval(D, nl$vals, left.open=TRUE) + 1]
  pmin(1, 2*pmin(le, 1-lt))
}

# D of the pairs (rows of ij) for the samples n of at most nine cases; NA when there are ties
spearmanD <- function(X, M, ij, n, pairwise) {
  if(!pairwise) { # all the cases are complete: vectorized
    tied <- apply(X, 2, anyDuplicated)>0
    Rk <- apply(X, 2, rank)
    G <- crossprod(Rk)
    D <- diag(G)[ij[,1]] + diag(G)[ij[,2]] - 2*G[ij]
    D[tied[ij[,1]] | tied[ij[,2]]] <- NA
    return(round(D))
  }
  vapply(seq_len(nrow(ij)), function(k) {
    idx <- M[,ij[k,1]] & M[,ij[k,2]]
    x <- X[idx, ij[k,1]]
    y <- X[idx, ij[k,2]]
    if(anyDuplicated(x) || anyDuplicated(y)) NA_real_ else sum((rank(x)-rank(y))^2)
  }, numeric(1))
}

# Tie sums of each variable of the pairs: t(t-1), t(t-1)(2t+5) and t(t-1)(t-2) over the
# groups of tied values, as the rows of a matrix with a column per pair
kendallTies <- function(X, M, ij, pairwise) {
  ts <- function(x) {
    t <- tabulate(match(x, unique(x)))
    c(sum(t*(t-1)), sum(t*(t-1)*(2*t+5)), sum(t*(t-1)*(t-2)))
  }
  if(!pairwise) {
    C <- apply(X, 2, ts)
    return(list(x=C[, ij[,2], drop=FALSE], y=C[, ij[,1], drop=FALSE]))
  }
  both <- vapply(seq_len(nrow(ij)), function(k) {
    idx <- M[,ij[k,1]] & M[,ij[k,2]]
    c(ts(X[idx, ij[k,2]]), ts(X[idx, ij[k,1]]))
  }, numeric(6))
  list(x=both[1:3, , drop=FALSE], y=both[4:6, , drop=FALSE])
}

# Program to apply evolving nets to correlations.
## d_netCorr ----
d_netCorr <- function(variables, nodes= NULL, weight=NULL, 
                      pairwise=FALSE,
                      frequency=pairwise, means=TRUE,
                      method="pearson", criteria="value", Bonferroni=FALSE,
                      minL=0, maxL=Inf,
                      sortL=NULL, decreasingL=TRUE, 
                      factorial=c("null", "pc", "nf", "vf", "of"), 
                      components=TRUE, backcomponents=FALSE, sequence=seq(.20, 1, .01), 
                      textFilter=c(1, .99), speed=50, dir=NULL, ...)
{
  arguments <- list(variables= variables, nodes=nodes, weight=weight,
                    pairwise=pairwise,
                    frequency=frequency, means=means,
                    method=method, criteria=criteria, Bonferroni=Bonferroni,
                    minL=minL, maxL=maxL, sortL=sortL, decreasingL=decreasingL, ...)
  if(exists("layout", arguments) & !identical(factorial, c("null", "pc","nf","rf","of"))) warning("Argument factorial is incompatible with layout")
  if(exists("layout", arguments)) C <- arguments$layout else C <- layoutFact(variables, factorial)
  G <- do.call (netCorr, arguments)
  g <- list()
  g$lineplots <- c("component", "Degree", "closeness", "betweenness", "eigen", "ratio" )
  g$mode <- "frame"
  g$speed <- speed
  if(!is.null(dir)){
    g$dir <- dir
  }
  ch=0
  if(length(textFilter)<2) textFilter[2]<- .99
  for(I in  sequence) {
    G$links$text <- ifelse(abs(G$links$value)>=min(textFilter[1], I+textFilter[2]), sprintf("%.2f", G$links$value), "")
    # bipolar links are filtered by their strength: the absolute value
    addArgs <- list(G, linkFilter = paste0(if(isTRUE(G$options$linkBipolar)) "abs(value)" else "value", ">", I), lwidth="value",
                    ltext="text", layout=C, main=paste0("Correlation: ", I))
    if(means) addArgs$size <- "mean" # the means are the sizes of the nodes
    H <- do.call(addNetCoin, addArgs)
    if(exists("hidden", H$links)) H$links <- H$links[!H$links$hidden, ]
    else H$links$hidden <- FALSE
    enodes <- unique(c(H$links$Source,H$links$Target))
    H$nodes <- H$nodes[H$nodes[[H$options$nodeName]] %in% enodes, , drop=FALSE]
    comps <- igraph::components(toIgraph(H))
    cc <- sum(comps$csize>1)
    if(((cc > ch | !components) | (cc < ch &  components & backcomponents)) & nrow(H$links)>0) {
      central <- calCentr(H)
      H$nodes <- cbind(H$nodes, central$nodes)
      H$nodes <- cbind(H$nodes, compon(comps$membership))
      names(H$nodes)[ncol(H$nodes)] <- "component"
      names(H$nodes)[names(H$nodes)=="degree"] <- "Degree"
      H$links$ratio <- abs(H$links$value)/I
      g[[paste0("C",I)]]  <- H
    }
    ch <- cc
  }
  multi <- do.call(multigraphCreate, g)
  return(multi)
}

corr <- function (a, b = a, weight = NULL )
{
  if (is.null(weight)) weight= rep(1/nrow(a), nrow(a))
  s<-complete.cases(cbind(a,b,weight))
  a<-as.matrix(a[s,]);b<-as.matrix(b[s,])
  # normalize weights
  weight <- weight[s] / sum(weight[s])
  
  # center matrices
  a <- sweep(a, 2, apply((a * weight),2,sum))
  b <- sweep(b, 2, apply((b * weight),2,sum))
  
  # compute weighted correlatio
  t(a*weight) %*% b / sqrt(apply((a**2 *weight),2,sum) %*% t(apply((b**2 *weight),2,sum))) 
  
}

corrp <- function (a, b = a, weight = NULL ){
  count <- 0
  CC <- matrix( NA , nrow=(ncol(a)), ncol=ncol(b))
  rownames(CC) <- colnames(a)
  colnames(CC) <- colnames(b)
  for(i in colnames(a)) {
    for(j in colnames(b)) {
      if(i==j) CC[i,j] <- 1
      else {
        CC[i,j] <- corr(b[j],a[i],weight=weight)
      }
    }
    count <- count +1
  }
  return(CC)
}

# For d_netCorr: Labelling the graph components 
compon  <- function(comps) {
  s <- table(comps)
  o <- sort(table(comps), decreasing=T)
  xdig <- paste0("%",max(c(nchar(trunc(abs(as.numeric(names(s))))))),"d")
  j <- 1
  for(n in as.numeric(names(o))) {
    d <- as.numeric(names(s[j]))
    comps <- ifelse(comps==n, paste0("C-",sprintf(xdig, d)), comps)
    j <- j+1
  }
  return(comps)
}

