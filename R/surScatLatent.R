## surScat and latent models ----
# Latent profile analysis (LPA) groups cases on quantitative variables, latent class analysis
# (LCA) on qualitative ones, so surScat runs the one its type calls for: LPA for pca, LCA for
# mca. Both work, as looking4clusters does, on the cases one row each and on the active
# variables as they come, and their groups go into the scattergram the same way (see addL4C).

latentModel   <- function(type) if(type=="pca") "LPA" else "LCA"
latentPackage <- c(LPA="mclust", LCA="poLCA")

# A latent class model starts from random parameters and easily stops at a local maximum,
# so it is estimated this many times over and the best of those runs is kept.
lcaNrep <- 5


## latentClusters ----
# Which latent model surScat runs, if any, out of the methods l4cClusters read (TRUE stands
# there as "latent"), now that type is known. Asked for with TRUE, the model is run only if
# its package is installed, and silently left out otherwise; asked for by name, a missing
# package or a model that does not suit the active variables is reported.
latentClusters <- function(methods, type) {
  model <- latentModel(type)
  if("latent" %in% methods)
    return(if(requireNamespace(latentPackage[[model]], quietly=TRUE)) model else character(0))
  named <- toupper(intersect(methods, c("lpa", "lca")))
  if(!length(named)) return(character(0))
  wrong <- setdiff(named, model)
  if(length(wrong))
    warning(wrong, " left out: ", if(model=="LPA") "LCA needs qualitative active variables (type=\"mca\")"
            else "LPA needs quantitative active variables (type=\"pca\")", call.=FALSE)
  if(!(model %in% named)) return(character(0))
  if(!requireNamespace(latentPackage[[model]], quietly=TRUE)) {
    warning(model, " left out: install '", latentPackage[[model]], "' to compute it", call.=FALSE)
    return(character(0))
  }
  model
}


## latentRun ----
# Estimate the latent model for every number of groups in ks, case by case, and hand back
# the class of every case under each of them, named by that number of groups, and the BIC of
# each, the lower the better (mclust reports it with the opposite sign, which is turned round
# here). A number of groups that cannot be estimated is reported and left out, with an NA
# BIC. seed is set again before every number of groups, so that each result is the same
# whatever else has been asked for alongside it.
# data: the active variables as a data frame, numeric for LPA and factors for LCA.
latentRun <- function(data, model, ks, seed=NULL) {
  reseed <- function() if(!is.null(seed)) set.seed(seed)
  if(model=="LPA") {
    x <- as.matrix(data)
    fit <- function(k) { # equal variances and no covariances: the usual LPA, "EEI" in mclust
      # Mclust calls mclustBIC from the frame of its caller, where it is not found unless
      # mclust is attached, so it is bound here
      mclustBIC <- mclust::mclustBIC
      m <- mclust::Mclust(x, G=k, modelNames="EEI", verbose=FALSE)
      if(is.null(m)) stop("mclust found no solution")
      list(class=as.vector(m$classification), bic=-m$bic)
    }
  } else {
    # poLCA takes categories coded 1, 2... and a formula of plain names
    x <- as.data.frame(lapply(data, function(v) as.integer(as.factor(v))))
    names(x) <- paste0("V", seq_along(x))
    f <- stats::as.formula(paste0("cbind(", paste(names(x), collapse=","), ")~1"))
    fit <- function(k) {
      m <- NULL
      utils::capture.output(m <- poLCA::poLCA(f, x, nclass=k, nrep=lcaNrep, verbose=FALSE,
                                                calc.se=FALSE))
      if(k>1 && m$npar >= prod(vapply(x, max, numeric(1))))
        stop("there are more parameters than response patterns")
      list(class=as.vector(m$predclass), bic=m$bic)
    }
  }
  clusters <- list()
  bic <- stats::setNames(rep(NA_real_, length(ks)), ks)
  for(j in seq_along(ks)) {
    k <- ks[j]
    reseed()
    r <- tryCatch(fit(k), error=function(e) {
      warning(model, " could not be computed for ", k, " groups: ", conditionMessage(e), call.=FALSE)
      NULL
    })
    if(is.null(r)) next
    clusters[[as.character(k)]] <- r$class
    bic[j] <- r$bic
  }
  list(clusters=clusters, bic=bic)
}
