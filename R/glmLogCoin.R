## glmCoin ----
glmCoin <- function(formulas, data, weights=NULL, pmax=.05,
        robust=FALSE, showArrows=TRUE,
        color="variable", lwidth="s.value", lcolor="estimate",
        language=c("en","es","ca"),
        igraph=FALSE, ...) {
  vcov <- stats::vcov
  if(robust){
    vcov <- "HC1"
  }

  cleanVariables <- function(x){
    x <- iconv(x,to="ASCII//TRANSLIT")
    x <- gsub(" ",".",x)
    return(gsub("[^a-zA-Z0-9_]",".",x))
  }
  originalNames <- colnames(data)
  colnames(data) <- cleanVariables(colnames(data))
  cleanedNamesIndex <- which(originalNames!=colnames(data))
  aux <- originalNames[cleanedNamesIndex]
  cleanedNamesIndex <- cleanedNamesIndex[order(nchar(aux),decreasing=TRUE)] # reorder to check first longest strings
  for(i in cleanedNamesIndex){
    formulas <- gsub(originalNames[i],colnames(data)[i],formulas,fixed=TRUE)
  }
  if (is.character(weights)){
    weights <- cleanVariables(weights)
    if(!weights %in% colnames(data)) stop("Weights variable '",weights,"' not found in data")
    weights<-data[[weights]]
  }

  prenet <- contr.gw(formulas, data=data, weights=weights, vcov=vcov, .keep_models=TRUE)
  equations <- eqSummary(prenet$models, formulas)
  arguments <- list(nodes = prenet$Nodes, links = prenet$Links,
    showArrows=showArrows, color=color, linkFilter=paste0("p.value<",pmax),
    lwidth=lwidth, lcolor=lcolor, language=language, ...)
  if(is.null(arguments$linkBipolar)){
    arguments$linkBipolar <- TRUE
  }
  if(length(cleanedNamesIndex)){
    cleanedToOriginal <- originalNames[cleanedNamesIndex]
    names(cleanedToOriginal) <- colnames(data)[cleanedNamesIndex]
    restoreNames <- function(x){ # labels are variable[*variable][:categories]
      vapply(as.character(x), function(s){
        sep <- regexpr(":", s, fixed=TRUE)
        vars <- if(sep>0) substr(s, 1, sep-1) else s
        rest <- if(sep>0) substr(s, sep, nchar(s)) else ""
        tokens <- strsplit(vars, "*", fixed=TRUE)[[1]]
        found <- match(tokens, names(cleanedToOriginal))
        tokens[!is.na(found)] <- cleanedToOriginal[found[!is.na(found)]]
        paste0(paste(tokens, collapse="*"), rest)
      }, character(1), USE.NAMES=FALSE)
    }
    restoreFormulas <- function(x){ # replace whole name tokens only
      for(i in cleanedNamesIndex){
        pattern <- paste0("(?<![a-zA-Z0-9._])", gsub(".","\\.",colnames(data)[i],fixed=TRUE), "(?![a-zA-Z0-9._])")
        x <- gsub(pattern, originalNames[i], x, perl=TRUE)
      }
      return(x)
    }
    if(!is.null(arguments$nodes)) arguments$nodes[,'name'] <- restoreNames(arguments$nodes[,'name'])
    if(!is.null(arguments$links) && nrow(arguments$links)){
      arguments$links[,'Source'] <- restoreNames(arguments$links[,'Source'])
      arguments$links[,'Target'] <- restoreNames(arguments$links[,'Target'])
      arguments$links[,'Model'] <- restoreFormulas(arguments$links[,'Model'])
    }
    if(!is.null(equations)) equations$model <- restoreFormulas(equations$model)
  }
  if(arguments$language[1]!="en"){
    colnames(arguments$nodes)[colnames(arguments$nodes)=="name"] <- getByLanguage(nameList,arguments$language[1])
  }
  net <- do.call(netCoin, arguments)
  net$equations <- equations
  if(igraph) net <- igraph::set_graph_attr(toIgraph(net), "equations", equations)
  return(net)
}

## See file contrast.gw.R for more

## Fit measures per equation for glmCoin ----
eqSummary <- function(models, spec){
  declared <- vapply(.parse_specs(spec), function(e){
    if(is.null(e$family)) "gaussian"
    else if(is.character(e$family)) e$family
    else e$family$family
  }, character(1))
  modelLab <- vapply(names(models), function(s) .split_top_level(s)$lhs, character(1), USE.NAMES=FALSE)
  rows <- lapply(seq_along(models), function(i){
    m <- models[[i]]
    mf <- stats::model.frame(m)
    y <- stats::model.response(mf)
    w <- stats::model.weights(mf)
    if(is.null(w)){ # multinom's model.frame drops the weights column
      w <- tryCatch(eval(m$call$weights, environment(stats::formula(m))), error=function(e) NULL)
      if(!is.null(w) && length(w)!=nrow(mf)) w <- NULL
    }
    nobs <- tryCatch(stats::nobs(m), error=function(e) nrow(mf))
    if(is.null(nobs) || is.na(nobs)) nobs <- nrow(mf)
    n <- if(is.null(w)) nobs else sum(w) # weighted N: deviances and logLik already scale with it
    nullArgs <- list(formula = y ~ 1)
    if(!is.null(w)) nullArgs$weights <- w
    m0 <- tryCatch({
      if(inherits(m,"multinom")) do.call(nnet::multinom, c(nullArgs, list(trace=FALSE)))
      else if(inherits(m,"glm")) do.call(stats::glm, c(nullArgs, list(family=m$family)))
      else do.call(stats::lm, nullArgs)
    }, error=function(e) NULL)
    dev  <- suppressWarnings(tryCatch(stats::deviance(m), error=function(e) NA_real_))
    dev0 <- if(is.null(m0)) NA_real_ else suppressWarnings(tryCatch(stats::deviance(m0), error=function(e) NA_real_))
    ll   <- suppressWarnings(tryCatch(as.numeric(stats::logLik(m)), error=function(e) NA_real_))
    ll0  <- if(is.null(m0)) NA_real_ else suppressWarnings(tryCatch(as.numeric(stats::logLik(m0)), error=function(e) NA_real_))
    r2 <- 1-dev/dev0
    mcfadden <- 1-ll/ll0
    if(!is.na(ll) && !is.na(ll0)){ # likelihood-based Cox-Snell; deviance-based for quasi families
      nagelkerke <- (1-exp(2*(ll0-ll)/n))/(1-exp(2*ll0/n))
    } else {
      nagelkerke <- (1-exp((dev-dev0)/n))/(1-exp(-dev0/n))
    }
    aic <- suppressWarnings(tryCatch(stats::AIC(m), error=function(e) NA_real_))
    bic <- suppressWarnings(tryCatch(stats::BIC(m), error=function(e) NA_real_))
    dfres <- tryCatch(stats::df.residual(m), error=function(e) NA_real_)
    if(is.null(dfres)) dfres <- NA_real_
    data.frame(model=modelLab[i], family=declared[i], n=round(n), n_obs=nobs,
      r2=round(r2,3), r2.mcfadden=round(mcfadden,3), r2.nagelkerke=round(nagelkerke,3),
      logLik=round(ll,2), aic=round(aic,2), bic=round(bic,2),
      deviance=round(dev,2), null.deviance=round(dev0,2), df.residual=dfres,
      stringsAsFactors=FALSE, row.names=NULL)
  })
  do.call(rbind, rows)
}


## logCoin() ----
logCoin<-function(data, variables=names(data), exogenous=NULL, noFirstCat=NULL, weight=NULL, 
                  order=2, pairwise=FALSE, twotails=FALSE, pmax=.05,
                  frequency=FALSE, percentage=FALSE, 
                  directed=FALSE, igraph=FALSE, ...) {
  arguments <- list(...)
  names(data) <- gsub(" ","_", names(data))
  variables   <- gsub(" ","_", variables)
  if(!is.null(exogenous))   exogenous  <- gsub(" ","_", exogenous)
  if (!is.null(noFirstCat)) noFirstCat <- gsub(" ","_", noFirstCat)
  if(!is.null(weight))      weight     <- gsub(" ","_", weight)
  variables <- union(setdiff(union(variables, noFirstCat), exogenous), exogenous)
  
  if(!is.null(weight)) {
    if(inherits(weight,"character")){
      variables <- setdiff(variables, weight)
      weight<-data.frame(weight=data[[weight]])
      data<-cbind(weight,data[,variables])
    }
    else{
      if(length(weight)!=dim(data)[1]) stop("Weights have not the correct dimensions")
      if (pairwise) data <- cbind(data[,variables],weight)[,1:length(data[,variables])]
      else data <- na.omit(cbind(data[,variables],weight))[,1:length(data[,variables])]
    }
  }
  else data<-cbind(data[,variables], data.frame(weight=rep(1, nrow(data))))
  
  
  formula   <- comb(setdiff(variables,exogenous), exogenous, "weight", order)
  varOrder  <- variables # To order variables later before coin
  #Check methods. No necessary because edgeList call these routines.
  #procedures<-i.method(c_method(procedures))
  #criteria<-i.method(c_method(criteria))
  #procedures<-union(procedures,unlist(arguments[c("lwidth","lweight","lcolor","ltext")]))
  
  #Names  
  if(!("language" %in% names(arguments))) arguments$language <- "en"
  nodes <- arguments$nodes
  if (inherits(nodes,"tbl_df")) nodes<-as.data.frame(nodes)
  name <- arguments$name <- nameByLanguage(arguments$name,arguments$language,arguments$nodes)
  if (!("level" %in% names(arguments))) level<-.95 else level <-arguments$level
  
  #Data.frame  
  if (all(inherits(data,c("tbl_df","tbl","data.frame"),TRUE))) data<-as.data.frame(data) # convert haven objects
  if (inherits(weight,"character")) variables <- setdiff(variables,weight)
  pivots <- setdiff(variables, noFirstCat)
  if (length(pivots)>5) stop("This function doesn't support more than 5 variables with first categories included")
  
  if (!pairwise & inherits(weight,"character")) {
    if (!is.null(weight)) weight <- data[rowSums(is.na(data[,variables]))<1,weight]
    data <- data[complete.cases(data[,variables]),]
  }
  
  arguments$scenarios <- sum(rowSums(!is.na(data))>0)
  
  data <- as_factor(data)
  dt <- aggregate(formula, data=data, FUN="sum")
  dt$weight <- round(dt$weight) # To estimate Poisson (without decimals)
  
  sc <- sum(rowSums(!is.na(dt))>0)
  
  dtab <- xtabs(weight ~ ., dt)
  fm <- loglm(comb(setdiff(variables,exogenous), exogenous, "", order), dtab)  # numerals as names.
  
  
  nam <- eti <- var <- cat <- varn <- labn <- NULL
  nvar <-1
  ncell <- 1
  for (x in variables) {
    dt[[x]] <- as_factor(dt[[x]])
    nval <- length(levels(dt[[x]]))
    nam <- c(nam, paste0(x, levels(dt[[x]])))
    eti <- c(eti, paste0(x, ":", levels(dt[[x]])))
    var <- c(var, rep(x, nval))
    cat <- c(cat, levels(dt[[x]]))
    varn <- c(varn, rep(nvar, nval))
    labn <- c(labn, 1:nval)
    contrasts(dt[[x]]) <- contr.first(dt[[x]])
    nvar <- nvar+1
    ncell <- ncell*nval
  }
  informa <- data.frame(eti, var, cat, varn, labn, row.names = nam)
  
  arguments$note <- paste0("<p>L2= ", sprintf("%3.2f",fm$lrt), "; d.f.= ", fm$df, "; p(>X^2)= ", 
                           sprintf("%3.3g",1-pchisq(fm$lrt, fm$df)), ".</p><p>Covariance structures: ", 
                           sprintf("%1.0f", sc)," (",sprintf("%1.0f", ncell),").</p><p>",
                           Reduce(paste, deparse(formula)),"</p>")
  
  coefs <- summary(glm(formula, family=poisson, data=dt))$coefficients
  
  if (length(pivots)==0) coefs <- list(data=dt, coefs=coefs)
  
  if (length(pivots) >0) {
    coefs <- addCoefs(c(1,1), coefs, formula, dt, pivots)
  }
  
  if (length(pivots) >1) {
    coefs <- addCoefs(c(1,2), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(0,1), coefs$coef, formula, coefs$data, pivots)
  }
  
  if (length(pivots) >2) {
    coefs <- addCoefs(c(1,3), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(1,1), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(0,2), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(0,1), coefs$coef, formula, coefs$data, pivots)
  }
  
  if (length(pivots) >3) {
    coefs <- addCoefs(c(1,4), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(1,1), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(1,2), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(0,1), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(0,3), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(1,1), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(0,2), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(0,1), coefs$coef, formula, coefs$data, pivots)
  }
  
  if (length(pivots) >4) {
    coefs <- addCoefs(c(1,5), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(1,1), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(1,2), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(0,1), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(1,3), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(1,1), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(0,2), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(0,1), coefs$coef, formula, coefs$data, pivots)
    
    coefs <- addCoefs(c(0,4), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(1,1), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(1,2), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(0,1), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(0,3), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(1,1), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(0,2), coefs$coef, formula, coefs$data, pivots)
    coefs <- addCoefs(c(0,1), coefs$coef, formula, coefs$data, pivots)
  }
  
  coefs <- as.data.frame(coefs$coef[sort(row.names(coefs$coef)),])
  nodes <- coefs[!grepl(":",rownames(coefs)), ][-1,1, drop=F]
  
  # nodes
  nodes$name <- row.names(nodes)
  nodes$Estimate <- round(nodes$Estimate, 3)
  nodes <- merge(nodes, informa [1:3], by="row.names")[,-1]
  if(frequency | percentage) nodes <- merge(nodes, margins(dtab, frequency, percentage), by="name", all.x=T)
  ordervars  <- factor(nodes$var, levels=varOrder)
  ordernodes <- nodes[order(ordervars), "name"]
  # coefs
  pmax<-abs(ifelse(twotails,qnorm(pmax/2,FALSE),qnorm(pmax,FALSE)))
  if(!twotails) {
    coefs <- coefs[coefs$`z value`>pmax & grepl(":",rownames(coefs)), ]
    coefs$`Pr(>|z|)` <- coefs$`Pr(>|z|)` /2
    names(coefs)[names(coefs)=="Pr(>|z|)"] <- "Pr(>z)"
  }
  else {
    coefs <- coefs[abs(coefs$`z value`)>pmax & grepl(":",rownames(coefs)), ]
    arguments$linkBipolar <- TRUE
    arguments$lcolor <- "Estimate"
  }
  coefs$order <- nchar(gsub("[^:]","",row.names(coefs)))+1
  links <- coefs[coefs$order==2,]
  links$Source <- sub(":.*","", row.names(links))
  links$Target <- sub(".*:","", row.names(links))
  if(nrow(links)){
    links[,c("Source", "Target")] <- olinks(links$Source, links$Target, ordernodes)
  }
  
  ulinks <- coefs[coefs$order>2,]
  if (nrow(ulinks)>0) {
    L <- unlist(strsplit(row.names(ulinks),":"))
    CC <- paste0("u", sprintf(paste0("%0",ceiling((nrow(ulinks)+1)/10),"d"),1:nrow(ulinks)))
    ulinks$Source <- CC
    M <- data.frame(Source=rep(CC, times=ulinks$order), Target=L)
    unodes <- data.frame(name=CC, Estimate=floor(min(nodes$Estimate)), eti="", var=paste0(".I", ulinks$order), cat="")
    if(frequency)  unodes$freq <- sapply(unodes$name, FUN=ofreq, x=M, table=dtab, index=informa)
    if(percentage) unodes$prop <- sapply(unodes$name, FUN=ofreq, x=M, table=dtab, index=informa)/marginSums(dtab)
    nodes  <- rbind(nodes, unodes)
    ulinks <- merge(M, ulinks)
    links  <- rbind(links, ulinks)
  }
  if(exists("language", arguments)) {
    names(nodes)[names(nodes)=="name"] <- getByLanguage(nameList,  arguments$language)
    names(nodes)[names(nodes)=="eti"]   <- getByLanguage(labelList, arguments$language)
  }
  arguments$nodes <- nodes
  if(nrow(links)){
    arguments$links <- links
    if (!"degreeFilter" %in% names(arguments)) arguments$degreeFilter <- 1
  }
  if (!"size" %in% arguments) {
    arguments$size= "Estimate"
    if(frequency)  arguments$size <- "freq"
    if(percentage) arguments$size <- "prop"
  }
  if (!"color" %in% names(arguments)) arguments$color <- "var"
  if (!"lwidth" %in% names(arguments)) arguments$lwidth <- "Estimate"
  if (!"label" %in% names(arguments)) arguments$label <- getByLanguage(labelList, arguments$language)
  G <- do.call(netCoin, arguments)
  return(G)
}


# Atencion a lwidth y en lenguaje espanol o catalan
glmCoin2 <- function(formulas, data, weights=NULL, pmax=.05, twotail=FALSE, showArrows=TRUE,
                     frequency = FALSE, percentage = TRUE, 
                     color="variable", lwidth="z.value", circle= NA, language=c("en","es","ca"),
                     igraph=FALSE, ...){
  cleanVariables <- function(x){
    x <- iconv(x,to="ASCII//TRANSLIT")
    x <- gsub(" ",".",x)
    return(gsub("[^a-zA-Z0-9_]",".",x))
  }
  originalNames <- colnames(data)
  colnames(data) <- cleanVariables(colnames(data))
  cleanedNamesIndex <- which(originalNames!=colnames(data))
  aux <- originalNames[cleanedNamesIndex]
  cleanedNamesIndex <- cleanedNamesIndex[order(nchar(aux),decreasing=TRUE)] # reorder to check first longest strings
  for(i in cleanedNamesIndex){
    formulas <- gsub(originalNames[i],colnames(data)[i],formulas,fixed=TRUE)
  }
  if (is.character(weights)){
    weights <- cleanVariables(weights)
    weights<-data[[weights]]
  }
  Links <- data.frame(A=NA,B=NA,C=NA,D=NA, E=NA, G=NA, H=NA)[-1,]
  names(Links)  <- headreg(language)
  if (lwidth=="z.value" & language[1]!="en") lwidth<-"val.z"
  Formulas  <- as.list(gsub("[[:space:]]","",unlist(strsplit(formulas,"\n"))))
  formulas  <- sapply(Formulas,function(X){substr(X,start=1,stop=as.numeric(gregexpr(",",X))-1)})
  familias  <- sapply(Formulas,family)
  # dependent <- sapply(Formulas,function(X){substr(X,start=1,stop=as.numeric(gregexpr("\u7E",X))-1)})
  # dependent <- gsub("\u60","",dependent)
  variables<-extract(formulas, data)
  
  for(instance in 1:length(formulas)) {
    m <- net_lm(formulas[instance], familias[instance], data, weights, pmax, twotail)
    if (nrow(m)>0) {
      names(m) <- headreg(language)[-7]
      formula <- gsub("`","",formulas[instance])
      for(i in cleanedNamesIndex){
        m[,'Source'] <- gsub(colnames(data)[i],originalNames[i],m[,'Source'],fixed=TRUE)
        m[,'Target'] <- gsub(colnames(data)[i],originalNames[i],m[,'Target'],fixed=TRUE)
        formula <- gsub(colnames(data)[i],originalNames[i],formula,fixed=TRUE)
      }
      m[[headreg(language)[7]]] <- formula
      Links <-rbind(Links,m)
    }
  }
  
  if (!("nodes" %in% names(list(...)))) {
    Nodes<-data.frame(name=iconv(union(Links$Source,Links$Target),to="UTF-8"),
                      variable=gsub(":.*","",iconv(union(Links$Source,Links$Target),to="UTF-8"))
                      ,stringsAsFactors = FALSE)
    row.names(Nodes)<-Nodes$name
    arguments<-list(nodes=Nodes, links=Links, showArrows = showArrows, color = color, lwidth = lwidth, language = language, ...)
  }
  else arguments<-list(links=Links, showArrows = showArrows, color = color, lwidth = lwidth, language = language, ...)
  
  if (frequency | percentage) arguments$nodes<-meanPer(data, variables, arguments$nodes, names(arguments$nodes)[1] , frequency, percentage, weights)
  
  #ADD to N percentages/means. vid extract.R and as.nodes(surCoin)
  
  if (!"name" %in% names(arguments)) {
    arguments$name <- nameByLanguage(arguments$name,language,arguments$nodes)
    names(arguments$nodes)[1]<-nameList[language[1]]
  }
  
  if (!is.na(circle)) arguments$layout <- layoutCircle(arguments$nodes, variables$D, circle)
  
  if(twotail){
    arguments$linkBipolar <- TRUE
    arguments$lcolor <- "Estimate"
  }
  
  if(nrow(arguments$nodes)+nrow(arguments$links)>0) {
    xNx <- do.call(netCoin,arguments)
  }
  else(stop("No nodes, no relations"))
  # xNx$nodes$name<-iconv(xNx$nodes$name,to="UTF-8")
  # xNx$links[[c("Source","Target")]]<-iconv(xNx$links[[c("Source","Target")]],to="UTF-8")
  if (igraph) return(toIgraph(xNx))
  else return(xNx)
}


## subfunctions ----
headreg<-function(language){
  switch (language[1],
          en = c('Source','Target','Estimate','Std.error', 'z.value', 'Pr(>|z|)','Equation'),
          es = c('Source','Target','Estimador','Err.t\uEDp.', 'val.z', 'Pr(>|z|)','Ecuaci\uF3n'),
          ca = c('Source','Target','Estimador','Err.t\uEDp.', 'val.z', 'Pr(>|z|)', 'Ecuaci\uF3')
  )
}

extract <- function(formulas, data) {
  formulas<-gsub("`","",formulas)
  dependent<-variables<-factors<-NULL
  for (formula in formulas) {
    dependent <- union(dependent, gsub("[[:space:]]","",unlist(strsplit(strsplit(formula,"\\~")[[1]][1],"\\+"))))
    variables <- union(variables, gsub("[[:space:]]","",unlist(strsplit(strsplit(formula,"\\~")[[1]][2],"\\+"))))
    factors   <- union(factors, variables[vapply(data[variables],inherits,TRUE,what="factor")])
  }
  independent<-setdiff(variables,factors)
  variables<-list(D=dependent, I=independent, F=factors)
  return(variables)
}

meanPer<-function(data, variables, frame, name=names(frame[1]), frequency= FALSE, percentage= TRUE, weights = NULL){
  if (is.null(weights)) weights <- rep(1, nrow(data))
  columns<-setdiff(names(frame),c(name, "n.","%"))
  l.frame<-length(frame)
  row.names(frame)<-frame[[name]]
  frame.order<-row.names(frame)
  quantitatives<-c(variables$D,variables$I)
  data<-na.omit(cbind(as.data.frame(data)[,unlist(variables)],weights))
  for (varF in variables$F) {
    beginF <- length(data)+1
    data<-dichotomize(data, varF, "", add=TRUE)
    names(data)[beginF:length(data)] <- paste0(varF,":",names(data)[beginF:length(data)])
  }
  setvariables <- setdiff(names(data),union(weights,variables$F))
  sta<-data.frame(names=setvariables)
  row.names(sta)<-sta$names
  sta$N. <- 0
  if (frequency & !is.null(weights)) sta$N. <-round(apply(data[,setvariables]*data$weights, 2, sum),0)
  sta$N. <- ifelse(sta$names %in% quantitatives, nrow(data), sta$N.)
  if (percentage) {
    means<-apply(data[,setvariables], 2, weighted.mean, data$weights)
    maxs <-apply(data[,setvariables], 2, max)
    mins <-apply(data[,setvariables], 2, min)
    sta$M.<-(means-mins)/(maxs-mins)*100
  }
  frame<-merge(frame,sta[,-1,drop=F], by="row.names", all.x = TRUE)[,-1]
  row.names(frame)<-frame[[name]]
  adds<-c()
  if (frequency)  {
    frame$n. <- frame$N.
    adds <- "n."
  }
  if (percentage) {
    frame$`%`<- frame$M.
    adds <- c(adds, "%")
  }
  out<- grep('^[KLMN]\\.$',names(frame))
  frame<-frame[,-c(out)]
  return(frame[frame.order,c(name,adds,columns)])
}


# Funciones previas
contr.first <- function(variable){
  if(!inherits(variable, "factor")) stop
  nlevels <- length(levels(variable))
  f <- matrix(c(rep(-1, nlevels-1),diag(nlevels-1)), nrow=nlevels, byrow=TRUE, dimnames=list(levels(variable),levels(variable)[-1]))
  return(f)
}

contr.last <- function(variable) {
  if(!inherits(variable, "factor")) stop
  nlevels <- length(levels(variable))
  l <- matrix(c(diag(nlevels-1), rep(-1, nlevels-1)), nrow=nlevels, byrow=TRUE, dimnames=list(levels(variable),levels(variable)[-nlevels]))
  return(l)
}

addCoefs <- function(vector=c(0,1), coefficients, formula, data, variables, family="poisson") {
  if(vector[1]==1) contrasts(data[[variables[[vector[2]]]]]) <- contr.last(data[[variables[[vector[2]]]]])
  else             contrasts(data[[variables[[vector[2]]]]]) <- contr.first(data[[variables[[vector[2]]]]])
  coefs <- rbind(coefficients,
                 summary(glm(formula, family=family, data=data))$coefficients)
  return(list(data=data, coef=coefs[unique(rownames(coefs)),]))
}


comb <-function (endogenous, exogenous=NULL, frequency="x", order=2) {
  endogenous <- setdiff(endogenous, exogenous)
  combin <-  function(x, order) {apply(combn(x, order), 2, paste0, collapse="*")}
  formula <- paste0(frequency,"~")
  ordEndo <- min(order, length(endogenous))
  l <- combin(endogenous, ordEndo)
  if (is.null(exogenous) | length(endogenous)>=order ) for (nl in l) formula <- paste0(formula,nl,"+")
  if (!is.null(exogenous)) {
    for (exo in exogenous) {
      if(order<length(endogenous)+length(exogenous)) l<- combin(endogenous, min(ordEndo,order-1))
      nexo <- ifelse(l[1]=="","","*")
      for (nl in l) formula <- paste0(formula,nl,nexo,exo, "+")
    }
  }
  return(as.formula(sub("\\+$","",formula)))
}

margins <- function(table, freq=T, prop=F) {
  names <- names(attr(table,"dimnames"))
  vect  <- c()
  n <- ifelse(prop, marginSums(table),1)
  for (x in names) {
    marg <- marginSums(table, x)
    names(marg) <- paste0(x, names(marg))
    vect <- c(vect, marg)
  }
  dtfrm <- data.frame(name=names(vect), freq=vect)
  if (prop) dtfrm$prop <- dtfrm$freq/n
  if (!freq) dtfrm$freq <- NULL
  return(dtfrm)
}

olinks <- function(source, target, order) {
  for (i in 1:length(source)) {
    if (which(source[[i]]==order) < which(target[[i]]==order)) {
      s <- target[[i]]
      target[[i]] <- source[[i]]
      source[[i]] <- s
    }
  }
  return(cbind(source, target))
}

ofreq <- function (u, x, table, index) {
  vector <- x[x$Source==u, "Target"]
  return(marginSums(table, index[vector, 4])[t(index[vector,5])])
}

net_lm<-function(formula, family=gaussian, data, weights=NULL, pmax=.05, twotail=FALSE, ...){
  arguments <- list(formula=formula, family=family, data= data, weights= weights, ...)
  glm       <- do.call(glm, arguments)
  links     <- linkregress(remodel(glm),pmax=pmax,twotails=twotail)
  return(links)
}

remodel<-function(model){
  DD <- model$model
  FO <- model$formula
  FF <- model$family$family
  if (!is.null(model$model$`(weights)`)) WW <- model$model$`(weights)`
  else WW<-NULL
  C  <- model$coefficients
  ## VV<-row.names(attr(terms(F),"factors"))[-1]
  VV <- gsub("[[:space:]]","",unlist(strsplit(strsplit(FO,"\\~")[[1]][2],"\\+")))
  V <- VV[vapply(DD[VV],inherits,TRUE,what="factor")]
  if (length(V)==0) return(model)
  Vq <- paste0("^",V)
  y  <- lapply(Vq,grep,names(model$coefficients))
  names(y) <- Vq
  if (is.null(names(y))) colnames(y)<-gsub("\\^","",colnames(y))
  else names(y)<-gsub("\\^","",names(y))
  
  c<-f<-m<-z<-list()
  for (i in names(y)) {
    z[i]<-list(C[y[[i]]])
    m[i]<-which.min(unlist(z[i]))
    f[i]<-ifelse(z[[i]][m[[i]]]<0,m[[i]],0)
    if(f[i]>0) DD[[i]] <- relevel(DD[[i]],sub(names(z[i]),"",names(z[[i]])[f[[i]]]))
    c[[i]]<-matriz(DD[[i]])
  }
  if (is.null(WW)) neomodel<-glm(FO, contrasts= c, data=DD, family=FF)
  else neomodel <- glm(FO, contrasts= c, data=cbind(DD,WW), family=FF, weights = WW)
  return(neomodel)
}


linkregress<-function(model, pmax=0.05, twotails=FALSE){
  pmax<-abs(ifelse(twotails,qnorm(pmax/2,FALSE),qnorm(pmax,FALSE)))
  q<-summary(model)$coefficients
  y<-names(model$xlevels) #names(model$model)[-1]
  for (i in y) {
    rownames(q)<-gsub(paste0("(^",i,")"),"\\1:",rownames(q))
  }
  m<-data.frame(Source=rownames(q),Target=names(model$model)[1], q, stringsAsFactors = FALSE)
  rownames(m)<-NULL
  names(m)[-c(1:2)]<-colnames(q)
  m<-m[-1,]
  crit<-names(m)[grep(" value",names(m))]
  if (twotails) m <- m[abs(m[[crit]])>pmax,]
  else m <- m[m[[crit]]>pmax,]
  row.names(m)<-NULL
  return(m)
}

matriz<-function(factor){
  T<-t(table(factor)/length(factor))
  M<-matrix(rep(T,length(T)),nrow=length(T),
            dimnames=list(colnames(T),colnames(T)))
  M[,2:dim(M)[1]]<--M[,2:dim(M)[1]]
  diag(M)[2:dim(M)[1]]<-1-M[2:dim(M)[1],1]
  m<-t(matrix(M[,2:dim(M)[2]],nrow=dim(M)[2],dimnames=list(colnames(T),colnames(T)[-1])))
  m.t<-rbind(constant=1/ncol(m),m)
  matriz<-matrix(solve(m.t)[,-1], nrow=nrow(m.t), dimnames=list(colnames(m.t),rownames(m.t)[-1]))
  return(matriz)
}

family<-function(formula) {
  FAM <- ifelse(as.numeric(gregexpr(",",formula))==-1,
                "",
                substr(formula,start=as.numeric(gregexpr(",",formula))+1,1E6))
  family<-c_family(FAM)
  return(family)
}
c_family<-function(method=NULL) {
  if (is.null(method) | method=="") family<-"GAU"
  else {
    family<-substr(toupper(method),1,3)  
    if (family=="QUA" & nchar(method)>5) family<-substr(toupper(method),1,6)
  }
  families<-matrix(c("gaussian","binomial","Gamma","inverse.gaussian","poisson","quasi","quasibinomial","quasipoisson"),
                   nrow=1, dimnames=list("Family", c("GAU","BIN","GAM","INV","POI","QUA","QUASIB","QUASIP")))
  return(families[,family])
}
