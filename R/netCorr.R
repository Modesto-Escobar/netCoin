# Program to apply nets to correlations
## netCorr ----
netCorr<-function(variables, weight=NULL, pairwise=FALSE,
                  minimum=-Inf, maximum=Inf, sort=FALSE, decreasing=TRUE,
                  frequency=FALSE, means=TRUE, 
                  method=c("pearson", "kendall", "spearman"), criteria="p", Bonferroni=FALSE,
                  minL=0, maxL=Inf,
                  sortL=NULL, decreasingL=TRUE,
                  igraph=FALSE, ...)
{
  arguments <- list(...)
  arguments$name <- nameByLanguage(arguments$name,arguments$language,arguments$nodes)
  method <- r.method(method[1])
  if(!("size" %in% names(arguments))) arguments$size <- "mean"
  if(!("lwidth" %in% names(arguments))) arguments$lwidth <- "value"
  if(!("lweight" %in% names(arguments))) arguments$lweight <- "value"
  if(!pairwise) variables<-na.omit(variables)
  if (inherits(weight,"character")) variables <- variables[, setdiff(names(variables), weight)]
  cases<-nrow(variables)
  if(exists("linkBipolar", arguments) && arguments$linkBipolar) xl <- 2
  else xl <-1
  if (criteria=="p" & maxL==Inf)  maxL<-.5*xl
  if (criteria=="p" & Bonferroni) maxL<-maxL/choose(cases,2)
  if (criteria %in% c("pearson", "kendall", "spearman")) criteria <- "value"
  statistics <-data.frame(name=colnames(variables),
                          mean=round(apply(variables,2,mean, na.rm=TRUE),2),
                          std=round(sqrt(apply(variables,2,var, na.rm=TRUE)),2),
                          min=apply(variables,2,min, na.rm=TRUE),
                          max=apply(variables,2,max, na.rm=TRUE))
  colnames(statistics)[1] <- arguments$name
  if(!is.null(arguments$nodes)) arguments$nodes <- merge(statistics, arguments$nodes, by=arguments$name, all.x=TRUE, sort=FALSE)
  else arguments$nodes <- statistics
  if (pairwise) use <- "pairwise.complete.obs"
  else use <- "complete.obs"
  R<-cor(variables[,arguments$nodes[,2]>=minimum & arguments$nodes[,2]<=maximum],method=method[1], use=use)
  E<-edgeList(R, "shape", min=-1, max=1, directed=FALSE, diagonal=FALSE)
  E$z<-E$value/sqrt(1-E$value^2)*sqrt(cases-2)
  if(exists("linkBipolar", arguments) && arguments$linkBipolar) E$p <- (1-pt(abs(E$z),40000-2))*2
  else xl <- E$p<-1-pt(E$z,cases-2)
  E<-E[E[[criteria]]>=minL & E[[criteria]]<=maxL,]
  if (!is.null(sortL)) E<-E[order((-1*decreasingL+!decreasingL)*E[[sortL]]),]
  arguments$links <- E
  if(exists("layout", arguments) && is.character(arguments$layout) && tolower(substr(arguments$layout,1,2))=="pc") arguments$layout <- layoutPCA(R)
  xNx <- do.call(netCoin,arguments)
  if (igraph) return(toIgraph(xNx))
  else return(xNx)
}

# Program to apply evolving nets to correlations.
## d_netCorr ----
d_netCorr <- function(variables, nodes= NULL, weight=NULL, 
                      pairwise=FALSE, minimum=-Inf, maximum=Inf, 
                      frequency=FALSE, means=TRUE, 
                      method=c("pearson", "kendall", "spearman"), criteria="value", Bonferroni=FALSE,
                      minL=0, maxL=Inf,
                      sortL=NULL, decreasingL=TRUE, 
                      factorial=c("null", "pc", "nf", "vf", "of"), 
                      components=TRUE, backcomponents=FALSE, sequence=seq(.20, 1, .01), 
                      textFilter=c(1, .99), speed=50, dir=NULL, ...)
{
  arguments <- list(variables= variables, nodes=nodes, weight=weight,
                    pairwise=pairwise, minimum=minimum, maximum=maximum, 
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
    H <- addNetCoin(G, linkFilter = paste0("value>",I), lwidth="value",
                    ltext="text", size="mean", layout=C,
                    main=paste0("Correlation: ", I))
    if(exists("hidden", H$links)) H$links <- H$links[!H$links$hidden, ]
    else H$links$hidden <- FALSE
    enodes <- unique(c(H$links$Source,H$links$Target))
    H$nodes <- H$nodes[H$nodes[[H$options$nodeName]] %in% enodes,]
    comps <- igraph::components(toIgraph(H))
    cc <- sum(comps$csize>1)
    if(((cc > ch | !components) | (cc < ch &  components & backcomponents)) & nrow(H$links)>0) {
      central <- calCentr(H)
      H$nodes <- cbind(H$nodes, central$nodes)
      H$nodes <- cbind(H$nodes, compon(comps$membership))
      names(H$nodes)[ncol(H$nodes)] <- "component"
      names(H$nodes)[names(H$nodes)=="degree"] <- "Degree"
      H$links$ratio <- H$links$value/I
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

r.method <- function(method) {
  if(!toupper(substr(method,1,1)) %in% c("P", "K", "S")) method <- "pearson"
  rs <- matrix(c("pearson","kendall", "spearman"), nrow=1, dimnames=list("method", c("P", "K", "S")))
  method <- rs[,toupper(substr(method, 1, 1))]
  return(method)
}