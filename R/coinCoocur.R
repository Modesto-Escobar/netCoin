## coin ----
# List of coincidences

coin<-function(incidences,minimum=1, maximum=nrow(incidences), sort=FALSE, decreasing=TRUE, 
               total=FALSE, subsample=FALSE, weight=NULL, pairwise=FALSE) {
  if (pairwise){
    n <- sum(rowSums(!is.na(incidences))>0)
    nomiss     <- ifelse(is.na(incidences),0,1)
    incidences <- replace(incidences,is.na(incidences),0)
  }
  else {
    if (!is.null(weight)) weight <- weight[rowSums(is.na(incidences))<1]
    incidences<-na.omit(incidences)
  }
  if (subsample){
    vector<-apply(incidences,1,sum)
    incidences<-incidences[vector>0,]
  }
  if (total & is.null(weight)) incidences<-data.frame(Total=1,incidences)
  if (all(is.na(incidences) | incidences==0 | incidences==1)) {
    if (!pairwise) n <- nrow(incidences)
    names(n)<-"n"
    if (is.null(weight)) f<-crossprod(as.matrix(incidences))
    else {
      if (length(weight)!=dim(incidences)[1]) warning("weight has not the appropiate length!")
      f<-crossprod(t(crossprod(as.matrix(incidences),diag(weight,length(weight)))),as.matrix(incidences))
      n<-maximum<-sum(weight)
    }    
    if (is.null(colnames(f))) dimnames(f)<-list(paste("X",1:ncol(f),sep=""),paste("X",1:ncol(f),sep=""))
    d<-diag(f)
    if (sort) d<-sort(d,decreasing=decreasing)
    S<-names(d[(d>=minimum &  d<=maximum)])
    if (total & is.null(weight)) S<-c("Total",S)
    if (total & !is.null(weight)) warning("total cannot be applied in weighted tables")
    if (length(S)>0) {
      if (!pairwise) structure(f[S,S], n=n, class=c("coin"))
      else {
        colnames(nomiss) <- colnames(incidences)
        nomiss <- nomiss[,S]
        incidences <- incidences[,S]
        if (is.null(weight)) {
          m<-crossprod(as.matrix(nomiss))
          x<-crossprod(1-as.matrix(nomiss),as.matrix(incidences))
        }
        else {
          if (length(weight)!=dim(nomiss)[1]) warning("weight has not the appropiate length!")
          m<-crossprod(t(crossprod(as.matrix(nomiss),diag(weight,length(weight)))),as.matrix(nomiss))
          x<-crossprod(t(crossprod(1-as.matrix(nomiss),diag(weight,length(weight)))),as.matrix(incidences))
        }
        structure(f[S,S], n=n, m=m[S,S], x=x[S,S], class=c("coin"))
      }
    }
    else cat("No variables left")
  }
  else warning("All data in incidence matrix has to be dichotomous.")
}

## coocur ----

coocur<-function (ocurrences, minimum = 1, maximum = Inf, sort = FALSE, decreasing=TRUE) 
{
  result <- matrix(nrow = ncol(ocurrences), ncol = ncol(ocurrences), dimnames = list(colnames(ocurrences), 
                                                                                     colnames(ocurrences)))
  for (val in c(1:nrow(result))) {
    for (cal in c(val:ncol(result))) {
      result[cal, val] <- sum(pmin(ocurrences[, cal], ocurrences[, val]))
      if (val != cal) 
        result[val, cal] = result[cal, val]
    }
  }
  d<-diag(result)
  if (sort) d<-sort(d,decreasing=decreasing)
  S<-names(d[(d>=minimum &  d<=maximum)])
  if (length(S)>=1) {
    result<-result[S,S]
    if (length(S)==1) names(result)<-S
    n <- sum(ocurrences[,S])
    m <- sum(apply(as.matrix(ocurrences[,S]), 1, max))
    attr(result, "n") <- n
    attr(result, "m") <- m
    structure(result, class = "cooc")
  }
  else cat("No variables left")
}

## Methods ----
print.coin<-function(x, ...) {
  cat("n= ",attr(x,"n"),"\n",sep="")
  print(lower(x,0))
}

print.cooc<-function(x, ...) {
  cat("n= ",attr(x,"n"),"; m= ", attr(x,"m"),"\n",sep="")
  print(lower(x,0))
}

summary.coin <- function(object, ...){
  cat(attr(object,"n"),"scenarios and", dim(object[,])[1], "events\n")
  diag(object[,])/attr(object,"n")
}

plot.coin <- function(x, dir=tempDir(), language=c("en","es","ca"), ...){
  N <- asNodes(x, language = language)
  colnames(N)[2] <- "incidences"
  E <- edgeList(x,c("Frequencies","Expected"))
  bar <- barplot_rd3(N, E, name = names(N)[1], coincidences = "coincidences", incidences = "incidences", expected = "expected", cex = 1, language = language)
  pie <- pieCoin(x, cex = 1, language = language)
  rd3_multigraph(barplot=bar, piechart=pie, dir = dir)
}


propCoin<-function(x, margin= 0, decimals=1) {
  if (!inherits(x,"coin")) stop("Error: input must be a coin object (see coin function)")
  if ("m" %in% names(attributes(x))) n <- attr(x, "m")
  else n <- attr(x,"n")
  x <- x[,]
  switch(format(margin),
         "0" = round(100 * x / n, decimals),
         "1" = round(100 * x / diag(x), decimals),
         "2" = round(t(100 * x / diag(x)), decimals))
}

# Transform a coin object into a data frame with name and frequency
asNodes<-function(C, frequency = TRUE, percentages = FALSE, language = c("en","es","ca")){
  nodes <- NULL
  if (inherits(C,"coin")) {
    if ("m" %in% names(attributes(C))) divider <- diag(attr(C,"m"))
    else divider <- attr(C,"n")
    if (!percentages & frequency) nodes<-data.frame(name=as.character(colnames(C)),frequency=diag(C))
    else if (!frequency & percentages) nodes<-data.frame(name=as.character(colnames(C)),"%"=diag(C)/divider*100,check.names=FALSE)
    else if (percentages & frequency)nodes<-data.frame(name=as.character(colnames(C)),frequency=diag(C), "%"=diag(C)/divider*100,check.names=FALSE)
    else nodes<-data.frame(name=as.character(colnames(C)),check.names=FALSE)
    if(language[1]!="en"){
      colnames(nodes)[colnames(nodes)=="frequency"] <- getByLanguage(frequencyList,language)
      colnames(nodes)[colnames(nodes)=="name"] <- getByLanguage(nameList,language)
    }     
  }
  else if (min(c("Source", "Target") %in% names(C))) nodes<-data.frame(name=sort(union(C$Source,C$Target)))
  else warning("Is neither a coin object or an edge data frame")
  return(nodes)
}

## Subfunctions ----
summaryNet <- function(x){
  cat(dim(x$nodes)[1], "nodes and", dim(x$links)[1], "links.\n")
  freq <- frequencyList[frequencyList %in% names(x$nodes)]
  if(length(freq)==1) {
    cat(freq," distribution of nodes:","\n", sep="")
    print(summary(x$nodes[[freq]]))
  }
  lwidth <- NULL
  if(!is.null(x$options$linkWidth))
    lwidth <- x$options$linkWidth
  else if(length(x$links)>2)
    lwidth <- names(x$links)[3]
  if(length(lwidth)==1){
    cat(lwidth, "'s distribution:","\n",sep="")
    print(summary(x$links[[lwidth]]))
  }
}


tempDir <- function(){
  dir.create("temp", showWarnings = FALSE)
  return(paste("temp",round(as.numeric(Sys.time())),sep="/"))
}