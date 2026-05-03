## Path and cobCoin
# Elaborate a netCoin object from a lavaan object.
## pathCoin ----
pathCoin<-function(model, estimates=c("b","se","z","pvalue","beta"), fitMeasures=c("chisq", "cfi", "rmsea"), ...){
  arguments <- list(...)
  if(!("language" %in% names(arguments))) arguments$language <- "en"
  if(!("linkBipolar" %in% names(arguments))) arguments$linkBipolar <- TRUE
  if(!("showArrows" %in% names(arguments))) arguments$showArrows <- TRUE
  arguments$name <- nameByLanguage(arguments$name,arguments$language,arguments$nodes)
  M<-pathParameter(model,estimates=estimates)
  names(M$nodes)[1]<-arguments$name
  
  
  if(!is.null(arguments$note)) arguments$note<-paste0(catFit(model,fitMeasures),arguments$note)
  else if(!is.null(fitMeasures)) arguments$note<-catFit(model,fitMeasures)
  
  if("nodes" %in% names(arguments)) {
    vvnodes<-setdiff(names(arguments$nodes),arguments$name)
    arguments$nodes<-merge(arguments$nodes,M$nodes,by.x=arguments$name,by.y=arguments$name,sort=F)
    arguments$nodes<-arguments$nodes[,c(names(M$nodes),vvnodes)]
  }
  else arguments$nodes<-M$nodes
  arguments$links<-M$links
  
  do.call(netCoin,arguments)
}

# produces a netCoin object to graph a CoWeb graphic (Upton 2000).
## cobCoin ----
cobCoin<- function(data, variables=names(data), degree=0, significance=.05, ...) {
  arguments <- list(...)
  if(!exists("color", arguments)) arguments$color<-"variable"
  if(!exists("group", arguments)) arguments$group<-"variable"
  if(!exists("lwidth", arguments)) arguments$lwidth<-"Haberman"
  if(!exists("lcolor", arguments)) arguments$lcolor<-"Haberman"
  if(!exists("linkBipolar", arguments)) arguments$linkBipolar <- TRUE
  arguments$groupText<-TRUE
  arguments$data <- data[,variables]
  arguments$maxL <- 1
  arguments$commonlabel <- ""
  N <- do.call(surCoin, arguments)
  control <- N$links[,1:2]
  control[[1]] <- gsub(":.*","",control[[1]])
  control[[2]] <- gsub(":.*","",control[[2]])
  N$links <- N$links[control[[1]]!=control[[2]],]
  N$links <- N$links[N$links$`p(Z)` < significance/2 | N$links$`p(Z)` > (1-significance/2),]
  N <- addNetCoin(N, layout=layoutCircle(N$nodes, deg=degree, name=N$options$nodeName))
  return(N)
}

#subFunctions ----
pathParameter<-function(model,estimates=c("b","se","z","pvalue","beta")){
  if(inherits(model,"lavaan")){
    links<-lavaan::parameterEstimates(model,standardized = T)
    names(links)<-gsub("^est$","b",names(links))
    names(links)<-gsub("^std.all$","beta",names(links))
    links<-links[links$op=="~",c("rhs","lhs",estimates)]
    names(links)[1:2]<-c("Source","Target")
    if(length(intersect(union(links$Source,links$Target),model@Data@ov$name))>0) {
      nodes<-as.data.frame(model@Data@ov,stringsAsFactors=F)[,intersect(names(model@Data@ov),c("name","mean","var"))]
      nodes$stdev<-sqrt(nodes$var)
      row.names(nodes)<-nodes$name
    }
    else {
      nodes<-data.frame(name=union(links$Source,links$Target))
    }
    # nodes<-nodes[,c("name","mean","stdev")]
    nodes$name<-iconv(nodes$name,"","UTF-8")
    links<-links[,c("Source","Target",estimates)]
    links$Source<-iconv(links$Source,"","UTF-8")
    links$Target<-iconv(links$Target,"","UTF-8")  
    structure(list(links = links, nodes = nodes))
  }
  else stop("Model has to be a lavaan object")
}

catFit<-function(model,fitMeasures){
  if("chisq" %in% fitMeasures) fitMeasures<-union(fitMeasures,c("df","pvalue"))
  fit<-lavaan::fitMeasures(model,fitMeasures)
  text<-NULL
  if("chisq" %in% names(fit)) text<-paste(text, paste0("Chi2=",format(fit["chisq"],digits=4)," (",fit["df"]," df)",", pvalue=", format(fit["pvalue"],digits=2)), sep=". ")
  if("cfi" %in% names(fit)) text<-paste(text, paste0("CFI= ",format(fit["cfi"],digits=3)), sep=". ") 
  if("rmsea" %in% names(fit)) text<-paste(text, paste0("RMSEA= ",format(fit["rmsea"],digits=3)), sep=". ")
  paste0(gsub("^. ","<p>",text),".</p>")
}

