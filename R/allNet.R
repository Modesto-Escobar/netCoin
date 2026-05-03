# Complete netCoin from an incidences matrix

allNet<-function(incidences, weight = NULL, subsample = FALSE, pairwise = FALSE,
                 minimum=1, maximum = nrow(incidences), sort = FALSE, decreasing = TRUE,
                 frequency = FALSE, percentages = TRUE, 
                 procedures = "Haberman", criteria = "Z", Bonferroni = FALSE,
                 support = -Inf, minL = -Inf, maxL = Inf,
                 directed = FALSE, diagonal = FALSE, sortL = NULL, decreasingL = TRUE,
                 igraph = FALSE, dir=NULL, ...)
{
  arguments <- list(...)
  arguments$dir<-dir
  if((criteria=="Z" | criteria=="hyp") & maxL==Inf) maxL=.5
  if(!("language" %in% names(arguments))) arguments$language <- "en"
  arguments$name <- nameByLanguage(arguments$name,arguments$language,arguments$nodes)
  if (!("size" %in% names(arguments)))
    if(percentages)
      arguments$size <- "%"
  if (!("level" %in% names(arguments))) level<-.95 else level <-arguments$level
  if (!pairwise) incidences<-na.omit(incidences)
  if (inherits(weight, "character")) incidences <- incidences[, setdiff(names(incidences), weight)]
  incidences <- incidences[,colSums(incidences)>0]
  if (all(is.na(incidences) | incidences==0 | incidences==1)) {
    C<-coin(incidences, minimum, maximum, sort, decreasing, weight=weight, subsample=subsample, pairwise = pairwise)
    if(exists("size",arguments))if(arguments$size=="frequency")frequency=TRUE
    O<-asNodes(C,frequency,percentages,arguments$language)
    names(O)[1]<-arguments$name
    if (is.null(arguments$nodes)){
      if(any(sapply(incidences,function(X) {"label" %in% names(attributes(X))}))) {
        label <- "label"
        if(arguments$language %in% c("es","ca")){
          label <- "etiqueta"
        }
        O[[label]] <- "NULL"
        O[[label]] <- ifelse(sapply(incidences, attr, "label")=="NULL", O[[arguments$name]], sapply(incidences, attr, "label"))
        arguments$label <- label
      }
      arguments$nodes<-O
    }else{
      nodesOrder<-as.character(arguments$nodes[[arguments$name]])
      arguments$nodes<-merge(O,arguments$nodes[,setdiff(names(arguments$nodes),frequencyList),drop=FALSE],by.x=arguments$name,by.y=arguments$name,all.y=TRUE, sort=FALSE)
      row.names(arguments$nodes)<-arguments$nodes[[arguments$name]]
      arguments$nodes<-arguments$nodes[nodesOrder,]
    }
    procedures<-union(procedures,unlist(arguments[c("lwidth","lweight","lcolor","ltext")]))
    arguments$links<-edgeList(C, procedures, criteria, level, Bonferroni, minL, maxL, support, 
                              directed, diagonal, sortL, decreasingL)
    for(lattr in c("lwidth","lweight","lcolor","ltext"))
      if(!is.null(arguments[[lattr]])) arguments[[lattr]]<-i.method(c_method(arguments[[lattr]]))
    if(is.character(arguments$layout)){
      if(tolower(substr(arguments$layout,1,2))=="mc")arguments$layout<-layoutMCA(incidences)
      else if(tolower(substr(arguments$layout,1,2))=="pc")arguments$layout<-layoutPCA(C)
    }
    arguments$scenarios <- attr(C,"n")
    xNx <- do.call(netCoin,arguments)
    if (igraph) return(toIgraph(xNx))
    else return(xNx)
  }
  else warning("Input is not a dichotomous matrix of incidences")
}

