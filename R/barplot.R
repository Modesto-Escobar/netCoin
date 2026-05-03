barCoin<-function(data, variables = colnames(data), commonlabel = NULL,
        dichotomies = c("_all","_none"), valueDicho = 1, weight = NULL,
        subsample = FALSE, sort = NULL, decreasing = TRUE, nodes = NULL,
        name = NULL, select = NULL, scalebar = FALSE, note = NULL, label = NULL,
        text = NULL, color = NULL, defaultColor = "#1f77b4", expected = FALSE,
        confidence = FALSE, level = .95, significance = FALSE, minimum = 1 ,
        maximum = nrow(data), percentages = FALSE, criteria = c("Z","hyp"),
        Bonferroni = FALSE, support = 1, minL = -Inf, maxL = 1,
        language = c("en","es","ca"), cex = 1.0, dir = NULL)
{
  language <- language[1]
  name <- nameByLanguage(name = name, language =language, nodes = nodes)
  dicho<-function(input,variables,value) {
    datum<-as.data.frame(ifelse(input[,variables]==value,1,0))
    if (all(inherits(input,c("tbl_df","tbl","data.frame"),TRUE))) {
      # L<-sapply(datum[,variables],attr,"label")
      M<-sapply(sapply(datum,attr,"label"),function(X) ifelse(is.null(X),NA,X))
      L<-ifelse(is.na(M),variables,M)
      names(datum)<-L
    }
    return(datum)
  }

# model of bar
  if (confidence) procedures <- c("Frequencies","Expected","Confidence")
  else procedures <- c("Frequencies","Expected")
  criteria <- criteria[1]
  if(criteria!="hyp") criteria <- "Z"
  procedures<-c(procedures,criteria)
  criteriacolname <- c(Z="p(Z)",hyp="p(Fisher)")[criteria]
  

# names  
  if (inherits(nodes,"tbl_df")) nodes<-as.data.frame(nodes)


# classification of variables
  if("_all"  %in% dichotomies) dichotomies<-variables
  if("_none" %in% dichotomies) {
    allvar<-variables  
    dichotomies <- NULL
  }
    else allvar <- union(dichotomies,variables)
#  if(!is.null(nodes)) allvar<-intersect(unlist(nodes[name]),allvar)
  
# treatment of weight in the data frame. Omission of na data.
    if(!is.null(weight)) {
    if(inherits(weight,"character")){
      data<-na.omit(data[,allvar])
      allvar<-setdiff(allvar,weight)
      variables<-setdiff(variables,weight)
      weight<-data[,weight]
    }
    else{
      if(length(weight)!=dim(data)[1]) stop("Weights have not the correct dimensions")
      data<-na.omit(cbind(data[,allvar],weight))[,1:length(data[,allvar])]
    }
  }
  else data<-na.omit(data[,allvar])

# set dichotomies    
  if(length(dichotomies)>0) {
    dichos<-dicho(data,dichotomies,valueDicho)
    variables<-setdiff(variables,dichotomies)
  }

# data.frame setting
  data[,variables]<-as_factor(data[,variables])
  if (all(inherits(data,c("tbl_df","tbl","data.frame"),TRUE))) data<-as.data.frame(data) # convert haven objects
  
# dichotomizing factor variables
  if (length(variables)>0){
    incidences<-dichotomize(data, variables, "", min=minimum, length=0, values=NULL, sparse=FALSE, add=FALSE, sort=TRUE)
    if(exists("dichos")) incidences<-cbind(dichos,incidences)
  } 
  else if(exists("dichos")) incidences<-dichos
  
  
# nodes filter  
  if (!is.null(nodes)) {
    nonAmong<-setdiff(as.character(nodes[[name]]),names(incidences))
    nodeList<-setdiff(as.character(nodes[[name]]),nonAmong)
    incidences<-incidences[,nodeList]
    if(length(nonAmong)>0)
      warning(paste0(toString(nonAmong)," is/are not present among incidences."))
  }
  
# coincidences elaboration
  if(!exists("incidences")) stop("There are no qualitative variables. Try netCorr.") 
  incidences<-na.omit(incidences)
  if(!is.null(nodes)) incidences <- incidences[,intersect(unlist(nodes[name]),colnames(incidences))]
  if (all(incidences==0 | incidences==1)) {
    C<-coin(incidences, minimum, maximum, sort=TRUE, decreasing=TRUE, weight=weight, subsample=subsample)

# nodes data.frame elaboration
    O<-asNodes(C, !percentages, percentages, language = language) # Attention to !percentages
    if(name!="name") names(O)[1]<-name
    names(O)[2] <- "incidences"
    if(!is.null(nodes)) {
      O<-merge(O,nodes[,setdiff(names(nodes),frequencyList),drop=FALSE],by.x=name,by.y=name,all.x=TRUE)
    }else {
      if (!is.null(commonlabel)) { # Preserve the prename (variable) of a node if specified in commonlabel
        label<-getByLanguage(labelList,language)
        provlabels<-as.character(O[[name]])
        O[[label]]<-ifelse(substr(O[[name]],1,regexpr('\\:',O[[name]])-1) %in% commonlabel,provlabels,substr(O[[name]],regexpr('\\:',O[[name]])+1,1000000L))
      }
    }

# making edgeList
    level <- checkLevel(level)
    E<-edgeList(C, procedures, criteria, level, Bonferroni, minL, maxL, support, 
                directed=FALSE, diagonal= FALSE, sort= NULL)

# definition of parameters
    if (expected || confidence){
      expected <- "expected"
    }else{
      expected <- NULL
    }
    if (confidence){
      if(significance){
        confidence <- "confidence"
        E[,c("conf.L","conf.U")] <- NULL
      }else{
        confidence <- c("conf.L","conf.U")
        E[,"confidence"] <- NULL
      }
    }else{
      confidence <- NULL
    }
    if(significance){
      significance <- criteriacolname
    }else{
      significance <- NULL
      E[,criteriacolname] <- NULL
    }

# convertion to percentages
    if (percentages) {
      E[,intersect(names(E),c("coincidences","expected","confidence","conf.L","conf.U"))]<-
      E[,intersect(names(E),c("coincidences","expected","confidence","conf.L","conf.U"))]/attr(C,"n")*100
    }
        
# preparing bar graph
    bar <- barplot_rd3(O, E, name = name, select = select,
        source = "Source", target = "Target",
        label = label, text = text, color = color,
        incidences = "incidences", coincidences = "coincidences",
        expected = expected, confidence = confidence, level = level, significance = significance,
        sort = sort, decreasing = decreasing,
        scalebar = scalebar, defaultColor = defaultColor, note = note, cex = cex,
        language = language, dir = dir)
    class(bar) <- c("barCoin",class(bar))
    return(bar)
  }
  else warning("Input is not a dichotomous matrix of incidences")
}

## Method ----
summary.barCoin <- function(object, ...){
  summaryNet(object)
}
