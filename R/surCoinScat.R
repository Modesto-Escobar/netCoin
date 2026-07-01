# surCoin is a wrapper to build a netCoin object from an original non-dichotomized data.frame. See below dichotomize()
## surCoin ----
surCoin<-function(data,variables=names(data), commonlabel=NULL, 
                  dichotomies=NULL, valueDicho=1, metric=NULL, exogenous=NULL,
                  weight=NULL, subsample=FALSE, pairwise=FALSE,
                  minimum=1, maximum=nrow(data), sort=FALSE, decreasing=TRUE,
                  frequency=FALSE, percentages=TRUE,
                  procedures="Haberman", criteria="Z", Bonferroni=FALSE,
                  support=-Inf, minL=-Inf, maxL=Inf,
                  directed=FALSE, diagonal=FALSE, sortL=NULL, decreasingL=TRUE,
                  igraph=FALSE, coin=FALSE, dir=NULL, ...)
{
  arguments <- list(...)
  if((criteria=="Z" | criteria=="hyp") & maxL==Inf) maxL=.5
  varOrder  <- variables # To order variables later before coin
  #Check methods. No necessary because edgeList call these routines.
  #procedures<-i.method(c_method(procedures))
  #criteria<-i.method(c_method(criteria))
  procedures<-union(procedures,unlist(arguments[c("lwidth","lweight","lcolor","ltext")]))
  
  #Metrics 
  if(!is.null(metric)) {
    variables<-setdiff(variables,metric)
    procedures<-intersect(union(procedures,"Z"),c("Pearson","Haberman","Z"))
    criteria<-intersect(union(criteria,"Z"),c("Pearson","Haberman","Z"))[1]
  }
  
  #Names  
  if(!("language" %in% names(arguments))) arguments$language <- "en"
  nodes <- arguments$nodes
  if (inherits(nodes,"tbl_df")) nodes<-as.data.frame(nodes)
  name <- arguments$name <- nameByLanguage(arguments$name,arguments$language,arguments$nodes)
  if (!("level" %in% names(arguments))) level<-.95 else level <-arguments$level
  
  #Data.frame  
  if (all(inherits(data,c("tbl_df","tbl","data.frame"),TRUE))) data<-as.data.frame(data) # convert haven objects
  if (inherits(weight,"character")) variables <- setdiff(variables,weight)
  
  #ignore constant dichotomies
  if(!is.null(dichotomies)){
    cons <- character()
    for(i in dichotomies){
      if(length(unique(data[,i]))==1){
        cons <- c(cons,i)
      }
    }
    dichotomies <- setdiff(dichotomies,cons)
    if(!length(dichotomies)){
      dichotomies <- NULL
    }
  }
  
  allvar<-union(union(metric,dichotomies),variables)
  
  if (!pairwise & inherits(weight,"character")) {
    if (!is.null(weight)) weight <- data[rowSums(is.na(data[,allvar]))<1,weight]
    data <- data[complete.cases(data[,allvar]),]
  }
  
  if(!is.null(weight)) {
    if(inherits(weight,"character")){
      weight<-data[,weight]
      data<-data[,allvar]
    }
    else{
      if(length(weight)!=dim(data)[1]) stop("Weights have not the correct dimensions")
      if (pairwise) data <- cbind(data[,allvar],weight)[,1:length(data[,allvar])]
      else data <- na.omit(cbind(data[,allvar],weight))[,1:length(data[,allvar])]
    }
  }
  else data<-data[,allvar]
  
  data[,variables]<-as_factor(data[,variables])
  arguments$scenarios <- sum(rowSums(!is.na(data))>0) # Number of scenarios
  
  #Size 
  if(!("size" %in% names(arguments)))
    if(percentages)
      arguments$size <- "%"
  
  #Dichotomies    
  if(!is.null(dichotomies)){
    dichos<-dicho(data, dichotomies, valueDicho, newlabel = FALSE)
    variables<-setdiff(variables,dichotomies)
  }
  
  #Dichotomize
  if (length(variables)>0){
    incidences<-dichotomize(data, variables, "", min=minimum, length=0, values=NULL, sparse=FALSE, add=FALSE, sort=sort, nas=NULL)
    if(!is.null(dichotomies)) incidences<-cbind(dichos,incidences)
  } 
  else if(exists("dichos")) incidences<-dichos
  
  #Nodes filter  
  if (!is.null(nodes)) {
    nonAmong  <-setdiff(as.character(nodes[[name]]),names(incidences))
    nodeList  <-setdiff(as.character(nodes[[name]]),nonAmong)
    incidences<-incidences[,nodeList, drop=FALSE]
    nonAmongM <- setdiff(metric, as.character(nodes[[name]]))
    metric    <- setdiff(metric,nonAmongM)
    if (length(metric)==0) metric <-NULL
    if (length(nonAmong)>0)
      warning(paste0(toString(nonAmong)," is/are not present in the data frame."))
    # nodes <- nodes[as.character(nodes[[name]]) %in% union(names(incidences),metric),]
  }
  
  #Nodes elaboration
  if (!exists("incidences") | ncol(incidences)<2) stop("There are no more than 1 qualitative category. Try netCorr.")
  if (all(is.na(incidences) | incidences==0 | incidences==1)) {
    incidences <- incidences[,names(incidences)[order(match(sub(":.*","",names(incidences)),varOrder))], drop=FALSE]
    C<-coin(incidences, minimum, maximum, sort, decreasing, weight=weight, subsample=subsample, pairwise = pairwise)
    if(coin) return(C)
    O<-asNodes(C,frequency,percentages,arguments$language)
    names(O)[1]<-name
    O$variable <- sub(":.*","",O[,name])
    if(!is.null(nodes)) {
      O<-merge(O,nodes[,setdiff(names(nodes),frequencyList),drop=FALSE],by.x=name,by.y=name,all.x=TRUE)
    }else {
      if (!is.null(commonlabel)) { # Preserve the prename (variable) of a node if specified in commonlabel
        arguments$label<-getByLanguage(labelList,arguments$language)
        provlabels<-as.character(O[[name]])
        O[[arguments$label]]<-ifelse(substr(O[[name]],1,regexpr('\\:',O[[name]])-1) %in% commonlabel,provlabels,substr(O[[name]],regexpr('\\:',O[[name]])+1,1000000L))
      }
    }
    
    #Links elaboration
    E<-edgeList(C, procedures, criteria, level, Bonferroni, minL, maxL, support, 
                directed, diagonal, sortL, decreasingL, pairwise)
    for(lattr in c("lwidth","lweight","lcolor","ltext"))
      if(!is.null(arguments[[lattr]])) arguments[[lattr]]<-i.method(c_method(arguments[[lattr]]))
    
    if(!is.null(arguments$layout)) {
      layout2 <- layout <- arguments$layout
      if (inherits(layout,"matrix") && is.null(metric)){
        if (!is.null(nodes)){
          if(nrow(layout)==nrow(nodes)){
            Oxy <- matrix(NA,nrow(O),2)
            rownames(Oxy) <- as.character(O[,name])
            rownames(layout) <- as.character(nodes[,name])
            layoutnames <- intersect(rownames(Oxy),rownames(layout))
            Oxy[layoutnames,] <- layout[layoutnames,]
            arguments$layout <- Oxy
          } else warning("layout must have a coordinate per node")
        } else warning("layout must be applied to the nodes variable")
      } else {
        if(is.character(layout)){
          if(tolower(substr(layout,1,2))=="mc")arguments$layout<-layoutMCA(incidences)
          else if(tolower(substr(layout,1,2))=="pc")arguments$layout<-layoutPCA(C)
        }
        else if(!is.null(metric)) arguments$layout<-NULL # There is metric information and not MCA or PCA
      }
    }
    
    if(!is.null(metric)) {
      #Metric nodes elaboration
      if(percentages) O$mean<-O$`%`/100
      O$min<-0
      O$max<-1
      means<-sapply(na.omit(data[,metric, drop=F]),mean)
      mins<-sapply(na.omit(data[,metric, drop=F]),min)
      maxs<-sapply(na.omit(data[,metric, drop=F]),max)
      P<-(means-mins)/(maxs-mins)*100
      O2<-data.frame(name=names(means),mean=means,min=mins,max=maxs,P=P,variable=names(means))
      colnames(O2)[1] <- name
      colnames(O2)[5] <- "%"
      if(!is.null(nodes)){
        O2 <- O2[as.character(O2[[name]]) %in% as.character(nodes[[name]]),] #nodes filter 
        for(col in as.character(O2[[name]]))
          O2[as.character(O2[[name]])==col,colnames(nodes)] <- nodes[as.character(nodes[[name]])==col,]
      }
      O<-rbind_all_columns(O,O2)
      
      #Metric links elaboration
      methods<-union(procedures,criteria)
      if (pairwise) R <- corrp(data[,metric, drop=F], cbind(incidences, data[,metric, drop=F]), weight=weight)
      else R <- corr(data[,metric, drop=F], cbind(incidences,data[,metric, drop=F]), weight=weight)      
      if (nrow(R)==1) row.names(R)<-metric
      allvar<-union(as.character(nodes[[name]]),c(names(incidences),metric))
      order1<-intersect(allvar,rownames(R))
      order2<-intersect(allvar,colnames(R))
      R<-R[order1,order2, drop=F]
      Pearson<-mats2edges(R)
      colnames(Pearson)[3]<-"Pearson"
      D<-as.data.frame(Pearson)
      if("Haberman" %in% methods){
        H<-R*sqrt(nrow(incidences))
        Haberman<-mats2edges(H)
        D<-cbind(D,Haberman[,3]); colnames(D)[length(D)]<-"Haberman"
      }
      if("Z" %in% methods) {
        t<-R/sqrt((1-pmin(1L,R))/(nrow(incidences)-2))
        Z <- mats2edges(1-pt(t,nrow(incidences)-2))
        D<-cbind(D,Z[,3]); colnames(D)[length(D)]<-"Z"
      }
      D<-D[,c("Source","Target",methods)]
      D<-D[D[criteria] > minL & D[criteria] < maxL,]
      colnames(D)<-sub("^Z$","p(Z)",colnames(D))
      if(is.null(E))E<-D
      else E<-rbind_all_columns(E,D)
      
      #Layout
      if (inherits(layout,"matrix")){
        if (!is.null(nodes)){
          if(nrow(layout2)==nrow(nodes)){
            Oxy <- matrix(NA,nrow(O),2)
            rownames(Oxy) <- as.character(O[,name])
            rownames(layout2) <- as.character(nodes[,name])
            layoutnames <- intersect(rownames(Oxy),rownames(layout2))
            Oxy[layoutnames,] <- layout2[layoutnames,]
            layout2 <- Oxy
          } else warning("layout must have a coordinate per node")
        } else warning("layout must be applied to the nodes variable")
        arguments$layout <- layout2
      }
    }
    if (!is.null(exogenous)) {
      exogenous2<-intersect(exogenous,c(metric,dichotomies))
      E$chaine<-ifelse(((substr(E$Source,1,regexpr("\\:",E$Source)-1) %in% exogenous) |
                          (E$Source %in% exogenous2))  &
                         ((substr(E$Target,1,regexpr("\\:",E$Target)-1) %in% exogenous) |
                            (E$Target %in% exogenous2)),"No","Yes")
      arguments$linkFilter<-paste(ifelse(is.null(arguments$linkFilter),"",paste(arguments$linkFilter,"&")),"chaine=='Yes'")
    }
    if("showArrows" %in% names(arguments$options) & exists("nodes")) E<-orderEdges(E,nodes[[name]])
    if(exists("ltext",arguments)) {
      if(toupper(arguments$ltext) == "Z") arguments$ltext <- "p(Z)"
      if(arguments$ltext =="Fisher") arguments$ltext <- "p(Fisher)"
    }
    
    if(!is.null(dir)){
      arguments$dir <- dir
    }
    arguments$nodes <- O
    arguments$links <- E
    xNx <- do.call(netCoin,arguments)
    if (igraph) {
      return(toIgraph(xNx))
    } else {
      return(xNx)
    }
  } else warning("Input is not a dichotomous matrix of incidences")
}

## surScat ----
# surScat is a wrapper to build a netCoin object from an original non-dichotomized data.frame and see frequencies.
surScat <- function(data, variables=names(data), active=variables, type= c("mca", "pca"), nclusters=2, maxN=2000, ...) {
  if(type[1]=="mca") {
    B <- as.data.frame(droplevels(as_factor(na.omit(data[,variables]))))
    b <- B[,active]
    m <- as.matrix(dichotomize(b,variables=names(b), sort=F, add=F, nas=NULL))
    cc <- layoutMca(m, rows=T)
  }
  else {
    B <- na.omit(data[,variables])
    b <- as.data.frame(sapply(B[,active], as.numeric))
    factors <- setdiff(variables, active)
    B[, factors] <- as.data.frame(droplevels(as_factor(B[,factors])))
    B[, active]  <- b
    m  <- prcomp(b, center = TRUE, scale. = TRUE)
    cc <- m$x[,1:2]
  }
  for(i in nclusters) {
    G <- stats::kmeans(cc, centers=i)
    g <- paste0("Grupos(",sprintf(paste0("%0",length(levels(G$cluster)),"d"),i),")")
    B[[g]] <- paste0("Grupo: ",sprintf(paste0("%0",length(levels(G$cluster)),"d"),G$cluster))
  }
  arguments <- list(...)
  arguments$name <- nameByLanguage(NULL,arguments$language,NULL)
  if(is.character(rownames(B)))  B[[arguments$name]] <- rownames(B)
  else B[[arguments$name]] <- sprintf(paste0("%0",nchar(nrow(B)),"d"),as.numeric(rownames(B)))
  B <- B[, c(active, setdiff(names(B), active))]
  if(nrow(B)>maxN) {
    set.seed(2020)
    rcases <- sample(1:nrow(B), maxN)
    B  <- B[rcases,]
    cc <- cc[rcases,]
  }
  arguments$nodes <- B
  arguments$layout <- cc
  arguments$color <- g
  arguments$frequencies <- TRUE
  arguments$showAxes <- TRUE
  arguments$showCoordinates <- TRUE
  arguments$statistics <- TRUE
  arguments$degreeFilter <- NULL
  if(is.null(arguments$label)) arguments$label <- ""
  if(is.null(arguments$controls)) arguments$controls <- c(1,2,4)  
  return(do.call(netCoin, arguments))
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


## There is also a layoutMCA for other graphs
layoutMca<-function(matrix, nfactors=2, rows=FALSE){ # Correspondencias simples clasicas aplicadas a dicotómicas.
  P=matrix/nrow(matrix)
  column.masses<-colSums(P)
  row.masses=rowSums(P)
  E=row.masses %o% column.masses
  R=P-E
  I=R/E
  Z=R/sqrt(E) # Corrected
  SVD=svd(Z)
  rownames(SVD$v)=colnames(P)
  CC <-list()
  CC$standard.coordinates.rows = sweep(SVD$u, 1, sqrt(row.masses), "/")
  CC$principal.coordinates.rows = sweep(CC$standard.coordinates.rows, 2, SVD$d, "*")
  CC$standard.coordinates.columns = sweep(SVD$v, 1, sqrt(column.masses), "/")
  CC$principal.coordinates.columns = sweep(CC$standard.coordinates.columns, 2, SVD$d, "*")
  CC <- lapply(CC, function(X){X<-X[,2:(nfactors+1)];colnames(X) <- paste0("F",1:nfactors); return(X)})
  if(rows) return(CC$principal.coordinates.rows)
  else return(CC$principal.coordinates.columns)
}

dicho<-function(input,variables,value,newlabel=TRUE) {
  if(length(value)!=length(variables)){
    value <- rep(value[1],length(variables))
  }
  
  for(i in seq_along(variables)){
    dicho <- variables[i]
    vector0 <- rep(0,nrow(input))
    vector0[as.character(input[,dicho])==value[i]] <- 1
    input[,dicho] <- vector0
  }
  
  datum <- as.data.frame(input[, variables, drop=FALSE])
  j=0
  for (i in variables) {
    j=j+1
    if (!is.null(attributes(input[[i]]))) {
      if (newlabel) names(datum)[j] <- ifelse(exists("label",attributes(input[[i]])),attr(input[[i]],"label"),i)
    }
  }
  return(datum)
}


orderEdges<-function(links,nodes){ #Used in surCoin to order arrows
  A<-unlist(sapply(paste0("^",links[,"Source"],"$"),grep,x=nodes))
  B<-unlist(sapply(paste0("^",links[,"Target"],"$"),grep,x=nodes))
  links[A>B,c("Source","Target")]<-links[A>B,c("Target","Source")]
  links<-links[!is.na(links$Source) & !is.na(links$Target),]
  return(links)
}

rbind_all_columns <- function(x, y) {
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  
  x[, c(as.character(y.diff))] <- NA
  y[, c(as.character(x.diff))] <- NA
  
  return(rbind(x, y))
}
