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
surScat <- function(data, variables=names(data), active=variables, weight=NULL, patterns=FALSE, jitter=0,
                     type=c("mca", "pca"), nclusters=2, clusterOn=c("factors", "variables"), nfactors=2, critFactors=0,
                     seed=2020, maxN=2000, ...) {
  clusterOn <- match.arg(clusterOn)
  if(nfactors<1) stop("nfactors must be at least 1")
  nExtract <- max(nfactors,2) # the scattergram always needs two dimensions, regardless of nfactors
  if(clusterOn=="variables" && (!missing(nfactors) || !missing(critFactors)))
    warning("nfactors and critFactors are ignored when clusterOn=\"variables\"")
  if(!is.null(seed)) {
    if(exists(".Random.seed", envir = globalenv())) {
      oldseed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", oldseed, envir = globalenv()), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }
  autoType <- missing(type)
  type <- match.arg(type)
  if(autoType && all(sapply(data[,active, drop=FALSE], is.numeric))) type <- "pca"
  if(inherits(weight,"character")) {
    if(!(weight %in% names(data))) stop("weight column not found in data")
    variables <- setdiff(variables, weight)
    active <- setdiff(active, weight)
    weight <- data[[weight]]
  }
  if(!is.null(weight) && length(weight)!=nrow(data)) stop("Weights have not the correct dimensions")
  keep <- complete.cases(data[,variables, drop=FALSE])
  if(!is.null(weight)) {
    keep <- keep & !is.na(weight)
    weight <- weight[keep]
  }
  D <- data[keep, variables, drop=FALSE]
  if(patterns) { # collapse cases into unique response patterns; weight = cases (or sum of weights) per pattern
    key <- do.call(paste, c(D, sep="\r"))
    counts <- rowsum(if(is.null(weight)) rep(1,length(key)) else weight, key)
    D <- D[!duplicated(key), , drop=FALSE]
    weight <- as.vector(counts[key[!duplicated(key)],])
    rownames(D) <- sprintf(paste0("%0",nchar(nrow(D)),"d"),seq_len(nrow(D)))
  }
  if(type=="mca") {
    B <- as.data.frame(droplevels(as_factor(D)))
    b <- B[,active, drop=FALSE]
    m <- as.matrix(dichotomize(b,variables=names(b), sort=F, add=F, nas=NULL))
    ff <- layoutMca(m, nfactors=nExtract, rows=T, weight=weight)
    vm <- m # dichotomized active variables, for clustering on clusterOn="variables"
  }
  else {
    if(length(active)<2) {
      if(autoType) stop("Only one quantitative active variable was found, so type=\"pca\" was chosen automatically. Include at least one more active variable, or set type=\"mca\" explicitly.")
      else stop("type=\"pca\" requires at least two active variables")
    }
    B <- D
    b <- as.data.frame(lapply(B[,active, drop=FALSE], as.numeric))
    factors <- setdiff(variables, active)
    B[, factors] <- as.data.frame(droplevels(as_factor(B[,factors, drop=FALSE])))
    B[, active]  <- b
    if(is.null(weight)) {
      m  <- prcomp(b, center = TRUE, scale. = TRUE)
      ff <- m$x[,1:nExtract, drop=FALSE]
      attr(ff,"eigenvalues") <- m$sdev^2 # full spectrum, for factor-selection criteria
      vm <- scale(as.matrix(b)) # standardized active variables, for clustering on clusterOn="variables"
    } else { # weighted PCA: weighted standardization and eigenvectors of the weighted correlation matrix
      w <- weight/sum(weight)
      ctr <- sweep(as.matrix(b), 2, colSums(as.matrix(b)*w))
      std <- sweep(ctr, 2, sqrt(colSums(ctr^2*w)), "/")
      eigres <- eigen(corr(b, weight=weight), symmetric=TRUE)
      ff <- (std %*% eigres$vectors)[,1:nExtract, drop=FALSE]
      attr(ff,"eigenvalues") <- eigres$values # full spectrum, for factor-selection criteria
      vm <- std
    }
  }
  # factors used for clustering: those among the first nfactors whose eigenvalue exceeds
  # critFactors times the mean eigenvalue of the full spectrum (generalizes Kaiser's rule to MCA)
  eig <- attr(ff, "eigenvalues")
  nUse <- max(1L, sum(eig[seq_len(nfactors)] > critFactors*mean(eig)))
  fm <- ff[, seq_len(nUse), drop=FALSE]
  cc <- ff[, 1:2, drop=FALSE] # the scattergram itself is always two-dimensional
  cm <- if(clusterOn=="variables") vm else fm
  arguments <- list(...)
  if(patterns) {
    vars <- names(B)
    nName <- make.unique(c(names(B),"n"))[ncol(B)+1]
    B[[nName]] <- weight
    if(is.null(arguments$size)) arguments$size <- nName
    txt <- rep("", nrow(B)) # HTML profile of each pattern: variable values plus number of cases
    for(v in vars) txt <- paste0(txt, "<b>", v, "</b>: ", as.character(B[[v]]), "<br/>")
    txt <- paste0(txt, "<b>", getByLanguage(casesList, arguments$language), "</b>: ", round(B[[nName]]))
    tName <- make.unique(c(names(B),"ntext"))[ncol(B)+1]
    B[[tName]] <- txt
    if(is.null(arguments$ntext)) arguments$ntext <- tName
  }
  groupsWord <- getByLanguage(groupsList, arguments$language)
  groupWord  <- getByLanguage(groupList,  arguments$language)
  for(i in nclusters) {
    G <- stats::kmeans(cm, centers=i)
    g <- paste0(groupsWord,"(",sprintf(paste0("%0",nchar(max(nclusters)),"d"),i),")")
    labels <- paste0(groupWord,": ",sprintf(paste0("%0",nchar(i),"d"),seq_len(i)))
    B[[g]] <- factor(labels[G$cluster], levels=labels, ordered=TRUE)
  }
  arguments$name <- nameByLanguage(NULL,arguments$language,NULL)
  if(is.character(rownames(B)))  B[[arguments$name]] <- rownames(B)
  else B[[arguments$name]] <- sprintf(paste0("%0",nchar(nrow(B)),"d"),as.numeric(rownames(B)))
  B <- B[, c(active, setdiff(names(B), active))]
  if(nrow(B)>maxN) {
    rcases <- sample(1:nrow(B), maxN)
    B  <- B[rcases,]
    cc <- cc[rcases,]
  }
  if(isTRUE(jitter)) jitter <- .02
  if(jitter>0) # visual spread only: clusters were computed on the exact coordinates
    cc <- cc + sweep(matrix(stats::runif(2*nrow(cc),-.5,.5), ncol=2), 2,
                     jitter*apply(cc, 2, function(z) diff(range(z))), "*")
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
layoutMca<-function(matrix, nfactors=2, rows=FALSE, weight=NULL){ # Correspondencias simples clasicas aplicadas a dicotómicas.
  if(is.null(weight)) P=matrix/nrow(matrix)
  else P=matrix*(weight/sum(weight)) # weighted row masses
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
  out <- if(rows) CC$principal.coordinates.rows else CC$principal.coordinates.columns
  attr(out, "eigenvalues") <- SVD$d[-1]^2 # full non-trivial spectrum (drops the null 1st dimension), for factor-selection criteria
  return(out)
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
