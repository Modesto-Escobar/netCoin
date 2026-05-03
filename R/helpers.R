## Modesto Escobar
## shared helpers of netCoin

languages <- c("en","es","ca")

nameList <- c('name','nombre','nom')
names(nameList) <- languages

frequencyList <- c("frequency","frecuencia","freq\uFC\uE8ncia","%")
names(frequencyList) <- languages

labelList <- c('label','etiqueta','etiqueta')
names(labelList) <- languages

## allnet, coin, layout, netcoin-core, netcorr, path, regression, surcoin
getByLanguage <- function(varlist,language){
  if(!is.null(language) && language[1] %in% names(varlist))
    language <- language[1]
  else
    language <- "en"
  return(unname(varlist[language]))
}

## allnet, netcoin-core, netcorr, path, regression, surcoin
nameByLanguage <- function(name,language,nodes){
  if(is.null(name)){
    name <- getByLanguage(nameList,language)
  }
  if(!is.null(nodes)){
    if(!(name %in% colnames(nodes)))
      warning(paste0("name: '",name,"' column missing in nodes data frame"))
    else if(sum(duplicated(nodes[[name]])))
      warning(paste0("name: '",name,"' column values must be unique"))
  }
  return(name)
}

### allnet, coin, edgelist-engine, netcorr, surcoin
c_method<-function(method) {
  if(is.null(method))return(NULL)
  method<-toupper(method)
  if ("ALL"==method[1]) method<-c("M","T","G","S","B","J","D","A","O","K","N","Y","P","V","C","R","E","L","H","Z","f","F","X","I","0","Q","U")
  method<-sub("FIS","W",method)
  method<-sub("HYP","W",method)
  method<-sub("HAM","NNH",method)
  method<-sub("CC" ,"UC" ,method)
  method<-sub("CP" ,"QC" ,method)  
  method<-sub("OD" ,"CD" ,method)
  method<-sub("RO" ,"TA" ,method)
  method<-sub("AND","BER",method)
  method<-sub("TET","VET",method)
  method<-sub("CON","LCO",method)
  method<-sub("II" ,"0"  ,method)
  method<-sub("COI","F", method)
  method<-substr(method,1,1)
  return(method)
}

## allnet, coin, edgelist-engine, netcorr, surcoin
i.method<-function(method) {
  similarities<-matrix(c("matching","Rogers","Gower","Sneath", "Anderberg",
                         "Jaccard","dice", "antiDice","Ochiai","Kulczynski",
                         "Hamann", "Yule", "Pearson", "odds", "Russell", "expected", "Haberman", "confidence", "Z",
                         "coincidences", "relative", "sConditional","tConditional", "c.conditional","c.probable","tetrachoric","Fisher"), 
                       nrow=1, dimnames=list("Similarity", c("M","T","G","S","B","J","D","A","O","K","N","Y","P","C","R","E","H","L","Z","F","X","I","0","U","Q","V","W")))
  similarities<-similarities[,method]
  if("L" %in% method) {
    cade <- c("conf.L","confidence","conf.U")
    posi <- match("confidence",similarities)
    if (posi ==1) similarities <- c(cade,similarities[-1])
    else if (posi == length(similarities)) similarities<-c(similarities[-length(similarities)], cade)
    else similarities <- c(similarities[c(1:max(1,posi-1))], cade,
                           similarities[(match("confidence",similarities)+1):length(similarities)]) 
  }
  return(similarities)
}

## allnet, coin, edgelist-engine, netcorr, surcoin
checkLevel <- function(level){
  if (level >=1 & level < 100)
    level <- level/100
  if (level <=0 | level >=100) {
    level <- .95
    warning("Not valid level")
  }
  return(level)
}

## allnet, surcoin
layoutMCA<-function(matrix) { # Correspondencias simples clasicas aplicadas a dicotomicas.
  matrix<-cbind(matrix,1-matrix)
  n<-sum(matrix)
  P=matrix/n
  column.masses<-colSums(P)
  row.masses=rowSums(P)
  E=row.masses %o% column.masses
  R=P-E
  I=R/E
  Z=R/sqrt(E) # Corrected
  SVD=svd(Z)
  rownames(SVD$v)=colnames(P)
  standard.coordinates.columns = sweep(SVD$v[1:(ncol(Z)/2),1:2], 1, sqrt(column.masses[1:(ncol(Z)/2)]), "/")
  principal.coordinates.columns = sweep(standard.coordinates.columns, 2, SVD$d[1:2], "*")
  colnames(principal.coordinates.columns)<-c("F1","F2")
  return(principal.coordinates.columns)
}


## allnet, netcorr, surcoin
layoutPCA<-function(coin) { # Coordenadas a partir de Pearson: Haberman/raiz(n)
  if(inherits(coin, "coin")) A<-eigen(sim(coin,"P"))
  else A <- eigen(coin)
  C<-sweep(A$vectors[,1:2],2,sqrt(A$values[1:2]),"*")
  rownames(C)<-rownames(coin)
  colnames(C)<-c("F1","F2")
  return(C)
}

# For d_netCorr Factorial coordinates.
layoutFact <- function(data, method=c("pc","nf","vf", "of")) {
  difference <- setdiff(method, c("null", "pc","nf","vf","of"))
  if(length(difference)>0) warning(paste0(difference, " method is not implemented. "))
  R <- list()
  if(method[1]=="null") return(NULL)
  if("pc" %in% method) R$`Principal Components`  <- princomp(na.omit(data), cor=T)$loadings[,1:2]
  if("nf" %in% method) R$`Non-rotated Factorial` <- factanal(na.omit(data), factors=2, rotation="none")$loadings[,1:2]
  if("vf" %in% method) R$`Varimax factorial`     <- factanal(na.omit(data), factors=2, rotation="varimax")$loadings[,1:2]
  if("of" %in% method) R$`Oblimin factorial`     <- factanal(na.omit(data), factors=2, rotation="oblimin")$loadings[,1:2]
  if(length(R)>1) return(R) else return(R[[1]])
}

## allnet, coin, edgelist-engine, netcorr, surcoin
mats2edges<-function(data,list=NULL,criteria=1,min=-Inf,max=Inf,support=-Inf,directed=FALSE,diagonal=FALSE){
  # Input control
  if (!is.null(list)) {
    if (!identical(dim(data),dim(list[[1]]))) {
      warning("data & list must have the same dimensions")
      return()
    }
    if (criteria!=1 & !(criteria %in% names(list))) warning("criteria are not in the matrices list")
  }
  if (is.null(rownames(data))) rownames(data)<-as.character(1:nrow(data))
  if (is.null(colnames(data))) colnames(data)<-as.character(1:ncol(data))
  
  #  type of edgelist
  if (nrow(data)<ncol(data)) { # For asymmetric matrices. Improve later for directed cases.
    l<-data>-Inf
    l[,rownames(l)]<-lower.tri(l[,rownames(l)])
  }
  else { # For symmetric matrices
    if (directed) l <- as.vector(lower.tri(data,diag=diagonal) | upper.tri(data))
    else l <- as.vector(lower.tri(data,diag=diagonal))
  }
  # data.frame building
  sources<-rep(colnames(data),each=dim(data)[1])
  targets<-rep(rownames(data),dim(data)[2])
  Mat <- data.frame(Source=sources,Target=targets)
  value <- as.vector(data)
  
  # data.frame for matrices list
  if (!is.null(list)) {
    c <- as.vector(list[[criteria]])
    a<-as.data.frame(lapply(list,as.vector))
    Mat<-cbind(Mat,a)[l==TRUE & c<=max & c>=min & value>=support,]
  }
  # data.frame for alone matrix
  else Mat<-{
    if (criteria!=1) warning("Criteria don't apply for only one matrix. Use min and max")
    cbind(Mat,value)[l==TRUE & value<=max & value>=min,]
  }
  # return
  if (nrow(Mat)>0) {
    row.names(Mat)<-NULL
    return(Mat)
  }
  else return(NULL)
}

rescale <- function(x) {
  to <- c(0, 1)
  from <- range(x, na.rm = TRUE, finite = TRUE)
  return((x - from[1]) / diff(from) * diff(to) + to[1])
}

toColorScale <- function(items){
  if(is.numeric(items)){
    return(hsv(1,1,rescale(items)))
  }else{
    colors <- c(
      "#1f77b4", # blue
      "#2ca02c", # green
      "#d62728", # red
      "#9467bd", # purple
      "#ff7f0e", # orange
      "#8c564b", # brown
      "#e377c2", # pink
      "#7f7f7f", # grey
      "#bcbd22", # lime
      "#17becf", # cyan
      "#aec7e8", # light blue
      "#98df8a", # light green
      "#ff9896", # light red
      "#c5b0d5", # light purple
      "#ffbb78", # light orange
      "#c49c94", # light brown
      "#f7b6d2", # light pink
      "#c7c7c7", # light grey
      "#dbdb8d", # light lime
      "#9edae5" # light cyan
    )
    items <- as.numeric(as.factor(items))
    items <- ((items-1) %% length(colors))+1
    return(colors[items])
  }
}
