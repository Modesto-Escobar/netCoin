## edgeList ----
# Links. See below funcs="shape"
edgeList <- function(data, procedures="Haberman", criteria="Z", level = .95, Bonferroni=FALSE, min=-Inf, max=Inf, support=-Inf, 
                     directed=FALSE, diagonal=FALSE, sort=NULL, decreasing=TRUE, pairwise=FALSE) {
  level <- checkLevel(level)
  if (tolower(substr(criteria,1,2))%in%c("z","hy") & substr(tolower(procedures[1]),1,2)!="sh") {
    if (max==Inf) max<-.50
    if (Bonferroni ) max<-max/choose(nrow(data),2) # Changes of Z max criterium (Bonferroni)
  }
  if (substr(tolower(procedures)[1],1,2)!="sh") { # For coin objects
    if (!inherits(data,"coin")) stop("Error: input must be a coin object (see coin function)")
    funcs<-c_method(procedures)
    if(!is.null(sort)) funcs<-union(c_method(sort),funcs)
    criteria<-c_method(criteria)
    todas<-union(funcs,criteria)
    matrices<-sim(data,todas,level=level, pairwise=pairwise)
    funcs<-i.method(funcs)
    criteria<-i.method(criteria)
    if (length(union(funcs,criteria))==1) {
      M<-new.env()
      M[[funcs]]<-matrices
      matrices<-as.list(M)
    }
    matrices<-matrices[i.method(todas)]
    Mat<-mats2edges(data[,],matrices,criteria=criteria,min=min,max=max,support=support,directed=directed,diagonal=diagonal)
  }
  else {
    if (!inherits(data,"matrix") && !inherits(data,"data.frame"))
      stop("Error: input must be a matrix (shape) or a data.frame (tree)")    
    if (inherits(data,"matrix")){
      if(min==-Inf)min<-1    
      #funcs="value"
      #M<-new.env()
      #M[[criteria]]<-M[[funcs]]<-data
      #matrices<-as.list(M)
      #data<-list(f=M[[funcs]],n=NA)
      Mat<-mats2edges(data,min=min,max=max,directed=directed,diagonal=diagonal)
    }
    if (inherits(data,"data.frame")) {
      lines<-sapply(data,as.character)
      lines<-rbind(c(lines[1,1],rep(NA,ncol(lines)-1)),lines) # Add one blank case in order to avoid mXm problem.
      lines<-ifelse(lines=="",NA,lines)
      adjlist<-split(lines,seq(nrow(lines))) # splits the character strings into list with different vector for each line
      adjlist<-sapply(adjlist,na.omit)
      Source=unlist(lapply(adjlist,function(x) rep(x[1],length(x)-1))) # establish first column of edgelist by replicating the 1st element (=ID number) by the length of the line minus 1 (itself)
      Target=unlist(lapply(adjlist,"[",-1)) # the second line I actually don't fully understand this command, but it takes the rest of the ID numbers in the character string and transposes it to list vertically
      return(as.data.frame(cbind(Source,Target),stringsAsFactors = FALSE,row.names=FALSE))
    }
  }  
  
  # Last transformations: c.Conditional c.Probable and sort
  
  if(length(Mat)>0) {
    if (!is.null(Mat$c.conditional)) 
      Mat$c.conditional<-factor(Mat$c.conditional,levels=c(0:8),
                                labels=c("Null","Mere","Conditional","Significant","Quite significant","Very significant","Subtotal","Suptotal","Total"))
    if (!is.null(Mat$c.probable)) 
      Mat$c.probable<-factor(Mat$c.probable,levels=c(0:8),
                             labels=c("Null","Mere","Probable","Significant","Quite significant","Very significant","Subtotal","Suptotal","Total"))
    if (!is.null(sort)) {
      if (substr(tolower(procedures)[1],1,2)!="sh") Mat<-Mat[order(Mat[[i.method(c_method(sort))]],decreasing = decreasing),]
      else Mat<-Mat[order(Mat$value,decreasing=decreasing),]
    }
  }
  else {
    warning("Check max and min values")
    return(NULL)
  }
  names(Mat)[names(Mat)=="Z"]<-"p(Z)"
  names(Mat)[names(Mat)=="Fisher"]<-"p(Fisher)"
  return(Mat)
}



## sim ----
# Similatities
sim<-function (input, procedures="Jaccard", level=.95, distance=FALSE, minimum=1, maximum=Inf, sort=FALSE, decreasing=FALSE, 
               weight=NULL, pairwise=FALSE) {
  level <- checkLevel(level)
  method<-c_method(procedures)
  if (is.matrix(input) && !inherits(input,"coin")) {
    if (is.null(colnames(input))) dimnames(input)<-list(NULL,paste("X",1:ncol(input),sep=""))
    input<-as.data.frame(input)
  }
  if (is.data.frame(input)) {
    C<-coin(input, minimum, maximum, sort, decreasing, weight=weight, pairwise=pairwise)
    a<-C[,]
  }
  else if (inherits(input,"coin")) {
    C <- input
    a <- input[(diag(input)>=minimum & diag(input)<=maximum),(diag(input)>=minimum &diag(input)<=maximum)]
  }
  else stop("Error: 1st parameter has to be a data frame or a coin object (see coin function)")
  
  if (pairwise) {
    N <- attr(C, "m")
    NN <- max(C)
    X <- attr(C, "x")
  }
  else {
    N <- attr(C,"n")
    X <- matrix(0, nrow=dim(a)[1], ncol = dim(a)[2])
  }
  
  if (sort==TRUE | decreasing==TRUE) {
    orderM <- order(diag(a),decreasing=decreasing)
    a <- a[orderM , orderM]
    if (pairwise) N <- N[orderM , orderM]
    if (pairwise) X <- X[orderM , orderM]
  }
  
  
  b<--(a-diag(a))-t(X)
  c<-t(b)
  d=N-a-b-c
  if (any(c("G","B","O","K","Y","P","V") %in% method)) m<-ifelse(a+d==N,1,ifelse(b+c==N,-1,0)) #Special values of distances
  s<-new.env()
  
  if ("M" %in% method) s$matching <- distant((a + d)/(a + b + c + d),distance)
  if ("T" %in% method) s$Rogers <- distant((a + d)/(a + 2 * (b + c) + d), distance)
  if ("G" %in% method) {
    s$Gower <- distant(a * d/sqrt((a + b) * (a + c) * (d + b) * (d + c)),distance)
    s$Gower <- ifelse(is.na(s$Gower),distant(pmax(m,0),distance),s$Gower)
  }
  if ("S" %in% method) s$Sneath <- distant(2*(a+d)/(2*(a+d)+(b+c)),distance)
  if ("B" %in% method) {
    s$Anderberg <- distant((a/(a+b)+a/(a+c)+d/(c+d)+d/(b+d))/4,distance)
    s$Anderberg <- ifelse(is.na(s$Anderberg),distant(pmax(m,0),distance),s$Anderberg)
  }
  if ("J" %in% method) s$Jaccard <- distant(a/(a + b + c),distance)
  if ("D" %in% method) s$dice <- distant(2 * a/(2 * a + b + c), distance)
  if ("A" %in% method) s$antiDice <- distant(a/(a + 2 * (b + c)),distance)
  if ("O" %in% method) {
    s$Ochiai <- distant(a/sqrt((a + b) * (a + c)),distance)
    s$Ochiai <- ifelse(is.na(s$Ochiai),distant(pmax(m,0),distance),s$Ochiai)    
  }
  if ("K" %in% method) {
    s$Kulczynski <- distant((a/(a+b)+a/(a+c))/2, distance)
    s$Kulczynski <- ifelse(is.na(s$Kulczynski),distant(pmax(m,0),distance),s$Kulczynski)
  }
  if ("N" %in% method) s$Hamann <- distant((a - (b + c) + d)/(a + b + c + d), distance)
  if ("Y" %in% method) {
    s$Yule <- distant((a*d-b*c)/(a*d+b*c))
    s$Yule <- ifelse(is.na(s$Yule),distant(m,distance),s$Yule)
  }
  if ("P" %in% method) {
    s$Pearson <- distant((a * d - b * c)/sqrt((a + b) * (a + c) * (b + d) *  (d + c)),distance)
    s$Pearson <- ifelse(is.na(s$Pearson),distant(m,distance),s$Pearson)
  }
  if ("V" %in% method) {
    s$tetrachoric<-((a*d/(b*c))^(pi/4)-1)/((a*d/(b*c))^(pi/4)+1)
    s$tetrachoric <- ifelse(is.na(s$tetrachoric),distant(m,distance),s$tetrachoric)
  }
  if ("C" %in% method) {
    s$odds <- (pmax(a,.5)*pmax(d,.5))/(pmax(b,.5)*pmax(c,.5))
    if (distance) s$odds<--s$odds
    diag(s$odds)<-ifelse(distance,-Inf,Inf)
  }
  if ("R" %in% method) {
    s$Russell <- distant(a/(a + b + c + d),distance)
    if (!distance) diag(s$Russell) <- 1
  }
  if ("E" %in% method) s$expected <- (a+b)*(a+c)/N
  if ("L" %in% method) {
    s$'conf.L' <- pmax(a-qt(level+(1-level)/2, N-1)*sqrt(((a+b)*(a+c)/N)*((1-(a+b)/N)*(1-(a+c)/N))),0)
    signo<-2*(((a+b)*(a+c)/N)<a)-1
    s$confidence <- pmax((a+b)*(a+c)/N+signo*qt(level,N-1)*sqrt(((a+b)*(a+c)/N)*((1-(a+b)/N)*(1-(a+c)/N))),0)
    diag(s$confidence) <- diag(a)
    s$'conf.U' <- pmin(a+qt(level+(1-level)/2, N-1)*sqrt(((a+b)*(a+c)/N)*((1-(a+b)/N)*(1-(a+c)/N))),N)
  }
  if ("H" %in% method) {
    s$Haberman <- sqrt(N) * (a * d - b * c)/sqrt((a + b) * (a + c) * (b + d) *  (d + c))
    if (pairwise) s$Haberman[is.na(s$Haberman)]<-sqrt(NN)
    else s$Haberman[is.na(s$Haberman)]<-sqrt(N)[is.na(s$Haberman)]
    if (distance) s$Haberman<-(N+s$Haberman)/(2*N)
  }
  if ("Z" %in% method) {
    s$Z <- 1-pt(sqrt(N) * (a * d - b * c)/sqrt((a + b) * (a + c) * (b + d) *  (d + c)),N)
    s$Z[is.na(s$Z)]<-0
  }
  if ("W" %in% method) s$Fisher<-1-phyper(a-1,pmin((a+b),(a+c)),N-pmin((a+b),(a+c)),pmax((a+b),(a+c)))
  if ("x" %in% method) {
    s$Fisher<-matrix(NA,nrow=nrow(a),ncol=ncol(a))
    for (Ro in c(1:nrow(a))) {
      for (Co in c(Ro:ncol(a))) {
        inMatrix=matrix(c(a[Ro,Co],b[Ro,Co],c[Ro,Co],d[Ro,Co]),nrow=2)
        (s$Fisher[Ro,Co]<-fisher.test(inMatrix,alternative="greater")$p.value)
        s$Fisher[Co,Ro] = s$Fisher[Ro,Co]
      }
    }
    # diag(s$Fisher)<-0
    rownames(s$Fisher)<-colnames(s$Fisher)<-rownames(a)
  }
  if ("F" %in% method) s$coincidences <- a
  if ("X" %in% method) s$relative <- a/N*100
  if ("I" %in% method) s$sConditional <-a/(a+c)*100
  if ("0" %in% method) s$tConditional <-a/(a+b)*100
  if ("U" %in% method) {
    Z <- 1-pt(sqrt(N) * (a * d - b * c)/sqrt((a + b) * (a + c) * (b + d) *  (d + c)),N)
    s$c.conditional<-matrix(ifelse(b+c==0, 8,
                                   ifelse(c==0,  7,                        
                                          ifelse(b==0,  6,
                                                 ifelse(Z<.001,5,
                                                        ifelse(Z<.01, 4,
                                                               ifelse(Z<.05, 3,
                                                                      ifelse(Z<.50, 2,       
                                                                             ifelse(a>0, 1, 0)))))))),nrow=nrow(a),dimnames=dimnames(a))
  }
  if ("Q" %in% method) {
    Z <- 1-pt((a/(a+c)-.50)/(1/(2*sqrt(a+c))),(a+c))
    s$c.probable<-matrix(ifelse(b+c==0, 8,
                                ifelse(c==0,   7,                        
                                       #                                  ifelse(b==0,  "Subtotal",
                                       ifelse(Z<.001, 5,
                                              ifelse(Z<.01,  4,
                                                     ifelse(Z<.05,  3,
                                                            ifelse(Z<.50,  2,      
                                                                   ifelse(a>0, 1, 0))))))),nrow=nrow(a),dimnames=dimnames(a))
  }  
  if (length(method)==1) return(s[[names(s)[1]]])
  else return(as.list(s,sorted=TRUE))
}

## expectedList ----
expectedList<- function(data, names=NULL, min=1, confidence=FALSE) {
  if (!inherits(data,"coin")) stop("Error: input must be a coin object")
  if (!is.null(names)) colnames(data[,])<-rownames(data[,])<-names
  a<-data[,]
  b<--(a-diag(a))
  c<--t((t(a)-diag(a)))
  d=attr(data,"n")-a-b-c
  attr(data,"e")<-(a+b)*(a+c)/(a+b+c+d)
  E<-edgeList(attr(data,"e"),"shape",min=0,max=Inf)
  F<-edgeList(data[,],"shape",min=0,max=Inf)
  if (!confidence) {
    dataL<-cbind(F,E[,3])[F[,3]>=min,]
    colnames(dataL)[3:4]<-c("coincidences","expected")
  }
  else {
    N<-a+b+c+d
    signo<-2*(((a+b)*(a+c)/N)<a)-1
    attr(data,"l") <- pmax((a+b)*(a+c)/N+signo*1.64*sqrt(((a+b)*(a+c)/N)*((1-(a+b)/N)*(1-(a+c)/N))),0)
    diag(attr(data,"l")) <- diag(a)
    L<-edgeList(attr(data,"l"),"shape",min=-Inf,max=Inf)
    dataL<-cbind(F,E[,3],L[,3])[F[,3]>=min,]
    colnames(dataL)[3:5]<-c("coincidences","expected","confidence")
  }
  return(dataL)
}


## subfunctions ----
# Convert similatities into dissimilarities
distant<-function(s,t=FALSE) {
  if (t==TRUE) s<-as.dist(1-s)
  return(s)
}
# http://pbil.univ-lyon1.fr/ade4/ade4-html/dist.binary.html


# Print lower matrices

lower<-function(matrix,decimals=3) { # Add an option to hiden diagonal
  m<-as.matrix(matrix)
  form=paste("%1.",decimals,"f",sep="")
  lower<-apply(m,1,function(x) sprintf(form,x))
  lower[upper.tri(lower)]<-""
  lower<-as.data.frame(lower, stringsAsFactors=FALSE)
  if (ncol(m)==1) rownames(lower)<-colnames(lower)<-names(matrix)
  rownames(lower)<-names(lower)
  return(lower)
}

