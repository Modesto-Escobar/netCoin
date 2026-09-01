caring_read_file <- function(filepath){
  ext <- gsub("^.*\\.","",filepath)

  if(ext=="sav" || ext=="dta"){
    if(ext=="sav")
      data <- haven::read_spss(file=filepath)
    if(ext=="dta")
      data <- haven::read_stata(file=filepath)
    data <- as.data.frame(data)
    colnames(data) <- make.unique(sapply(colnames(data),function(col){
      label <- attr(data[[col]],"label")
      if(!is.null(label) && length(label)==1){
        return(label)
      }
      return(col)
    }), sep = "_")
    for(col in colnames(data)){
      if(inherits(data[[col]],"haven_labelled")){
        data[[col]] <- factor(haven::as_factor(data[[col]]))
      }
    }
  }else if(ext=="xlsx"){
    data <- openxlsx::read.xlsx(filepath)
    for(col in colnames(data))
      if(inherits(data[[col]],"character"))
        data[[col]] <- factor(data[[col]])
  }else{
    enc <- as.data.frame(readr::guess_encoding(filepath))[1,1]
    if(enc!="UTF-8")
      enc <- "Latin-1"
    data <- data.table::fread(file=filepath,data.table=FALSE,encoding=enc)
    if(enc!="UTF-8")
      colnames(data) <- enc2utf8(colnames(data))
    for(i in seq_along(data)){
      if(is.character(data[,i])){
        if(enc!="UTF-8")
          data[,i] <- enc2utf8(data[,i])
        data[,i] <- as.factor(data[,i])
      }
    }
  }

  colnames(data) <- gsub('"',"'",colnames(data))

  for(col in colnames(data)){
    if(is.factor(data[[col]])){
      if(all(!is.na(suppressWarnings(as.numeric(levels(data[[col]])))))){
        data[[col]] <- as.numeric(data[[col]])
      }
    }
  }

  return(data)
}

caring_create_graphs <- function(data, arguments){
  initialvariables <- arguments[['variables']]

  if(!is.null(arguments[['weight']])){
    arguments[['weight']] <- arguments[['weight']][1]
    if((!arguments[['weight']] %in% colnames(data)) || !is.numeric(data[,arguments[['weight']]])){
      arguments[['weight']] <- NULL
    }else{
      arguments[['variables']] <- setdiff(arguments[['variables']],arguments[['weight']])
    }
  }

  pmax <- NULL
  if(is.null(arguments[['maxL']]) || !is.finite(arguments[['maxL']])){
    if(is.null(arguments[['criteria']]) || (arguments[['criteria']]=="Z" || arguments[['criteria']]=="hyp")){
      arguments[['maxL']] <- 0.5
    }
  }else{
    pmax <- arguments[['maxL']]
  }

  if(is.null(arguments[['dichotomies']])){
    dichotomies <- character(0)
    for(col in arguments[['variables']]){
      if(is.numeric(data[,col])){
        val <- unique(data[,col])
        if(!length(setdiff(val,c(0,1,NA))))
          dichotomies <- c(dichotomies,col)
      }
    }
    if(length(dichotomies)){
      arguments[['dichotomies']] <- dichotomies
      arguments[['valueDicho']] <- 1
    }
  }

  variables <- arguments[['variables']]
  if(!is.null(arguments[['dichotomies']])){
    variables <- setdiff(variables,arguments[['dichotomies']])
    if(!is.null(arguments[['metric']])){
      arguments[['metric']] <- setdiff(arguments[['metric']],arguments[['dichotomies']])
    }
  }
  if(!is.null(arguments[['metric']])){
    arguments[['variables']] <- setdiff(arguments[['variables']],arguments[['metric']])
  }

  arguments[['frequency']] <- TRUE
  arguments[['percentage']] <- TRUE
  arguments[['color']] <- "variable"
  arguments[['data']] <- data

  plots <- arguments[['plot']]
  arguments[['plot']] <- NULL

  multiArgs <- list(dir = arguments[['dir']])
  arguments[['dir']] <- NULL

  if("surCoin" %in% plots){
    if(!is.null(arguments[['metric']]) && length(setdiff(arguments[['variables']],arguments[['metric']]))==0){
      corrArgs <- arguments
      corrArgs[['variables']] <- arguments[['data']][,arguments[['metric']]]
      keys <- union(names(formals(netCorr)),names(formals(netCoin)))
      corrArgs <- corrArgs[intersect(names(corrArgs),keys)]
      net1 <- do.call(netCorr,corrArgs)
    }else{
      keys <- union(names(formals(surCoin)),names(formals(netCoin)))
      surArgs <- arguments[intersect(names(arguments),keys)]
      net1 <- do.call(surCoin,surArgs)
    }
    multiArgs[[plots[which(plots=="surCoin")+1]]] = net1
  }
  if("surScat" %in% plots){
    scatArgs <- arguments[intersect(names(arguments),union(formalArgs("surScat"),formalArgs("netCoin")))]
    if(!is.null(scatArgs[['patterns']])){
      scatArgs[['vPatterns']] <- scatArgs[['patterns']]
      scatArgs[['patterns']] <- TRUE
    }
    if(!length(scatArgs[['nclusters']])){
      scatArgs[['nclusters']] <- min(nrow(data)-1,6):2
    }
    if(length(scatArgs[['type']]) && scatArgs[['type']]==""){
      scatArgs[['type']] <- NULL
    }
    scatArgs[['degreeFilter']] <- NULL
    net2 <- do.call(surScat,scatArgs)
    multiArgs[[plots[which(plots=="surScat")+1]]] = net2
  }
  logs <- intersect(c("logCoin2","logCoin3","logCoin4","logCoin5"),plots)
  if(length(logs)){
    for(log in logs){
      logArgs <- arguments[intersect(names(arguments),union(formalArgs("logCoin"),formalArgs("netCoin")))]
      if(!is.null(arguments[['dichotomies']])){
        noFirstCat <- character(0)
        for(i in seq_along(arguments[['dichotomies']])){
          col <- arguments[['dichotomies']][i]
          valueDicho <- arguments[['valueDicho']][i]
          datacol <- as.character(logArgs[['data']][,col])
          logArgs[['data']][,col] <- factor(datacol,levels=union(valueDicho,setdiff(unique(datacol),valueDicho)))
          noFirstCat <- c(noFirstCat,col)
        }
        if(length(noFirstCat)){
          logArgs[['noFirstCat']] <- noFirstCat
        }
      }
      if(!is.null(pmax)){
        logArgs[['pmax']] <- pmax
      }
      logArgs[['color']] <- "var"
      logArgs[['order']] <- as.numeric(sub("logCoin","",log,fixed=TRUE))
      net3 <- do.call(logCoin,logArgs)
      multiArgs[[plots[which(plots==log)+1]]] = net3
    }
  }
  if("allNet" %in% plots){
    inciArgs <- arguments[intersect(names(arguments),union(formalArgs("allNet"),formalArgs("netCoin")))]
    inciArgs[['incidences']] <- data[,initialvariables]
    if(!inciArgs[['color']] %in% colnames(inciArgs[['incidences']])){
      inciArgs[['color']] <- NULL
    }
    for(n in colnames(inciArgs[['incidences']])){
      inciArgs[['incidences']][,n] <- as.numeric(as.logical(inciArgs[['incidences']][,n]))
    }
    net4 <- do.call(allNet,inciArgs)
    multiArgs[[plots[which(plots=="allNet")+1]]] = net4
  }
  if("glmCoin" %in% plots){
    exogenous <- arguments[['exogenous']]
    if(!length(exogenous)){
      stop("missing independent variables")
    }

    glmArgs <- arguments[intersect(names(arguments),union(formalArgs("glmCoin"),formalArgs("netCoin")))]
    glmArgs[['data']] <- data[,initialvariables]

    if(!is.null(pmax)){
      glmArgs[['pmax']] <- pmax
    }

    chaine <- arguments[['chaine']]
    family <- arguments[['family']]

    metric <- arguments[['metric']]
    dichotomies <- arguments[['dichotomies']]
    textvariables <- setdiff(initialvariables,unique(c(metric,dichotomies,exogenous)))
    finalvariables <- initialvariables

    if(!is.null(dichotomies)){
      for(i in seq_along(dichotomies)){
        dic <- dichotomies[i]
        value <- arguments[['valueDicho']][i]
        newvar <- gsub(" ","",paste0(dic,":",value))
        glmArgs[['data']][,newvar] <- ifelse(glmArgs[['data']][,dic]==value, 1, 0)
        glmArgs[['data']][,dic] <- NULL
        finalvariables[initialvariables==dic] <- newvar
      }
    }

    initialchaine <- chaine
    if(!is.null(textvariables)){
      textfamilies <- vapply(textvariables,function(dic){
        return(family[initialvariables==dic])
      },character(1))
      for(i in seq_along(textvariables)){
        dic <- textvariables[i]
        values <- unique(glmArgs[['data']][,dic])
        for(value in values){
          newvar <- gsub(" ","",paste0(dic,":",value))
          glmArgs[['data']][,newvar] <- ifelse(glmArgs[['data']][,dic]==value, 1, 0)
          chaine <- c(chaine,initialchaine[initialvariables==dic])
          family <- c(family,textfamilies[i])
          finalvariables <- c(finalvariables,newvar)
        }
        keep <- finalvariables!=dic
        chaine <- chaine[keep]
        family <- family[keep]
        finalvariables <- finalvariables[keep]
        glmArgs[['data']][,dic] <- NULL
      }
    }

    chaine[finalvariables %in% exogenous] <- 0
    family[chaine==0] <- NA

    # data.frame para la elaboración de la fórmula (ecuación)
    A <- data.frame(v=finalvariables,
                n=chaine,
                m=family,
                stringsAsFactors = FALSE)

    A <- A[order(A$n, as.numeric(rownames(A)), decreasing=TRUE),]
    formulas <- ""

    # Bucle sobre cada fila del dataframe
    for (i in seq_len(nrow(A))) {
      subvariables <- A$v[A$n < A$n[i]]
      if (length(subvariables) > 0) {
        ecuacion <- paste(A$v[i], "~", paste(subvariables, collapse = "+"), ",", A$m[i])
        if (formulas != "") {
          formulas <- paste(formulas, ecuacion, sep = "\n")
        } else {
          formulas <- ecuacion
        }
      }
    }

    glmArgs[['formulas']] <- formulas
    glmArgs[['note']] <- gsub("\n","</br>",formulas)

    net5 <- do.call(glmCoin,glmArgs)
    multiArgs[[plots[which(plots=="glmCoin")+1]]] = net5
  }

  return(do.call(multigraphCreate,multiArgs))
}

caring <- function(filepath,arguments){

  data <- caring_read_file(filepath)

  graphs <- caring_create_graphs(data,arguments)

  return(graphs)
}
