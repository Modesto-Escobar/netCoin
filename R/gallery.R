gallery <- gallery_rd3

gallery2 <- function(...){
  warning("Using 'gallery2' function is deprecated. Use 'exhibit' instead.")
  do.call(gallery2_rd3,list(...))
}
exhibit <- gallery2_rd3

gallery3 <- gallery3_rd3

netGalleryWrapper <- function(treeGallery_function, tree, deep, initialType, tableformat, ...){
  args <- list(...)
  args$tree <- tree
  if(deep){
    args$deep <- TRUE
  }
  args$initialType <- initialType
  args$tableformat <- tableformat
  dir <- NULL
  if(!is.null(args$dir)){
    dir <- args$dir
    args$dir <- NULL
  }
  net <- do.call(treeGallery_function,args)
  if(!is.null(net$options$nodeTypes) && !is.null(net$options$nodeText)){
    net$nodes[[net$options$nodeText]] <- sapply(seq_len(nrow(net$nodes)),function(i){
      ntext <- net$nodes[i,net$options$nodeText]
      if(!is.na(ntext) && grepl('class="info-template"',ntext,fixed=TRUE)){
        name <- net$nodes[i,net$options$nodeName]
        relatives <- list()
        for(type in net$options$nodeTypes){
          aux <- character(0)
          if(length(net$nodes_relatives)){
            aux <- net$nodes_relatives[[i]][net$nodes_relativesTypes2[[i]]==type] 
          }else{
            if(name %in% net$tree[,1]){
              aux <- net$tree[net$tree[,1]==name & net$tree[,4]==type,2]
            }
            if(name %in% net$tree[,2]){
              aux <- net$tree[net$tree[,2]==name & net$tree[,3]==type,1]
              if(!length(aux)){
                parent <- net$tree[net$tree[,2]==name,1]
                for(p in parent){
                  aux <- c(aux,net$tree[net$tree[,1]==p & net$tree[,4]==type,2])
                }
              }
            }
          }
          if(length(aux)){
            aux <- aux[!is.na(aux)]
            aux <- unique(aux)
            aux <- aux[aux!=name]
            if(length(aux)){
              label <- aux
              if(!is.null(net$options$nodeLabel) && (net$options$nodeLabel!=net$options$nodeName)){
                labels <- net$nodes[[net$options$nodeLabel]]
                names(labels) <- net$nodes[[net$options$nodeName]]
                label <- labels[aux]
              }
              indices <- seq_len(nrow(net$nodes))-1
              names(indices) <- net$nodes[[net$options$nodeName]]
              aux <- indices[aux]
              aux <- paste0('<span nodename="',aux,'">',label,'</span>')
              aux <- aux[order(label)]
              relatives[[type]] <- paste0(aux,collapse="; ")
            }
          }
        }
        if(length(relatives)){
          relatives <- paste0(paste0('<b>',names(relatives),':</b> ',relatives),collapse="<br/>")
          ntext <- sub('<p class="template-text">',paste0('<div class="tree-relatives">',relatives,'</div><p class="template-text">'),ntext,fixed=TRUE)
        }
      }
      return(ntext)
    })
  }
  class(net) <- c("netGallery",class(net))
  if(!is.null(dir)){
    rD3plot:::galleryCreate(net,dir)
  }
  return(net)
}

netGallery <- function(tree, deep = FALSE, initialType = NULL, tableformat = FALSE, ...){
  return(netGalleryWrapper(treeGallery_rd3, tree, deep, initialType, tableformat, ...))
}

netGallery2 <- function(tree, initialType = NULL, tableformat = FALSE, ...){
  warning("Using 'netGallery2' function is deprecated. Use 'netExhibit' instead.")
  return(netExhibit(tree, initialType, tableformat, ...))
}

netExhibit <- function(tree, initialType = NULL, tableformat = FALSE, ...){
  return(netGalleryWrapper(treeGallery2_rd3, tree, FALSE, initialType, tableformat, ...))
}

asGallery <- function(net){
  if(inherits(net,"network_rd3")){
    options <- net$options
    gallery <- gallery_rd3(nodes = net$nodes, name = options$nodeName, label = options$nodeLabel,
      color = options$nodeColor, border = options$nodeBorder, ntext = options$nodeText, info = options$nodeInfo, image = options$imageItems,
      zoom = options$zoom, main = options$main, note = options$note,
      showLegend = options$showLegend, frequencies = options$frequencies,
      help = options$help, helpOn = options$helpOn, cex = options$cex, language = options$language)
    return(gallery)
  }else{
    stop("net: must be a network object")
  }
}
