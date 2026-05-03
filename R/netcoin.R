## Programs to apply net coin analysis
# Image is a files vector with length and order equal to nrow(nodes). Place as nodes field
# Batch

## netCoin ----
netCoin <- function(nodes = NULL, links = NULL, tree = NULL,
        community = NULL, layout = NULL,
        name = NULL, label = NULL, group = NULL, groupText = FALSE,
        labelSize = NULL, size = NULL, color = NULL, shape = NULL,
        border = NULL, legend = NULL, sort = NULL, decreasing = FALSE,
        ntext = NULL, info = NULL, image = NULL, imageNames = NULL,
        centrality = NULL,
        nodeBipolar = FALSE, nodeScaleLimits = NULL, nodeFilter = NULL, degreeFilter = NULL,
        lwidth = NULL, lweight = NULL, lcolor = NULL, ltext = NULL,
        intensity = NULL, linkBipolar = FALSE, linkScaleLimits = NULL, linkFilter = NULL,
        repulsion = 25, distance = 10, zoom = 1,
        fixed = showCoordinates, limits = NULL,
        main = NULL, note = NULL, showCoordinates = FALSE, showArrows = FALSE,
        showLegend = TRUE, frequencies = FALSE, showAxes = FALSE,
        axesLabels = NULL, scenarios = NULL, help = NULL, helpOn = FALSE,
        mode = c("network","heatmap"), roundedItems = FALSE, controls = 1:8,
        cex = 1, background = NULL, defaultColor = "#1f77b4",
        language = c("en","es","ca"), dir = NULL)
{
  if(is.null(links) &&  is.null(nodes)){
    stop("You must explicit a nodes or links data frame.")
  }

  if(inherits(nodes, 'netCoin')){
    stop("Using netCoin function to change netCoin object attributes is deprecated. Use addNetCoin instead.")
  }

  name <- nameByLanguage(name,language,nodes)
  if(!is.null(nodes)){
      if (all(inherits(nodes,c("tbl_df","tbl","data.frame"),TRUE))) nodes<-as.data.frame(nodes) # convert haven objects
  }
  if(!is.null(links)){
      if (all(inherits(links,c("tbl_df","tbl","data.frame"),TRUE))) links<-as.data.frame(links) # convert haven objects
  }

  color <- setAttrByValueKey("color",color,nodes)
  shape <- setAttrByValueKey("shape",shape,nodes)
  lcolor <- setAttrByValueKey("lcolor",lcolor,links)

  net <- network_rd3(nodes = nodes, links = links, tree = tree,
        community = community, layout = layout,
        name = name, label = label, group = group, groupText = groupText,
        labelSize = labelSize, size = size, color = color, shape = shape,
        border = border, legend = legend,
        sort = sort, decreasing = decreasing, ntext = ntext, info = info,
        image = image, imageNames = imageNames,
        nodeBipolar = nodeBipolar, nodeScaleLimits = nodeScaleLimits, nodeFilter = nodeFilter, degreeFilter = degreeFilter,
        source = "Source", target = "Target",
        lwidth = lwidth, lweight = lweight, lcolor = lcolor, ltext = ltext,
        intensity = intensity, linkBipolar = linkBipolar, linkScaleLimits = linkScaleLimits, linkFilter = linkFilter,
        repulsion = repulsion, distance = distance, zoom = zoom,
        fixed = fixed, limits = limits,
        main = main, note = note, showCoordinates = showCoordinates, showArrows = showArrows,
        showLegend = showLegend, frequencies = frequencies, showAxes = showAxes,
        axesLabels = axesLabels, scenarios = scenarios, help = help, helpOn = helpOn,
        mode = mode, roundedItems = roundedItems, controls = controls, cex = cex,
        background = background, defaultColor = defaultColor,
        language = language, dir = dir)
  class(net) <- c("netCoin",class(net))

  if(!is.null(centrality)){
    columns <- calCentr(net, centrality)$nodes
    for(col in setdiff(colnames(columns),c("nodes","degree"))){
      net$nodes[[col]] <- columns[[col]]
    }
  }

  return(net)
}

summary.netCoin <- function(object, ...){
  summaryNet(object)
}


setAttrByValueKey <- function(name,item,items){
    if(is.list(item) && !is.data.frame(item)){
      checkedlist <- list()
      for(k in names(item)){
        if(!k %in% colnames(items) || !(is.character(items[[k]]) || is.factor(items[[k]]))){
          warning(paste0(name,": the names in the list must match character columns of the items, but '",k,"' doesn't"))
        }else{
          if(!is.character(item[[k]]) || is.null(names(item[[k]]))){
            warning(paste0(name,": each item in the list must be a named character vector describing value-",name,", but '",k,"' doesn't"))
          }else{
            checkedlist[[k]] <- unname(item[[k]][items[[k]]])
          }
        }
      }
      if(length(checkedlist)){
        item <- as.data.frame(checkedlist)
      }else{
        item <- NULL
      }
    }
    return(item)
}

##addNetCoin ----
addNetCoin <- function(x, ...){
    arguments <- list(...)

    for(n in c("nodes","links","tree")){
      if(!(n %in% names(arguments))){
        arguments[[n]] <- x[[n]]
      }
    }

    options <- x$options

    getOpt <- function(opt,item=opt){
      if(item %in% names(arguments)){
        return(arguments[[item]])
      }else{
        if(!is.null(options[[opt]])){
          return(options[[opt]])
        }else{
          return(NULL)
        }
      }
    }

    attributes <- c(
      "name" = "nodeName",
      "cex" = "cex",
      "distance" = "distance",
      "repulsion" = "repulsion",
      "zoom" = "zoom",
      "scenarios" = "scenarios",
      "limits" = "limits",
      "main" = "main",
      "note" = "note",
      "help" = "help",
      "background" = "background",
      "language" = "language",
      "nodeBipolar" = "nodeBipolar",
      "linkBipolar" = "linkBipolar",
      "helpOn" = "helpOn",
      "frequencies" = "frequencies",
      "defaultColor" = "defaultColor",
      "controls" = "controls",
      "mode" = "mode",
      "axesLabels" = "axesLabels",
      "fixed" = "fixed",
      "showCoordinates" = "showCoordinates",
      "showArrows" = "showArrows",
      "showLegend" = "showLegend",
      "showAxes" = "showAxes",
      "roundedItems" = "roundedItems",
      "label" = "nodeLabel",
      "labelSize" = "nodeLabelSize",
      "group" = "nodeGroup",
      "groupText" = "groupText",
      "size" = "nodeSize",
      "color" = "nodeColor",
      "shape" = "nodeShape",
      "border" = "nodeBorder",
      "legend" = "nodeLegend",
      "ntext" = "nodeText",
      "info" = "nodeInfo",
      "sort" = "nodeOrder",
      "decreasing" = "decreasing",
      "image" = "imageItems",
      "imageNames" = "imageNames",
      "lwidth" = "linkWidth",
      "lweight" = "linkWeight",
      "lcolor" = "linkColor",
      "ltext" = "linkText",
      "intensity" = "linkIntensity"
    )

    for(item in names(attributes)){
      arguments[[item]] <- getOpt(attributes[[item]],item)
    }

    return(do.call(netCoin,arguments))
}

##savePajek ----
savePajek<-function(net, file="file.net", arcs=NULL, edges=NULL, partitions= NULL, vectors=NULL){
  if(length(setdiff(partitions,names(net[["nodes"]])))>0) stop("At least one partition is not amongst ",paste(names(net$nodes),collapse=", "),".")
  if(length(setdiff(vectors,names(net[["nodes"]])))>0) stop("At least one vector is not amongst ",paste(names(net$nodes),collapse=", "),".")
  if(length(setdiff(arcs,names(net[["links"]])))>0) stop("At least one arc is not amongst ",paste(names(net$links),collapse=", "),".")
  if(length(setdiff(edges,names(net[["links"]])))>0) stop("At least one edge is not amongst ",paste(names(net$links),collapse=", "),".")
  
  if(!grepl("\\.",file))file<-paste0(file,".net")
  if(!is.null(vectors) | !is.null(partitions)) file<-gsub(".net",".paj",file)
  connec<-file(file,"w")
  writeLines(paste0("*Network ",net[["options"]]$main),con=connec)
  close(connec)
  connec<-file(file,"a")
  writeLines(paste0("*Vertices ",as.character(nrow(net[["nodes"]]))),con=connec)
  writeLines(paste0(seq(1:nrow(net[["nodes"]])),' "',net[["nodes"]]$name,'" '), con=connec)
  N<-cbind(n=seq(1:nrow(net[["nodes"]])),net[["nodes"]][1])
  L<-cbind(N[unlist(net[["links"]]$Source),1],N[unlist(net[["links"]]$Target),1])
  
  if(!is.null(arcs)) {
    cont=1
    for(weights in arcs) {
      writeLines(paste0("*Arcs : ",cont,' "',weights,'"'), con=connec)
      writeLines(paste(L[,1],L[,2],net[["links"]][[weights]]), con=connec)
      cont=cont+1
    }
  }
  if(!is.null(edges)) {
    ifelse(exists("cont"),cont<-cont,cont<-1)
    for(weights in edges) {
      writeLines(paste0("*Edges : ",cont,' "',weights,'"'), con=connec)
      writeLines(paste(L[,1],L[,2],net[["links"]][[weights]]), con=connec)
      cont=cont+1
    }
  }
  if(!is.null(partitions)){
    for(partition in partitions) {
      writeLines(paste0("*Partition ", partition), con=connec)
      writeLines(paste0("*Vertices ", nrow(net$nodes)), con=connec)
      writeLines(as.character(as.numeric(as.factor(net[["nodes"]][[partition]]))),con=connec)
    }
  }
  if(!is.null(vectors)){
    for(vector in vectors) {
      writeLines(paste0("*Vector ", vector), con=connec)
      writeLines(paste0("*Vertices ", nrow(net$nodes)), con=connec)
      Line<-as.character(net[["nodes"]][[vector]])
      Line[is.na(Line)]<-"0"
      writeLines(Line, con=connec)
    }   
  }
  close(connec)
}

##saveGhml----
saveGhml <- function(net, file="netCoin.graphml"){
  if(!inherits(net, "netCoin")) stop("This program only works with netCoin objects")
  if(!grepl("\\.",file))file<-paste0(file,".graphml")
  graph <- toIgraph(net)
  write_graph(graph, file=file, format="graphml")
}

##shinyCoin----
shinyCoin <- function(x){
  shiny_rd3(x)
}








