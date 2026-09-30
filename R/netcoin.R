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
        showLegend = TRUE, frequencies = FALSE, statistics = FALSE, showAxes = FALSE,
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

  if(!is.null(nodes)){
    nodes <- as.data.frame(nodes)
  }
  if(!is.null(links)){
    links <- as.data.frame(links)
  }
  name <- nameByLanguage(name,language,nodes)

  color <- setAttrByValueKey("color",color,nodes)
  shape <- setAttrByValueKey("shape",shape,nodes)
  lcolor <- setAttrByValueKey("lcolor",lcolor,links)
  layout <- layoutCompute(layout, nodes, links, name, lweight)

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
        showLegend = showLegend, frequencies = frequencies, statistics = statistics, showAxes = showAxes,
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

## layoutCompute ----
# rD3plot computes a single layout from its name, but draws several planes when it is given
# them as a list of coordinate matrices. So several names, as in layout=c("fo","fr","ka"), or
# a list mixing names and matrices, are turned here into that list: each name is computed by
# rD3plot itself on a bare network of the same nodes and links, so that "fo" weighs the links
# by lweight as it would on its own, and the nodes come out in the order the final network
# keeps. A single name or a single matrix goes through untouched.
# Each plane is named as it was given, the layout name when it had none.
layoutCompute <- function(layout, nodes, links, name, lweight) {
  if(is.character(layout) && length(layout) > 1) layout <- as.list(layout)
  if(!is.list(layout) || is.data.frame(layout)) return(layout)
  nm <- names(layout)
  if(is.null(nm)) nm <- rep("", length(layout))
  for(i in seq_along(layout)) {
    x <- layout[[i]]
    if(is.character(x)) {
      if(length(x) != 1 || is.na(x))
        stop("each element of layout must be a single layout name or a coordinate matrix")
      if(nm[i] == "") nm[i] <- x
      # An unknown name makes rD3plot stop with a message about none of this
      net <- tryCatch(network_rd3(nodes = nodes, links = links, name = name,
                                  lweight = lweight, layout = x),
                      error = function(e) NULL)
      if(is.null(net$nodes$fx))
        stop("\"", x, "\" is not a layout netCoin can compute", call. = FALSE)
      layout[[i]] <- cbind(net$nodes$fx, net$nodes$fy)
    } else if(nm[i] == "") nm[i] <- paste0("layout", i)
    # rD3plot takes the axis labels of each plane from its column names, and the page stops
    # loading when a plane has none, so a plane without them gets blank ones
    if(is.matrix(layout[[i]]) && is.null(colnames(layout[[i]])))
      colnames(layout[[i]]) <- c("", "")
  }
  names(layout) <- nm
  layout
}

## layoutSpecial ----
# Replace in layout the names only the calling function can compute, such as "pc" for the
# principal components of the correlations in netCorr, by their coordinates. special is a
# list of functions keyed by the two-letter code, so that only the ones asked for are run.
# A single name comes back as a matrix, as it always did; among several, the matrix takes
# the place of the name and keeps it as the name of the plane, and the names left over are
# computed later by layoutCompute.
layoutSpecial <- function(layout, special) {
  if(!is.character(layout) && !(is.list(layout) && !is.data.frame(layout))) return(layout)
  single <- is.character(layout) && length(layout) == 1
  L <- as.list(layout)
  nm <- names(L)
  if(is.null(nm)) nm <- rep("", length(L))
  for(i in seq_along(L)) {
    x <- L[[i]]
    if(is.character(x) && length(x) == 1 && !is.na(x)) {
      code <- tolower(substr(x, 1, 2))
      if(code %in% names(special)) {
        if(nm[i] == "") nm[i] <- x
        L[[i]] <- special[[code]]()
      }
    }
  }
  if(single) return(L[[1]])
  names(L) <- nm
  L
}

summary.netCoin <- function(object, ...){
  summaryNet(object)
}


## print.netCoin ----
# A netCoin object printed as-is (there being no method for it before this one) dumped its
# whole list structure: every row of the node table, the ntext column -one long repetitive HTML
# string per row, built for the tooltip of every node- the whole options list, and an
# "attr(,...)" footer of the bookkeeping attributes the object carries for itself
# (clusterColumns, caseToPattern...). None of it reads as a summary. This replaces that dump
# with one: how many nodes and links, which columns the node table holds (its ntext and
# clusterization columns left out, the latter named on their own line instead), which
# clusterizations \link{addClusters} added to \link{surScat}'s own, and which planes
# \link{addAxes} added to the one it drew. \code{summary} (see \link{summaryNet}) covers
# different ground: how nodes and links are distributed over their frequency/width column.
print.netCoin <- function(x, ...) {
  cat("<netCoin object>\n")

  clusterCols <- intersect(attr(x, "clusterColumns"), names(x$nodes))

  if(is.data.frame(x$nodes)) {
    shown <- setdiff(names(x$nodes), c(x$options$nodeText, clusterCols))
    cat("Nodes: ", nrow(x$nodes), "  (", paste(shown, collapse=", "), ")\n", sep="")
  }

  if(is.data.frame(x$links))
    cat("Links: ", nrow(x$links), "  (", paste(names(x$links), collapse=", "), ")\n", sep="")
  else
    cat("Links: 0\n")

  if(length(clusterCols))
    cat("Clusterizations: ", paste(clusterCols, collapse=", "), "\n", sep="")

  # currentLayouts stops on an object with no coordinates of its own, which a plain netCoin()
  # object (as opposed to one from surScat) usually is
  if(is.data.frame(x$nodes) && all(c("fx","fy") %in% names(x$nodes)))
    cat("Planes: ", paste(names(currentLayouts(x)), collapse=", "), "\n", sep="")

  invisible(x)
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
      "statistics" = "statistics",
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

    # An absent nodeLabel means the nodes are drawn with no label at all, which is what
    # surScat asks for through label="", whereas an absent label argument makes netCoin
    # fall back to the name column. The two NULLs stand for opposite things, and the loop
    # above turns the first into the second, since assigning NULL drops the argument and
    # lets the default back in. The emptied label is thus restated as such.
    if(!("label" %in% names(arguments)) && !is.null(options) && is.null(options[["nodeLabel"]]))
      arguments$label <- ""

    net <- do.call(netCoin,arguments)

    # netCoin builds the object anew out of its arguments, so whatever is not one of them
    # has to be carried over: the planes addAxes added, which live in $layouts, the cluster
    # centroids on each of them, which live in $clusters, and the attributes surScat attaches
    # to say how the nodes stand for the cases (caseToPattern, caseWeight, sampledNodes,
    # sampledFrom) or which columns hold clusterizations. Without this, a call meant to change
    # a colour silently cost the object its second plane, its cluster centroids, its ability
    # to collapse a case-level clusterization, and the guards maxN relies on.
    # A caller replacing the nodes themselves is taken at its word: none of that survives a
    # node table of another length, so it is not carried over.
    if(identical(nrow(net$nodes), nrow(x$nodes))) {
      if(is.null(arguments$layout)) {
        if(length(x$layouts) && !length(net$layouts))
          net$layouts <- x$layouts
        if(length(x$clusters) && !length(net$clusters))
          net$clusters <- x$clusters
      }
      for(a in setdiff(names(attributes(x)), names(attributes(net))))
        attr(net, a) <- attr(x, a)
    }

    return(net)
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
  # graphml only understands numeric, character or logical attributes; a factor column of
  # the node or link table (an ordered k-means group, a categorical variable...) survives
  # toIgraph as a factor, which write_graph then rejects with "Attribute not numeric.
  # Invalid value" instead of writing it out as the text it displays.
  for(a in igraph::vertex_attr_names(graph)) {
    v <- igraph::vertex_attr(graph, a)
    if(is.factor(v)) igraph::vertex_attr(graph, a) <- as.character(v)
  }
  for(a in igraph::edge_attr_names(graph)) {
    v <- igraph::edge_attr(graph, a)
    if(is.factor(v)) igraph::edge_attr(graph, a) <- as.character(v)
  }
  write_graph(graph, file=file, format="graphml")
}

##shinyCoin----
shinyCoin <- function(x){
  shiny_rd3(x)
}








