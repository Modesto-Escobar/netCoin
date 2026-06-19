\name{netCoin}
\alias{netCoin}
\title{Networked coincidences.}
\description{
\code{netCoin} produces a netCoin object of coincidences. Its input has to be two data.frames: one of attributes of events or nodes, and the other of attributes of the edges or links. 
}
\value{
This function returns a \code{netCoin} object.
If the 'dir' attribute is specified, the function creates a folder in the computer with an HTML document named index.html which contains the produced graph. This file can be directly opened with your browser and sent to a web server to work properly.
}

\usage{
netCoin(nodes = NULL, links = NULL, tree = NULL,
        community = NULL, layout = NULL,
        name = NULL, label = NULL, group = NULL, groupText = FALSE,
        labelSize = NULL, size = NULL, color = NULL, shape = NULL,
        border = NULL, legend = NULL, sort = NULL, decreasing = FALSE,
        ntext = NULL, info = NULL, image = NULL, imageNames = NULL,
        centrality = NULL,
        nodeBipolar = FALSE, nodeScaleLimits = NULL, nodeFilter = NULL,
        degreeFilter = NULL,
        lwidth = NULL, lweight = NULL, lcolor = NULL, ltext = NULL,
        intensity = NULL, linkBipolar = FALSE, linkScaleLimits = NULL,
        linkFilter = NULL,
        repulsion = 25, distance = 10, zoom = 1,
        fixed = showCoordinates, limits = NULL,
        main = NULL, note = NULL, showCoordinates = FALSE, showArrows = FALSE,
        showLegend = TRUE, frequencies = FALSE, statistics = FALSE,
        showAxes = FALSE, axesLabels = NULL, scenarios = NULL, help = NULL,
        helpOn = FALSE, mode = c("network","heatmap"), roundedItems = FALSE,
        controls = 1:8, cex = 1, background = NULL, defaultColor = "#1f77b4",
        language = c("en","es","ca"), dir = NULL)
}
\arguments{
\item{nodes}{a data frame with at least one vector of names.}
\item{links}{a data frame with at least two vectors with source and target, including names of nodes.}
\item{tree}{a data frame with two vectors: source and target, describing relationships between nodes.}
\item{name}{name of the vector with names in the nodes data frame. By default, if language="en", name is "name".}
\item{label}{name of the vector with labels in the nodes data frame.}
\item{group}{name of the vector with groups in the nodes data frame.}
\item{groupText}{show names of the groups.}
\item{community}{algorithm to make communities: edge_betweenness("ed"), fast_greedy("fa"), label_prop("la"), leiden_eigen("le"), louvain("lo"), optimal("op"), spinglass("sp"), walktrap("wa")}
\item{centrality}{calculates the centrality measures of a network. See \link{calCentr}.}
\item{labelSize}{name of the vector with label size in the nodes data frame.}
\item{size}{name of the vector with size in the nodes data frame.}
\item{color}{name of the vector with color variable in the nodes data frame.}
\item{shape}{name of the vector with shape variable in the nodes data frame.}
\item{border}{name of the column with border width in the nodes data frame.}
\item{legend}{name of the vector with the variable to represent as a legend in the nodes data frame.}
\item{ntext}{name of the vector with html text in the nodes data frame.}
\item{info}{name of the vector with information to display in a panel in the nodes data frame.}
\item{sort}{name of the vector with node order in the nodes data frame (only for heatmap).}
\item{decreasing}{decreasing or increasing sort of the nodes (only for heatmap).}
\item{intensity}{name of the vector with intensity variable in the links data frame (only for heatmap).}
\item{lwidth}{name of the vector with width variable in the links data frame.}
\item{lweight}{name of the vector with weight variable in the links data frame.}
\item{lcolor}{name of the vector with color variable in the links data frame.}
\item{ltext}{name of the vector with labels in the links data frame.}
\item{nodeFilter}{condition for filtering nodes.}
\item{linkFilter}{condition for filtering links.}
\item{degreeFilter}{numeric vector to filter the resulting network by degree.}
\item{nodeBipolar}{a logical value that polarizes negative and positive node values in the graphical representation. Default = FALSE.}
\item{nodeScaleLimits}{vector indicating the node scale limits, must be a numeric vector of length 2 or 3.}
\item{linkBipolar}{a logical value that polarizes negative and positive link values in the graphical representation. Default = FALSE.}
\item{linkScaleLimits}{vector indicating the link scale limits, must be a numeric vector of length 2 or 3.}
\item{defaultColor}{a character vector giving a valid html color.}
\item{repulsion}{a percentage for repulsion between nodes.}
\item{distance}{a percentage for distance of links.}
\item{zoom}{a number between 0.1 and 10 to start displaying zoom.}
\item{fixed}{prevent nodes from being dragged.}
\item{scenarios}{a note showing number of scenarios.}
\item{main}{upper title of the graph.}
\item{note}{lower title of the graph.}
\item{frequencies}{a logical value true if the frequencies can be shown in barplots. Default = FALSE.}
\item{statistics}{a logical value true if the possibility of calculating statistics with the nodes will be added to the final graph.}
\item{help}{help text of the graph.}
\item{helpOn}{Should the help be shown at the beginning?}
\item{background}{background color or image of the graph.}
\item{layout}{a matrix with two columns or an algorithm to elaborate the coordinates: davidson.harel
drl("da"), circle("ci"), Force-Atlas-2("fo"), fruchterman.reingold("fr"), gem("ge"), grid("gr"), kamada.kawai("ka"), lgl("lg"), mds("md"), random("ra"), reingold.tilford("re"), star("sta"), sugiyama("sug")}
\item{limits}{vector indicating size references to display layout, must be a numeric vector of length 4: x1, y1, x2, y2.}
\item{cex}{number indicating the amount by which plotting text should be scaled relative to the default. Default = 1.}
\item{roundedItems}{Display items with rounded borders.}
\item{controls}{a numeric vector indicating which controls will be shown. 1 = sidebar, 2 = selection buttons, 3 = export buttons, 4 = nodes table, 5 = links table, 6 = search, 7 = zoom, 8 = legend. NULL hide all controls, negative values deny each control and 0 deny all.}
\item{mode}{a character vector indicating the graph mode allowed: network, heatmap or both (both by default).}
\item{showCoordinates}{a logical value true if the coordinates are to be shown in tables and axes. Default = FALSE.}
\item{showArrows}{a logical value true if the directional arrows are to be shown. Default = FALSE.}
\item{showLegend}{a logical value true if the legend is to be shown. Default = TRUE.}
\item{showAxes}{a logical value true if the axes are to be shown. Default = FALSE.}
\item{axesLabels}{a character vector giving the axes names.}
\item{language}{a character vector (es=spanish; en=english; ca=catalan).}
\item{image}{name of the vector with image files in the nodes data frame.}
\item{imageNames}{name of the vector with names for image files in the nodes data frame.}
\item{dir}{a "character" string representing the directory where the web files will be saved.}
}
\references{
  Escobar, M. and Martinez-Uribe, L. (2020)
  Network Coincidence Analysis: The \code{netCoin} \code{R} Package.
  \emph{Journal of Statistical Software}, \bold{93}, 1-32.
  \doi{10.18637/jss.v093.i11}.
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
\examples{
# A character column (with separator)
frame <- data.frame(A = c("Man; Women", "Women; Women",
                        "Man; Man", "Undet.; Women; Man"))
data <- dichotomize(frame, "A", sep = "; ")[2:4]
C <- coin(data) # coincidence matrix
N <- asNodes(C) # node data frame
E <- edgeList(C) # edge data frame
netCoin(N, E) # netCoin object
}

