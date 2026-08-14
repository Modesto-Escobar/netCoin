\name{addAxes}
\alias{addAxes}
\title{Add axes to a netCoin object from surScat.}
\description{
\code{addAxes} adds one or several planes of coordinates to a netCoin object from surScat,
which the viewer then offers in a selector next to the graph.
}
\value{
A 'netCoin' object.
}

\usage{
addAxes(scatObj, axes, name=NULL, which=NULL, weight=NULL)
}
\arguments{
\item{scatObj}{A 'netCoin' object.}
\item{axes}{a matrix or data.frame with an even number of columns of coordinates, a pair per
plane; a list of such matrices; or a 'looking4clusters' object, whose dimensionality reductions
are taken as planes.}
\item{name}{name of the new plane. When several planes are added at once, each one keeps its
own name, and 'name' becomes a prefix.}
\item{which}{reduction to take when 'axes' is a 'looking4clusters' object, given as a name or
an index. When it is left empty, every reduction held by the object is added as its own plane.}
\item{weight}{a vector of weights, one per case, used when collapsing cases into patterns. NULL
(default) takes the weights surScat was called with, so that the plane being added is collapsed
just as the one already drawn, and collapses by the plain arithmetic mean when there were none.
A vector of ones asks for the plain mean of a scattergram that does have weights.}
}
\details{
A scattergram drawn by surScat holds a single pair of coordinates, the one stated by its
\code{xaxis} and \code{yaxis} arguments. \code{addAxes} adds further pairs to the same
scattergram, so that the very same nodes, with their groups and their profiles, can be looked
at from several projections without building a graph for each. The plane that was already
drawn becomes the first entry of the selector, named after the axes it was drawn with.

Each pair of columns of \code{axes} makes up a plane: the first column is drawn horizontally
and the second vertically. Planes are named after their two columns, as in "PC1-PC2", unless
\code{name} states otherwise. A 'looking4clusters' object names them after its reductions
("pca", "tsne", "cmds", "nmf", "umap", plus any added through its own \code{addreduction}
function). The 'looking4clusters' package needs not be installed: the object is read through
its list structure.

The coordinates must come from the same data, in the same order, as the one given to surScat:
only their length is checked, so keeping the order is up to the caller. Coordinates given case
by case are collapsed into patterns, by their (weighted) mean, when surScat was called with
\code{vPatterns}. The weights are the ones surScat itself collapsed with, which it keeps, unless
\code{weight} states otherwise: a plane added to a weighted scattergram is therefore comparable
with the one already drawn without having to repeat the weights.

They must therefore hold one row per node, or one per case when the nodes are patterns. Anything
else is an error, which is also what happens when surScat capped the nodes through its
\code{maxN} argument: the node table then holds a random sample of the cases, or of the patterns,
so coordinates computed on the whole data no longer line up with it and cannot be matched to it.
Compute the coordinates on the node table of the object, which is that sample, or raise
\code{maxN} so that every case is drawn.

The labels of the axes belong to the graph and not to each plane, so they keep describing the
plane that surScat drew. Naming the added planes after their axes keeps the reader oriented.
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
\examples{
\dontrun{
data(iris)
S <- surScat(iris, variables=names(iris)[1:4], type="pca", nclusters=3)

# a second plane from further principal components
P <- prcomp(iris[,1:4], scale.=TRUE)$x
S <- addAxes(S, P[,3:4])

# every reduction of a looking4clusters object
# S <- addAxes(S, l4cObject)
plot(S)
}
}
