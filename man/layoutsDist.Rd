\name{layoutsDist}
\alias{layoutsDist}
\title{Distances between all the planes of a scattergram.}
\description{
\code{layoutsDist} measures how far apart every pair of planes of a netCoin object from
surScat lies, and hands the result back as the dist object hclust, cmdscale and the like
expect.
}
\value{
An object of class 'dist' holding one distance per pair of planes, with the alignment recorded
as its "method". It carries as attributes the rows compared per pair ("rows"), the weights
actually applied ("weights") and where they came from ("wSource") and, with
\code{align="procrustes"}, the rotation in degrees ("angle") and whether there was a
reflection ("reflection"), both as square matrices.
}

\usage{
layoutsDist(x, align=c("cloud", "procrustes", "canvas", "none"), weights=NULL,
            which=NULL, normalize=TRUE, who="layoutsDist")
}
\arguments{
\item{x}{a 'netCoin' object from surScat, or a list of planes, each a matrix or data.frame of
two columns of coordinates holding one row per node.}
\item{align}{how the two planes are put on an equal footing before measuring, and with it what
is being asked. See \code{\link{layoutDist}}.}
\item{weights}{what every node counts for. The name of a numeric node column, or a vector
holding one value per node or per case. NULL (default) weighs every node by the cases it
stands for when the nodes are patterns, and counts one per node otherwise. \code{weights=1}
forces the plain count of nodes.}
\item{which}{the planes to compare, as names or as indices. NULL (default) takes every one of
them.}
\item{normalize}{whether to divide by the number of nodes compared, or by the total of the
weights where there are any, so that every entry is how far any one node travels. TRUE by
default.}
\item{who}{the name errors and messages are issued under.}
}
\details{
A scattergram may hold several planes: the one it is drawn with, the reductions of a
'looking4clusters' object and any pair of variables added through \code{\link{addAxes}}. All
of them hold one row per node, in the same order, so node \emph{i} is the same case in every
one of them and it makes sense to ask how far it travels from one to another.

The distances are computed pair by pair with \code{\link{layoutDist}}, which is where the
alignments and the weighting are described. Being a plain dist, the planes can be clustered
with \code{hclust} or placed against one another with \code{cmdscale} with no further ado.

Two readings are worth keeping apart. Under \code{align="cloud"} a zero means two planes hold
the same configuration of points; under \code{align="canvas"} it means the viewer shows them
as the same picture, which is a weaker statement, as rD3plot stretches the range of each axis
to the width or the height of the canvas and so draws as one two planes differing in the scale
of their axes. A distance collapsing under \code{align="procrustes"} points to two reductions
agreeing but for their orientation.

The rotations are handed back as a full square matrix rather than inside the dist because they
are antisymmetric: rotating A onto B is undoing the rotation of B onto A. Reflections are
symmetric, and come as a square matrix alongside them for the sake of one shape.
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
\seealso{
\code{\link{layoutDist}}, \code{\link{addAxes}}, \code{\link{surScat}}
}
\examples{
\dontrun{
data(iris)
S <- surScat(iris, variables=names(iris)[1:4], type="pca", nclusters=3)
P <- prcomp(iris[,1:4], scale.=TRUE)$x
S <- addAxes(S, P[,3:4])
S <- addAxes(S, as.matrix(iris[,1:2]), name="sepal")

D <- layoutsDist(S)
D
as.matrix(D)
plot(hclust(D))

# what the reader is actually shown
layoutsDist(S, align="canvas")

# only some of the planes, and the rotations behind the shape
layoutsDist(S, which=c("sepal", "PC3-PC4"))
Q <- layoutsDist(S, align="procrustes")
round(attr(Q, "angle"))
attr(Q, "reflection")
}
}
