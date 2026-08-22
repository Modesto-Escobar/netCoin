\name{layoutDist}
\alias{layoutDist}
\title{Distance between two planes of a scattergram.}
\description{
\code{layoutDist} measures how far one and the same node travels from one plane of coordinates
to another, added up over the nodes and divided by how many they are.
}
\value{
A single number, the distance between the two planes. It carries as attributes the rows
compared ("rows") and those dropped ("dropped"), the alignment applied ("align"), the divisor
applied to each plane ("scales", one per plane under "cloud" and "procrustes", one per axis
under "canvas") and, with \code{align="procrustes"} and planes of two columns, the rotation in
degrees ("angle") and whether there was a reflection ("reflection").
}

\usage{
layoutDist(a, b, align=c("cloud", "procrustes", "canvas", "none"), weights=NULL,
           normalize=TRUE, who="layoutDist")
}
\arguments{
\item{a, b}{two planes, each a matrix or data.frame of coordinates holding one row per node,
in the same order.}
\item{align}{how the two planes are put on an equal footing before measuring, and with it what
is being asked. See the details.}
\item{weights}{a vector holding what every node counts for, applied to the sum and to the
centring, the scaling and the rotation alike. NULL (default) counts every node once; where the
nodes are patterns it is worth passing the cases each one stands for, so that a rare pattern
does not move as much as a common one.}
\item{normalize}{whether to divide the sum by the number of nodes compared, or by the total of
the weights where there are any. TRUE by default, which makes the figure how far any one node
travels: it no longer depends on the size of the survey, so two different objects are
comparable, and under \code{align="canvas"} it falls between 0 and sqrt(2), the diagonal of the
canvas. FALSE returns the plain sum.}
\item{who}{the name errors are issued under.}
}
\details{
Every plane comes in its own units: the coordinates of a tsne reduction reach 60 where those of
a classical one stop at 3, and nmf or a pair of observed variables live besides in the positive
quadrant. Adding up the distances as they come measures little else than that difference of
units, so the two planes are put on an equal footing first. How they are levelled decides what
is being asked:

\code{"cloud"} centres every plane on its centroid and divides it by its root mean square
radius. This strips translation and size away; what is left is a difference of configuration,
orientation included. The divisor is one for the whole plane rather than one per column:
standardising column by column would stretch the second axis until it matched the first, which
undoes what a reduction states, that its first direction carries more variation than its
second.

\code{"procrustes"} does the above and rotates besides the first plane, reflecting it where
that helps, until it lies as close to the second as it can. It measures shape alone: two
reductions ordering the cases alike but coming out rotated or mirrored end up at nearly no
distance, which is the reasonable reading wherever the orientation of the axes means nothing,
as in tsne or umap. No scaling factor of the fit is added, though the classical Procrustes
statistic carries one: the sizes have been levelled already, and that factor would break the
symmetry a matrix of distances between planes needs.

\code{"canvas"} brings every column to [0, 1] by its minimum and its maximum, which is what
rD3plot does when it sets the domain of each axis to the range of the coordinates and its own
range to the width or the height of the canvas. It measures the picture the reader is shown
rather than the cloud behind it, and two planes differing only in the scale of their axes come
out at zero. The canvas is taken square, as the shape of the window is not known from R.

\code{"none"} leaves the coordinates as they are, which fits where the units are already
comparable, as between two pairs of variables measured from 0 to 100.

Rows holding NA in either plane are dropped, and the alignment is computed on the rows about to
be compared rather than on the whole of them.
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
\seealso{
\code{\link{layoutsDist}}, \code{\link{addAxes}}, \code{\link{surScat}}
}
\examples{
\dontrun{
data(iris)
S <- surScat(iris, variables=names(iris)[1:4], type="pca", nclusters=3)
P <- prcomp(iris[,1:4], scale.=TRUE)$x
S <- addAxes(S, P[,3:4])
L <- S$layouts

layoutDist(L[[1]], L[[2]])                        # the cloud
layoutDist(L[[1]], L[[2]], align="procrustes")    # shape alone
layoutDist(L[[1]], L[[2]], align="canvas")        # the picture

z <- layoutDist(L[[1]], L[[2]], align="procrustes")
attr(z, "angle")
attr(z, "reflection")
}
}
