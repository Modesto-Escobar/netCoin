\name{layoutsAlign}
\alias{layoutsAlign}
\title{Planes of a scattergram laid over one another}
\description{
\code{layoutsAlign} brings the planes of a netCoin object from surScat to a common footing, as
\code{\link{layoutsDist}} does to measure how far apart they are, and hands back the transformed
coordinates themselves.
}
\value{
When \code{x} is a 'netCoin' object, the same object with the aligned planes in place of the
original ones in its \code{$layouts} (or each one next to its original, when \code{suffix} is
given), the
first plane also as its fx/fy columns, and the cluster
centroids in \code{$clusters} computed again on the new coordinates. The details of the
alignment go in its attribute "alignment".

When \code{x} is a list of planes, a named list holding the aligned ones, with those details as
attributes of the list.

The details are: the method (\code{align}), the reference plane (\code{to}), whether
\code{gpa} was used, the rows the fit was estimated on (\code{rows}) and those left out
(\code{dropped}), the center, scale and rotation matrix applied to every plane, and, with
\code{align="procrustes"} and two columns, the angle of each rotation in degrees and whether it
involved a reflection. With \code{gpa=TRUE} there are also the \code{consensus} plane (the mean
of all the aligned ones), the number of \code{iterations} and whether the fit
\code{converged}.
}

\usage{
layoutsAlign(x, align=c("procrustes", "cloud", "canvas", "none"), to=1, gpa=FALSE,
             weights=NULL, which=NULL, suffix=NULL, tol=1e-10, maxit=100,
             who="layoutsAlign")
}
\arguments{
\item{x}{a 'netCoin' object from surScat, or a list of planes, each a matrix or data.frame of
two columns of coordinates holding one row per node.}
\item{align}{how the planes are put on an equal footing: \code{"procrustes"} (default) centres
them, brings them to the same size and rotates (or reflects) each one onto the reference;
\code{"cloud"} does the same without rotating; \code{"canvas"} stretches every axis to the same
range, as the viewer does; \code{"none"} leaves them as they are. See
\code{\link{layoutDist}}.}
\item{to}{the reference plane, as its name or its index among all the planes of \code{x}. The
first one by default, which in a scattergram from surScat is the plane it was drawn with (for
instance "PC1-PC2"). The names are those printed under "Planes" by \code{print(x)}, or
\code{names(x$layouts)}: a reduction taken from a 'looking4clusters' object by
\code{\link{addAxes}} keeps its own name, such as "pca" or "tsne", unless \code{name} was given
there.}
\item{gpa}{with \code{align="procrustes"}, rotate every plane towards the mean of them all
(generalized Procrustes analysis) instead of towards \code{to}. The reference then only sets the
orientation and the units in which the whole set is handed back.}
\item{weights}{what every node counts for in the fit, as in \code{\link{layoutsDist}}.}
\item{which}{the planes to align, as names or as indices. NULL (default) takes every one of
them. The reference need not be among them; planes left out are kept as they were.}
\item{suffix}{only for a 'netCoin' object. NULL (default) replaces its planes by the aligned
ones. A string keeps the original planes and places each aligned one right after its original,
so that the selector of the viewer offers them side by side, named as the original followed by
\code{suffix} (with \code{suffix="'"}, "sepal" gives "sepal'"). The
reference plane is not added again, as it is left unchanged, except with \code{gpa=TRUE}.}
\item{tol, maxit}{convergence tolerance and maximum number of iterations of the generalized
Procrustes fit.}
\item{who}{the name errors and messages are issued under.}
}
\details{
The fit is estimated on the rows that are complete in every plane involved, the same ones for
all of them, and then applied to every row of each plane: a node missing from one plane keeps
its place in the others.

The result is expressed in the units of the reference plane, which is therefore left just as it
was (except under \code{gpa=TRUE}, where every plane moves towards the consensus). A plane
compared with it under \code{align="procrustes"} is thus turned to be read on the same axes, and
the distance \code{\link{layoutsDist}} reports between them is the mean distance between the two
sets of coordinates divided by the scale of the reference.

Under \code{gpa=TRUE} no plane is favoured over the others. The mean has no orientation of its
own, so the whole set is finally turned, as one, to lie as close as possible to the reference.
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
\seealso{
\code{\link{layoutsDist}}, \code{\link{layoutDist}}, \code{\link{addAxes}}, \code{\link{surScat}}
}
\examples{
\dontrun{
data(iris)
S <- surScat(iris, variables=names(iris)[1:4], type="pca", nclusters=3)
P <- prcomp(iris[,1:4], scale.=TRUE)$x
S <- addAxes(S, P[,3:4])
S <- addAxes(S, as.matrix(iris[,1:2]), name="sepal")

# every plane turned onto PC1-PC2, the one drawn first
A <- layoutsAlign(S)
A
attr(A, "alignment")$angle

# onto another plane, by its name
A <- layoutsAlign(S, to="sepal")

# the original planes kept, the aligned ones added as "PC3-PC4'" and "sepal'"
A <- layoutsAlign(S, suffix="'")
names(A$layouts)

# towards the mean of them all
G <- layoutsAlign(S, gpa=TRUE)
head(attr(G, "alignment")$consensus)

# the coordinates alone, from a list of planes
L <- layoutsAlign(S$layouts, which=c("sepal", "PC3-PC4"), to=1)
}
}
