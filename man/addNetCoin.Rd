\name{addNetCoin}
\alias{addNetCoin}
\title{Networked coincidences.}
\description{
\code{addNetCoin} changes netCoin object attributes.
}
\value{
This function returns a \code{netCoin} object.
If the 'dir' attribute is specified, the function creates a folder in the computer with an HTML document named index.html which contains the produced graph. This file can be directly opened with your browser and sent to a web server to work properly.
}

\usage{
addNetCoin(x, ...)
}
\arguments{
\item{x}{A 'netCoin' object.}
\item{...}{Any \link{netCoin} argument.}
}
\details{
The object is rebuilt out of its own options, with the ones stated here replacing them, so
anything not among the arguments of \link{netCoin} would be lost in the rebuilding and is
carried over instead: the planes \code{\link{addAxes}} added, and what a scattergram records
about how its nodes stand for the cases, which \code{\link{addClusters}},
\code{\link{addCommunities}} and \code{replaceClusters} read afterwards. Stating
\code{layout} replaces the planes, as asked; replacing \code{nodes} with a table of another
length carries none of it over, since none of it would describe those nodes.

Nodes drawn with no label at all, as \code{\link{surScat}} draws them, keep it that way. A
label is restored by naming its column, as in \code{label="name"}.
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
net <- netCoin(N, E) # netCoin object
net <- addNetCoin(net, color="frequency") # netCoin object with color applied to nodes
}

