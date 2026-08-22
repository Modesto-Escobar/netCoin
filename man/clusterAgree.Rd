\name{clusterAgree}
\alias{clusterAgree}
\title{Agreement between the classifications of a scattergram.}
\description{
\code{clusterAgree} measures the agreement between a model classification and the other
classifications held by a netCoin object from surScat through ARI, NMI and VI, which match no
categories and take therefore classifications with any number of groups.
}
\value{
A data.frame with one row per classification, sorted as \code{sort} states. It carries as
attributes the contingency table of every classification ("tables"), the number of categories
of the model ("kModel"), the unit of analysis ("unit"), the weights actually applied
("weights") and, where there were any, the patterns spanning several groups ("spanning").
}

\usage{
clusterAgree(x, model, columns=NULL, weights=NULL, modal=TRUE, sort=NULL, digits=3)
}
\arguments{
\item{x}{a 'netCoin' object from surScat, or the data.frame of its nodes.}
\item{model}{the classification everything else is compared with: the name of a node column,
or a vector holding one value per node. Where the nodes are patterns, a vector holding one
value per case is taken as well, and the comparison then goes case by case.}
\item{columns}{the classifications to compare. NULL (default) takes those surScat recorded in
the "clusterColumns" attribute and, where there are none, every factor or character column of
the node table.}
\item{weights}{what every unit counts for. The name of a numeric node column, or a vector
holding one value per node or per case. NULL (default) weighs every node by the cases it
stands for when the nodes are patterns, and counts one per unit otherwise. \code{weights=1}
forces the plain count of nodes.}
\item{modal}{where the nodes are patterns, whether to keep the modal group of each pattern.
See \code{\link{clusterMatch}}.}
\item{sort}{the column to sort the result by, from the largest down; a vector of names sorts
by all of them at once. The prefix "-" turns one column around. NULL (default) keeps the order
the classifications came in.}
\item{digits}{decimal places the numeric columns are rounded to.}
}
\details{
Splitting finer never harms precision --- in the limit, one group per case leaves it at 1 ---
so the measures of a confusion matrix cannot arbitrate between classifications with different
numbers of groups, and \code{\link{clusterMatch}} refuses to compare them at all. These three
indices are there for that: they read the two classifications as partitions and need no
matching of categories.

\code{ARI}, the adjusted Rand index, walks through the pairs of cases, counting whether the two
classifications put them together or keep them apart, and corrects for chance: 0 is what
independent partitions yield, 1 an exact match. It may go slightly negative.

\code{NMI}, the normalised mutual information, tells how much one classification cuts down the
uncertainty about the other, from 0 to 1. It is not corrected for chance, so it tends to rise
with the number of categories and does not arbitrate well between solutions of different size.

\code{VI}, the variation of information, is a distance in nats: 0 is an exact match and it
grows with the discrepancy. \code{VI.n} divides it by log(N) to keep it between 0 and 1. Note
when sorting that in these two the small value is the good one, so \code{sort="-VI"} is
usually what is wanted.

Where the nodes are patterns, weighting and modal groups work just as in
\code{\link{clusterMatch}}. Weights are taken to be on a scale of cases, which is what the
default weighting yields: the count of pairs behind ARI is extended continuously, and weights
given as proportions would make it meaningless.
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
\seealso{
\code{\link{clusterMatch}}, \code{\link{surScat}}, \code{\link{addClusters}}
}
\examples{
\dontrun{
data(iris)
S <- surScat(iris, variables=names(iris)[1:4], type="pca", nclusters=2:4)

# every classification, whatever its number of groups
A <- clusterAgree(S, "kmeans(3)", sort="ARI")
A
attr(A, "tables")[[1]]

# the most concordant on top, VI being a distance
clusterAgree(S, "kmeans(3)", sort="-VI")
}
}
