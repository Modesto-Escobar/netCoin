\name{addCommunities}
\alias{addCommunities}
\title{Add communities of cases to a netCoin object from surScat.}
\description{
\code{addCommunities} adds one cluster column per community detection method to a netCoin object
from surScat, taken from the network of correlations between its cases.
}
\value{
A 'netCoin' object.
}

\usage{
addCommunities(scatObj, data, variables=names(data), community="lo",
  scale=TRUE, method=c("pearson", "kendall", "spearman"),
  criteria="value", minL=0.95, maxL=Inf, pairwise=FALSE,
  minSize=0.01, name=NULL, sort=TRUE, weight=NULL,
  maxGroups=25, maxNodes=NULL)
}
\arguments{
\item{scatObj}{A 'netCoin' object returned by \code{\link{surScat}}.}
\item{data}{the data.frame the scattergram was drawn from.}
\item{variables}{the variables the cases are compared over. All of them by default.}
\item{community}{one or several community detection methods, each given by its two-letter code
or by its name: edge_betweenness ("ed"), fast_greedy ("fa"), label_prop ("la"),
leading_eigen ("le"), louvain ("lo"), optimal ("op"), spinglass ("sp"), walktrap ("wa").
One cluster column is added per method.}
\item{scale}{standardize each variable before transposing.}
\item{method}{correlation coefficient: "pearson", "kendall" or "spearman".}
\item{criteria}{which statistic the links are filtered on, "value" for the correlation itself
or "p" for its significance.}
\item{minL}{lowest value of \code{criteria} a pair of cases needs to be linked.}
\item{maxL}{highest value of \code{criteria} a pair of cases may have to be linked.}
\item{pairwise}{use pairwise complete observations for the correlations.}
\item{minSize}{communities holding fewer cases than this are gathered into a residual group,
read as a proportion of the cases when below 1 and as a count of them otherwise. 1 keeps every
community apart.}
\item{name}{name of the new cluster column. Each column keeps the name of its method by
default, and 'name' becomes a prefix when several methods are asked for at once.}
\item{sort}{sort the groups by their mean coordinate on the first axis. The residual group is
pinned last whatever its coordinate.}
\item{weight}{a vector of weights for the mean coordinate the groups are sorted by. It does
not weigh the correlations themselves.}
\item{maxGroups}{clusterizations with more groups than this are skipped with a warning. See
\code{\link{addClusters}}.}
\item{maxNodes}{how large a network a method may be asked about, overriding the ceiling each
one carries. NULL keeps those.}
}
\details{
The cases are correlated with one another over the variables, which is the transposition of the
usual reading of a data.frame: each case becomes a column and each variable a row, so that a link
joins two cases whose profiles run parallel. The communities of that network are a grouping of
cases, of the same shape as any other clusterization, and reach the node table through
\code{\link{addClusters}}, with the naming and the labelling described there. A louvain of five
groups is a column named "louvain(5)", labelled "Louvain: 1" ... "Louvain: 5".

A correlation between two cases centres each of them on its own mean across the variables, so
what it compares is the shape of their profiles. Standardizing the variables first puts them all
on one footing; without it, the variable of widest range writes the profile on its own, which is
why \code{scale} is on by default. Each correlation rests on as many points as there are
variables, so comparing cases over very few of them says very little, and warns.

A network of cases leaves isolated ones behind, those reaching no other case above \code{minL},
and igraph returns each of them as a community of its own. \code{minSize} gathers those and any
other small community into one residual group, which is pinned after the real ones and labelled
with a word instead of a number, so that it reads apart from them in a legend. Without it a
scattergram of three readable groups reports dozens and runs into \code{maxGroups}.

The whole network is built before the links are filtered, so its cost grows with the square of
the number of cases: some thousands of them are a matter of seconds and a few gigabytes, and
well beyond that it stops being feasible. The methods themselves do not all bear the same size
either, and each one carries the ceiling it can meet: cluster_optimal is exponential and refuses
above a hundred nodes, whereas cluster_louvain answers on thousands.

When surScat collapsed the cases into patterns through \code{vPatterns}, the communities are
found case by case and then reduced to the most frequent one of each pattern, as any other
clusterization is. When surScat capped the nodes through \code{maxN}, the network is built on the
cases that were drawn, which the object records, so that the communities line up with the nodes
that are actually plotted. The two together cannot be resolved: a sample of patterns cannot be
traced back to the cases behind it, and that is an error.
}
\seealso{
\code{\link{addClusters}}, \code{\link{surScat}}, \code{\link{netCorr}}.
}
\examples{
\dontrun{
data <- data.frame(a=rnorm(300), b=rnorm(300), c=rnorm(300), d=rnorm(300))
S <- surScat(data, maxN=Inf)
S <- addCommunities(S, data, community=c("louvain", "walktrap"))
plot(S)
}
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
