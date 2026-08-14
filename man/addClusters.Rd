\name{addClusters}
\alias{addClusters}
\title{Add cluster columns to a netCoin object from surScat.}
\description{
\code{addClusters} adds one or several cluster columns to a netCoin object from surScat.
}
\value{
A 'netCoin' object.
}

\usage{
addClusters(scatObj, clusters, name=NULL, sort=TRUE,
  weight=NULL, sourceData=NULL, which=NULL, maxGroups=25)
}
\arguments{
\item{scatObj}{A 'netCoin' object.}
\item{clusters}{vector, factor, data.frame holding one clusterization per column, or
clustering object (kmeans, poLCA, tidyLPA, looking4clusters).}
\item{name}{name of the new cluster column. 'clusters' by default. When several
clusterizations are added at once, each column keeps its own name, and 'name' becomes a prefix.}
\item{sort}{sort clusters if requested (by mean coordinate on first layout axis).}
\item{weight}{a vector of weights.}
\item{sourceData}{optional data.frame for validation that clusters are in correct order.}
\item{which}{clusterization to take when 'clusters' is a 'looking4clusters' object, given as
a name or an index. Methods run for several numbers of groups yield their optimal solution,
unless a specific number is asked for with the "name:k" form, as in \code{which="kmeans:4"}.
When it is left empty, every clusterization held by the object is added as its own column,
one per number of groups, named after the method and that number, as in "kmeans(4)".}
\item{maxGroups}{clusterizations with more groups than this are skipped with a warning, which
keeps a continuous variable passed in by mistake from becoming a column of hundreds of groups.
NULL adds them all.}
}
\details{
Each column is an ordered factor whose levels follow the group number. Their labels are prefixed
with the name of the clusterization they come from, capitalized and without its trailing number of
groups, so a column named "hclust(10)" is labelled "Hclust: 01" ... "Hclust: 10". Group numbers are
padded with zeros to the width of the largest one, which keeps them in numeric order wherever they
are sorted as text.

Each clusterization is named after where it comes from and how many groups it holds, as in
"kmeans(4)" for a method of a 'looking4clusters' object. The solutions of a 'tidyLPA' object are
named "LPA(4)", and labelled "LPA: 1" ... "LPA: 4", since a latent profile analysis estimated for
several numbers of profiles is one analysis. A 'tidyLPA' object estimated for several model
specifications at once is no longer a single analysis, so its solutions carry the number tidyLPA
gives each specification, as in "LPA1(4)" and "LPA6(4)".

Several clusterizations can be added in one call, either as the columns of a data.frame or as a
'looking4clusters' object. They must all come from the same data, in the same order, as the one
given to surScat: only their length is checked, so keeping the order is up to the caller.

A clusterization must hold one assignment per node, or one per case when surScat collapsed the
cases into patterns through \code{vPatterns}, in which case the cases of each pattern are reduced
to their most frequent group. Anything else is an error, which is also what happens when surScat
capped the nodes through its \code{maxN} argument: the node table then holds a random sample of
the cases, or of the patterns, so a clusterization computed on the whole data no longer lines up
with it and cannot be matched to it. Compute the clusters on the node table of the object, which
is that sample, or raise \code{maxN} so that every case is drawn.

A 'looking4clusters' object usually holds many clusterizations at once: those added through
its own \code{addcluster} function, plus one solution per number of groups for each method it
has run. Leaving \code{which} empty brings all of them in, whereas \code{which} picks a single one.
The 'looking4clusters' package needs not be installed: the object is read through its list
structure.
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
