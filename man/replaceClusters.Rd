\name{replaceClusters}
\alias{replaceClusters}
\title{Replace cluster columns in a netCoin object from surScat.}
\description{
\code{replaceClusters} replaces cluster columns in a netCoin object from surScat.
Removes all existing cluster columns and adds the new one.
}
\value{
A 'netCoin' object.
}

\usage{
replaceClusters(scatObj, clusters, name=NULL, sort=TRUE,
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
The columns to remove are the ones surScat and \code{\link{addClusters}} recorded as holding
clusterizations, so variables of the study are never removed, whatever they are named. Objects
built before this was recorded fall back to the name of the column surScat creates.
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
