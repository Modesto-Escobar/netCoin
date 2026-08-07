\name{addClusters}
\alias{addClusters}
\title{Add a cluster column to a netCoin object from surScat.}
\description{
\code{addClusters} adds a cluster column to a netCoin object from surScat.
}
\value{
A 'netCoin' object.
}

\usage{
addClusters(scatObj, clusters, name=NULL, sort=TRUE,
  weight=NULL, sourceData=NULL)
}
\arguments{
\item{scatObj}{A 'netCoin' object.}
\item{clusters}{vector, factor, or clustering object (kmeans, poLCA, tidyLPA).}
\item{name}{name of the new cluster column. 'clusters' by default.}
\item{sort}{sort clusters if requested (by mean coordinate on first layout axis).}
\item{weight}{a vector of weights.}
\item{sourceData}{optional data.frame for validation that clusters are in correct order.}
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
