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
replaceClusters(scatObj, clusters,
  name=NULL, sort=TRUE, weight=NULL, sourceData=NULL)
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
