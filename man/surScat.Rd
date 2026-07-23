\name{surScat}
\alias{surScat}
\title{Interactive scattergram of the cases of a data frame.}
\description{
\code{surScat} produces an interactive scattergram of the cases (rows) of a data frame. Its axes are the two first dimensions of a multiple correspondence analysis (MCA) of the active qualitative variables, or of a principal component analysis (PCA) of the active quantitative variables. Cases are colored by groups obtained by k-means clustering of these coordinates.
}
\value{
This function creates a netCoin object and, if stated, a folder in the computer with an HTML document named index.html which contains the produced graph. This file can be directly opened with your browser and sent to a web server to work properly.
}

\usage{
surScat(data, variables=names(data), active=variables, weight=NULL,
        patterns=FALSE, jitter=0, type=c("mca", "pca"), nclusters=2,
        clusterOn=c("factors", "variables"), nfactors=2, critFactors=0,
        seed=2020, maxN=2000, ...)
}
\arguments{
  \item{data}{a data frame.}
  \item{variables}{a vector of variables included in the previous data frame.}
  \item{active}{a vector of active variables used to compute the coordinates.}
  \item{weight}{a numeric vector of case weights, or the name of a column of the data frame containing them (this column is then excluded from variables). Weights are used to compute the MCA or PCA coordinates (weighted category masses in MCA, weighted correlations in PCA), so the axes represent the weighted structure of the data. The k-means clustering is then computed on these weighted coordinates without further weighting. NULL (default) performs an unweighted analysis.}
  \item{patterns}{if TRUE, cases sharing the same response pattern in all the variables are collapsed into a single node, whose size represents the number of supporting cases (or the sum of their weights): a column n is added to the nodes and used as size argument of \link{netCoin} (unless another size is stated). In addition, a column ntext with an HTML description of each pattern (the values of all its variables plus its number of cases) is generated and passed as ntext argument of \link{netCoin}, unless another ntext is stated. Coordinates and clusters are identical to those of the case by case analysis. Recommended for MCA with few qualitative variables, where cases pile up on a small number of distinct positions.}
  \item{jitter}{amount of random noise added to the coordinates before plotting, expressed as a fraction of the range of each axis (e.g. .02 = 2 percent; TRUE is a shortcut for .02). It only spreads the points visually: clusters and coordinates are computed on the exact positions. It can be combined with patterns to separate patterns differing only in non active variables.}
  \item{type}{Factorial type: mca for qualitative active variables, pca for quantitative active variables. If not stated, it is chosen automatically: pca if all the active variables are quantitative, mca otherwise.}
  \item{nclusters}{number of clusters for k-means. It can be a vector, in which case one group variable is added per value.}
  \item{clusterOn}{what k-means clusters on: "factors" (default) uses a selection of the MCA/PCA factorial coordinates (see nfactors and critFactors). "variables" uses the active variables instead (their dichotomized indicators for mca, their standardized values for pca), so groups reflect the full active variable space rather than only the extracted factors. nfactors and critFactors are ignored in this case, with a warning if they are explicitly set.}
  \item{nfactors}{when clusterOn="factors", the maximum number of MCA/PCA factors considered for clustering (critFactors then selects how many of these are actually used, at least 1). It only affects clustering: the scattergram plot itself always uses the first two factors, regardless of clusterOn, nfactors or critFactors. Ignored (with a warning if set explicitly) when clusterOn="variables".}
  \item{critFactors}{when clusterOn="factors", a factor among the first nfactors is used for clustering only if its eigenvalue exceeds critFactors times the mean eigenvalue of the full factorial spectrum. This generalizes Kaiser's rule (eigenvalue above average) to MCA, whose eigenvalues are not on the same scale as PCA's: critFactors=1 reproduces the classic Kaiser criterion in (unweighted or weighted) PCA, since the mean eigenvalue of a correlation matrix is always 1, while in MCA it adapts to the actual inertia scale of the analysis. The default, 0, applies no filtering, so all nfactors factors are used. At least one factor is always used. Ignored (with a warning if set explicitly) when clusterOn="variables".}
  \item{seed}{random seed for reproducible k-means clustering, row sampling and jitter. The RNG state of the session is restored on exit. Use NULL to avoid setting a seed.}
  \item{maxN}{Maximum number or rows.}
  \item{...}{Any \link{netCoin} argument.}
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
# A data frame with two variables Gender and Opinion
frame<-data.frame(Gender=c(rep("Man",3),rep("Woman",3)),
                  Opinion=c("Yes","Yes","No","No","No","Yes"))
surScat(frame) # network object

# PCA example: quantitative active variables, clustered on those variables
# (not on the two plotted PCA axes), shaped by a non-active variable
data(iris)
surScat(iris, active=c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
        shape="Species", nclusters=2:3, clusterOn="variables")
}

