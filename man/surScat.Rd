\name{surScat}
\alias{surScat}
\title{Interactive scattergram of the cases of a data frame.}
\description{
\code{surScat} produces an interactive scattergram of the cases (rows) of a data frame. Its axes are, by default, the two first dimensions of a multiple correspondence analysis (MCA) of the active qualitative variables, or of a principal component analysis (PCA) of the active quantitative variables. Cases are colored by groups obtained by k-means clustering of these coordinates.

Either axis, or both, may be changed with xaxis and yaxis: a different factorial dimension, or an observed variable of the data frame. Clustering is not affected by this choice, so a plot of an observed variable against a factor still colors its cases according to the whole structure found in the active variables.

A plotted factorial axis is labelled with the percentage of variance (PCA) or of inertia (MCA) it accounts for, computed on the full factorial spectrum. In MCA the raw inertia percentages are known to be strongly pessimistic, so Benzecri's adjustment is applied: only the factors whose eigenvalue exceeds 1/Q (Q being the number of active variables) are retained, their eigenvalues are transformed into (Q/(Q-1))^2 (eigenvalue - 1/Q)^2, and the percentages are computed on this adjusted inertia. The adjustment is skipped, and the raw inertia percentages are kept, when the second factor does not itself exceed the 1/Q threshold, since it would then be labelled as 0 percent. An axis holding an observed variable is labelled with the label of that variable, or with its name if it has none. These labels are passed to \link{netCoin} as its axesLabels argument, unless another axesLabels is stated.

Every node carries an HTML description, in a column named ntext which is passed to \link{netCoin} as its ntext argument, unless another ntext is stated. When cases are collapsed it describes the pattern (see patterns). Otherwise it describes the case itself: its value on each of the two axes, followed in brackets by the mean of the k-means group it belongs to on that same axis, and the number of that group (the last k-means solution, when nclusters states several). A factorial axis is described by the coordinate which is plotted, whereas an observed variable is described in its own units and with its own labels, even if its axis has been standardized.
}
\value{
This function creates a netCoin object and, if stated, a folder in the computer with an HTML document named index.html which contains the produced graph. This file can be directly opened with your browser and sent to a web server to work properly.
}

\usage{
surScat(data, variables=names(data), active=variables, weight=NULL,
        patterns=FALSE, vPatterns=variables, jitter=0, type=c("mca", "pca"),
        xaxis=NULL, yaxis=NULL, scaleAxes=TRUE, nclusters=2,
        clusterOn=c("factors", "variables"), scaleClusters=NULL, nstart=25,
        nfactors=2, critFactors=0, sortClusters=TRUE, columns=NULL,
        seed=2020, maxN=2000, ...)
}
\arguments{
  \item{data}{a data frame.}
  \item{variables}{a vector of variables included in the previous data frame.}
  \item{active}{a vector of active variables used to compute the coordinates.}
  \item{weight}{a numeric vector of case weights, or the name of a column of the data frame containing them (this column is then excluded from variables). Weights are used to compute the MCA or PCA coordinates (weighted category masses in MCA, weighted correlations in PCA), so the axes represent the weighted structure of the data. The k-means clustering is then computed on these weighted coordinates without further weighting. When patterns collapses cases (see patterns and vPatterns), weights are also used for the mean/mode summaries and for the case count of each pattern. NULL (default) performs an unweighted analysis.}
  \item{patterns}{if TRUE, cases sharing the same response pattern in vPatterns are collapsed into a single node, whose size represents the number of supporting cases (or the sum of their weights): a column n is added to the nodes and used as size argument of \link{netCoin} (unless another size is stated). In addition, a column ntext with an HTML description of each pattern (the values of its vPatterns variables plus its number of cases) is generated and passed as ntext argument of \link{netCoin}, unless another ntext is stated. Coordinates, k-means clustering and the eigenvalue-based factor selection are always computed case by case, exactly as when patterns=FALSE; collapsing only affects how the resulting node table is built afterwards: each pattern is plotted at the (weighted) mean position of its cases, variables outside vPatterns take their (weighted) mean if quantitative or their (weighted) mode if qualitative, and the k-means group columns list every distinct group reached by the pattern's cases, separated by "|" (a single value if all its cases share the same group). Recommended for MCA with few qualitative variables, where cases pile up on a small number of distinct positions. Automatically set to TRUE if vPatterns is stated.}
  \item{vPatterns}{a subset of variables defining the response pattern used to collapse cases when patterns=TRUE (or is implied by stating vPatterns). Defaults to variables, i.e. cases are collapsed only if they agree on every variable. Stating a subset collapses cases sharing a pattern on just those variables: the remaining variables (including active ones) and the k-means group columns are then summarized as described under patterns, and the ntext description is restricted to just the vPatterns variables.}
  \item{jitter}{amount of random noise added to the coordinates before plotting, expressed as a fraction of the range of each axis (e.g. .02 = 2 percent; TRUE is a shortcut for .02). It only spreads the points visually: clusters and coordinates are computed on the exact positions. It can be combined with patterns to separate patterns differing only in non active variables, and it is particularly useful when an axis holds a discrete observed variable (see xaxis and yaxis), whose few distinct values pile the cases up into columns.}
  \item{type}{Factorial type: mca for qualitative active variables, pca for quantitative active variables. If not stated, it is chosen automatically: pca if all the active variables are quantitative, mca otherwise.}
  \item{xaxis, yaxis}{what is plotted on the horizontal and on the vertical axis. NULL (default) plots the first and the second factor respectively. A factorial dimension can be stated as a number or as a name such as "F3" or "PC2" (any of these prefixes is accepted whatever the type of analysis), so that yaxis=3 plots the first factor against the third one. The name of a variable of the data frame plots that observed variable instead: it must be quantitative or an ordered factor, since placing a nominal variable on an axis would impose an arbitrary order on its categories. Should a column of the data frame be named like a factorial axis, the column prevails. A variable stated here does not need to be included in variables, as it is added to them: it thus takes part in the complete-cases filter and appears in the node table as an illustrative variable. It may also be an active one, although a variable is then plotted against factors it has itself helped to build. Optionally, either argument accepts a second element, TRUE or FALSE, stating whether that axis is to be standardized (see scaleAxes), as in xaxis=c("age", TRUE) or yaxis=c("F1", FALSE); c() turns such a pair into a character vector, which is read back correctly, but list("age", TRUE) is accepted as well. Note that combining an observed axis with patterns=TRUE plots each pattern at the (weighted) mean of that variable over its cases, unless the variable is itself part of vPatterns.}
  \item{scaleAxes}{if TRUE (default), the plotted factorial axes are each expressed in standard-deviation units (mean 0, sd 1), equivalent to dividing by the square root of their eigenvalue. Both MCA and PCA coordinates already have mean 0, so this only rescales their spread, making the two axes comparable and independent of how much variance/inertia they capture. Points are not bounded to any fixed range: their spread only follows the usual properties of a standardized variable. FALSE keeps the raw MCA/PCA coordinates, whose scale grows with the eigenvalue of each axis. An axis holding an observed variable (see xaxis and yaxis) is not affected: it keeps its own units, which would otherwise be lost. Both defaults can be overridden axis by axis with the second element of xaxis and yaxis.}
  \item{nclusters}{number of clusters for k-means. It can be a vector, in which case one group variable is added per value.}
  \item{clusterOn}{what k-means clusters on: "factors" (default) uses a selection of the MCA/PCA factorial coordinates (see nfactors and critFactors). "variables" uses the active variables instead (their dichotomized indicators for mca, their standardized values for pca), so groups reflect the full active variable space rather than only the extracted factors. nfactors and critFactors are ignored in this case, with a warning if they are explicitly set.}
  \item{scaleClusters}{whether the columns k-means clusters on are each standardized (mean 0, sd 1) beforehand, so that every one of them has the same say in the distances. NULL (default) decides it by what is being clustered: TRUE for factorial coordinates, and also for the active variables of a pca, which come standardized anyway; FALSE for the dichotomized indicators of an mca (clusterOn="variables"). TRUE and FALSE force it either way.

Standardizing matters because k-means weighs each column by its own spread. On raw factorial coordinates a dominant first factor cuts the groups along itself alone, and they then read as bands on a plot whose axes are standardized (see scaleAxes), which is why TRUE suits them.

Dichotomized indicators are the opposite case. Standardizing one held by a proportion p of the cases turns a "yes" into sqrt((1-p)/p) standard deviations, so a category held by 2\% of the sample weighs seven times more than a balanced one, and squared distances add up over categories: the groups end up singling out the few cases holding a rare category. Left raw, the squared distance between two cases is the number of categories they differ in, which is the usual metric for qualitative data.

Note that FALSE on factorial coordinates reproduces the behaviour of versions previous to this argument.}
  \item{nstart}{number of random starts of the k-means algorithm, 25 by default. k-means only converges to a local optimum, the one its starting centres lead to, so it is run this many times and the best of those runs is kept. A single start (nstart=1, the behaviour of versions previous to this argument) reaches a clearly worse grouping in about one out of five runs, without any sign of it. Raising it costs proportionally more computing time, which is nonetheless a small part of the whole call.}
  \item{nfactors}{when clusterOn="factors", the maximum number of MCA/PCA factors considered for clustering (critFactors then selects how many of these are actually used, at least 1). It only affects clustering: the scattergram plot itself always uses the factors stated by xaxis and yaxis (the two first ones by default), regardless of clusterOn, nfactors or critFactors. Ignored (with a warning if set explicitly) when clusterOn="variables".}
  \item{critFactors}{when clusterOn="factors", a factor among the first nfactors is used for clustering only if its eigenvalue exceeds critFactors times the mean eigenvalue of the full factorial spectrum. This generalizes Kaiser's rule (eigenvalue above average) to MCA, whose eigenvalues are not on the same scale as PCA's: critFactors=1 reproduces the classic Kaiser criterion in (unweighted or weighted) PCA, since the mean eigenvalue of a correlation matrix is always 1, while in MCA it adapts to the actual inertia scale of the analysis. The default, 0, applies no filtering, so all nfactors factors are used. At least one factor is always used. Ignored (with a warning if set explicitly) when clusterOn="variables".}
  \item{sortClusters}{if TRUE (default), the clusters of each k-means solution are renumbered so that group 1 is the one whose (weighted) mean on the horizontal axis is the lowest, group 2 the next one, and so on. k-means numbers its clusters after its own initialization, which is arbitrary, so this is what gives a meaning to the order of the group columns, which are ordered factors. The axis it sorts on is the one actually plotted, that is, the first factor unless xaxis states an observed variable, so that the numbering always follows what the reader sees. It applies as well when clusterOn="variables". Only the numbering changes: the partition of the cases is exactly the same. FALSE keeps the numbering produced by k-means, as previous versions did for a given seed.}
  \item{columns}{order of the columns of the resulting node data frame. NULL (default) orders them as the vPatterns variables (only when vPatterns has been explicitly stated, since it defaults to every variable), the n column of pattern counts, the active variables, the remaining ones (n is placed right after the active variables instead when the vPatterns block is empty), the k-means group columns, the ntext column and the name column, after which \link{netCoin} appends the fx and fy coordinates of the scattergram. A character vector places the columns it names first, in the stated order, and keeps every remaining column behind them in the default order, so no column is ever dropped. Such a vector may also contain the reserved names _vpattern, _actives, _others and _nclusters, which are expanded, where they are stated, into the vPatterns variables, the active variables, the variables which are neither of them, and the k-means group columns respectively. A column named twice (by itself and within a block) is kept at its first position. Names which are not a column of the node data frame are ignored with a warning.}
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

# An observed variable on the horizontal axis, read against the second component,
# and the first component against the third one
surScat(iris, active=c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
        xaxis="Petal.Width", shape="Species")
surScat(iris, active=c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
        yaxis="PC3", shape="Species")
}

