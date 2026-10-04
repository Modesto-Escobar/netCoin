\name{surScat}
\alias{surScat}
\title{Interactive scattergram of the cases of a data frame.}
\description{
\code{surScat} produces an interactive scattergram of the cases (rows) of a data frame. Its axes are, by default, the two first dimensions of a multiple correspondence analysis (MCA) of the active qualitative variables, or of a principal component analysis (PCA) of the active quantitative variables. Cases are colored by groups obtained by k-means clustering of these coordinates.

Either axis, or both, may be changed with xaxis and yaxis: a different factorial dimension, or an observed variable of the data frame. Clustering is not affected by this choice, so a plot of an observed variable against a factor still colors its cases according to the whole structure found in the active variables.

A plotted factorial axis is labelled with the percentage of variance (PCA) or of inertia (MCA) it accounts for, computed on the full factorial spectrum. In MCA the raw inertia percentages are known to be strongly pessimistic, so Benzecri's adjustment is applied: only the factors whose eigenvalue exceeds 1/Q (Q being the number of active variables) are retained, their eigenvalues are transformed into (Q/(Q-1))^2 (eigenvalue - 1/Q)^2, and the percentages are computed on this adjusted inertia. The adjustment is skipped, and the raw inertia percentages are kept, when the second factor does not itself exceed the 1/Q threshold, since it would then be labelled as 0 percent. An axis holding an observed variable is labelled with the label of that variable, or with its name if it has none. These labels are passed to \link{netCoin} as its axesLabels argument, unless another axesLabels is stated.

Every node carries an HTML description, in a column named ntext which is passed to \link{netCoin} as its ntext argument, unless another ntext is stated. When cases are collapsed it describes the pattern (see patterns). Otherwise it describes the case itself: its value on each of the two axes, followed in brackets by the mean of the k-means group it belongs to on that same axis, and the number of that group (the last k-means solution, when nclusters states several). A factorial axis is described by the coordinate which is plotted, whereas an observed variable is described in its own units and with its own labels, even if its axis has been standardized.

The resulting object also holds \code{$clusters}, the (weighted) centroid of every group of every k-means solution, on the plane just drawn: a flat list with one matrix per clusterization, its name joining the clusterization and the plane with a dot, as in \code{S$clusters[["Groups(2).PC1-PC2"]]}, and one row per group. When the nodes are patterns whose cases fall in several groups (see vPatterns), each pattern takes part in the centroid of every group it spans with the (weighted) cases it holds of that group, as recorded in the attribute "clusterShares" (see \link{addClusters}). \link{addAxes} extends it with an entry for every clusterization on whichever plane it adds, and \link{addClusters} with an entry on every plane the object already holds for whichever clusterization it adds, so that \code{$clusters} always holds one entry per clusterization-plane pair. Printing it, on its own or as part of the object, only states how many clusterizations and planes it covers.

It also holds \code{$clusterMeans}, the (weighted) mean of every group of every clusterization on every numeric column of the node table but fx and fy: a list with one data frame per clusterization, as in \code{S$clusterMeans[["Groups(2)"]]}, with one row per group, in the same order as the rows of its entries in \code{$clusters}, and one column per variable. These are the centroids of the groups when the graph places a variable on an axis rather than a coordinate of a plane: a centroid being the mean on each axis separately, any two of its columns give the centroid on the plane they span, and either of them can also be paired with a column of an entry of \code{$clusters}. They are computed just as \code{$clusters}, patterns spanning several groups included, and kept up to date by \link{addClusters} and \link{replaceClusters}.

Further planes and clusterizations can be obtained from \link{looking4clusters} with the arguments layouts and clusters. They are added to the plane and to the k-means groups that surScat draws, never put in their place: the main plane is still the one stated by type, xaxis and yaxis, and the selector of the graph offers the others next to it. looking4clusters works on the same cases as the k-means of surScat, one row each, but on the active variables as they come: the quantitative ones unstandardized (so its "pca" differs from the PCA of surScat, which standardizes them) and the qualitative ones as a matrix of 0/1 indicators. When the nodes are patterns, its results are collapsed as \link{addAxes} and \link{addClusters} collapse them; when maxN has drawn a sample, only the cases of the nodes drawn are kept. Case weights are not taken into account by looking4clusters, whose methods accept none. The same holds for the latent models clusters may also ask for: a latent profile analysis (LPA) when the active variables are quantitative, and a latent class analysis (LCA) when they are qualitative, which take the active variables as they come too, and whose classes are added in the same way. Once more than one plane is held, the axes are labelled X axis and Y axis, unless axesLabels is stated, since a single pair of labels has to serve for every plane.
}
\value{
This function creates a netCoin object and, if stated, a folder in the computer with an HTML document named index.html which contains the produced graph. This file can be directly opened with your browser and sent to a web server to work properly.

When nclusters holds "auto", the object carries an attribute "clusterIndex": a data frame with one row per number of groups weighed, holding that number (k), the value of the index (in a column named after clusterIndex) and whether it was the one chosen (optimal).
}

\usage{
surScat(data, variables=names(data), active=variables, weight=NULL,
        patterns=FALSE, vPatterns=variables, jitter=0, type=c("mca", "pca"),
        xaxis=NULL, yaxis=NULL, scaleAxes=TRUE, nclusters=2,
        clusterOn=c("factors", "variables"), scaleClusters=NULL, nstart=25,
        nfactors=2, critFactors=0, sortClusters=TRUE, columns=NULL,
        seed=2020, maxN=2000, layouts=NULL, clusters=NULL, suffix=NULL,
        force_execution=FALSE, clusterIndex=c("ch", "silhouette", "bic"), ...)
}
\arguments{
  \item{data}{a data frame.}
  \item{variables}{a vector of variables included in the previous data frame. Only the cases complete over all of them (and over weight) are analysed, and it defaults to every column of the data frame, which on a survey with filtered questions may leave few cases or none: it is therefore best stated. When less than half the cases are complete, a warning names the variables missing the most; when fewer than two are, it is an error.}
  \item{active}{a vector of active variables used to compute the coordinates.}
  \item{weight}{a numeric vector of case weights, or the name of a column of the data frame containing them (this column is then excluded from variables). Weights are used to compute the MCA or PCA coordinates (weighted category masses in MCA, weighted correlations in PCA), so the axes represent the weighted structure of the data. The k-means clustering is then computed on these weighted coordinates without further weighting. When patterns collapses cases (see patterns and vPatterns), weights are also used for the mean/mode summaries and for the case count of each pattern. NULL (default) performs an unweighted analysis.}
  \item{patterns}{if TRUE, cases sharing the same response pattern in vPatterns are collapsed into a single node, whose size represents the number of supporting cases (or the sum of their weights): a column n is added to the nodes and used as size argument of \link{netCoin} (unless another size is stated). In addition, a column ntext with an HTML description of each pattern (the values of its vPatterns variables plus its number of cases) is generated and passed as ntext argument of \link{netCoin}, unless another ntext is stated. Coordinates, k-means clustering and the eigenvalue-based factor selection are always computed case by case, exactly as when patterns=FALSE; collapsing only affects how the resulting node table is built afterwards: each pattern is plotted at the (weighted) mean position of its cases, variables outside vPatterns take their (weighted) mean if quantitative or their (weighted) mode if qualitative, and the k-means group columns list every distinct group reached by the pattern's cases, separated by "|" (a single value if all its cases share the same group). Recommended for MCA with few qualitative variables, where cases pile up on a small number of distinct positions. Automatically set to TRUE if vPatterns is stated.}
  \item{vPatterns}{a subset of variables defining the response pattern used to collapse cases when patterns=TRUE (or is implied by stating vPatterns). Defaults to variables, i.e. cases are collapsed only if they agree on every variable. Stating a subset collapses cases sharing a pattern on just those variables: the remaining variables (including active ones) and the k-means group columns are then summarized as described under patterns, and the ntext description is restricted to just the vPatterns variables.}
  \item{jitter}{amount of random noise added to the coordinates before plotting, expressed as a fraction of the range of each axis (e.g. .02 = 2 percent; TRUE is a shortcut for .02). It only spreads the points visually: clusters and coordinates are computed on the exact positions. It can be combined with patterns to separate patterns differing only in non active variables, and it is particularly useful when an axis holds a discrete observed variable (see xaxis and yaxis), whose few distinct values pile the cases up into columns.}
  \item{type}{Factorial type: mca for qualitative active variables, pca for quantitative active variables. If not stated, it is chosen automatically: pca if all the active variables are quantitative, mca otherwise.}
  \item{xaxis, yaxis}{what is plotted on the horizontal and on the vertical axis. NULL (default) plots the first and the second factor respectively. A factorial dimension can be stated as a number or as a name such as "F3" or "PC2" (any of these prefixes is accepted whatever the type of analysis), so that yaxis=3 plots the first factor against the third one. The name of a variable of the data frame plots that observed variable instead: it must be quantitative or an ordered factor, since placing a nominal variable on an axis would impose an arbitrary order on its categories. Should a column of the data frame be named like a factorial axis, the column prevails. A variable stated here does not need to be included in variables, as it is added to them: it thus takes part in the complete-cases filter and appears in the node table as an illustrative variable. It may also be an active one, although a variable is then plotted against factors it has itself helped to build. Optionally, either argument accepts a second element, TRUE or FALSE, stating whether that axis is to be standardized (see scaleAxes), as in xaxis=c("age", TRUE) or yaxis=c("F1", FALSE); c() turns such a pair into a character vector, which is read back correctly, but list("age", TRUE) is accepted as well. Note that combining an observed axis with patterns=TRUE plots each pattern at the (weighted) mean of that variable over its cases, unless the variable is itself part of vPatterns.}
  \item{scaleAxes}{if TRUE (default), the plotted factorial axes are each expressed in standard-deviation units (mean 0, sd 1), equivalent to dividing by the square root of their eigenvalue. Both MCA and PCA coordinates already have mean 0, so this only rescales their spread, making the two axes comparable and independent of how much variance/inertia they capture. Points are not bounded to any fixed range: their spread only follows the usual properties of a standardized variable. FALSE keeps the raw MCA/PCA coordinates, whose scale grows with the eigenvalue of each axis. An axis holding an observed variable (see xaxis and yaxis) is not affected: it keeps its own units, which would otherwise be lost. Both defaults can be overridden axis by axis with the second element of xaxis and yaxis.}
  \item{nclusters}{number of clusters for k-means. It can be a vector, in which case one group variable is added per value, the last one being the main solution: the one that colors the graph and describes each case. It may also hold "auto", for the number of groups to be chosen by clusterIndex. On its own, nclusters="auto" weighs from 2 to 10 groups and keeps only the chosen solution; next to a range, as in nclusters=c("auto", 2:6), it weighs the numbers of that range, keeps a group variable for each of them and places the chosen one last, so that it becomes the main solution. Numbers of groups reaching the number of distinct cases are left out. The values of the index for every number of groups weighed are kept in the attribute "clusterIndex" of the resulting object (see Value).}
  \item{clusterIndex}{the index by which nclusters="auto" chooses the number of groups: the one with the highest value of "ch" or "silhouette", or with the lowest value of "bic". "ch" (default) is the Calinski-Harabasz index, the between-group variance over the within-group variance, each divided by its degrees of freedom; it comes straight out of k-means, so it adds no computing time. "silhouette" is the mean silhouette width of the cases (package cluster), which compares the distance of every case to its own group with its distance to the nearest other group; it needs every distance between cases, so beyond 2000 cases it is computed on a random sample of 2000 of them. 

Both judge the groups on all the active variables, whatever space k-means clusters on (see clusterOn): standardized for type="pca", and as the 0/1 indicators of their categories for type="mca", where the squared distance between two cases is the number of categories in which they differ. That is, k-means still makes its groups on the factors by default, but the number of groups is the one that best separates the cases on everything that was measured, not only on the plane. Judged on the factors themselves, the indices would mislead: a few factors make a low-dimensional space which, short of clear-cut groups, can be cut ever finer for an ever better index, so that the largest number of groups weighed would nearly always be chosen; with qualitative variables this is worse still, as the cases pile up on a few points of the plane. A cut that only splits the plane gains nothing on the whole of the active variables. Mind that this judge tends, on the contrary, to prefer few groups, often 2, when the data hold no clear structure. With clusterOn="variables" k-means clusters and the indices judge on the same variables. Neither index takes case weights into account, as k-means does not either.

 "bic" (or "BIC") is the Bayesian information criterion of the latent model that suits the active variables (see clusters): a latent profile analysis for type="pca", estimated by package mclust, or a latent class analysis for type="mca", estimated by package poLCA. It is the usual criterion for the number of classes of such models and it may also choose a single group. The number it chooses is the one of every method, the k-means of surScat included. Should the package be missing, or the model fail for every number of groups, Calinski-Harabasz is used instead, with a warning, being the index of k-means closest to it: both rest on the within-group variance, which governs the likelihood of a mixture of spherical normal distributions, of which k-means is the limiting case. Ignored unless nclusters holds "auto".}
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
  \item{seed}{random seed for reproducible k-means clustering, row sampling and jitter, and also for the methods of \link{looking4clusters} (see layouts and clusters), each of which starts from it. The RNG state of the session is restored on exit. Use NULL to avoid setting a seed.}
  \item{maxN}{maximum number of nodes to draw. A scattergram of more nodes than this, whether cases or patterns, is capped by drawing maxN of them at random, reproducibly for a given seed. The node table then holds that sample and nothing else, so a clusterization or a pair of coordinates computed on the whole data no longer lines up with it: \link{addClusters} and \link{addAxes} report the mismatch rather than guess, and the way to add either is to compute it on the node table of the resulting object. Raise maxN above the number of cases, or patterns, to keep them all.}
  \item{layouts}{dimensionality reductions computed by \link{looking4clusters} and added as further planes: any of "pca", "tsne", "cmds" (or "mds"), "nmf" and "umap", in the order they are to be offered. TRUE asks for all of them. NULL (default) adds none. surScat draws no links, so the network layouts accepted by \link{netCoin} ("fo", "la"...) do not apply to it: they are ignored with a warning. tsne, nmf and umap need the packages Rtsne, NMF and uwot respectively; nmf only accepts non-negative data, so it suits the 0/1 indicators of qualitative variables, but not quantitative variables holding negative values. Unlike the factorial axes, those of tsne and umap have no meaning of their own: only the proximity between points can be read.}
  \item{clusters}{clustering methods whose groups are added as further columns: any of "pam", "hclust" and "kmeans" of \link{looking4clusters}, and the latent models "lpa" and "lca". TRUE asks for all of the former and for the latent model that suits the active variables, provided its package is installed: a latent profile analysis (package mclust, with equal variances and no covariances, model "EEI") for type="pca", or a latent class analysis (package poLCA, estimated 5 times over from random starts, the best being kept) for type="mca". Stated by name, a latent model that does not suit the active variables, or whose package is missing, is left out with a warning; so is a number of groups for which it cannot be estimated. Their columns are named in capitals, as in "LPA(3)" or "LCA(3)". The k-means of looking4clusters is not a copy of surScat's own: it works on the active variables as they come (see the Description), not on the standardized factorial coordinates, and it is started once, whereas surScat starts its own nstart times. Only the numbers of groups stated by nclusters are kept, in the same order, and each column is named after its method and its number of groups, as in "pam(3)". pam and hclust need the packages fpc, cluster and dendextend. NULL (default) adds none.}
  \item{suffix}{a string. Every plane added by layouts is also offered aligned by Procrustes rotation to the main plane, under its own name followed by this suffix, so that both its original and its aligned version can be seen (see \link{layoutsAlign}). NULL (default) adds no aligned plane.}
  \item{force_execution}{looking4clusters leaves out pam, hclust and cmds when there are more than 5000 cases, since they need the whole matrix of distances between them, and says so with a warning. TRUE computes them anyway.}
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

# The number of groups chosen by the Calinski-Harabasz index among 2 to 6, keeping
# every solution, and the values of the index
S <- surScat(iris, active=c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
             shape="Species", nclusters=c("auto", 2:6))
S
attr(S, "clusterIndex")

# Further planes and clusterizations from looking4clusters: the PCA of the raw (not
# standardized) variables and a UMAP, each also aligned to the main plane (suffix),
# and the kmeans, pam and hclust groups of looking4clusters for the same numbers of
# groups as the k-means of surScat
\donttest{
if(requireNamespace("uwot", quietly=TRUE) && requireNamespace("fpc", quietly=TRUE) &&
   requireNamespace("cluster", quietly=TRUE) && requireNamespace("dendextend", quietly=TRUE)) {
  surScat(iris, active=c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
          shape="Species", nclusters=2:3, layouts=c("pca", "umap"), clusters=TRUE,
          suffix="'")
}

# Qualitative active variables: looking4clusters works on their 0/1 indicators,
# which suit nmf. clusters may also name only some of its methods.
data(ess)
if(requireNamespace("Rtsne", quietly=TRUE) && requireNamespace("NMF", quietly=TRUE) &&
   requireNamespace("fpc", quietly=TRUE) && requireNamespace("cluster", quietly=TRUE) &&
   requireNamespace("dendextend", quietly=TRUE)) {
  surScat(ess, variables=c("Gender", "Age", "Social participation",
                           "Political participation"),
          weight="cweight", patterns=TRUE, nclusters=3,
          layouts=c("tsne", "nmf", "mds"), clusters=c("hclust", "kmeans"), suffix="*")
}

# Latent profiles (LPA) of the quantitative active variables, next to the k-means
# groups, for 2 to 4 profiles
if(requireNamespace("mclust", quietly=TRUE)) {
  P <- surScat(iris, active=c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
               shape="Species", nclusters=2:4, clusters="lpa")
  P
  table(P$nodes$Species, P$nodes[["LPA(3)"]])
}

# Latent classes (LCA) of the qualitative active variables, their number chosen by the
# BIC of the latent class model among 1 to 4, which then sets that of k-means too
if(requireNamespace("poLCA", quietly=TRUE)) {
  L <- surScat(ess, variables=c("Gender", "Age", "Social participation",
                                "Political participation"),
               weight="cweight", patterns=TRUE, nclusters=c("auto", 1:4),
               clusterIndex="bic", clusters="lca")
  L
  attr(L, "clusterIndex")
}
}
}

