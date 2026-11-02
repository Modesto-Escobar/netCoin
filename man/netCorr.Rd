\name{netCorr}
\alias{netCorr}
\title{Networked correlations.}
\description{
\code{netCorr} produces a network object of correlations. Its input has to be at least one set of quantitative variables.
}
\value{
The function creates a netCoin object and eventually a folder in the computer with an HTML document named index.html which contains the produced graph. This file can be directly opened with your browser and sent to a web server to work properly.

The links hold the \code{value} of the association, the statistic \code{z} of its test and its probability \code{p}, all of them of the first method. With \code{pairwise=TRUE} and missing values they also hold \code{n}, the number of cases of each pair. The remaining methods add their coefficient (named after the method), \code{z(method)} and \code{p(method)}.
}
\usage{
netCorr(variables, weight=NULL, pairwise=FALSE,
        sort=FALSE, decreasing=TRUE,
        frequency=pairwise, means=TRUE,
        method="pearson",
        criteria="p", Bonferroni=FALSE, minL=0, maxL=Inf,
        sortL=NULL, decreasingL=TRUE,
        igraph=FALSE, ...)
}
\arguments{
\item{variables}{a data frame with at least two quantitative variables.}
\item{weight}{weights of the cases: the name of a column of \code{variables} (which is then left out of the analysis) or a numeric vector with a weight per case. Pearson, Spearman and partial correlations are weighted; the tests use the effective sample size of Kish, \eqn{(\sum w)^2 / \sum w^2}. Kendall and dcor do not admit weights.}
\item{pairwise}{Pairwise mode of handling missing values if TRUE. Listwise by default. The tests of each pair rest on the cases of that pair. Partial correlations always use the cases complete in all the variables.}
\item{sort}{if TRUE, the variables (nodes, and so the rows and columns of the correlation matrix) are sorted by frequency, that is, by their number of valid cases, and by mean when the frequencies are equal. Default=FALSE.}
\item{decreasing}{decreasing or increasing sort of the variables. Default=TRUE.}
\item{frequency}{a logical value true if the frequency (number of valid cases) of each variable is to be shown in the nodes. By default it is TRUE only if \code{pairwise=TRUE}, the only case in which it can differ among the variables.}
\item{means}{a logical value true if the mean, standard deviation, minimum and maximum of each variable are to be shown in the nodes. Default=TRUE.}
\item{method}{one or several measures of association: \code{"pearson"}, \code{"spearman"}, \code{"kendall"}, \code{"partial"} (partial correlation of each pair controlling for all the other variables) and \code{"dcor"} (bias-corrected distance correlation of Szekely and Rizzo, which also detects non-linear dependence). The first one gives the \code{value}, \code{z} and \code{p} of the links; the others are added as further columns. See details.}
\item{criteria}{statistic to be used for selection criteria: \code{"p"}, \code{"z"}, \code{"value"}, the name of a method (its coefficient), or \code{"z(method)"} and \code{"p(method)"}. \code{"p"}, \code{"z"} and \code{"value"} refer to the first method.}
\item{Bonferroni}{Bonferroni criterium of the signification test: when \code{criteria} is a probability, \code{maxL} is divided by the number of contrasts, \code{choose(k, 2)} for \code{k} variables.}
\item{minL}{minimum value of the statistic to include the edge in the list. With bipolar links (the default, see \code{...}) the absolute value of \code{value}, \code{z} and of the coefficients is compared, so \code{minL=.3} keeps the correlations with |r| >= .3.}
\item{maxL}{maximum value of the statistic to include the edge in the list. When the criterium is a probability and it is not set, it is 0.5 (1 for bipolar links).}
\item{sortL}{sort the list according to the values of a statistic. See below}
\item{decreasingL}{order in a decreasing way.}
\item{igraph}{Produces an igraph object instead of a netCoin object if TRUE}
\item{...}{Any \link{netCoin} argument. Unless they are given, \code{linkBipolar=TRUE} (links of negative correlation are kept, and tested two-sided), \code{lcolor="value"}, \code{lwidth="value"} and \code{lweight="value"}.}
}
\details{
The test of each association is two-sided by default (\code{linkBipolar=TRUE}), and one-sided (positive association) with \code{linkBipolar=FALSE}, when the links of negative association are also left out by the default \code{minL=0}.

\itemize{
\item Pearson and partial: \eqn{t = r \sqrt{df/(1-r^2)}} with \eqn{n-2} degrees of freedom (\eqn{n-k} for partial correlations of \eqn{k} variables), \eqn{n} being the number of cases of the pair.
\item Spearman: the same \eqn{t} (\code{z}) with \eqn{n-2} degrees of freedom, as \code{\link[stats]{cor.test}} does with \code{exact=FALSE}; the exact probability is used for \eqn{n \le 9} without ties.
\item Kendall: \code{z} is S divided by its standard deviation, corrected for ties, and \code{p} comes from the normal distribution; the exact probability is used for \eqn{n < 50} without ties.
\item dcor: the t-test of independence of Szekely and Rizzo (2013).
}
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
\examples{
# A character column (with separator)
data(iris)
netCorr(iris[,1:4],ltext="value",
 main="Correlations between measurements of Iris Species",
 note="Anderson, Edgar (1935) y Fisher, R. A. (1936)") # network object

# Several measures at once, selecting by the probability of Kendall's tau
netCorr(iris[,1:4], method=c("pearson", "kendall", "dcor"), criteria="p(kendall)", maxL=.01)

# Links of absolute correlation of .7 or more
netCorr(iris[,1:4], criteria="value", minL=.7)

# Only positive correlations, tested one-sided
netCorr(iris[,1:4], linkBipolar=FALSE)
}

