\name{d_netCorr}
\alias{d_netCorr}
\title{Dynamic networked correlations.}
\description{
\code{netCorr} produces a network object of dynamic correlations. Its input has to be at least one set of quantitative variables. 
}
\value{
The function creates a list object to create a multigraph and eventually a folder in the computer with an HTML document named index.html which contains the produced dynamic graphs. This file can be directly opened with your browser and sent to a web server to work properly.
}
\usage{
d_netCorr(variables, nodes = NULL, weight=NULL, 
        pairwise=FALSE,
        frequency=pairwise, means=TRUE,
        method="pearson", criteria="value", Bonferroni=FALSE,
        minL=0, maxL=Inf,
        sortL=NULL, decreasingL=TRUE,
        factorial=c("null", "pc", "nf", "vf", "of"),
        components=TRUE, backcomponents=FALSE,
        sequence=seq(.20, 1, .01), textFilter=c(1, .99), speed=50,
        dir=NULL, ...)
}
\arguments{
\item{variables}{a data frame with at least two quantitative variables.}
\item{nodes}{a data frame with at least one vector of names and other information from the nodes}
\item{weight}{weights of the cases: the name of a column of \code{variables} or a numeric vector. See \link{netCorr}.}
\item{pairwise}{Pairwise mode of handling missing values if TRUE. Listwise by default.}
\item{frequency}{a logical value true if the frequency (number of valid cases) of each variable is to be shown in the nodes. TRUE by default only if \code{pairwise=TRUE}.}
\item{means}{a logical value true if the mean, standard deviation, minimum and maximum of each variable are to be shown in the nodes, and the means give the size of the nodes. Default=TRUE.}
\item{method}{measures of association. Pearson correlation by default; see \link{netCorr} for the others.}
\item{criteria}{statistic to be use for selection criteria.}
\item{Bonferroni}{Bonferroni criterium of the signification test.}
\item{minL}{minimum value of the statistic to include the edge in the list.}
\item{maxL}{maximum value of the statistic to include the edge in the list.}
\item{sortL}{sort the list according to the values of a statistic. See below}
\item{decreasingL}{order in a decreasing way.}
\item{factorial}{factorial layout: Principal components (pc), factorial (nf), factorial with varimax rotation (vf), and factorial with oblimin rotation (of) are possible}
\item{components}{display following graph in sequence only when the components are unequal}
\item{backcomponents}{display following graph in sequence even if there are less components}
\item{sequence}{evolution of the dinamic graphs. c(first threshold, last threshold, step)}
\item{textFilter}{limits for showing the correlations as text on the graph. This limit is the minimun value of the first (absolute), and the second (threshold plus its value)}
\item{speed}{speed of the dynamic evolution from 0 to 100}
\item{dir}{a "character" string representing the directory where the web files will be saved.}
\item{...}{any \link{netCoin} argument.}
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
\examples{
# A character column (with separator)
data(iris)
d_netCorr(iris[,1:4],ltext="value", sequence=seq(.80, 1, .05), zoom=3, textFilter=.85,
 main="Correlations between measurements of Iris Species", components=FALSE,
 note="Anderson, Edgar (1935) y Fisher, R. A. (1936)") # network object
}

