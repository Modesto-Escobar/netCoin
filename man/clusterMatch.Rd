\name{clusterMatch}
\alias{clusterMatch}
\title{Compare the classifications of a scattergram, category by category.}
\description{
\code{clusterMatch} compares a model classification with the other classifications held by a
netCoin object from surScat, matching their categories one to one and reporting sensitivity,
precision, F1 and specificity.
}
\value{
A data.frame with one row per classification, sorted as \code{sort} states. It carries as
attributes the confusion matrix of every classification ("tables"), its measures category by
category ("byClass"), the matching of categories ("mapping"), the unit of analysis ("unit"),
the weights actually applied ("weights") and, where there were any, the patterns spanning
several groups ("spanning").
}

\usage{
clusterMatch(x, model, columns=NULL, weights=NULL, align=TRUE, weighted=TRUE,
             modal=TRUE, sort=NULL, digits=3)
}
\arguments{
\item{x}{a 'netCoin' object from surScat, or the data.frame of its nodes.}
\item{model}{the classification everything else is compared with: the name of a node column,
or a vector holding one value per node. Where the nodes are patterns, a vector holding one
value per case is taken as well, and the comparison then goes case by case.}
\item{columns}{the classifications to compare. NULL (default) takes those surScat recorded in
the "clusterColumns" attribute and, where there are none, every factor or character column of
the node table.}
\item{weights}{what every unit counts for. The name of a numeric node column, or a vector
holding one value per node or per case. NULL (default) weighs every node by the cases it
stands for when the nodes are patterns, and counts one per unit otherwise. \code{weights=1}
forces the plain count of nodes.}
\item{align}{whether to realign the categories of every classification with those of the model
before comparing them. TRUE by default; FALSE compares them as they come, which only makes
sense where the categories already mean the same thing in both.}
\item{weighted}{whether to add, next to the plain means of F1 and specificity, the means
weighted by the size of every category (.w).}
\item{modal}{where the nodes are patterns, whether to keep the modal group of each pattern.
See the details.}
\item{sort}{the column to sort the result by, from the largest down; a vector of names sorts
by all of them at once, the first as the main criterion. The prefix "-" turns one column
around: \code{c("-k", "F1")}. NULL (default) keeps the order the classifications came in.}
\item{digits}{decimal places the numeric columns are rounded to.}
}
\details{
The groups of a classification carry no name of their own: group 1 of one method need not be
group 1 of another, as the numbering is arbitrary. Before comparing, the categories of every
classification are therefore realigned with those of the model through the optimal assignment
(Hungarian algorithm), which maximizes the diagonal of the contingency table. That demands a
one-to-one correspondence, and \code{clusterMatch} compares therefore only the classifications
holding as many categories as the model: without a bijection these measures are not defined.
\code{\link{clusterAgree}} is there for the rest.

The measures follow the definitions of caret::confusionMatrix (one against the rest,
prediction down the rows and reference along the columns) but are computed here so as not to
depend on caret. Every one of them is averaged over the categories in two ways: the plain mean
treats them alike, and a small group badly classified sinks it; the one weighted by the size of
each group measures what happens to a case drawn at random. Sensitivity is given no weighted
mean, as it would coincide with \code{Accuracy}.

Where surScat has collapsed the cases into patterns, every node stands for several cases and
counting it once would give the rare pattern the same say as the common one. The nodes are
therefore weighed by the cases they stand for (the "caseToPattern" attribute, and "caseWeight"
where surScat was given weights), so that the tables come out in cases.

Pattern nodes bring a limitation of the object itself: surScat keeps in every cluster column
all the groups the pattern spans, as text ordered from the commonest down ("Group: 2|Group:
1"), and loses which group each case went to. The comparison is then made with the modal group
of every pattern, and how many patterns span more than one is reported: where they are many,
the figures are optimistic and it is worth comparing case by case, passing the model as a
vector of cases and classifications that have not been collapsed.
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
\seealso{
\code{\link{clusterAgree}}, \code{\link{surScat}}, \code{\link{addClusters}}
}
\examples{
\dontrun{
data(iris)
S <- surScat(iris, variables=names(iris)[1:4], type="pca", nclusters=2:4)

R <- clusterMatch(S, "kmeans(3)")
R
attr(R, "byClass")[[1]]    # measures category by category
attr(R, "tables")[[1]]     # the realigned confusion matrix
attr(R, "mapping")[[1]]    # which category was matched with which

# sorted by F1, and by the simplest classifications first
clusterMatch(S, "kmeans(3)", sort="F1")
clusterMatch(S, "kmeans(3)", sort=c("-k", "F1"))

# with an outside vector rather than a node column
clusterMatch(S, iris$Species)
}
}
