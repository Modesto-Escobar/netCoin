\name{marginsTable}
\alias{marginsTable}
\alias{print.marginsTable}
\title{Marginal Effects Table from Regression Graphs}
\description{
Builds a marginal-effects table from the object returned by \link{glmCoin}, using only the information stored in it: nodes (frequency and percentage of each category), links (average marginal effects, standard errors and p values) and the \code{equations} data frame (n and pseudo-R2 of each model). Supports several dependent variables (including multinomial responses), each producing its own block of columns.
}
\usage{
marginsTable(G, by = c("category", "variable"),
             criterion = c("p", "effect"), contrast = c("gw", "base"),
             columns = c("n", "marg", "effect", "se", "sig"),
             sig = c(0.05, 0.01, 0.001), note = FALSE,
             digits = 1, footer = c("n", "r2"), file = NULL, format = NULL)
}
\arguments{
  \item{G}{A netCoin object created by \link{glmCoin}.}
  \item{by}{"category" produces one row per category of the explanatory variables; "variable" collapses to one row per variable, using its representative category.}
  \item{criterion}{How the representative category is chosen when \code{by="variable"}: "p" (smallest p value) or "effect" (largest absolute effect).}
  \item{contrast}{"gw" (default) shows marginal effects relative to the weighted overall mean (GW contrasts); "base" shows them relative to the first (base) category of each qualitative variable, whose row displays \code{ref.} (slopes of numeric variables are the same under both). Requires the \code{estimate.base} columns stored in \code{G$links} by the current version of \link{glmCoin}.}
  \item{columns}{Which columns to display (all by default): "n" is the global frequency column, and "marg" (marginal percentage), "effect" (marginal effect), "se" (standard error) and "sig" (significance stars) form the block repeated for each dependent variable. At least one block column is required.}
  \item{sig}{Numeric vector with the significance thresholds of the stars, sorted internally in decreasing order: one star for p below the largest threshold, two for the next, and so on (any number of thresholds is allowed, e.g. \code{sig=c(.1,.05,.01,.001)} for a four-star scale).}
  \item{note}{If TRUE, a legend with the meaning of the stars (built from \code{sig}) is added below the table, both when printing and in the exported formats.}
  \item{digits}{Number of decimals for percentages and effects.}
  \item{footer}{Character vector with the columns of \code{G$equations} to show as bottom rows, one value per model (multinomial responses share a single value across their category blocks). Available columns: \code{n}, \code{n_obs}, \code{r2}, \code{r2.mcfadden}, \code{r2.nagelkerke}, \code{logLik}, \code{aic}, \code{bic}, \code{deviance}, \code{null.deviance}, \code{df.residual} and \code{family}. Use \code{NULL} to omit the footer.}
  \item{file}{Optional file name to export the table. The format is inferred from the extension (.txt, .tex or .csv) unless \code{format} is given.}
  \item{format}{Export format: "text", "latex" or "csv" (the latter opens directly in Excel). If \code{file} is missing, a default name is used.}
}
\value{
A character matrix of class \code{marginsTable}, printed with grouped headers (one block of columns per dependent variable: marginal percentage, marginal effect, standard error and significance stars). The bottom rows show the statistics selected in \code{footer} (by default the n and the pseudo-R2 of each equation), taken from \code{G$equations}. When exporting, the table is returned invisibly.
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
\seealso{\link{glmCoin}}
\examples{
data(birthwt, package="MASS")
birthwt$smoke <- factor(birthwt$smoke, labels=c("no","yes"))
birthwt$race  <- factor(birthwt$race,  labels=c("white","black","other"))
g <- glmCoin("low ~ smoke + race, binomial", birthwt)
marginsTable(g)
marginsTable(g, by="variable")

## Poisson model: Dobson (1990) Page 93, Randomized Controlled Trial
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
Dobson <- data.frame(counts=counts, outcome=outcome, treatment=treatment)
gp <- glmCoin("counts ~ outcome + treatment, poisson", Dobson)
marginsTable(gp, contrast="base", note=TRUE)
}
