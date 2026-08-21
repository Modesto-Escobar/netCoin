\title{Coefficient Table from Regression Graphs}
\name{coefTable}
\alias{coefTable}
\alias{print.coefTable}
\description{
Builds a side-by-side coefficient table from one or several objects returned by \link{glmCoin}, using only the information stored in them: links (estimates, standard errors, confidence intervals and p values) and the \code{equations} data frame (n and pseudo-R2 of each model). Each column block corresponds to one combination of object, equation and \code{Target}, so the same function covers hierarchical regressions (one response, nested sets of predictors), several responses in a single object, multinomial responses (one block per category) and any combination of these. It is the coefficient counterpart of \link{marginsTable}, which tabulates marginal effects.
}
\usage{
coefTable(G, labels = NULL, columns = c("b", "sig", "se"),
          blocks = c("entry", "variable", "none"),
          sig = c(0.05, 0.01, 0.001), digits = 3, digits.se = digits,
          footer = c("n_obs", "r2", "dr2"), note = FALSE,
          file = NULL, format = NULL)
}
\arguments{
  \item{G}{A netCoin object created by \link{glmCoin}, or a list of them. A named list uses its names as column labels.}
  \item{labels}{Optional character vector of column-block labels, one per model/target combination. If missing, labels are built automatically from the components that actually vary (model, response and, for multinomial responses, category), so that no redundant text is shown.}
  \item{columns}{Which columns form the block repeated for each model: "b" (estimate), "sig" (significance stars), "se" (standard error in parentheses), "ci" (confidence interval) and "p" (p value). At least one is required.}
  \item{blocks}{How predictor rows are grouped. "entry" (default) groups them by the equation in which each predictor first appears, which reproduces the theoretical blocks of a hierarchical regression; "variable" groups the categories of each explanatory variable, whose name heads the block; "none" shows a flat list. A named list of predictor names may also be given to define the groups explicitly.}
  \item{sig}{Numeric vector with the significance thresholds of the stars, sorted internally in decreasing order: one star for p below the largest threshold, two for the next, and so on (any number of thresholds is allowed, e.g. \code{sig=c(.1,.05,.01,.001)} for a four-star scale).}
  \item{digits}{Number of decimals for the coefficients.}
  \item{digits.se}{Number of decimals for standard errors and confidence limits. Defaults to \code{digits}.}
  \item{footer}{Character vector with the statistics to show as bottom rows. Besides any column of \code{G$equations} (\code{n}, \code{n_obs}, \code{r2}, \code{r2.mcfadden}, \code{r2.nagelkerke}, \code{logLik}, \code{aic}, \code{bic}, \code{deviance}, \code{null.deviance}, \code{df.residual}), three are computed here: \code{k} (number of predictors), \code{dr2} (increment in R2 over the reference model) and \code{ref} (the reference model itself). The reference is found automatically as the largest model strictly nested in each equation and fitted on the same response. Categories of a multinomial response share a single equation, so its statistics are shown once, under the first of its blocks. Use \code{NULL} to omit the footer.}
  \item{note}{If TRUE, a legend with the meaning of the stars (built from \code{sig}) is added below the table, both when printing and in the exported formats.}
  \item{file}{Optional file name to export the table. The format is inferred from the extension (.txt, .tex or .csv) unless \code{format} is given.}
  \item{format}{Export format: "text", "latex" or "csv" (the latter opens directly in Excel). If \code{file} is missing, a default name is used.}
}
\details{
Rows are the union of the predictors of all the equations, in order of first appearance; an empty cell means that the predictor was not included in that equation, which is not the same as a zero coefficient.

Because R2 and its increments are only comparable when the equations are fitted on the same cases, a warning is issued whenever the models differ in their number of observations, and the fact is recorded in the note when \code{note=TRUE}. This typically happens when a predictor added in a later step has missing values.
}
\value{
A character matrix of class \code{coefTable}, printed with grouped headers (one block of columns per model, response or category). The bottom rows show the statistics selected in \code{footer}. When exporting, the table is returned invisibly.
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
\seealso{\link{glmCoin}, \link{marginsTable}}
\examples{
data(birthwt, package="MASS")
birthwt$smoke <- factor(birthwt$smoke, labels=c("no","yes"))
birthwt$race  <- factor(birthwt$race,  labels=c("white","black","other"))

## Hierarchical regression: one response, nested sets of predictors.
## Blocks are deduced from the equation in which each predictor enters.
m1 <- glmCoin("bwt ~ age + lwt", birthwt)
m2 <- glmCoin("bwt ~ age + lwt + smoke", birthwt)
m3 <- glmCoin("bwt ~ age + lwt + smoke + race", birthwt)
coefTable(list(m1, m2, m3), note=TRUE)

## Confidence intervals instead of standard errors, with fewer decimals
coefTable(list(m1, m2, m3), columns=c("b","sig","ci"), digits=2, digits.se=1)

## Several dependent variables in a single glmCoin object
g <- glmCoin(c("bwt ~ age + lwt", "low ~ age + lwt, binomial"), birthwt)
coefTable(g, footer=c("n_obs","r2"))

## Multinomial response: one block of columns per category
data(housing, package="MASS")
h <- glmCoin("Sat ~ Infl + Type + Cont, multinomial", housing, weights="Freq")
coefTable(h, blocks="variable", columns=c("b","sig"))
}
