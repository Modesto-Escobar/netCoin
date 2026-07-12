\name{glmCoin}
\alias{glmCoin}
\title{Regression Graphs}
\description{
produces a netCoin object from a set of glm regressions.
}
\usage{
glmCoin(formulas, data, weights=NULL, pmax=.05,
        robust=FALSE, showArrows=TRUE,
        color="variable", lwidth="s.value", lcolor="estimate",
        language=c("en","es","ca"),
        igraph=FALSE, ...)
}
\arguments{
  \item{formulas}{
A set of formulas separated, folowed by the family and a return.
For example:
model <- "counts ~ outcome + treatment, poisson
          counts ~ outcome, poisson"
The possible families are \emph{gaussian} (default), \emph{poisson}, \emph{binomial} and \emph{multinomial}.}
  \item{data}{Data frame containing the variables in the model.}
  \item{weights}{Optional vector of weights to be used in the fitting process.}
  \item{pmax}{Selection of links with Pr(>|z|) less than p (one-tail by default).}
  \item{robust}{Robust estimation of standard errors if TRUE.}
  \item{showArrows}{a logical value true if the directional arrows are to be shown. Default = FALSE.}
  \item{color}{Nodes' attribute to be used for expressing color ("variable" by default).}
  \item{lwidth}{Nodes' attribute to be used for widht of arrows ("s.value" by default).}
  \item{lcolor}{Nodes' attribute to be used for color of arrows ("estimate" by default).}
  \item{language}{Language of the graph controls.}
  \item{igraph}{Produces an igraph object instead of a netCoin object if TRUE.}
  \item{...}{Any \link{netCoin} argument.}
}
\value{
This function creates a netCoin object (or igraph) and, if stated, a folder in the computer with an HTML document named index.html which contains the produced graph. This file can be directly opened with your browser and sent to a web server to work properly.

The links of the graph carry the average marginal effects as GW contrasts (\code{estimate}, \code{std.error}, \code{p.value}: effects relative to the weighted overall mean) and also relative to the first category of each qualitative variable (\code{estimate.base}, \code{std.error.base}, \code{p.value.base}; for the base category itself the estimate is 0 and its standard error and p value are NA). Both sets can be displayed with \link{marginsTable}.

The netCoin object also includes a data frame named \code{equations} with one row per declared equation and the following goodness-of-fit measures: \code{model} (formula), \code{family} (declared family), \code{n} (number of observations; the sum of weights in weighted models), \code{n_obs} (number of observed cases, regardless of weights), \code{r2} (deviance-based R2: 1-deviance/null deviance, which equals the classical R2 for gaussian models), \code{r2.mcfadden} (McFadden's pseudo-R2, NA for quasi families), \code{r2.nagelkerke} (Nagelkerke's pseudo-R2), \code{logLik}, \code{aic}, \code{bic} (NA for quasi families), \code{deviance}, \code{null.deviance} and \code{df.residual}. When \code{igraph=TRUE}, this data frame is attached as the graph attribute \code{equations}.
}

\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}

\examples{
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
Dobson <- data.frame(counts=counts, outcome=outcome, treatment=treatment)
model <- "counts ~ outcome + treatment, poisson"
glmCoin(model,Dobson)

## Multinomial model: housing satisfaction in Copenhagen (Madsen, 1976)
data(housing, package="MASS")
g <- glmCoin("Sat ~ Infl + Type + Cont, multinomial", housing, weights="Freq")
g$equations
}
