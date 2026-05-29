\name{netCoin-package}
\alias{netCoin-package}
\docType{package}
\title{The netCoin package.}
\description{
netCoin combines the data analysis capabilities of R with the interactive visualization libraries of JavaScript to create \bold{networks of coincidences, co-occurrences, correlations and regressions} that can be explored directly in a browser, embedded in HTML files or integrated into Shiny applications. It also generates \bold{networked HTML galleries} and other \bold{plots of coincidences}.
}
\details{
Coincidence analysis detects what events, characters, objects, attributes, or characteristics tend to occur together within certain limits.

These given limits are call scenarios (\eqn{S}) and are considered to be the units of analysis, and as such they have to be placed in the rows of a matrix or data.frame.

In each \eqn{i} scenario, a series of \eqn{J} events \eqn{X_j}, which are to be represented as dichotomous variables \eqn{X_{j}} in columns, may occur (1) or may not occur (0). Scenarios and events constitute an incidence matrix (\bold{I}).

\bold{Incidence matrix}

\tabular{lrrrrr}{
          \tab \eqn{X_1} \tab \eqn{X_2} \tab \eqn{X_3} \tab  \eqn{\;} \dots \tab \eqn{X_J} \cr
\eqn{S_1} \tab  0        \tab 1         \tab 0         \tab  \dots \tab 1\cr 
\eqn{S_2} \tab  1        \tab 0         \tab 1         \tab  \dots \tab 0\cr 
\dots     \tab  \dots    \tab \dots     \tab \dots     \tab  \dots \tab \dots\cr 
\eqn{S_n} \tab  1        \tab 1         \tab 0         \tab  \dots \tab 1\cr 
}
From this incidences matrix, a coincidence (\bold{C}) matrix can be obtained with the function \code{\link{coin}}. In this matrix the main diagonal represents frequencies of \eqn{X_j}, while the others elements are number of coincidences between two events.

\bold{Coincidence matrix}

\tabular{lrrrrr}{
          \tab \eqn{X_1} \tab \eqn{X_2} \tab \eqn{X_3} \tab  \eqn{\;} \tab \eqn{X_J} \cr
\eqn{X_1} \tab  2        \tab 1         \tab 1         \tab  \dots \tab 1\cr 
\eqn{X_2} \tab  1        \tab 2         \tab 0         \tab  \dots \tab 2\cr 
\eqn{X_3} \tab  1        \tab 0         \tab 1         \tab  \dots \tab 0\cr 
\dots     \tab  \dots    \tab \dots     \tab \dots     \tab  \dots \tab \dots\cr 
\eqn{X_J} \tab  1        \tab 2         \tab 0         \tab  \dots \tab 2\cr
}
Once there is a \code{coin} object, a similarity matrix can be obtained. Similarity matrices available in netCoin are:
\itemize{
\item Matching (m), Rogers & Tanimoto (t) Gower (g) Sneath (s) and Anderberg (and).
\item Jaccard (j), dice (d), antiDice (a), Ochiai (o) and Kulczynski (k).
\item Hamann (ham), Yule (y), Pearson (p), odds ratio (od) and Rusell (r).
}
Other measures that can be obtained from \code{coin} are:
\itemize{
\item Relative frequencies (x), conditional frequencies (i) coincidence degree (cc) and probable degree of coincidence (cp).
\item Haberman (h) and Z value of Haberman (z)
}
To obtain similarity and other measures matrices, the function \code{\link{sim}} elaborates a list of them.

\bold{Similarity matrix}

\tabular{lrrrrr}{
          \tab \eqn{X_1} \tab \eqn{X_2} \tab \eqn{X_3} \tab  \eqn{\;} \tab \eqn{X_J} \cr
\eqn{X_1} \tab  1.73     \tab -.87      \tab .87      \tab  \dots \tab -.87\cr 
\eqn{X_2} \tab  -.87     \tab 1.73      \tab -1.73     \tab  \dots \tab 1.73\cr 
\eqn{X_3} \tab  .87      \tab -1.73     \tab 1.73      \tab  \dots \tab -1.73\cr 
\dots     \tab  \dots    \tab \dots     \tab \dots     \tab  \dots \tab \dots\cr 
\eqn{X_J} \tab  -.87     \tab 1.73      \tab -1.73     \tab  \dots \tab 1.73\cr 
}
\code{\link{edgeList}} makes a collecion of edges composed by a list of similarity measures whenever a criterium (generally p(Z)<.50) is met.

\bold{Edge list}

\tabular{lrrrr}{
  \tab source \tab target \tab Haberman  \tab P(z)      \cr
1 \tab X1     \tab X3     \tab 0.8660254 \tab 0.22509243\cr
2 \tab X2     \tab X4     \tab 1.7320508 \tab 0.09084506\cr 
}
In order to make a graph, two data frames are needed: a nodes data frames with names and other nodes attributes (see \code{\link{asNodes}}) and an edge data frame (see \code{\link{edgeList}}). For more information go to \code{\link{netCoin}}.
}
% \section{Introduction}{Introduction example}
\section{Author}{Modesto Escobar, Department of Sociology and Communication, University of Salamanca.}
\references{
Escobar, M. (2009): "Redes Semanticas en Textos Periodisticos: Propuestas Tecnicas para su Representacion", en \emph{Empiria}, 17, 13-39.

Escobar, M.(2015): "Studying Coincidences with Network Analysis and Other Multivariate Tools", in \emph{The Stata Journal}, 15(4), 1118-1156.

Escobar, M. and J. Gomez Isla (2015): "The Expression of Identity through the Image: The Photographic Archives of Miguel de Unamuno and Joaquin Turina", en \emph{Revista Espanola de Investigaciones Sociologicas}, 152, 23-46.

Escobar, M. and L. Martínez-Uribe (2020): "Network Coincidence Analysis: The netCoin R Package", en \emph{Journal of Statistical Software}, 93(11), 1-31.

}

