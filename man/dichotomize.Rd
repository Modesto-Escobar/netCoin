\name{dichotomize}
\alias{dichotomize}
\title{
  Dichotomize
}
\description{
  Converts factor or character columns of a data frame into a set of dichotomous (dummy/indicator) columns. Their names will correspond to the unique labels or texts of every category discovered across the selected variables. If multiple columns are provided, categories are unified per row.
}
\usage{
  dichotomize (data, variables, sep = "", regex = FALSE, 
               text = c("none", "lower", "upper"),
               min = 1, length = 0, sort = TRUE,
               values = NULL, stopwords = NULL, nas = "None",
               sparse = FALSE, add = FALSE)
}
\arguments{
  \item{data}{A data frame containing the factor or textual columns to be converted. Cells can contain simple values or multiple events delimited by a separator.}
  \item{variables}{A character vector of column names to be converted into dichotomous vectors.}
  \item{sep}{A character string (or vector) used to divide columns with multiple events. If this separator is \code{""} (default), every unique cell is treated as a single category. Use \code{"C"} to concatenate multiple columns before dichotomization.}
  \item{regex}{Logical. If \code{FALSE} (default), \code{sep} is matched exactly. If \code{TRUE}, \code{sep} is treated as a regular expression (e.g., \code{"\\\s+"} to split by whitespace).}
  \item{text}{Character. Text transformation option. \code{"none"} (default) leaves text intact. \code{"lower"} or \code{"upper"} converts text to lower or upper case, respectively, and completely removes punctuation and numbers.}
  \item{min}{Integer or numeric. Converts to dichotomous vectors only the labels or texts that have a frequency greater than or equal to the value of this parameter. If the value is between 0 and 1, it is interpreted as a percentage of the total rows.}
  \item{length}{Integer. The strict maximum number of dichotomous columns generated, prioritizing the most frequent categories. \code{0} implies no limit.}
  \item{sort}{Logical. If \code{TRUE}, the new columns are ordered by their frequencies in descending order.}
  \item{values}{A character vector of specific labels or texts strictly selected for conversion into dichotomous columns.}
  \item{stopwords}{Character vector. An optional list of words to be removed from the text before dichotomization. Only applied if \code{text} is set to \code{"lower"} or \code{"upper"}.}
    \item{nas}{Character. Variable name used to account for \code{NA} (missing) values across the set of variables. Set to \code{NULL} to ignore missing values.}
  \item{sparse}{Logical. If \code{TRUE}, returns a sparse matrix (requires the \code{Matrix} package) instead of a data.frame.}
  \item{add}{Logical. If \code{TRUE}, appends the newly generated dichotomous columns to the original input data.frame.}
}
\value{
  A data frame composed of the original variables plus the added dichotomous columns, or a sparse matrix if \code{sparse = TRUE}.
}
\references{
  Escobar, M. and Martinez-Uribe, L. (2020)
  Network Coincidence Analysis: The \code{netCoin} \code{R} Package.
  \emph{Journal of Statistical Software}, \bold{93}, 1-32.
  \doi{10.18637/jss.v093.i11}.
}
\author{
  Modesto Escobar, Department of Sociology and Communication, University of Salamanca, and Luis Martinez Uribe, Fundacion Juan March.
}
\examples{
  # A character column
  frame1 <- data.frame(A = c("Man", "Women", "Man", "Undet."))
  dichotomize(frame1, "A", sep = "; ")
  
  # A character column (with separator)
  frame2 <- data.frame(A = c("Man; Women", "Women; Women",
                             "Man; Man", "Undet.; Women; Man"))
  dichotomize(frame2, "A", sep = "; ")
  
  # A character column and another factor column (same separator)
  frame3 <- data.frame(A = c("Man; Women", "Women; Women",
                             
                             "Man; Man", "Undet.; Women; Man"),
                       C = factor(c(1:4), labels = c("Paris", "New York",
                                                     "London; New York", "<NA>")))
  dichotomize(frame3, c("A", "C"), sep = "; ")
  
  # A set of simple character or factor (same levels) variables.
  # In this case, you must use "C" separator.
  frame4 <- data.frame(A = c("Man", "Women","Man", "Undet",NA),
                       B = c("Women","Women","Man","Women",NA),
                       C = c(NA,NA,NA,"Man",NA))
  dichotomize(frame4, c("A","B","C"), sep="C")
  
  # Text mining example: Splitting by whitespace, removing punctuation/numbers, 
  # transforming to lowercase and filtering stopwords
  frame5 <- data.frame(text_col = c("I love apples and 10 oranges!", 
                                    "Apples are great.", 
                                    "I hate bananas."))
  dichotomize(frame5, "text_col", sep = "\\\\s+", regex = TRUE, 
              text = "lower", stopwords = c("i", "and", "are"))
}
