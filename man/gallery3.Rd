\name{gallery3}
\alias{gallery3}
\title{Images in a grid gallery.}
\description{
\code{gallery3} produces an interactive image gallery (alternative display).
}
\value{
Object of class \code{gallery_rd3}.
}

\usage{
gallery3(nodes, name = NULL, label = NULL, subtitle = NULL, order = NULL,
    decreasing = FALSE, ntext = NULL,
    mainframeHeight = NULL, mainframeWidth = NULL, mainframeImage = 0,
    image = NULL, imageCopy = NULL, imageRatio = NULL, noFilterCols = NULL,
    zoom = NULL, scale = NULL, main = NULL, note = NULL,
    tableButton = FALSE, export = FALSE,
    search = TRUE, cex = 1, language = c("en", "es", "ca"), dir = NULL)
}
\arguments{
\item{nodes}{a data frame with at least one column with the names of the gallery's elements.}
\item{name}{name of the vector with names in the nodes data frame.}
\item{label}{column name with image labels in the nodes data frame.}
\item{subtitle}{column name with image subtitles in the nodes data frame.}
\item{order}{name of the column with node order in the nodes data frame.}
\item{decreasing}{decreasing or increasing sort of the nodes.}
\item{ntext}{column name with html text in the nodes data frame.}
\item{mainframeHeight}{a number to proportionally change mainframe height.}
\item{mainframeWidth}{a number to proportionally change mainframe width.}
\item{mainframeImage}{a number to display image in mainframe. 0 = auto (default), 1 = before text, 2 = top right}
\item{image}{column name which indicates the image paths in the nodes data frame.}
\item{imageCopy}{column name which indicates the image copyright and attribution in the nodes data frame.}
\item{imageRatio}{the aspect ratio of the images. Is the ratio of its width to its height. 3/4 by default.}
\item{noFilterCols}{disable filtering by especified column names.}
\item{zoom}{Deprecated. Replaced by argument 'scale'.}
\item{scale}{a number between 0.1 and 1 as initial displaying scale.}
\item{main}{upper title of the graph.}
\item{note}{lower title of the graph.}
\item{search}{Should gallery search be displayed? (default: TRUE)}
\item{tableButton}{A button will appear to show the data as a table.}
\item{export}{A button will appear to export gallery as excel.}
\item{cex}{number indicating the amount by which plotting text should be scaled relative to the default.}
\item{language}{a character string indicating the language of the graph (en=english (default); es=spanish; ca=catalan).}
\item{dir}{a character string representing the directory where the web files will be saved.}
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
\examples{
data("finches")
finches$species <- system.file("extdata", finches$species,
        package="rD3plot") # copy path to the species field
gallery <- gallery3(finches, image="species")
\dontrun{
plot(gallery)
}
}
