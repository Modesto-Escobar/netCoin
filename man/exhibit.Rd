\name{exhibit}
\alias{exhibit}
\title{Images in a grid gallery.}
\description{
\code{exhibit} produces an interactive image gallery (alternative display).
}
\value{
Object of class \code{gallery_rd3}.
}

\usage{
exhibit(nodes, name = NULL, label = NULL, subtitle = NULL, order = NULL,
    decreasing = FALSE, ntext = NULL,
    mainframeHeight = NULL, mainframeWidth = NULL, mainframeImage = 0,
    image = NULL, imageCopy = NULL, noFilterCols = NULL,
    zoom = NULL, scale = NULL, main = NULL, note = NULL,
    logo = TRUE, pages = c("bottom","dynamic","top"), mainposition = c("dynamic","fixed"),
    frequencies = FALSE, tutorial = FALSE, tableButton = FALSE, export = FALSE,
    search = TRUE, colorScheme = 0,  language = c("en", "es", "ca"), dir = NULL)
}
\arguments{
\item{nodes}{a data frame with at least one column with the names of the gallery's elements.}
\item{name}{name of the vector with names in the nodes data frame.}
\item{label}{column name with image labels in the nodes data frame.}
\item{subtitle}{column name with image subtitles in the nodes data frame.}
\item{order}{name of the column with node order in the nodes data frame.}
\item{decreasing}{decreasing or increasing sort of the nodes.}
\item{ntext}{column name with html text in the nodes data frame.}
\item{mainframeHeight}{a number between 0.1 and 1 to reduce mainframe height.}
\item{mainframeWidth}{a number to proportionally change mainframe width.}
\item{mainframeImage}{a number to display image in mainframe. 0 = auto (default), 1 = before text, 2 = top right}
\item{image}{column name which indicates the image paths in the nodes data frame.}
\item{imageCopy}{column name which indicates the image copyright and attribution in the nodes data frame.}
\item{noFilterCols}{disable filtering by especified column names.}
\item{zoom}{Deprecated. Replaced by argument 'scale'.}
\item{scale}{a number between 0.1 and 1 as initial displaying scale.}
\item{main}{upper title of the graph.}
\item{note}{lower title of the graph.}
\item{logo}{switch to FALSE if you want to hide netCoin logo.}
\item{pages}{position of pagination info. In the header (top), footer (bottom), or depending on the scroll (dynamic).}
\item{mainposition}{If 'fixed', the upper title will always be visible in the header, by default ('dynamic'), it will be hidden when scrolling.}
\item{frequencies}{a logical value true if barplots representing node attributes frequencies will be added to the final graph.}
\item{tutorial}{Should tutorial be displayed?}
\item{tableButton}{A button will appear to show the data as a table.}
\item{export}{A button will appear to export gallery as excel.}
\item{search}{Should gallery search be displayed? (default: TRUE)}
\item{colorScheme}{Select a color theme for display (0-11).}
\item{language}{a character string indicating the language of the graph (en=english (default); es=spanish; ca=catalan).}
\item{dir}{a character string representing the directory where the web files will be saved.}
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
\examples{
data("finches")
finches$species <- system.file("extdata", finches$species,
        package="netCoin") # copy path to the species field
exhibit(finches, image="species")
}
