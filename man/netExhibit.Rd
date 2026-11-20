\name{netExhibit}
\alias{netExhibit}
\title{Images in a grid gallery.}
\description{
\code{netExhibit} produces an interactive image gallery (alternative display). 
}
\value{
Object of class \code{treeGallery_rd3}.
}

\usage{
netExhibit(tree, initialType = NULL, tableformat = FALSE, ...)
}
\arguments{
\item{tree}{a data frame with two columns: source and target, describing relationships between nodes. It indicates a hierarchy between nodes which can be dynamically explored. Optionally, another two columns describing types can be passed.}
\item{initialType}{A character vector indicating which node type will be shown at start.}
\item{tableformat}{If the tree is especified as a table, one column per level.}
\item{...}{Any \link{exhibit} argument.}
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
\examples{
# Create some film data
tree <- data.frame(
  film = c("The Shawshank Redemption", "The Godfather", "The Dark Knight", "The Godfather Part II",
    "12 Angry Men", "Schindler's List", "The Lord of the Rings: The Return of the King",
    "Pulp Fiction", "The Lord of the Rings: The Fellowship of the Ring",
    "The Good, the Bad and the Ugly"),
  director = c("Frank Darabont", "Francis Ford Coppola", "Christopher Nolan",
    "Francis Ford Coppola", "Sidney Lumet", "Steven Spielberg", "Peter Jackson",
    "Quentin Tarantino", "Peter Jackson", "Sergio Leone"),
  actor = c("Tim Robbins|Morgan Freeman|Bob Gunton", "Marlon Brando|Al Pacino|James Caan",
    "Christian Bale|Heath Ledger|Aaron Eckhart", "Al Pacino|Robert De Niro|Robert Duvall",
    "Henry Fonda|Lee J. Cobb|Martin Balsam", "Liam Neeson|Ralph Fiennes|Ben Kingsley",
    "Elijah Wood|Viggo Mortensen|Ian McKellen", "John Travolta|Uma Thurman|Samuel L. Jackson",
    "Elijah Wood|Viggo Mortensen|Ian McKellen", "Clint Eastwood|Eli Wallach|Lee Van Cleef")
)

gallery <- netExhibit(tree, initialType="film", tableformat=TRUE)

\dontrun{
plot(gallery)
}
}
