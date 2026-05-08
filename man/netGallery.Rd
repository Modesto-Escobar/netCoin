\name{netGallery}
\alias{netGallery}
\title{Images in a grid gallery.}
\description{
\code{netGallery} produces an interactive image gallery. 
}
\value{
Object of class \code{netGallery}.
}

\usage{
netGallery(tree, deep = FALSE, initialType = NULL, tableformat = FALSE, ...)
}
\arguments{
\item{tree}{a data frame with two columns: source and target, describing relationships between nodes. It indicates a hierarchy between nodes which can be dynamically explored. Optionally, another two columns describing types can be passed.}
\item{deep}{The tree is especified in a structure that preserves paths. Each column will be a tree level, colnames as node types.}
\item{initialType}{A character vector indicating which node type will be shown at start. No effects in 'deep' mode.}
\item{tableformat}{If the tree is especified as a table, one column per level. No effects in 'deep' mode.}
  \item{...}{Any \link{gallery} argument.}
}
\author{
Modesto Escobar, Department of Sociology and Communication, University of Salamanca.
}
\examples{
# simply linked tree
tree <- data.frame(
  parent=c("a","a","b","b","b","ab","ab","ab","ab","ba","ba","baa","ba"),
  child=c("aa","ab","ba","bb","bc","aba","abb","abc","abd","baa","bab","baaa","aba")
)

gallery <- netGallery(tree)

\dontrun{
plot(gallery)
}

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

gallery <- netGallery(tree, initialType="film", tableformat=TRUE, color="type", scale=2)

\dontrun{
plot(gallery)
}

# example with path preservation
tree <- rbind(
c("a1","b1","c1|c2"),
c("a2","b2","c3|c4"),
c("a3","b2","c4|c5"),
c("a4","b3|b4","c6|c7|c8")
)

colnames(tree) <- c("a","b","c")

gallery<-netGallery(tree,deep=TRUE,color="type")
\dontrun{
plot(gallery)
}
}
