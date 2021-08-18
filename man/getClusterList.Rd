\name{getClusterList}
\alias{getClusterList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getCluterList
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Gets a list of all the clusters in the seurat object passed and returns it
}
\usage{
getClusterList(Data, PValues)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
%%     ~~Describe \code{Data} here~~
Seurat object containing the clusters
}
  \item{PValues}{
%%     ~~Describe \code{PValues} here~~
PValues returned from getPValues
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list containing the clusters in the Seurat object
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Matt Heffernan, University of Illinois at Chicago
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

getClusterList(SData, Ps)

## The function is currently defined as
function (Data, PValues) 
{
    nClusters = nlevels(seqwell@meta.data$seurat_clusters)
    ClusterList <- list()
    i = 0
    while (i < nClusters) {
        ClusterList[[paste("C", as.character(i), sep = "")]] = PValues[PValues$cluster == 
            i, ]
        i = i + 1
    }
    return(ClusterList)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
