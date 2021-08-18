\name{getDataFiles}
\alias{getDataFiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getDataFiles
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calls many helper functions to compute the cluster lists, functional module, and edge weights and then sends the data to CSV files so that it can be output by Cytoscape.
}
\usage{
getDataFiles(SData, Ps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SData}{
%%     ~~Describe \code{SData} here~~
A Seurat object containing the scppi data 
}
  \item{Ps}{
%%     ~~Describe \code{Ps} here~~
Pvalues containing the Pvalues for each gene in the clusters
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
Does not return any variable.  Creates 2 CSV files for every cluster.  One for the nodes and another for the edges in Cytoscape.  
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

getDataFiles(SData = seqwell, Ps = PValues)

## The function is currently defined as
function (SData, Ps) 
{
    countTable <- SData[["RNA"]]@data
    CList <- getClusterList(SData, Ps)
    FMList <- list()
    i = 0
    for (c in CList) {
        pValues <- c$p_val
        names(pValues) <- c$gene
        FMList[[paste("FM", as.character(i), sep = "")]] <- FM(pValues)
        i = i + 1
    }
    networkList <- list()
    i = 0
    for (f in FMList) {
        networkList[[paste("N", as.character(i), sep = "")]] <- getEdgeWeights(f, 
            countTable)
        i = i + 1
    }
    i = 1
    nClusters <- length(CList)
    while (i <= nClusters) {
        SendToCSV(FMList[[i]], networkList[[i]], CList[[i]], 
            (i - 1))
        i = i + 1
    }
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
