\name{getEdgeWeightsG}
\alias{getEdgeWeightsG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getEdgeWeightsG
%%  ~~function to do ... ~~
}
\description{
This functions takes in a functional module from scPPIN and a countTable for gene reads as a data frame and computes the edge weight for each connection in the functional module and returns a network data frame with edge weights included. Used for creating Cytoscape graphs. 
}
\usage{
getEdgeWeightsG(FM, countTable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FM}{A scPPIN object that is returned from the function getFM.  Used here to generate the network}
  \item{countTable}{A data frame containing counts of gene expression for different cell barcodes.  Used here to generate edge weights}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns a functional Modeule (igraph) with weight attribute set for edges based on spearman's correlation
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Matt Heffernan, University of Illinois at Chicago
}
\note{
Called by getCyGraphs
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

getEdgeWeightsG(FM=FM0, countTable=countTable)

## The function is currently defined as
function(FM, countTable)
{
  network <- as.data.frame(get.edgelist(FM))
  
  nEdges <- nrow(network)
  i=1
  while(i <= nEdges)
  {
    G1 = network[i,1]
    G2 = network[i,2]
    G1C <- countTable[G1,]
    G2C <- countTable[G2,]
    corr <- unname(cor.test(G1C,G2C,method="spearman", exact=FALSE)$estimate)
    FM <- set_edge_attr(FM, "weight", index = E(FM, P = c(G1,G2),directed = FALSE), value = corr)
    i = i+1
  }
  return(FM)
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
