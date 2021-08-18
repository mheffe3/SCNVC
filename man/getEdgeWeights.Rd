\name{getEdgeWeights}
\alias{getEdgeWeights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getEdgeWeights
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This functions takes in a functional module from scPPIN and a countTable for gene reads as a data frame and computes the edge weight for each connection in the functional module and returns a network data frame with edge weights included.  
}
\usage{
getEdgeWeights(functionalModule, countTable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{functionalModule}{
%%     ~~Describe \code{functionalModule} here~~
A scPPIN object that is returned from the function FM.  Used here to generate the network
}
  \item{countTable}{
%%     ~~Describe \code{countTable} here~~
A data frame containing counts of gene expression for different cell barcodes.  Used here to generate edge weights
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
returns a data frame containing the network and the edge weights associated with each connection
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

getEdgeWeights(f,countTable)

## The function is currently defined as
function (functionalModule, countTable) 
{
    network <- as.data.frame(get.edgelist(functionalModule))
    colnames(network) <- c("TF", "GENE")
    Rho <- vector()
    nEdges <- nrow(network)
    i = 1
    while (i <= nEdges) {
        G1 = network[i, 1]
        G2 = network[i, 2]
        G1C <- countTable[G1, ]
        G2C <- countTable[G2, ]
        corrData <- cor.test(G1C, G2C, method = "spearman", exact = FALSE)
        corr <- unname(corrData$estimate)
        Rho[i] <- corr
        i = i + 1
    }
    network$Rho <- Rho
    return(network)
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
