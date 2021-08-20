\name{getCyGraphs}
\alias{getCyGraphs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getCyGraphs
%%  ~~function to do ... ~~
}
\description{
Calls all the functions necessary to get the functional modules from the Seurat data and PValues and then graphs in Cytoscape.  
}
\usage{
getCyGraphs(SData, Ps,NSLB = 100, EWLB = 0.05, cnames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SData}{A Seurat object containing preprocessed data }
  \item{Ps}{Pvalues containing the Pvalues for each gene in the clusters}
  \item{NSLB}{Lower bound for Node Score filtering}
  \item{EWLB}{Lower Bound For Edge Weight filtering}
  \item{cnames}{A list containing names for each of the clusters.  Defaults to NULL and uses C1, C2... instaed}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
No return value instead creates a filtered and unfiltered network in Cytoscape for each cluster
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Matt Heffernan, University of Illinois at Chicago
}
\note{
Requires Cytoscape to be opened and installed.  Calls getClusterList, getFM, getEdgeWeightsG, GraphInCytoscape, and CyFilter
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

getCyGraphs(SData=SData, Ps=Ps,NSLB = 100, EWLB = 0.05, cnames = NULL)

## The function is currently defined as
function(SData, Ps,NSLB = 100, EWLB = 0.05, cnames = NULL)
{
  countTable <- SData[["RNA"]]@data
  CList <- getClusterList(SData, Ps, cnames)
  FMList <- list()
  i = 0
  for(c in CList)
  {
    pValues <- c$p_val
    names(pValues) <- c$gene
    FMList[[paste("FM", as.character(i),sep="")]] <- getFM(pValues)
    i = i+1
    
  }
  networkList <- list()
  i = 0
  for(f in FMList)
  {
    networkList[[paste("N", as.character(i),sep="")]] <- getEdgeWeightsG(f,countTable)
    i = i+1
    
  }
  
  i = 1
  nClusters <- length(CList)
  hasnames = TRUE
  if(is.null(cnames) || length(cnames) != nClusters)
  {
    hasnames = FALSE
  }
  while(i <= nClusters)
  {
    if(hasnames)
    {
      GraphInCytoscape(FMList[[i]], name=cnames[i])
      CyFilter(FMList[i], name=cnames[i], NSLB = NSLB, EWLB = EWLB)
      
    }
    else
    {
      GraphInCytoscape(FMList[[i]], name=paste("C", as.character(i-1), sep=""))
      CyFilter(FM, name=paste("C", as.character(i-1), sep=""), NSLB = NSLB, EWLB = EWLB)
    }
    
    i = i+1
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
