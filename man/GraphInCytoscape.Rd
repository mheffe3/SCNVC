\name{GraphInCytoscape}
\alias{GraphInCytoscape}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{GraphInCytoscape
%%  ~~function to do ... ~~
}
\description{
Graphs a functional Module (igraph) from scPPIN in Cytoscape.  The more red the node the more highly expressed it is.  If node is an ellipse it is a TF else it's not.  Edge thickness corresponds to the correlation between the two nodes the thicker the more highly correlated.  Edge color corresponds to correlation, the more red the more positively correlated the more blue the more negatively correlated.
}
\usage{
GraphInCytoscape(FM, name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FM}{functional Module from scPPIN to be graphed
}
\item{name}{name of cluster being graphed
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
No return value simply graphs the cluster as a network in Cytoscape
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Matt Heffernan, University of Illinois at Chicago
}
\note{
to use Cytoscape must be installed and opened.  Called by getCyGraph
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

GraphInCytoscape <- function(FM=FM0, name="endothelial")

## The function is currently defined as
function(FM, name)
{
  createNetworkFromIgraph(FM, title= name)
  
  setNodeColorMapping('nodeScore', c(min(V(FM)$nodeScore), mean(V(FM)$nodeScore), max(V(FM)$nodeScore)), c('#F5EDDD', '#F59777', '#F55333'))
  TFs <- readLines("TFs.txt")
  table <- as_data_frame(FM,what=c("vertices"))
  filtered <- table[table$name %in% TFs,]
  TF <- filtered[,'name']
  setNodeShapeBypass(node.names=TF, new.shapes="ELLIPSE")
  
  
  setEdgeLineWidthMapping('weight',c(min(E(FM)$weight), mean(E(FM)$weight), max(E(FM)$weight)), widths = c(1,10,20))
  setEdgeColorMapping('weight',c(min(E(FM)$weight), mean(E(FM)$weight), max(E(FM)$weight)), colors= c("#00008B", "#FFFFFF", "#EE4B2B"))
  
  
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
