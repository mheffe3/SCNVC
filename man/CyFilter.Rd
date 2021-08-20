\name{CyFilter}
\alias{CyFilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CyFilter
%%  ~~function to do ... ~~
}
\description{
Filters the Cytoscape graph to include only NodeScores > NSLB OR EdgeWeights > |0.05| Then creates a subnetwork from it
}
\usage{
CyFilter(FM, name, NSLB = 100, EWLB = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FM}{functional module from scPPIN}
  \item{name}{name of FM cluster}
  \item{NSLB}{Lower bound for Node Score filtering}
  \item{EWLB}{Lower Bound For Edge Weight filtering}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
No return value simply outputs a subnetwork to Cytoscape
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Matt Heffernan, University of Illinois at Chicago
}
\note{
Requires Cytoscape to be opened and installed.  Called by getCyGraphs
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

CyFilter <- function(FM=FM0, name="endothelial", NSLB = 100, EWLB = 0.05)

## The function is currently defined as
function(FM, name, NSLB = 100, EWLB = 0.05)
{
  createColumnFilter(filter.name="Node Score Filter", column = "nodeScore", c(NSLB, max(V(FM0)$nodeScore)), predicate = "BETWEEN")
  createColumnFilter(filter.name = "Edge Weight Filter", column = "weight", c(-EWLB, EWLB), predicate = "IS_NOT_BETWEEN", type = "edges")
  createCompositeFilter(filter.name = "Edge Weight & Node Score Filter", c("Node Score Filter", "Edge Weight Filter"), "ANY")
  createSubnetwork(subnetwork.name = paste(name, "Filtered Network", sep=""))
  
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
