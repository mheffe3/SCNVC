\name{getFM}
\alias{getFM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getFM
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calls scPPIN functions to compute a steiner tree that finds the best solution to the connections in a cluster from the single cell data and returns it in functional module object.  Called by getDataFiles
}
\usage{
getFM(pvalues)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvalues}{
%%     ~~Describe \code{pvalues} here~~
A vector containing the named pvalues for each node in a cluster
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
Returns a functional module containing the solution to the connections of the ppi network solved by dapcstp
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Matt Heffernan, University of Illinois at Chicago
}
\note{
Called by getDataFiles and getCyGraphs
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

getFM(pValues)

## The function is currently defined as
function (pvalues) 
{
    ppin <- loadPPIN()
    checksum = sum(pvalues < 0) + sum(pvalues > 1)
    if (checksum > 0) {
        warning("p-values are not in (0,1]!")
    }
    FDR <- 10^{
        -2
    }
    pvalues[pvalues < 10^-308] <- 10^-308
    functionalModule <- detectFunctionalModule(ppin, pvalues, 
        FDR)
    return(functionalModule)
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
