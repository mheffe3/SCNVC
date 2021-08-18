\name{getPValues}
\alias{getPValues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getPValues
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calls the FindAllMarkers function from the Seurat package with a Seurat object and returns the PValues associated with that object
}
\usage{
getPValues(SData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SData}{
%%     ~~Describe \code{SData} here~~
Seurat object containing single cell RNA seq data
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
Returns a data frame containing the pvalues for each gene in the data corresponding to the genes expression level 
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

PValues <- getPValues(SData)

## The function is currently defined as
function (SData) 
{
    PValues <- FindAllMarkers(SData, return.thresh = 1, logfc.threshold = 0)
    return(PValues)
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
