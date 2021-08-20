\name{SendToCSV}
\alias{SendToCSV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{SendToCSV
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Takes in a functional Module from scPPIN, a network dataframe, gene pvalues for expression, and the cluster number, and writes two csvs that can be used by Cytoscape for visualization of the network.
}
\usage{
SendToCSV(functionalModule, network, Ps, C)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{functionalModule}{
%%     ~~Describe \code{functionalModule} here~~
A object containing the solution to the cluster network computed by scPPIN
}
  \item{network}{
%%     ~~Describe \code{network} here~~
A dataframe containing the network and edge weights for the cluster
}
  \item{Ps}{
%%     ~~Describe \code{Ps} here~~
A dataframe containg the PValues for each genes expression
}
  \item{C}{
%%     ~~Describe \code{C} here~~
The cluster number
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
Does not return any varibale, but instead creates two csv files in the current working directory
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Matt Heffernan, University of Illinois at Chicago
}
\note{
called by getDataFiles
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

SendToCSV(FM, network, Ps, 0)

## The function is currently defined as
function (functionalModule, network, Ps, C) 
{
    filename = paste("UterusNetwork", C, ".csv", sep = "")
    writeCSV(network, filename, row.names = FALSE)
    filename = paste("UterusTableNS", C, ".csv", sep = "")
    table2 <- as_data_frame(functionalModule, what = c("vertices"))
    writeCSV(table2, filename, row.names = FALSE)
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
