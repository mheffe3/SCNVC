\name{getFile}
\alias{getFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getFile
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Takes in a filename/path and checks if the file exists or is empty.  Then it checks that it is either a .csv or a .txt file.  If it is empty or does not exist it stops and reports this. If it is not a .csv or .txt it stops and reports this.  If it is a .csv or .txt it will create a dataframe from the data and return it filtered for cell barcodes with colsums >= 500
}
\usage{
getFile(filepath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filepath}{
%%     ~~Describe \code{filepath} here~~
The file name or path as a string
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
returns a dataframe of the scppin form the file if it is in the proper format
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

getFile("GSM4008666_Adult-Uterus1_dge.txt")

## The function is currently defined as
function (filepath) 
{
    if ((file.exists(filepath)) & (file.size(filepath) > 0)) {
        if (grepl(".csv", filepath, fixed = T)) {
            Data <- read.csv(filepath, header = T, row.names = 1)
        }
        else if (grepl(".txt", filepath, fixed = T)) {
            Data <- read.delim(filepath, header = T, row.names = 1)
        }
        else {
            stop("FILE NOT IN PROPER FORMAT (.csv OR .txt)")
        }
    }
    else {
        stop("FILE DOES NOT EXIST OR IS EMPTY")
    }
    Data <- Data[, colSums(Data) >= 500]
    return(Data)
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
