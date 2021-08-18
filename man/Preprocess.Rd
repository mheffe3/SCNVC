\name{Preprocess}
\alias{Preprocess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Preprocess
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Takes in a dataframe containing a scppin and calls Seurat functions on it to get it in the proper format to use scPPIN functions with it.  
}
\usage{
Preprocess(Data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
%%     ~~Describe \code{Data} here~~
A dataframe containing scppin data
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
returns a Seurat object with Seurat functions called on it
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

Preprocess(uterusData)

## The function is currently defined as
function (Data) 
{
    SData <- CreateSeuratObject(counts = Data, min.cells = 3, 
        min.features = 300, names.delim = "\\.")
    SData[["percent.mito"]] <- PercentageFeatureSet(object = SData, 
        pattern = "^MT-")
    SData <- subset(SData, subset = nFeature_RNA > 300 & nFeature_RNA < 
        2500 & percent.mito < 20)
    SData <- NormalizeData(object = SData, normalization.method = "LogNormalize", 
        scale.factor = 10000)
    SData <- ScaleData(SData, vars.to.regress = c("nCount_RNA", 
        "percent.mito"))
    SData <- FindVariableFeatures(SData, selection.method = "disp")
    length(SData@assays$RNA@var.features)
    var.gene <- SData@assays$RNA@var.features
    var.gene <- var.gene[!grepl(pattern = "*RPS", x = var.gene)]
    var.gene <- var.gene[!grepl(pattern = "*RPL", x = var.gene)]
    var.gene <- var.gene[!grepl(pattern = "*MT", x = var.gene)]
    length(var.gene)
    SData <- RunPCA(object = SData, features = var.gene, npcs = 50, 
        do.print = TRUE, pcs.print = 1:5, nfeatures.print = 5)
    SData <- FindNeighbors(object = SData, reduction = "pca", 
        dims = 1:15)
    SData <- FindClusters(object = SData, reduction = "pca", 
        dims = 1:15, resolution = c(0.6))
    return(SData)
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
