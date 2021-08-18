
# library(Seurat)
# library(scPPIN)
# library(data.table)
# library(kutils)
# library(RCy3)

#' @import Seurat
#' @import scPPIN
#' @import RCy3
#' @importFrom kutils writeCSV
#' 
#' @export
getFile <- function(filepath)
{
  if((file.exists(filepath)) & (file.size(filepath) > 0))
  {
    if(grepl(".csv", filepath, fixed = T))
    {
      Data <- read.csv(filepath, header = T, row.names = 1)
    }
    else if(grepl(".txt", filepath, fixed = T))
    {
      Data <- read.delim(filepath, header = T, row.names = 1)
    }
    else
    {
      stop("FILE NOT IN PROPER FORMAT (.csv OR .txt)")
    }
    
  }
  else
  {
    stop("FILE DOES NOT EXIST OR IS EMPTY")
  }
  
  Data <- Data[,colSums(Data)>=500]
  return(Data)
  
}

#' @export
Preprocess <- function(Data)
{
  SData <- CreateSeuratObject(counts = Data, min.cells = 3,min.features = 300,names.delim = "\\.")
  SData[["percent.mito"]] <- PercentageFeatureSet(object = SData, pattern = "^MT-")
  SData <- subset(SData , subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mito < 20)
  SData <- NormalizeData(object = SData, normalization.method = "LogNormalize",scale.factor = 10000)
  SData <- ScaleData(SData, vars.to.regress=c("nCount_RNA", "percent.mito"))
  
  SData <- FindVariableFeatures(SData,selection.method = "disp")
  
  length(SData@assays$RNA@var.features)
  var.gene<-SData@assays$RNA@var.features
  var.gene<-var.gene[!grepl(pattern = "*RPS",x=var.gene)]
  var.gene<-var.gene[!grepl(pattern = "*RPL",x=var.gene)]
  var.gene<-var.gene[!grepl(pattern = "*MT",x=var.gene)]
  length(var.gene)
  SData <- RunPCA(object = SData,  features = var.gene, npcs = 50, do.print = TRUE,pcs.print = 1:5, nfeatures.print = 5)
  
  SData <- FindNeighbors(object = SData, reduction = "pca", dims = 1:15)
  
  SData <- FindClusters(object = SData, reduction = "pca", dims = 1:15,resolution =c(0.6))
  
  return(SData)
}


#' @export
getPValues <- function(SData)
{
  PValues <- FindAllMarkers(SData, return.thresh=1, logfc.threshold = 0.0)
  return(PValues)
}



#scPPIN


#' @export
getClusterList <- function(Data, PValues)
{
  nClusters = nlevels(seqwell@meta.data$seurat_clusters)
  ClusterList <- list()
  i = 0
  while(i < nClusters)
  {
    ClusterList[[paste("C",as.character(i),sep="")]] = PValues[PValues$cluster==i,]
    i = i+1
  }
  return(ClusterList)
  
}


#' @export
FM <- function(pvalues)
{
  ppin <- loadPPIN()
  checksum = sum(pvalues<0) + sum(pvalues>1)
  if (checksum>0){
    warning('p-values are not in (0,1]!')
  }
  FDR <- 10^{-2}
  pvalues[pvalues<10^-308] <- 10^-308
  functionalModule <-detectFunctionalModule(ppin,pvalues,FDR)
  #./dapcstp DAPCSTP_temp.stp --type pcstp -o DAPCSTP_temp.sol
  return(functionalModule)
}


#' @export
getEdgeWeights <- function(functionalModule, countTable)
{
  network <- as.data.frame(get.edgelist(functionalModule))
  colnames(network) <- c('TF','GENE')
  Rho <- vector()
  
  nEdges <- nrow(network)
  i=1
  while(i <= nEdges)
  {
    G1 = network[i,1]
    G2 = network[i,2]
    G1C <- countTable[G1,]
    G2C <- countTable[G2,]
    corrData <- cor.test(G1C,G2C,method="spearman", exact=FALSE)
    corr <- unname(corrData$estimate)
    Rho[i] <- corr
    i = i+1
  }
  network$Rho <- Rho
  return(network)
}


#' @export
SendToCSV <- function(functionalModule, network, Ps, name)
{
  filename = paste(name, "N.csv", sep="")
  writeCSV(network, filename, row.names = FALSE)
  
  filename = paste(name, "E.csv", sep="")
  table <- as_data_frame(functionalModule,what=c("vertices"))
  writeCSV(table2, filename, row.names = FALSE)
  
}


#' @export
getDataFiles <- function(SData, Ps)
{
  countTable <- SData[["RNA"]]@data
  CList <- getClusterList(SData, Ps)
  FMList <- list()
  i = 0
  for(c in CList)
  {
    pValues <- c$p_val
    names(pValues) <- c$gene
    FMList[[paste("FM", as.character(i),sep="")]] <- FM(pValues)
    i = i+1
    
  }
  networkList <- list()
  i = 0
  for(f in FMList)
  {
    networkList[[paste("N", as.character(i),sep="")]] <- getEdgeWeights(f,countTable)
    i = i+1
    
  }
  
  i = 1
  nClusters <- length(CList)
  while(i <= nClusters)
  {
    SendToCSV(FMList[[i]],networkList[[i]],CList[[i]],(i-1))
    i = i+1
  }
  
}


#Function Calls


# uterusData <- getFile("GSM4008666_Adult-Uterus1_dge.txt")
# seqwell <- Preprocess(uterusData)
# PValues <- getPValues(seqwell)
# getDataFiles(SData = seqwell, Ps = PValues)
