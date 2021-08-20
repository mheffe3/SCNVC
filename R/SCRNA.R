
# library(Seurat)
# library(scPPIN)
# library(kutils)
# library(RCy3)







#' @import Seurat
#' @import scPPIN
#' @import RCy3
#' @importFrom kutils writeCSV



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
getClusterList <- function(SData, PValues, names = NULL)
{
  nClusters = nlevels(SData@meta.data$seurat_clusters)
  ClusterList <- list()
  i = 0
  if(!is.null(names) && length(names) == nClusters)
  {
    while(i < nClusters)
    {
      ClusterList[[names[i]]] = PValues[PValues$cluster==i]
      i = i+1
      
    }
   
  }
  
  else
  {
    warning("EITHER NO NAMES WERE GIVEN OR CHARCTER VECTOR PROVIDED IS NOT THE SAME LENGTH AS THE NUMBER OF CLUSTERS. USING C-ith FOR NAMES")
    while(i < nClusters)
    {
      ClusterList[[paste("C",as.character(i),sep="")]] = PValues[PValues$cluster==i,]
      i = i+1
    }
  }
  
  return(ClusterList)
  
}





#' @export
getFM <- function(pvalues)
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
getEdgeWeightsF <- function(FM, countTable)
{
  network <- as.data.frame(get.edgelist(FM))
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
    corr <- unname(cor.test(G1C,G2C,method="spearman", exact=FALSE)$estimate)
    Rho[i] <- corr
    i = i+1
  }
  network$Rho <- Rho
  return(network)
}




#' @export
getEdgeWeightsG <- function(FM, countTable)
{
  network <- as.data.frame(get.edgelist(FM))
  
  nEdges <- nrow(network)
  i=1
  while(i <= nEdges)
  {
    G1 = network[i,1]
    G2 = network[i,2]
    G1C <- countTable[G1,]
    G2C <- countTable[G2,]
    corr <- unname(cor.test(G1C,G2C,method="spearman", exact=FALSE)$estimate)
    FM <- set_edge_attr(FM, "weight", index = E(FM, P = c(G1,G2),directed = FALSE), value = corr)
    i = i+1
  }
  return(FM)
}



#' @export
SendToCSV <- function(functionalModule, network, Ps, name = C)
{
    filename = paste(name, "N.csv", sep="")
    writeCSV(network, filename, row.names = FALSE)
    
    filename = paste(name, "E.csv", sep="")
    table <- as_data_frame(functionalModule,what=c("vertices"))
    writeCSV(table2, filename, row.names = FALSE)
  
}



#' @export
GraphInCytoscape <- function(FM, name)
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




#' @export
CyFilter <- function(FM, name, NSLB = 100, EWLB = 0.05)
{
  createColumnFilter(filter.name="Node Score Filter", column = "nodeScore", c(NSLB, max(V(FM0)$nodeScore)), predicate = "BETWEEN")
  createColumnFilter(filter.name = "Edge Weight Filter", column = "weight", c(-EWLB, EWLB), predicate = "IS_NOT_BETWEEN", type = "edges")
  createCompositeFilter(filter.name = "Edge Weight & Node Score Filter", c("Node Score Filter", "Edge Weight Filter"), "ANY")
  createSubnetwork(subnetwork.name = paste(name, "Filtered Network", sep=""))
  
}




#' @export
getDataFiles <- function(SData, Ps, cnames = NULL)
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
    networkList[[paste("N", as.character(i),sep="")]] <- getEdgeWeightsF(f,countTable)
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
      SendToCSV(FMList[[i]],networkList[[i]],CList[[i]], name= cnames[i])
      
    }
    else
    {
      SendToCSV(FMList[[i]],networkList[[i]],CList[[i]], name= paste("C", as.character(i-1), sep=""))
    }
    
    i = i+1
  }
  
}







#' @export
getCyGraphs <- function(SData, Ps,NSLB = 100, EWLB = 0.05, cnames = NULL)
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










#Function Calls


# Data <- getFile("GSM4008666_Adult-Uterus1_dge.txt")
# SData <- Preprocess(Data)
# PValues <- getPValues(SData)
# getDataFiles(SData, PValues)
# getCyGraphs(SData, PValues, 100, 0.05)
