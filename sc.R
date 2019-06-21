#install.packages(c("ROCR","ggplot2","Hmisc","reshape","gplots","stringr","NMF","mixtools",
#"lars","reshape2","vioplot","fastICA","tsne","Rtsne","fpc","ape","VGAM","gdata",
#"knitr","useful","igraph","rARPACK","igraph","Rmisc", "gridExtra","randomForest","Matrix",
#"pracma","RcppAnnoy","RSpectra","ggjoy", "ggridges","pbapply","gmodels","RANN", "shape", "pvclust",
#"MASS","cowplot","rgl","networkD3","irlba","fpc","VGAM","ape","ggrepel"))
#"XLConnect"

#source("https://bioconductor.org/biocLite.R")
#biocLite(c("preprocessCore","Biostrings","GraphAlignment","MAST","tweeDEseq"))
# sourcefiles
source(paste0(dirpath, "UpdateObject.R"))
#source("./UpdateObject.R")

# Base Packages
#library(grDevices)
list.of.packages <- c("ggplot2","reshape","gplots","stringr","reshape2","vioplot","tsne","fpc","Rtsne","ape","VGAM","gdata","knitr","useful","dplyr","rARPACK","igraph", "Rmisc","gridExtra","Matrix","Rcpp","pracma","RcppAnnoy","RSpectra","pbapply","shape","cowplot","irlba","randomForest","xgboost","BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# CRAN
#require(ROCR)
#require(scde)
require(ggplot2)
#require(Hmisc)
require(reshape)
require(gplots)
require(stringr)
#require(NMF)
#require(mixtools)
#require(lars)
#require(XLConnect)
require(reshape2)
require(vioplot)
#require(fastICA)
require(tsne)
require(fpc)
require(Rtsne)
require(ape)
require(VGAM)
require(gdata)
require(knitr)
require(useful)
require(dplyr)
require(rARPACK)
require(igraph)
require(Rmisc)
require(gridExtra)
#require(randomForest)
require(Matrix)
library(Rcpp)
library(pracma)
library(RcppAnnoy)
library(RSpectra)
#library(ggjoy)
#library(ggridges)
require(pbapply)
#library(gmodels)
#require(RANN)
#library(jackstraw)
library(shape)
#library(pvclust)
library(cowplot)
#library(MASS)
library(irlba)

#library(rgl)
#require(networkD3)

# Bioconductor pacakges
# require(Biostrings)
# library(GraphAlignment)
# require(preprocessCore)
# require(MAST)
# require(tweeDEseq)
# library(scran)
# require(rfFC)

#' scR class object

scR <- methods::setClass("scR", 
                      slots = c(
                        data = "ANY", 
                        count.data="ANY", 
                        regress.data="ANY",
                        scale.data="ANY",
                        var.genes="vector",
                        is.expr="numeric",
                        ident="vector",
                        data.info="data.frame",
                        project.name="character",
                        dr = "list",
                        de.list = "list",
                        de.union = "data.frame",
                        hvg.info = "data.frame",
                        imputed = "ANY",
                        cell.names = "vector",
                        cluster.tree = "list",
                        adjmat = "dgCMatrix",
                        calc.params = "list",
                        kmeans = "ANY",
                        spatial = "ANY",
                        misc = "ANY",
                        clust.avg.exp = "ANY",
                        clust.avg.perc = "ANY",
                        ident.fxn = "function",
                        version = "ANY",
                        data.ngene="vector",
                        pca.x="data.frame",
                        pca.rot="data.frame",
                        pca.var="numeric",
                        pca.obj="list",
                       tsne.rot="data.frame", 
                       ica.rot="data.frame", 
                       ica.x="data.frame",
                       ica.obj="list"
))

#' @include seurat.R
NULL
# Set up dim.reduction class

dim.reduction <- setClass(
  Class = "dim.reduction",
  slots = list(
    cell.embeddings = "matrix",
    gene.loadings = "matrix",
    gene.loadings.full = "matrix",
    dim.var = "numeric",
    eigvals = "numeric",
    key = "character",
    jackstraw="ANY",
    misc = "ANY"
  )
)

#' setting up a new object
setup <- function(
  object,
  project = "NewProject",
  min.cells = 0,
  min.genes = 500,
  max.genes = 10000,
  normalization.method = "logTPM",
  is.expr = 0,
  scale.factor=NULL,
  meta.data = NULL,
  pseudocount.use = 1,
  save.raw = TRUE,
  threshold.use = NULL,
  threshold.quantile = 0.995,
  version = "Oct2017",
  verbose=TRUE
){
  
  object@is.expr = is.expr
  object@version = version
  # Filter on the number of genes detected
  count.data <- object@count.data
  if (is.expr > 0) {
    suppressMessages(count.data[count.data < is.expr] <- 0)
  }
  
  num.genes = Matrix::colSums(count.data > object@is.expr)
  num.trans = Matrix::colSums(count.data)
  
  cells.use = names(num.genes[which(num.genes>min.genes & num.genes < max.genes)]) 
  object@count.data = object@count.data[,cells.use]
  object@data = count.data[,cells.use]
  
  # Filter genes on the number of cells expressing
  genes.use = rownames(object@data)
  if (min.cells > 0){
    num.cells=Matrix::rowSums(object@data > is.expr)
    genes.use = names(num.cells[which(num.cells>min.cells)])
    object@data <- object@data[genes.use,]
  }
  
  # Group identity of cells
  object@ident=factor(unlist(lapply(colnames(object@data),object@ident.fxn)))
  names(object@ident)=colnames(object@data)
  
  object@cell.names=names(object@ident)
  # if there are more than 100 idents, set all idents to project name
  ident.levels <- length(x = unique(x = object@ident))
  if ((ident.levels > 100 || ident.levels == 0) || ident.levels == length(x = object@ident)) {
    object@ident <- factor(rep(project, length(object@cell.names))); names(object@ident) = object@cell.names
  }
  
  object@data.ngene=num.genes[cells.use]
  object@data.info=data.frame(object@data.ngene); colnames(object@data.info)[1]="nGene"
  num.trans = num.trans[cells.use]
  object@data.info = cbind(object@data.info,num.trans)
  colnames(object@data.info)[2] = "nTranscripts"
  
  object@data.info[names(object@ident), "orig"] <- object@ident
  if (!is.null(meta.data)) {
    object <- AddMetaData(object=object, metadata=meta.data)
  }
  
  if (!is.null(normalization.method)){
    
    object <- NormalizeData(object=object,
                            normalization.method = normalization.method,
                            scale.factor=scale.factor,
                            pseudocount.use = pseudocount.use,
                            verbose = verbose)
  } else {
    if (verbose) print(paste0("No normalization. Log-transforming after adding ", pseudocount.use))
    object@data = log(object@data + pseudocount.use)
  }
  
  if (!is.null(threshold.use) | !is.null(threshold.quantile)){
    if (!is.null(threshold.use)){
      if (verbose) print(paste0("Thresholding all normalized values to ", threshold.use))
      valsToChange = Matrix::which(object@data > threshold.use)
      object@data[valsToChange] = threshold.use
    } else {
      pos.vals = object@data[object@data > 0]
      threshold.use = quantile(pos.vals, threshold.quantile)
      if (verbose) print(paste0("Thresholding all normalized values to ", 100*threshold.quantile, " percentile of the expressed vals =", threshold.use))
      valsToChange = Matrix::which(object@data > threshold.use)
      object@data[valsToChange] = threshold.use
    }
  }
  
  # spatial.obj <- new(
  #   Class = "spatial.info",
  #   mix.probs = data.frame(nGene)
  # )
  
  #object@spatial <- spatial.obj
  
  # to save memory downstream, especially for large objects if raw.data no
  # longer needed
  if (!(save.raw)) {
    object@count.data <- matrix()
  }
  
  # parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("CreateSeuratObject"))]
  # parameters.to.store$raw.data <- NULL
  # parameters.to.store$meta.data <- NULL
  # object <- SetCalcParams(object = object,
  #                         calculation = "CreateSeuratObject",
  #                         ... = parameters.to.store)
  
  return(object)
  
}

#' Normalize Assay Data
#'
#' Normalize data for a given assay
#'
#' @param object scR object
#' @param normalization.method Method for normalization. Default is
#' TPM
#' @param scale.factor Sets the scale factor for cell-level normalization. If null,median of all counts is used
#' @param display.progress display progress bar for scaling procedure.
#'
#' @return Returns object after normalization. Normalized data is stored in data
#' or scale.data slot, depending on the method
#'
#' @export
#'
#' @examples
#' pbmc_small
#' pmbc_small <- NormalizeData(object = pbmc_small)
#'
NormalizeData <- function(
  object,
  normalization.method = "logTPM",
  scale.factor = NULL,
  pseudocount.use = 1,
  verbose = TRUE
) {
  #parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("NormalizeData"))]
  # object <- SetCalcParams(
  #   object = object,
  #   calculation = "NormalizeData",
  #   ... = parameters.to.store
  # )
  
  if(is.null(normalization.method)) {
    cat("No normalization is performed \n", file=stderr())
    return(object)
  }
  if (normalization.method == "logTPM") {
    data <- object@data
    if (is.null(data)) {
      stop(paste("Data slot has not been set has not been set"))
    }
    
    if (is.null(scale.factor)){
      num.trans = Matrix::colSums(data)
      scale.factor = median(num.trans)
      if (verbose) print(paste0("Performing median normalization - normalizing each cell to transcript count ", scale.factor))
    } else {
      if (verbose) print(paste0("Performing median normalization - normalizing each cell to transcript count ", scale.factor))
      num.trans = Matrix::colSums(data)
    }
    
    normalized.data <-  scale.factor * t(t(data) / num.trans)
    rm(data)
    if (verbose) print(paste0("Log-transforming TPM values after adding ", pseudocount.use))
    normalized.data <- log(normalized.data + pseudocount.use)
  }
  # if (normalization.method == "genesCLR") {
  #   raw.data <- GetAssayData(
  #     object = object,
  #     assay.type = assay.type,
  #     slot = "raw.data"
  #   )
  #   if (is.null(x = raw.data)) {
  #     stop(paste("Raw data for", assay.type, "has not been set"))
  #   }
  #   normalized.data <- CustomNormalize(
  #     data = raw.data,
  #     custom_function = function(x) log1p((x)/(exp(sum(log1p((x)[x > 0]), na.rm=TRUE) / length(x+1)))),
  #     across = "genes"
  #   )
  #   object <- SetAssayData(
  #     object = object,
  #     assay.type = assay.type,
  #     slot = "data",
  #     new.data = normalized.data
  #   )
  # }
  
  object@data = normalized.data
  return(object)
}



calc.drop.prob=function(x,a,b) {
  return(exp(a+b*x)/(1+exp(a+b*x)))
}




weighted.euclidean=function(x,y,w) {
  v.dist=sum(sqrt(w*(x-y)^2))
  return(v.dist)
}

#from Jean Fan - thanks!!
custom.dist <- function(my.mat, my.function,...) {
  n <- ncol(my.mat)
  mat <- matrix(0, ncol = n, nrow = n)
  colnames(mat) <- rownames(mat) <- colnames(my.mat)
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      mat[i,j] <- my.function(my.mat[,i],my.mat[,j],...)
    }}
  return(as.dist(mat))
}

setGeneric("buildClusterTree", function(object, genes.use=NULL,regress.use=FALSE,pcs.use=NULL,do.plot=TRUE,do.reorder=FALSE,reorder.numeric=FALSE, stat.fxn=expMean, linkage.method="complete", dist.fun="complete", nboot.use=1000, do.bootstrap=FALSE, do.scale=FALSE) standardGeneric("buildClusterTree"))
setMethod("buildClusterTree","scR",
  function(object,genes.use=NULL,regress.use=FALSE,pcs.use=NULL,do.plot=TRUE,do.reorder=FALSE,reorder.numeric=FALSE, stat.fxn=expMean, linkage.method="complete", dist.fun="euclidean", nboot.use=1000, do.bootstrap=FALSE, do.scale=FALSE) {
            genes.use=set.ifnull(genes.use,object@var.genes)
            print(length(genes.use))
            print(regress.use)
            ident.names=as.character(unique(object@ident))
            
            if (is.null(pcs.use)) {
              if (!regress.use){
                genes.use=ainb(genes.use,rownames(object@data))
              } else {
                genes.use=ainb(genes.use,rownames(object@regress.data))
              }
              print(length(genes.use))
              # genes.use=ainb(genes.use,rownames(object@regress.data))
              
              data.avg=average.expression(object,genes.use = genes.use, fxn=stat.fxn, regress.use = regress.use)
              if (do.scale){
                print("Performing row scaling")
                data.avg = t(scale(t(data.avg), center=TRUE, scale=TRUE))
              } else{
                print("Not performing row scaling")
              }
              print(paste0("Using ", dist.fun, " distance function. Supply argument dist.fun if you want a different function"))
              if (dist.fun == "correlation"){
                data.dist = as.dist(1 - cor(data.avg[genes.use,]))
                
              } else {
                data.dist=dist(t(data.avg[genes.use,]), method=dist.fun)
              }
            }
            if (!is.null(pcs.use)) {
              print("Using PCs. Might be prone to errors, check and modify code")
              data.pca=average.pca(object)
              data.use=data.pca[pcs.use,]
              data.eigenval=(object@pca.obj[[1]]$sdev)^2
              data.weights=(data.eigenval/sum(data.eigenval))[pcs.use]; data.weights=data.weights/sum(data.weights)
              data.dist=custom.dist(data.pca[pcs.use,],weighted.euclidean,data.weights)
            }
            
            if (!do.bootstrap){
              print(paste0("Performing non-bootstrap hierarchical clustering using ", linkage.method, " linkage method. Supply argument linkage.method if you want a different function"))
              data.tree=as.phylo(hclust(data.dist, method=linkage.method))
              object@cluster.tree[[1]]=data.tree
            } else {
              if (is.null(pcs.use)){
                data.use = data.avg[genes.use,]
              } else {
                data.use = data.pca[pcs.use,]
              }
              print(paste0("Performing bootstrapped hierarchical clustering (pvclust) using ", dist.fun, " distance function and ",  
                           linkage.method, " linkage method. Supply arguments dist.fun and linkage.method if you want a different function"))
              print(paste0("Employing ", nboot.use, " bootstrap trials. Supply argument nboot if you wish to modify"))
              data.tree = pvclust(data.use, method.dist=dist.fun, method.hclust = linkage.method, nboot = nboot.use,r=seq(.5,1.4, by=.1), use.cor="all.obs", iseed=1234)
              object@cluster.tree[[1]]=data.tree
            }
            
            
            if (do.reorder) {
              print("Reordering tree. Currently incompatible with do.bootstrap=TRUE")
              old.ident.order=sort(unique(object@ident))
              data.tree=object@cluster.tree[[1]]
              all.desc=getDescendants(data.tree,data.tree$Nnode+2); all.desc=old.ident.order[all.desc[all.desc<=(data.tree$Nnode+1)]]
              object@ident=factor(object@ident,levels = all.desc,ordered = TRUE) 
              if (reorder.numeric) {
                object=set.ident(object,object@cell.names,as.integer(object@ident))
                object@data.info[object@cell.names,"tree"]=as.integer(object@ident)
              }
              object=buildClusterTree(object,genes.use,pcs.use,do.plot=FALSE,do.reorder=FALSE, stat.fxn=stat.fxn, linkage.method=linkage.method)
            }
            
            if (do.plot){
              if (!do.bootstrap){
                plotClusterTree(object)
              } else {
                plot(object@cluster.tree[[1]])
              }
            }
            return(object)
          }
)


setGeneric("plotClusterTree", function(object) standardGeneric("plotClusterTree"))
setMethod("plotClusterTree","scR",
          function(object) {
            data.tree=object@cluster.tree[[1]]
            plot(data.tree,direction="downwards")
            nodelabels()
          }
)





setGeneric("reorder.by.tree", function(object, genes.use=NULL,do.numeric=FALSE,...) standardGeneric("reorder.by.tree"))
setMethod("reorder.by.tree","scR",
          function(object,genes.use=NULL,do.numeric=FALSE...) {
            old.ident.order=sort(unique(object@ident))
            data.tree=object@cluster.tree[[1]]
            all.desc=getDescendants(data.tree,data.tree$Nnode+2); all.desc=old.ident.order[all.desc[all.desc<=(data.tree$Nnode+1)]]
            object@ident=factor(object@ident,levels = all.desc,ordered = TRUE) 
            if (do.numeric) {
              set.ident(object,object@cell.names,as.integer(object@ident))
            }
            return(object)
            
          }
)
setGeneric("plotNoiseModel",function(object, cells.use=NULL, col.use="black",lwd.use=2,do.new=TRUE,x.lim=10,...) standardGeneric("plotNoiseModel"))
setMethod("plotNoiseModel","scR",
          function(object, cells.use=NULL, col.use="black",lwd.use=2,do.new=TRUE,x.lim=10,...) {
            cells.use = set.ifnull(cells.use, rownames(object@drop.coefs))
            cells.use=cells.use[cells.use %in% rownames(object@drop.coefs)]
            cell.coefs=object@drop.coefs[cells.use,]
            if (length(cells.use) == 1){
              main = cells.use
              data=data.frame(object@data[object@trusted.genes,cells.use]); colnames(data) = cells.use;
              idents=data
              code_humpAvg=object@pop.avg
              code_humpAvg[code_humpAvg>9]=9
              code_humpAvg[is.na(code_humpAvg)]=0
              idents$code_humpAvg=code_humpAvg
              data[data>=object@drop.expr]=1
              data[data<object@drop.expr]=0
              data$bin=cut(code_humpAvg,20)
              data$avg=code_humpAvg
              rownames(idents)=rownames(data)
             getAB(cn=cells.use, data=data,data2=idents,code2="humpAvg",hasBin=TRUE,doPlot=TRUE)
            } else {
              main = "Single Cell FNR curves"
              if (do.new) plot(1,1,pch=16,type="n",xlab="Average expression",ylab="P(detection) or 1-FNR",xlim=c(0,x.lim),ylim=c(0,1), main=main)
              unlist(lapply(1:length(cells.use), function(y) {
                x.vals=seq(0,x.lim,0.05)
                y.vals=unlist(lapply(x.vals,calc.drop.prob,cell.coefs[y,1],cell.coefs[y,2]))
                lines(x.vals,y.vals,lwd=lwd.use,col=col.use,...)
              }))
            }
            
          }
)




#needs mod and RD

subsetData = function(
  object, 
  subset.name=NULL, 
  cells.use=NULL,
  ident.use = NULL,
  ident.remove = NULL,
  accept.low=-Inf, 
  accept.high=Inf,
  max.cells.per.ident=Inf,
  random.seed = 1,
  ...) {
  
  data.use <- NULL
  cells.use <- set.ifnull(cells.use, object@cell.names)
  if (!is.null(ident.use)){
    ident.use <-setdiff(ident.use, ident.remove)
    cells.use <- which.cells(object, ident.use)
  }
  
  if ((is.null(ident.use)) && ! is.null(ident.remove)) {
    ident.use <- setdiff(unique(object@ident), ident.remove)
    cells.use <- WhichCells(object, ident.use)
  }
  
  if (! is.null(subset.name)) {
    data.use <- fetch.data(object, subset.name, ...)
    if (length(x = data.use) == 0) {
      print("Error: Data of size 0 obtained")
      return(object)
    }
    subset.data <- data.use[, subset.name]
    pass.inds <- which(x = (subset.data > accept.low) & (subset.data < accept.high))
    cells.use <- rownames(data.use)[pass.inds]
  }
  
  cells.use <- intersect(cells.use,  object@cell.names)
  # To add downsample cells by uniformly sampling across ident
 # cells.use <-  WhichCells(
#    object = object,
#    cells.use = cells.use,
#    max.cells.per.ident = max.cells.per.ident,
#    random.seed = random.seed
#  )
  
  
  object@data=object@data[,cells.use]

  object@count.data=object@count.data[,cells.use]
  object@regress.data = object@regress.data[,cells.use]
  object@data.ngene=object@data.ngene[cells.use]
  object@ident=drop.levels(object@ident[cells.use])
  object@data.info = object@data.info[cells.use,]
  object@data.info$orig = drop.levels(object@data.info[cells.use,"orig"])
  object@tsne.rot=object@tsne.rot[cells.use,]
  object@pca.rot=object@pca.rot[cells.use,]
  object@ica.rot=object@ica.rot[cells.use,]
  object@cell.names=cells.use
  
  if (length(object@dr) > 0) {
    for (i in 1:length(object@dr)) {
      if(length(object@dr[[i]]@cell.embeddings) > 0){
        object@dr[[i]]@cell.embeddings <- object@dr[[i]]@cell.embeddings[cells.use, ,drop = FALSE]
      }
    }
  }

# Need to handle this smartly    
 # if (length(object@de.list) > 0) {
  #  for (i in 1:length(object@de.list)) {
   #   if(length(object@dr[[i]]@cell.embeddings) > 0){
  #      object@dr[[i]]@cell.embeddings <- object@dr[[i]]@cell.embeddings[cells.use, ,drop = FALSE]
  #    }
 #   }
 # }
  
 return(object)
          
}         



#needs mod and RD
setGeneric("combineObjects",  function(object1,object2, cells1.use=NULL, cells2.use=NULL, tag1=NULL, tag2=NULL,min.cells.use=NULL,min.genes.use=NULL, threshold.quantile=0.999, min.count.use=NULL, is.expr.use=NULL, do.renorm=TRUE,sep="_", add.tag.to.cell.names=TRUE, ident.transfer = "orig",current.ident.name="type",use.gene.union = TRUE, ...) standardGeneric("combineObjects"))
setMethod("combineObjects","scR",
          function(object1,object2, cells1.use=NULL, cells2.use=NULL, tag1=NULL, tag2=NULL, min.cells.use=NULL,min.genes.use=NULL,threshold.quantile=0.999, min.count.use=NULL, is.expr.use=NULL, do.renorm=TRUE,sep="_",add.tag.to.cell.names=TRUE, ident.transfer="orig",current.ident.name="type",use.gene.union = TRUE, ...) {
            tag1 = set.ifnull(tag1, "object1")
            tag2 = set.ifnull(tag2, "object2")
            cells1.use = set.ifnull(cells1.use, colnames(object1@count.data))
            cells2.use = set.ifnull(cells2.use, colnames(object2@count.data))
            if (use.gene.union){
              genes.use = union(rownames(object1@count.data), rownames(object2@count.data))
            } else {
              genes.use = intersect(rownames(object1@count.data), rownames(object2@count.data))
            }
            
            
            if (add.tag.to.cell.names){
              cells1.names = paste0(tag1,sep,cells1.use)
              cells2.names = paste0(tag2,sep,cells2.use)
            } else {
              cells1.names = cells1.use
              cells2.names = cells2.use
            }
                  
            
            
            min.cells.use=set.ifnull(min.cells.use, 20)
            min.genes.use=set.ifnull(min.genes.use, 0.8*min(c(object1@data.info[cells1.use,"nGene"],object2@data.info[cells2.use,"nGene"])))
            min.count.use=set.ifnull(min.count.use, 10)
            is.expr.use=set.ifnull(is.expr.use, min(object1@is.expr, object2@is.expr))
            
            Count.mat1 = object1@count.data[,cells1.use]
            Count.mat2 = object2@count.data[,cells2.use]
            #Count data (assumes both count matrices have same number of rows)
            genes.use1 = intersect(genes.use, rownames(Count.mat1))
            genes.use2 = intersect(genes.use, rownames(Count.mat2))
            
            genes.useNot1 = setdiff(genes.use, genes.use1)
            genes.useNot2 = setdiff(genes.use, genes.use2)

            if (length(genes.useNot1) > 0){
              print(1)
              shamMatrix1 = Matrix(0,nrow = length(genes.useNot1), ncol=ncol(Count.mat1))
              colnames(shamMatrix1) = colnames(Count.mat1)
              rownames(shamMatrix1) = genes.useNot1
              Count.mat1 = rbind(Count.mat1, shamMatrix1)
            }
            
            if (length(genes.useNot2) > 0){
              print(2)
              shamMatrix2 = Matrix(0,nrow = length(genes.useNot2), ncol=ncol(Count.mat2))
              colnames(shamMatrix2) = colnames(Count.mat2)
              rownames(shamMatrix2) = genes.useNot2
              Count.mat2 = rbind(Count.mat2, shamMatrix2)
            }
            
            Count.mat1 = Count.mat1[genes.use,]
            Count.mat2 = Count.mat2[genes.use,]
            Count.mat = cbind(Count.mat1[, cells1.use], Count.mat2[,cells2.use])
            colnames(Count.mat) = c(cells1.names, cells2.names)



            
            # Create new object
            object=scR(count.data=Count.mat,ident.fxn=getStat1)
            object =setup(object,project="fovea10X",min.cells = min.cells.use,min.genes = min.genes.use,is.expr=is.expr.use, threshold.quantile = threshold.quantile)
            
            if (!do.renorm){
              print("Reverting to non-normalized values")
              genes.use = rownames(object@data)[rownames(object@data) %in% rownames(object1@data)]
              genes.use = genes.use[genes.use %in% rownames(object2@data)]
              logTPM1 = object1@data[rownames(object@data), cells1.use]
              logTPM2 = object2@data[rownames(object@data), cells2.use]
              logTPM = cbind(logTPM1, logTPM2); rm(logTPM1); rm(logTPM2)
              colnames(logTPM) = c(cells1.names, cells2.names)
              object@data = logTPM
            } 
            # orig sample
            for (field in ident.transfer){
              orig_sample_id = c(paste0(tag1,"_", as.character(object1@data.info[cells1.use, field])),
                                 paste0(tag2,"_", as.character(object2@data.info[cells2.use, field])))
              names(orig_sample_id) = c(cells1.names, cells2.names)
              object@data.info[,paste0(field,"_old")] = factor(orig_sample_id)
            }
           
            
            # ident
            ident1 = paste0(tag1,"_", as.character(object1@ident))
            ident2 = paste0(tag2,"_", as.character(object2@ident))
            ident.new = c(ident1, ident2); names(ident.new) = c(cells1.names, cells2.names)
            ident.new = factor(ident.new)
            object@ident = ident.new
            object@data.info[,current.ident.name] = ident.new
            
            # Add identities
           

            return(object)
          }         
)

#needs mod and RD, probably remove
setGeneric("loadMetrics", function(object, metrics.file=NULL, col.names=NULL,sep.use="\t",row.add="_rsem",...) standardGeneric("loadMetrics"))
setMethod("loadMetrics","scR",
          function(object, metrics.file=NULL, col.names=NULL,sep.use="\t",row.add="_rsem",...) {
            metrics.file=set.ifnull(metrics.file,paste(object@project.dir,"summary/",object@project.name,".all.aln.metrics.txt",sep=""))
            col.names=set.ifnull(col.names,c("n.read","n.aln.read","pct.aln.rsem","pct.aln.ribo","pct.aln.spike","pct.aln.coli"))
            metrics.data=read.table(metrics.file,sep=sep.use,row.names="V1",...)
            colnames(metrics.data)=col.names
            rownames(metrics.data)=paste(sub.string(rownames(metrics.data),"-","_"),row.add,sep="")
            object@data.info=cbind(object@data.info,metrics.data[rownames(object@data.info),])
            return(object)
          }
)

setGeneric("project.pca", function(object,do.print=TRUE,pcs.print=5,pcs.store=20,genes.print=30,replace.pc=FALSE) standardGeneric("project.pca"))
setMethod("project.pca", "scR", 
          function(object,do.print=TRUE,pcs.print=5,pcs.store=20,genes.print=30,replace.pc=FALSE) {
            object@pca.x.full=data.frame(as.matrix(object@scale.data)%*%as.matrix(object@pca.rot))
            if (ncol(object@jackStraw.fakePC)>0) {
              object@jackStraw.empP.full=data.frame(sapply(1:ncol(object@jackStraw.fakePC),function(x)unlist(lapply(abs(object@pca.x.full[,x]),empP,abs(object@jackStraw.fakePC[,x])))))
              colnames(object@jackStraw.empP.full)=paste("PC",1:ncol(object@jackStraw.empP),sep="")
              rownames(object@jackStraw.empP.full)=rownames(object@scale.data)
            }
            
            if (replace.pc==TRUE) {
              object@jackStraw.empP=object@jackStraw.empP.full
              object@pca.x=object@pca.x.full
            }
            
            if (do.print) {
              pc_scores=object@pca.x.full
              for(i in 1:pcs.print) {
                code=paste("PC",i,sep="")
                sx=pc_scores[order(pc_scores[,code]),]
                print(code)
                print(rownames(sx[1:genes.print,]))
                print ("")
                
                print(rev(rownames(sx[(nrow(sx)-genes.print):nrow(sx),])))
                print ("")
                print ("")      
              }
            }
            return(object)
          }
)

run_tsne = function(object,cells.use=NULL,pcs.use=1:10,k.seed=1,do.fast=FALSE,add.iter=0, fast.method="BH",reduction.use="pca",reduction.key="tsne",reduction.name="tsne",...)  {
  cells.use=set.ifnull(cells.use,colnames(object@data))
  if (reduction.use=="pca") data.use=object@pca.rot[cells.use,pcs.use]
  if (reduction.use=="ica") data.use=object@ica.rot[cells.use,pcs.use]
  #data.dist=as.dist(mahalanobis.dist(data.use))
  
  if (!(do.fast)) {
    require(tsne)
    set.seed(k.seed); data.tsne=data.frame(tsne(data.use,whiten=FALSE,...))
  } else {
    
    if (fast.method == "BH"){
      require(Rtsne)
      set.seed(k.seed); data.tsne=Rtsne(as.matrix(data.use),pca=FALSE,verbose=TRUE,...)
      data.tsne=data.frame(data.tsne$Y)
    } else if (fast.method == "FFT"){
      source(paste0(dirpath,fast_tsne_script), chdir = T)
      dir.create("temp")
      data.tsne <- fftRtsne(as.matrix(data.use), ...);
    } else {
      stop("Error: fast.method is not valid. Choose either BH or FFT")
    }
    
    
  }
  if (add.iter > 0) {
    data.tsne=data.frame(tsne(data.use,initial_config = as.matrix(data.tsne),max_iter = add.iter,...))
  }
  
  colnames(data.tsne)=paste("tsne_",1:ncol(data.tsne),sep="")
  rownames(data.tsne)=cells.use
  object@tsne.rot=as.data.frame(data.tsne)
  
  tsne.obj <- new(
    Class = "dim.reduction",
    gene.loadings = matrix(rand(10),2,5),
    cell.embeddings = as.matrix(data.tsne),
    dim.var = c(1,2,3),
    key = reduction.key
  )
  
  eval(expr = parse(text = paste0("object@dr$", reduction.name, "<- tsne.obj")))
  
  return(object)
}




setGeneric("add_tsne", function(object,cells.use=NULL,pcs.use=1:10,do.plot=TRUE,k.seed=1,add.iter=1000,...) standardGeneric("add_tsne"))
setMethod("add_tsne", "scR", 
          function(object,cells.use=NULL,pcs.use=1:10,do.plot=TRUE,k.seed=1,add.iter=1000,...) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            data.use=object@pca.rot[cells.use,pcs.use]
            #data.dist=as.dist(mahalanobis.dist(data.use))
            set.seed(k.seed); data.tsne=data.frame(tsne(data.use,initial_config = as.matrix(object@tsne.rot[cells.use,]),max_iter = add.iter,...))
            colnames(data.tsne)=paste("tsne_",1:ncol(data.tsne),sep="")
            rownames(data.tsne)=cells.use
            object@tsne.rot=data.tsne
            return(object)
          }
)


setGeneric("ica", function(object,ic.genes=NULL,do.print=TRUE,ics.print=5,ics.store=25,genes.print=30,use.imputed=FALSE,use.scaled=TRUE,use.regress=FALSE,reduction.name="ica",reduction.key="IC",...) standardGeneric("ica"))
setMethod("ica", "scR", 
          function(object,ic.genes=NULL,do.print=TRUE,ics.print=5,ics.store=25,genes.print=30,use.imputed=FALSE,use.scaled=TRUE,use.regress=FALSE,reduction.name="ica",reduction.key="IC",...) {
            
            if (use.regress){
              if (object.size(object@regress.data) < 100){
                print("Warning : Regression has not been done. Using object@data")
                data.use = object@data
              } else {
                data.use = object@regress.data
              }
            } else {
              data.use = object@data
            }
            
            ic.genes=set.ifnull(ic.genes,object@var.genes)
            ic.genes = ic.genes[ic.genes%in%rownames(data.use)]
            if (use.scaled & !use.imputed){
              if (object.size(object@scale.data) < 1000){ 
                data.use=t(scale(t(object@data[ic.genes,]))) 
              } else {
                print("using scale.data that has already been computed")
                data.use = object@scale.data[ic.genes,]
              }
              
            }
            if (!use.scaled & !use.imputed) data.use = t(scale(t(object@data[ic.genes,]), center=TRUE, scale=FALSE))
            if (use.scaled & use.imputed) data.use=t(scale(t(object@imputed[ic.genes,])))
            if (!use.scaled & use.imputed) data.use=t(scale(t(object@imputed[ic.genes,]), center=TRUE, scale=FALSE))
            
            ic.genes.var = apply(data.use[ic.genes,],1,var)
            ic.data = data.use[ic.genes[ic.genes.var>0],]
            ica.obj = fastICA(t(ic.data),n.comp=ics.store,...)
            object@ica.obj=list(ica.obj)
            ics.store=min(ics.store,nrow(ic.data))
            ics.print=min(ics.print,nrow(ic.data))
            ic_scores=data.frame(ic.data%*%ica.obj$S)
            colnames(ic_scores)=paste("IC",1:ncol(ic_scores),sep="")
            object@ica.x=ic_scores
            object@ica.rot=data.frame(ica.obj$S[,1:ics.store])
            colnames(object@ica.rot)=paste("IC",1:ncol(object@ica.rot),sep="")
            print(head(object@ica.rot))
            rownames(object@ica.rot)=colnames(ic.data)
            
            ica.obj2 <- new(
              Class = "dim.reduction",
              gene.loadings = as.matrix(object@ica.x),
              cell.embeddings = as.matrix(object@ica.rot),
              dim.var = 1,
              eigvals = 1,
              key = reduction.key
            )
            
            eval(expr = parse(text = paste0("object@dr$", reduction.name, "<- ica.obj2")))
            
            if (do.print) {
              for(i in 1:ics.print) {
                code=paste("IC",i,sep="")
                sx=ic_scores[order(ic_scores[,code]),]
                print(code)
                print(rownames(sx[1:genes.print,]))
                print ("")
                
                print(rev(rownames(sx[(nrow(sx)-genes.print):nrow(sx),])))
                print ("")
                print ("")      
              }
            }
            return(object)
          }
)


#' weight.var - if TRUE, then scores are weighted by their singular values, if not not
pca2 = function(object,pc.genes=NULL,use.scaled=TRUE, use.col.scaled = FALSE,do.print=TRUE,pcs.print=5,pcs.store=40,genes.print=30,use.regress=TRUE,weight.var=TRUE,reduction.name="pca",reduction.key="PC",...) {
            
            if (use.regress){
              if (object.size(object@regress.data) < 100){
                print("Warning : Regression has not been done. Using object@data")
                data.use = object@data
              } else {
                data.use = object@regress.data
              }
            } else {
              data.use = object@data
            }
  
            pc.genes=set.ifnull(pc.genes,object@var.genes)
            pc.genes = pc.genes[pc.genes%in%rownames(data.use)]
            pc.genes.var = apply(data.use[pc.genes,],1,function(x) var(x))
            genes.use = pc.genes[pc.genes.var>0]
            
            print(paste0("Performing PCA based on ", length(genes.use), " variable genes"))
            if (use.scaled){
              if (object.size(object@scale.data) < 1000){ 
                if (!use.col.scaled){
                  data.use=t(scale(t(data.use[genes.use,]))) 
                } else {
                  print("Scaling cells instead of genes")
                  data.use = scale(data.use[genes.use,])
                }
                
              } else {
                print("Warning: using scale.data that has already been computed. set use.scaled=FALSE if you wish otherwise")
                data.use = object@scale.data[genes.use,]
                }
            }
            if (!use.scaled) data.use = t(scale(t(data.use[genes.use,]), center=TRUE, scale=FALSE))

            
            #pca.obj = irlba(t(data.use),n=pcs.store,center = FALSE, scale.=FALSE, ...)
            svd.obj = irlba(t(data.use),nv = pcs.store, nu = pcs.store, center=FALSE, scale = FALSE, ...)
            
            eigvals = svd.obj$d^2 # eigenvalues
            pca.rotation = svd.obj$v # the rotation matrix (genes x PCs)
            
            if (weight.var == TRUE){
              # Weight by eigenvalues
              pca.x = svd.obj$u %*% diag(svd.obj$d)
            } else {
              pca.x = svd.obj$u
            }
            
            rownames(pca.x) = colnames(data.use)
            colnames(pca.x) = paste0("PC", 1:pcs.store)
            rownames(pca.rotation) = rownames(data.use)
            colnames(pca.rotation) = paste0("PC", 1:pcs.store)
            
            pca.obj <- new(
              Class = "dim.reduction",
              gene.loadings = as.matrix(pca.rotation),
              cell.embeddings = as.matrix(pca.x),
              dim.var = 100*eigvals / sum(eigvals),
              eigvals = eigvals,
              key = reduction.key
            )
            
           # pca.obj$eigvals = eigvals # TO DO: include a line that computes variances if the user wishes so
            #pca.obj$x = pca.x
            #pca.obj$rotation = pca.rotation
            #object@pca.var = 100*pca.obj$eigvals / sum(pca.obj$eigvals)  
            
            pcs.store=min(pcs.store,ncol(pca.x))
            pcs.print=min(pcs.print,ncol(pca.x))
            object@pca.rot=data.frame(pca.x[,1:pcs.store])
            rownames(object@pca.rot) = colnames(object@data)
            object@pca.x=data.frame(pca.rotation[,1:pcs.store])
            rownames(object@pca.x) = genes.use
            
            eval(expr = parse(text = paste0("object@dr$", reduction.name, "<- pca.obj")))
            
            if (do.print) {
              pc_scores=object@pca.x
              for(i in 1:pcs.print) {
                code=paste("PC",i,sep="")
                sx=pc_scores[order(pc_scores[,code]),]
                print(code)
                print(rownames(sx[1:genes.print,]))
                print ("")
                
                print(rev(rownames(sx[(nrow(sx)-genes.print):nrow(sx),])))
                print ("")
                print ("")      
              }
            }
            return(object)
         
}



# Diffusion map
DiffMap <- function(object,dims.use=1:20, reduction.type="pca", cells.use=NULL, dim.store=20,reduction.key="DC", reduction.name="dmap", ...) {
            require(destiny)
            cells.use = set.ifnull(cells.use, colnames(object@data))
            data.use=GetCellEmbeddings(object, 
                                       reduction.type = reduction.type, 
                                       dims.use = dims.use, 
                                       cells.use = cells.use)
            data.use <- as.ExpressionSet(as.data.frame(data.use))
            dm <- DiffusionMap(data.use,n_eigs = dim.store, ...)
            cell.embeddings <- as.matrix(dm@eigenvectors)
            rownames(cell.embeddings) <- cells.use
            colnames(cell.embeddings) <- paste0(reduction.key,1:dim.store)
            diff.obj <- new(
              Class = "dim.reduction",
              gene.loadings = matrix(rand(20),2,5),
              cell.embeddings = cell.embeddings,
              dim.var = dm@eigenvalues,
              key = reduction.key
            )
            eval(expr = parse(text = paste0("object@dr$", reduction.name, "<- diff.obj")))
            return(object)
          
}


setGeneric("doMCL_clustering", function(object,clust.name=NULL,cells.use=NULL,pcs.use=1:10,inflation.range = seq(from=1.05,to=2,len=10),num.nn=100,num.reps=10,downsample.frac=0.9,
                                        scRPath=NULL, scRFxnsPath=NULL,...) standardGeneric("doMCL_clustering"))
setMethod("doMCL_clustering", "scR", 
          function(object,clust.name=NULL,cells.use=NULL,pcs.use=1:10,inflation.range = seq(from=1.05,to=2,len=10),num.nn=100,num.reps=10,downsample.frac=0.9,
                   scRPath=NULL, scRFxnsPath=NULL,...) {
            
            scRPath = set.ifnull(scRPath,'/ahg/regevdata/users/karthik/SOFTWARE/scR/scR.R')
            scRFxnsPath = set.ifnull(scRFxnsPath,'/ahg/regevdata/users/karthik/SOFTWARE/scR/scRFxns.R')
            cells.use=set.ifnull(cells.use,colnames(object@data))
            data.use=object@pca.rot[cells.use,pcs.use]
            clust.name=set.ifnull(clust.name, object@project.name)
            
            dir.name = paste0("MCL_output_",clust.name)
            if (file.exists(dir.name)){
              setwd(paste0("./", dir.name))
            } else {
              dir.create(dir.name)
              setwd(paste0("./", dir.name))
            }
            
            if (!file.exists("Data")) dir.create("Data")
            setwd("./Data")
            save("data.use", file="data_use.Robj")
            setwd("..")
            
            if (!file.exists("Rscripts")) dir.create("Rscripts") 
            if (!file.exists("nn_graphs")) dir.create("nn_graphs")
            if (!file.exists("clusters")) dir.create("clusters")
            if (!file.exists("bsub_logs")) dir.create("bsub_logs")
            
            run.filename=paste0("runMCL_",clust.name,"_downsample", downsample.frac,".sh")
            file.create(run.filename)
            write('#!/bin/bash', file=run.filename, append=TRUE)
            write('\n', file=run.filename, append=TRUE)
            
            for (i in inflation.range){
              job.filename=paste0("./Rscripts/MCL_",clust.name,"_inflation",round(i,2),"_downsample", downsample.frac,".R")
              file.create(job.filename)
              write(paste0('source(\'',scRPath,'\')'),file=job.filename, append=TRUE)
              write(paste0('source(\'',scRFxnsPath,'\')'), file=job.filename, append=TRUE)
              write('\n', file=job.filename, append=TRUE)
              write('load(\"./Data/data_use.Robj\")', file=job.filename, append=TRUE)
              write('clust.ident_all = c()', file=job.filename, append=TRUE)
              write('\n', file=job.filename, append=TRUE)
              
              write('#Actual run \n', file=job.filename, append=TRUE)
              filename = paste0("./nn_graphs/edges_",clust.name,"_inflation", round(i,2),"_NO_downsample_nn", num.nn,".tsv")
              filename.clusters = gsub("edges","clusters",filename)
              filename.clusters = gsub("nn_graphs","clusters",filename.clusters)
              filename.clusters = gsub("tsv","txt",filename.clusters)
              CMD = paste0("get_nn_graph(X=data.use, do.jaccard=TRUE, k=", num.nn, ", write.file = ","\"" , filename,"\"", ", downsample.frac = ", 1,")")
              write(CMD, file=job.filename, append=TRUE)
              CMD = paste0("run_MCL(nn_graph_file = ","\"" , filename,"\"" , ", out_file = ","\"" , filename.clusters,"\"" , ", inflation.use=",i,")")
              write(CMD, file=job.filename, append=TRUE)
              CMD = paste0("clust.ident = read_MCL_clusters(mcl_out_file = ", "\"", filename.clusters, "\"" , ", cells.return = rownames(data.use))")
              write(CMD, file = job.filename,append=TRUE)
              CMD = paste0("clust.ident_all = cbind(clust.ident_all, clust.ident)")
              write(CMD, file = job.filename, append=TRUE)
              CMD = paste0("colnames(clust.ident_all)[ncol(clust.ident_all)] = paste0(colnames(clust.ident_all)[ncol(clust.ident_all)], ncol(clust.ident_all)-1)")
              write(CMD, file = job.filename, append=TRUE)
              write("print(corner(clust.ident_all))", file = job.filename, append=TRUE)
              write('\n', file=job.filename, append=TRUE)
              
              write('#Downsampled runs \n', file=job.filename, append=TRUE)
              if (num.reps > 0){
                for (j in 1:num.reps){
                     filename = paste0("./nn_graphs/edges_",clust.name,"_inflation", round(i,2),"_downsample",downsample.frac, "_rep", j, "_nn", num.nn,".tsv")
                     filename.clusters = gsub("edges","clusters",filename)
                     filename.clusters = gsub("nn_graphs","clusters",filename.clusters)
                     filename.clusters = gsub("tsv","txt",filename.clusters)
                     
                     CMD = paste0("get_nn_graph(X=data.use, k=", num.nn, ", write.file = ","\"" , filename,"\"", ", downsample.frac = ", downsample.frac,")")
                     write(CMD, file=job.filename, append=TRUE)
                     
                     CMD = paste0("run_MCL(nn_graph_file = ","\"" , filename,"\"" , ", out_file = ","\"" , filename.clusters,"\"" , ", inflation.use=",i,")")
                     write(CMD, file=job.filename, append=TRUE)
                     
                     CMD = paste0("clust.ident = read_MCL_clusters(mcl_out_file = ", "\"", filename.clusters, "\"" , ", X = data.use)")
                     write(CMD, file = job.filename,append=TRUE)
                     
                     CMD = paste0("clust.ident_all = cbind(clust.ident_all, clust.ident)")
                     write(CMD, file = job.filename, append=TRUE)
                     write("print(corner(clust.ident_all))", file = job.filename, append=TRUE)
                     write('\n', file=job.filename, append=TRUE)
                }
                clust.filename = paste0("./clusters/AllClusters_",clust.name,"_inflation", round(i,2),".txt")
                CMD = paste0("write.table(clust.ident_all, file = ", "\"" , clust.filename,"\"", ", sep = \"\\t\", quote=FALSE)")
                write(CMD, file=job.filename, append=TRUE)
                
              }
              
              
              
              CMD1 = sprintf("bsub -P %s -q regevlab -e %s.err -o %s.out -R",clust.name, gsub("Rscripts","bsub_logs",job.filename), gsub("Rscripts","bsub_logs",job.filename)) 
              CMD2 = sprintf("R CMD BATCH --no-save --no-restore %s", job.filename)
              CMD = paste0(CMD1,' \"rusage[mem=50000]span[hosts=1]\" ', CMD2)     
              write(CMD, file=run.filename, append=TRUE)
              write('\n', file=run.filename, append=TRUE)
              
            }
            #Run job
            CMD = sprintf("chmod +x %s", run.filename)
            system(CMD)
            CMD = sprintf("./%s", run.filename)
            system(CMD)
            
            setwd("..")
            print("Running MCL jobs")
          }
)

setGeneric("NN_analysis", function(object,cells.use=NULL,pcs.use=1:10,nn.test=100,...) standardGeneric("NN_analysis"))
setMethod("NN_analysis", "scR", 
          function(object,cells.use=NULL,pcs.use=1:10,nn.test=100,...) {
            require(RANN)
            cells.use=set.ifnull(cells.use,colnames(object@data))
            data.use=object@pca.rot[cells.use,pcs.use]
            
            nearest=nn2(data.use,data.use,k=nn.test+1)
            nearest$nn.idx = nearest$nn.idx[,-1]
            nearest$nn.dists = nearest$nn.dists[,-1]
            thresh.use = quantile(as.numeric(nearest$nn.dists[,1:3]),0.9)
            
            num.nn = apply(nearest$nn.dists,1, function(x) sum(x <= thresh.use) )
            
            nnHist = c()
            for (i in 0:nn.test){
              
              nnHist = c(nnHist, sum(num.nn >= i))
            }
            
            plot(c(0:nn.test), nnHist / nnHist[1], xlab = "Number of nearest neighbors", ylab = "Fraction of cells")
            

          }
)

setGeneric("doLouvain_stability", function(object,n.iter=100, downsamp.frac=0.9, cells.use=NULL,pcs.use=1:10, num.nn=100,adj.full=FALSE,set.ident=TRUE, do.full.clust=FALSE,use.original=FALSE,do.jaccard=FALSE, do.noise=TRUE, do.prune=TRUE,...) standardGeneric("doLouvain_stability"))
setMethod("doLouvain_stability", "scR", function(object,n.iter=100,downsamp.frac = 0.9, cells.use=NULL,pcs.use=1:10,num.nn=100,adj.full=FALSE,set.ident=TRUE, do.full.clust=FALSE, use.original=FALSE, do.jaccard=FALSE, do.noise=TRUE,do.prune=TRUE,...) {
  
  cells.use=set.ifnull(cells.use,colnames(object@data))
  clust.assign.mat = matrix(0, nrow=n.iter, ncol=length(cells.use))
  colnames(clust.assign.mat) = cells.use
  
  inSampleMat = Matrix(0, nrow=length(cells.use), ncol=length(cells.use), sparse=TRUE)
  colnames(inSampleMat) = cells.use; rownames(inSampleMat) = cells.use
  inClusterMat = Matrix(0, nrow=length(cells.use), ncol=length(cells.use), sparse=TRUE)
  colnames(inClusterMat) = cells.use; rownames(inClusterMat) = cells.use
    
  for (i in 1:n.iter){
    print(paste0("Niter ", i))
    cells.use1 = sample(cells.use, round(downsamp.frac * length(cells.use)))
    b1 = combn(match(cells.use1, cells.use),2); b2 = b1[c(2,1),]
    inSampleMat[t(b1)] = inSampleMat[t(b1)] + 1;
    inSampleMat[t(b2)] = inSampleMat[t(b2)] + 1;
    clust.out = doLouvain_clustering(object, pcs.use=pcs.use, num.nn=sample(c(round(0.8*num.nn):round(1.2*num.nn)), 1), cells.use=cells.use1, do.jaccard=do.jaccard, return.clust=TRUE, do.noise=do.noise)
    k = length(clust.out$COM)
    clust.assign = clust.out$COM[[k]]
    names(clust.assign) = cells.use1
    clust.assign.mat[i, cells.use1] = clust.assign
    
    for (j in unique(clust.assign)){
      cells.in.j = names(which(clust.assign == j))
      if (length(cells.in.j) == 1) next
      b1 = combn(match(cells.in.j, cells.use),2); b2 = b1[c(2,1),]
      inClusterMat[t(b1)] = inClusterMat[t(b1)] + 1;
      inClusterMat[t(b2)] = inClusterMat[t(b2)] + 1;
    }
    
  }
  
  
  cell.ord = names(sort(object@ident))
  diag(inSampleMat)  = 1
  ConfusionMat = inClusterMat[cell.ord, cell.ord] / inSampleMat[cell.ord, cell.ord]

  annVal = c()
  for (i in sort(unique(object@ident))){
    annVal = c(annVal,rep(i, sum(object@ident==i)))
  }
  png("ConfusionMatrix_unfiltered_Louvain.png",h=800,w=800)
  print(aheatmap(as.matrix(ConfusionMat),Rowv=NA, Colv=NA, labRow=NULL, labCol=NULL, main="Confusion Matrix (unfiltered)"), annCol=annVal, annColors = list(rainbow(length(unique(annVal)))))
  dev.off()
  
  #Pruning
  #Compute the median stability of each cell within its cluster and within all other clusters. Remove cells that have a median stability < 0.7 in all clusters
  #If a cell has a higher stability in another cluster move it to the new cluster
  unstable_cells = rep(0, length(cell.ord));  names(unstable_cells) = cell.ord
  if (do.prune){
  clusts = sort(unique(object@ident))
  for (clust in clusts){
    cells = which.cells(object, clust)
    #Within a cluster compute the median stability of every cell within the cluster
    median_in_clust = apply(ConfusionMat[cells,cells],1, median)
    median_out_clust = matrix(0, nrow = length(clusts)-1, length(cells)); colnames(median_out_clust)=cells; rownames(median_out_clust) = setdiff(clusts,clust)
    for (clust2 in setdiff(clusts,clust)){
      cells2 = which.cells(object, clust2)
      median_out_clust[clust2,] = apply(ConfusionMat[cells,cells2],1, median)
    }
    
    max_val = apply(median_out_clust,2, max)
    max_clust = rownames(median_out_clust)[apply(median_out_clust,2, which.max)]
    names(max_clust) = colnames(median_out_clust)
    unstable = names(median_in_clust)[median_in_clust < 0.7]
    cells.to.move = names(median_in_clust)[median_in_clust < 0.7 & max_val > 0.7]
    unstable = setdiff(unstable, cells.to.move)
    if (length(cells.to.move) > 0){
      print(length(cells.to.move))
      object@ident[cells.to.move] = max_clust[cells.to.move]
    }
    
    unstable_cells[unstable] = 1 
    
    
  }
  object@data.info$m = NULL
  object@data.info$m = object@ident
  #Recompute confusion matrix
  cell.ord = names(sort(object@ident))
  ConfusionMat = inClusterMat[cell.ord, cell.ord] / inSampleMat[cell.ord, cell.ord]
  
  annVal = c()
  for (i in sort(unique(object@ident))){
    annVal = c(annVal,rep(i, sum(object@ident==i)))
  }
  png("ConfusionMatrix_pruned_louvain.png",h=800,w=800)
  print(aheatmap(as.matrix(ConfusionMat),Rowv=NA, Colv=NA, labRow=NULL, labCol=NULL, main="Confusion Matrix (unfiltered)"), annCol=annVal, annColors = list(rainbow(length(unique(annVal)))))
  dev.off()
  
  cluster_instability = table(object@ident[names(unstable_cells[unstable_cells==1])]) / table(object@ident);
  unstable_clusters = which(cluster_instability > 0.1, useNames = TRUE)
  p=ggplot(data.frame(cluster_instability)) + geom_bar(aes(x=Var1, y=Freq), stat="identity", fill="blue") + theme_bw() + xlab("Cluster") + ylab("Instability Score") + gg.xax() + gg.yax() + ylim(c(0,1))
  pdf("Cluster_instability_louvain.pdf",w=10,h=6)
  print(p)
  dev.off()
  
  
  
}

#Cluster Stability matrix
clusts = sort(unique(object@ident))
MedianConfusion = matrix(0, nrow=length(clusts), ncol=length(clusts))
colnames(MedianConfusion) = clusts; rownames(MedianConfusion) = clusts

for (clust in clusts){
  cells = which.cells(object, clust)
  MedianConfusion[clust,clust] = median(ConfusionMat[cells,cells])
  #Within a cluster compute the median stability of every cell within the cluster
    for (clust2 in setdiff(clusts,clust)){
    cells2 = which.cells(object, clust2)
    MedianConfusion[clust, clust2] = median(ConfusionMat[cells,cells2])
  }
}

pdf("ConfusionClust_louvain.pdf",h=8,w=8)
print(aheatmap(as.matrix(MedianConfusion),Rowv=NA, Colv=NA, labRow=NULL, labCol=NULL, main="Median Cluster Stability"))
dev.off()

return(object)
  
}
)

setGeneric("doGraphCluster_stability", function(object,n.iter=100, downsamp.frac=0.85, cells.use=NULL,pcs.use=1:10, num.nn=100,adj.full=FALSE,set.ident=TRUE, method="Louvain",do.full.clust=FALSE,use.original=FALSE,do.jaccard=FALSE, do.noise=TRUE, do.prune=TRUE,is.symm=TRUE,...) standardGeneric("doGraphCluster_stability"))
setMethod("doGraphCluster_stability", "scR", function(object,n.iter=100,downsamp.frac = 0.85, cells.use=NULL,pcs.use=1:10,num.nn=100,adj.full=FALSE,set.ident=TRUE, method="Louvain",do.full.clust=FALSE, use.original=FALSE, do.jaccard=FALSE, do.noise=TRUE,do.prune=TRUE,is.symm=TRUE,...) {
  
  cells.use=set.ifnull(cells.use,colnames(object@data))
  clust.assign.mat = matrix(0, nrow=n.iter, ncol=length(cells.use))
  colnames(clust.assign.mat) = cells.use
  
  for ( i in 1:n.iter){
    print(paste0("Niter ", i))
    cells.use1 = sample(cells.use, round(downsamp.frac * length(cells.use)))
    clust.out = doGraph_clustering(object, pcs.use=pcs.use, num.nn=sample(c(round(0.5*num.nn):round(3*num.nn)), 1), cells.use=cells.use1, do.jaccard=do.jaccard, return.clust=TRUE, do.noise=do.noise, method=method,is.symm=is.symm)
    clust.assign =clust.out
    clust.assign.mat[i, names(clust.assign)] = clust.assign    
    
  }
  
  return(clust.assign.mat)
  
}
)

doSIMLR_embedding = function(object,cells.use=NULL, genes.use=NULL, do.cluster=FALSE, cluster.name = "cSIMLR",
                             num.clusters = 20,num.nn=30, num.pc = 20, num.graphs=15, sigma.min=0.5, sigma.max=3, num.sigma=7, regress.use=TRUE,
                             reduction.key="SIMLR", reduction.name="SIMLR", replace.tsne=FALSE, set.ident=TRUE, cluster.method="Louvain", ...){
  genes.use = set.ifnull(genes.use, object@var.genes)
  cells.use = set.ifnull(cells.use, colnames(object@data))
  
  if (length(genes.use) < 10){
    stop("Specify valid var.genes / insufficient features in object@var.genes")
  }
  
  print("Performing SIMLR embedding")
  print(paste0("Using ", length(genes.use), " variable genes"))
  if (regress.use){
    print("Using batch corrected data")
    data.use = object@regress.data[genes.use, cells.use]
  } else {
    print("Using normalized uncorrected data")
    data.use = object@data[genes.use, cells.use]
  }
  
 
  print(paste0("Seeking ", num.clusters, " clusters"))
  print(paste0("First reducing the dimensionality to ", num.pc, " PCs"))
  print(paste0("Using ", num.graphs, " k-nearest neighbor graphs with average k = ",num.nn))
  print(paste0("Using ", num.sigma, " different bandwidths between sigma.min = ", sigma.min, " and sigma.max = ", sigma.max))
  print(paste0("Total ", num.sigma * num.graphs, " kernels"))
  res_large_scale = SIMLR_Large_Scale(X=as.matrix(object@data[genes.use,]), c=num.clusters, k=num.nn, 
                                      num.pc = num.pc, num.k = num.graphs, sigma.min=sigma.min, sigma.max=sigma.max, num.sigma=num.sigma)
  simlr.rot = data.frame(SIMLR_1 = res_large_scale$ydata[,1], SIMLR_2 = res_large_scale$ydata[,2])
  rownames(simlr.rot) = cells.use
  
  simlr.obj <- new(
    Class = "dim.reduction",
    gene.loadings = matrix(rand(10),2,5),
    cell.embeddings = as.matrix(simlr.rot),
    dim.var = c(1,2,3),
    key = reduction.key
  )
  
  eval(expr = parse(text = paste0("object@dr$", reduction.name, "<- simlr.obj")))
  
  if (replace.tsne){
    object@tsne.rot = as.data.frame(simlr.rot)
  }
  
  if (do.cluster){
    clusters = doGraph_clustering_withAdj(object, Adj = res_large_scale$Adj, method=cluster.method, return.clust = TRUE, ...)
  }
  
  object@data.info[names(clusters), cluster.name] = clusters
  
  if (set.ident){
    object = set.all.ident(object, id=cluster.name)
  }
  
  return(object)
  
}

setGeneric("doGraph_clustering", function(object,cells.use=NULL,pcs.use=1:10, num.nn=100,adj.full=FALSE,set.ident=TRUE, method="Louvain", do.full.clust=FALSE, return.clust=FALSE,use.original=FALSE,do.jaccard=FALSE, do.noise=FALSE,is.symm=TRUE, return.graph=FALSE,return.adj=FALSE,use.reduction="pca",...) standardGeneric("doGraph_clustering"))
setMethod("doGraph_clustering", "scR", function(object,cells.use=NULL,pcs.use=1:10,num.nn=100,adj.full=FALSE,set.ident=TRUE, method="Louvain",do.full.clust=FALSE, return.clust=FALSE, use.original=FALSE, do.jaccard=FALSE, do.noise=FALSE,is.symm=TRUE, return.graph=FALSE,return.adj=FALSE,use.reduction="pca",...) {
  
  print(1)
  cells.use=set.ifnull(cells.use,colnames(object@data))
  if (!use.original){
    if (use.reduction == "ica"){
      data.use=object@ica.rot[cells.use,pcs.use]
    } else {
      data.use=object@pca.rot[cells.use,pcs.use]
    }
  } else {
    data.use=as.data.frame(t(object@scale.data[,cells.use]))
  }
  
 
#   if (!do.full.clust){
#     Adj = getAdjMatrix(data.use,nn=num.nn,edge.weights=FALSE,do.jaccard=do.jaccard,do.sparse=TRUE,full.eval=adj.full)
#   } else {
#     Adj = as.matrix(data.use) %*% t(as.matrix(data.use)) 
#   }
  
  if (method=="Louvain"){
       if (!do.full.clust){
         data.use = as.matrix(data.use)
         Adj = getAdjMatrix(data.use,nn=num.nn,edge.weights=FALSE,do.jaccard=do.jaccard,do.sparse=TRUE,full.eval=adj.full)
        } else {
         Adj = as.matrix(data.use) %*% t(as.matrix(data.use)) 
        }
    if (do.noise){
      n = sum(Adj>0)
      #Remove and add 5% of the edges
      to_add = sample(which(Adj==0), round(0.025*n))
      to_remove = sample(which(Adj>0), round(0.025*n))
      Adj[to_add] = Adj[to_remove]
      Adj[to_remove] = 0
      mult.noise = exp(matrix(runif(n^2, min = -0.2, max = 0.2), nrow=n, ncol=n))
      rownames(mult.noise) = rownames(Adj); colnames(mult.noise) = colnames(Adj);
      #Multiplicative noise
      Adj = Adj * mult.noise
    }
    if (do.jaccard){weights=TRUE} else {weights=NULL}
    if (return.adj) return(Adj)
    g <- graph.adjacency(Adj, weighted=weights, mode="undirected")
    if (return.graph) return(g)
    graph.out = cluster_louvain(g)
    #clust.out =   cluster_jl(Adj, s=1, self=1, debug=1, ddebug=0, verbose=TRUE)
    #clust.out=clust.out[[1]]
    #k = length(clust.out$COM)
    #clust.assign = clust.out$COM[[k]]
    #data.names=names(object@ident[cells.use])
    #names(clust.assign) = data.names;
    
    clust.assign = factor(graph.out$membership, levels=sort(unique(graph.out$membership)))
    names(clust.assign) = graph.out$names
    k=order(table(clust.assign), decreasing = TRUE)
    new.levels = rep(1,length(unique(graph.out$membership)))
    new.levels[k] = 1:length(unique(graph.out$membership))
    levels(clust.assign) = new.levels
    clust.assign = factor(clust.assign, levels=1:length(unique(graph.out$membership)))
    
    if (length(colnames(object@data)) == length(cells.use)){
      object@data.info$m = NULL
      object@data.info[names(clust.assign),"m"]=clust.assign
      if (set.ident) {
        object@ident=clust.assign; names(object@ident)=names(clust.assign);               
      }
    }
  }
      
  
  if (method=="Infomap"){
    if (!do.full.clust){
      Adj = getAdjMatrix(data.use,nn=num.nn,edge.weights=FALSE,do.jaccard=do.jaccard,do.sparse=TRUE,full.eval=adj.full)
    } else {
      Adj = as.matrix(data.use) %*% t(as.matrix(data.use)) 
    }
    if (do.noise){
      n = sum(Adj>0)
      #Remove and add 5% of the edges
      to_add = sample(which(Adj==0), round(0.025*n))
      
      to_remove = sample(which(Adj>0), round(0.025*n))
      Adj[to_add] = Adj[to_remove]
      Adj[to_remove] = 0
      
      mult.noise = exp(runif(n, min = -0.2, max = 0.2))
      #Multiplicative noise
      edge.ind = which(Adj>0)
      Adj[edge.ind] = Adj[edge.ind] * mult.noise
    }
    if (do.jaccard){weights=TRUE} else {weights=NULL}
    if (return.adj) return(Adj)
    g <- graph.adjacency(Adj, weighted=weights, mode="undirected")
    if (return.graph) return(g)
    #set.seed(42)
    graph.out = cluster_infomap(g)
    #clust.out =   cluster_jl(Adj, s=1, self=1, debug=1, ddebug=0, verbose=TRUE)
    #clust.out=clust.out[[1]]
    #k = length(clust.out$COM)
    #clust.assign = clust.out$COM[[k]]
    #data.names=names(object@ident[cells.use])
    #names(clust.assign) = data.names;
    
    clust.assign = factor(graph.out$membership, levels=sort(unique(graph.out$membership)))
    names(clust.assign) = graph.out$names
    k=order(table(clust.assign), decreasing = TRUE)
    new.levels = rep(1,length(unique(graph.out$membership)))
    new.levels[k] = 1:length(unique(graph.out$membership))
    levels(clust.assign) = new.levels
    clust.assign = factor(clust.assign, levels=1:length(unique(graph.out$membership)))
    
    k=order(table(clust.assign), decreasing = TRUE)
    new.levels = rep(1,length(unique(graph.out$membership)))
    new.levels[k] = 1:length(unique(graph.out$membership))
    levels(clust.assign) = new.levels
    clust.assign = factor(clust.assign, levels=1:length(unique(graph.out$membership)))
                          
     if (length(colnames(object@data)) == length(cells.use)){
        object@data.info$m = NULL
        object@data.info[names(clust.assign),"m"]=clust.assign
         if (set.ident) {
            object@ident=clust.assign; names(object@ident)=names(clust.assign);               
            }
       }
    
}
#   if (method=="Infomap"){
#     #Init clust
#     if (!exists("Network", mode="function")){
#       source("/Users/karthik/Documents/SOFTWARE/Infomap/examples/R/load-infomap.R")
#       library(igraph)
#       #source("/home/unix/karthik/bin/Infomap/examples/R/load-infomap.R")
#     }
#     
#     
#     conf <- init("--two-level --silent")
#     network <- Network(conf);
#     
#     network = getInfomapNetwork(network,data.use,nn=num.nn,edge.weights=FALSE,do.jaccard=do.jaccard,do.sparse=TRUE, is.symm=is.symm, do.noise=do.noise)
#     network$finalizeAndCheckNetwork(TRUE, length(cells.use))
#     
#     #Infomap clustering
#     
#     tree <- HierarchicalNetwork(conf)
#     run(network, tree);
#     
#     clusterIndexLevel <- 1 # 1, 2, ... or -1 for top, second, ... or lowest cluster level
#     leafIt <- tree$leafIter(clusterIndexLevel)
#     modules <- integer(length = network$numNodes())
#     
#     while (!leafIt$isEnd()) {
#       modules[leafIt$originalLeafIndex + 1] = leafIt$clusterIndex() + 1
#       leafIt$stepForward()
#     }
#     
#     # Create igraph community data
#     #comm <- create.communities(modules, algorithm = 'Infomap')
#     #print(comm)
#     
#     clust.assign = modules
#     names(clust.assign) = cells.use
#     
#     
#     data.names=names(object@ident[cells.use])
#     if (length(colnames(object@data)) == length(cells.use)){
#       object@data.info$m = NULL
#       object@data.info[data.names,"m"]=factor(clust.assign, levels=sort(unique(clust.assign)))
#       if (set.ident) {
#         object@ident=factor(clust.assign,levels=sort(unique(clust.assign))); names(object@ident)=data.names;               
#       }
#     }
#     
#   }
  
  if (return.clust){
    return(clust.assign) 
  } else {
    return(object) 
  }
}
)


setGeneric("doGraph_clustering_withAdj", function(object,Adj=NULL,normalize=TRUE,set.ident=TRUE, method="Louvain", return.clust=FALSE, do.noise=FALSE,is.symm=TRUE,clust.id="cSIMLR",...) standardGeneric("doGraph_clustering_withAdj"))
setMethod("doGraph_clustering_withAdj", "scR", function(object,Adj=NULL,normalize=TRUE,set.ident=TRUE, method="Louvain", return.clust=FALSE, do.noise=FALSE,is.symm=TRUE,clust.id="cSIMLR",...) {
  
  print(1)
  cells.use=rownames(Adj)
  if (is.null(Adj)) stop("Error: Adjacency matrix is needed")
  if (normalize) Adj = Adj/max(Adj)
  
  if (method=="Louvain"){
    if (do.noise){
      n = sum(Adj>0)
      #Remove and add 5% of the edges
      to_add = sample(which(Adj==0), round(0.025*n))
      to_remove = sample(which(Adj>0), round(0.025*n))
      Adj[to_add] = Adj[to_remove]
      Adj[to_remove] = 0
      
      mult.noise = exp(matrix(runif(n^2, min = -0.2, max = 0.2), nrow=n, ncol=n))
      rownames(mult.noise) = rownames(Adj); colnames(mult.noise) = colnames(Adj);
      #Multiplicative noise
      Adj = Adj * mult.noise
    }
    g <- graph.adjacency(Adj, weighted=TRUE, mode="undirected")
    graph.out = cluster_louvain(g)
    clust.assign = factor(graph.out$membership, levels=sort(unique(graph.out$membership)))
    names(clust.assign) = graph.out$names
    k=order(table(clust.assign), decreasing = TRUE)
    new.levels = rep(1,length(unique(graph.out$membership)))
    new.levels[k] = 1:length(unique(graph.out$membership))
    levels(clust.assign) = new.levels
    clust.assign = factor(clust.assign, levels=1:length(unique(graph.out$membership)))
    
    if (length(colnames(object@data)) == length(cells.use)){
      object@data.info$m = NULL
      object@data.info[names(clust.assign),clust.id]=clust.assign
      if (set.ident) {
        object@ident=clust.assign; names(object@ident)=names(clust.assign);               
      }
    }
  }
  
  
  if (method=="Infomap"){
    if (do.noise){
      n = sum(Adj>0)
      #Remove and add 5% of the edges
      to_add = sample(which(Adj==0), round(0.025*n))
      to_remove = sample(which(Adj>0), round(0.025*n))
      Adj[to_add] = Adj[to_remove]
      Adj[to_remove] = 0
      
      mult.noise = exp(matrix(runif(n^2, min = -0.2, max = 0.2), nrow=n, ncol=n))
      rownames(mult.noise) = rownames(Adj); colnames(mult.noise) = colnames(Adj);
      #Multiplicative noise
      Adj = Adj * mult.noise
    }
    
    g <- graph.adjacency(Adj, weighted=TRUE, mode="undirected")
    set.seed(42)
    graph.out = cluster_infomap(g)
    
    clust.assign = factor(graph.out$membership, levels=sort(unique(graph.out$membership)))
    names(clust.assign) = graph.out$names
    k=order(table(clust.assign), decreasing = TRUE)
    new.levels = rep(1,length(unique(graph.out$membership)))
    new.levels[k] = 1:length(unique(graph.out$membership))
    levels(clust.assign) = new.levels
    clust.assign = factor(clust.assign, levels=1:length(unique(graph.out$membership)))
    
    if (length(colnames(object@data)) == length(cells.use)){
      object@data.info$m = NULL
      object@data.info[names(clust.assign),clust.id]=clust.assign
      if (set.ident) {
        object@ident=clust.assign; names(object@ident)=names(clust.assign);               
      }
    }
    
    # clust.assign = factor(graph.out$membership, levels=sort(unique(graph.out$membership)))
    # names(clust.assign) = graph.out$names
    # k=order(table(clust.assign), decreasing = TRUE)
    # new.levels = rep(1,length(unique(graph.out$membership)))
    # new.levels[k] = 1:length(unique(graph.out$membership))
    # levels(clust.assign) = new.levels
    # clust.assign = factor(clust.assign, levels=1:length(unique(graph.out$membership)))
    # 
    # k=order(table(clust.assign), decreasing = TRUE)
    # new.levels = rep(1,length(unique(graph.out$membership)))
    # new.levels[k] = 1:length(unique(graph.out$membership))
    # levels(clust.assign) = new.levels
    # clust.assign = factor(clust.assign, levels=1:length(unique(graph.out$membership)))
    # 
    # if (length(colnames(object@data)) == length(cells.use)){
    #   object@data.info$m = NULL
    #   object@data.info[names(clust.assign),"m"]=clust.assign
    #   if (set.ident) {
    #     object@ident=clust.assign; names(object@ident)=names(clust.assign);               
    #   }
    # }
    
  }
  
  if (return.clust){
    return(clust.assign) 
  } else {
    return(object) 
  }
}
)


setGeneric("vizNetwork", function(object,cells.use=NULL,pcs.use=1:10, num.nn=10,ident.use="m", do.jaccard=TRUE, is.symm=TRUE,node.title=NULL, randomSeed.use=123,do.freeze=TRUE,...) standardGeneric("vizNetwork"))
setMethod("vizNetwork", "scR", function(object,cells.use=NULL,pcs.use=1:10,num.nn=10,ident.use="m", do.jaccard=TRUE,is.symm=TRUE,node.title=NULL,randomSeed.use=123,do.freeze=TRUE,...) {
  
  
  cells.use=set.ifnull(cells.use,colnames(object@data))
  cells.ident = object@data.info[cells.use, ident.use]
  data.use=object@pca.rot[cells.use,pcs.use]
  Adj = getAdjMatrix(data.use,nn=num.nn,edge.weights=FALSE,do.jaccard=do.jaccard,do.sparse=TRUE,full.eval=FALSE)
  if (do.jaccard){weights=TRUE} else {weights=NULL}
  g <- graph.adjacency(Adj, weighted=weights, mode="undirected")
  
  nodes = igraph::as_data_frame(g, what="vertices")
  colnames(nodes)[1] = "id"
  nodes$shadow=TRUE
  if (is.null(node.title)){
    node.title = cells.use
  }
  nodes$title = node.title
  nodes$clust = object@data.info[cells.use, ident.use]
  #nodes$label <- nodes$type.label # Node label
  #nodes$size <- nodes$audience.size # Node size
  #nodes$borderWidth <- 2 # Node border width
  require(visNetwork)
  require('RColorBrewer')
  cols.use = brewer.pal(length(unique(cells.ident)), "Spectral")
  ident.index = match(cells.ident, unique(cells.ident))
  nodes$color.background <- cols.use[ident.index]  
  nodes$color.border <- cols.use[ident.index]# "lightgrey"
  nodes$color.highlight.background <- "orange"
  nodes$color.highlight.border <- "darkred"

  
  edges = igraph::as_data_frame(g, what="edges")
  if (ncol(edges) > 2) colnames(edges)[3] = "width"
  edges$width <- edges$width # line width
  #print(head(edges))
  edges$value = 1
  edges$color <- "gray"    # line color  
  #links$arrows <- "middle" # arrows: 'from', 'to', or 'middle'
  edges$smooth <- TRUE    # should the edges be curved?
  #links$shadow <- FALSE    # edge shadow
  if (!do.freeze){
    net <- visNetwork(nodes, edges) %>% visLayout(randomSeed=randomSeed.use)
  } else {
    net <- visNetwork(nodes, edges) %>% visInteraction(dragNodes=FALSE, dragView=FALSE, zoomView=FALSE) %>% visLayout(randomSeed=randomSeed.use)
  }
  return(net)

}
)

setGeneric("vizNetwork.feature", function(object,net, feature.use=NULL,use.raw=FALSE, use.count=FALSE,scale.range=NULL,cols.use=terrain.colors,randomSeed=123,do.export=FALSE,...) standardGeneric("vizNetwork.feature"))
setMethod("vizNetwork.feature", "scR", function(object,net, feature.use=NULL,use.raw=FALSE, use.count=FALSE,scale.range=NULL,cols.use=terrain.colors,randomSeed=123,do.export=FALSE,...) {
  
  require(visNetwork)
  require('RColorBrewer')
  
  cells.use=rownames(net[1]$x$nodes)
  feature.data=data.frame(t(fetch.data(object,feature.use,cells.use = cells.use, use.raw=use.raw, use.count=use.count)))
  data.gene = as.numeric(feature.data); names(data.gene) = colnames(feature.data)
  
  if (is.null(scale.range)){  
    data.range = c(min(data.gene), max(data.gene)) 
  } else {
    data.range = c(scale.range[1], scale.range[2])
  }
  
  data.gene[data.gene>=scale.range[2]] = scale.range[2]
  data.gene[data.gene<=scale.range[1]] = scale.range[1]

  if(data.range[2] > 5){
    b=1
  }else{
    b=(data.range[2] - data.range[1]) /5
  }
  breaks = seq(data.range[1], data.range[2], by=b)
  data.cut=as.numeric(as.factor(cut(data.gene,breaks = breaks)))
  data.cut[is.na(data.cut)]=1
  cols.use=rev(cols.use(length(breaks)))[data.cut]
  net$x$nodes$color.background <- cols.use
  net$x$nodes$color.border <- "lightgrey"
  #nodes$color.highlight.background <- "orange"
  #nodes$color.highlight.border <- "darkred"
  
  
  #edges$width <- edges$width # line width
  #edges$value = 1
  net$x$edges$color <- "lightgrey"    # line color  
  net$x$edges$smooth <- TRUE    # should the edges be curved?
  #edges$shadow <- FALSE    # edge shadow
  net=visNetwork(net$x$nodes, net$x$edges, main=feature.use) %>% visLayout(randomSeed=randomSeed)
  
  return(net)
}
)

setGeneric("doBackSPIN_clustering", function(object,cells.use=NULL, genes.use=NULL, return.clust=FALSE,max_splits=5,SPINiter=10,doHeatmap=FALSE,use.count=FALSE, min.cells = 10, B.iter=10, ...) standardGeneric("doBackSPIN_clustering"))
setMethod("doBackSPIN_clustering", "scR", function(object,cells.use=NULL, genes.use=NULL,return.clust=FALSE,max_splits=5,SPINiter=10,doHeatmap=FALSE, use.count=FALSE, min.cells=10, B.iter=10,...) {
  
  cells.use=set.ifnull(cells.use,colnames(object@data))
  genes.use=set.ifnull(genes.use, object@var.genes)
  
  if (use.count){
    data.use=object@count.data[genes.use,cells.use]
  } else {
    data.use=object@data[genes.use,cells.use]
  }

  
  bs.result = BackSPIN(X=data.use, max_splits=max_splits,SPINiter=SPINiter, min.cells=10, B.iter=B.iter)
  if (doHeatmap){
    doHeatMap(object,cells.use=cells.use,genes.use=bs.result$gene.order,disp.min=-2,disp.max=2,draw.line=TRUE,order.by.ident=TRUE,col.use=pyCols,do.Rowv=FALSE)
  }
  
  if (!return.clust){
    object@ident=factor(bs.result$cell.clust)
    object@gene.ident = factor(bs.result$gene.clust[genes.use])
    return(object)
  } else {
    to.return = list()
    to.return$cell.clust = factor(bs.result$cell.clust)
    to.return$gene.clust = factor(bs.result$gene.clust[genes.use])
    return(to.return)
    
  }
  

}
)


setGeneric("pairwise_markers", function(object,genes.use=NULL, min.de.genes = 10, diff.thres=log(2),diff.test='bimod', pval=0.01) standardGeneric("pairwise_markers"))
setMethod("pairwise_markers", "scR", 
          function(object,genes.use=NULL, min.de.genes=10, diff.thres=log(2),diff.test='bimod', pval=0.01) {
            genes.use = set.ifnull(genes.use, rownames(object@data))
            zz = file("FINAL_CLUSTER_PAIRWISE_MARKERS.txt",open="wt")
            
            pval.thresh=pval/dim(object@data)[1]
            num.clust=length(levels(object@ident))   
            pass.thresh=1e6*data.frame(diag(num.clust)); 
            
            for(i in 1:num.clust) {
              print(paste0("Testing Cluster ", i))
              for(j in ((i+1):num.clust)) {
                if (j>num.clust) break
                if (pass.thresh[i,j]==0) {
                  marker=find.markers(object,i,j,genes.use=genes.use,thresh.use=diff.thres,test.use= diff.test)
                  marker.pass=subset(marker,myP<pval.thresh)
                  sink(zz)
                  print(paste("Clusters ",i,j, "-# DE = ",nrow(marker.pass)))
                  print(marker.pass)
                  sink()
                  
                  
                  pass.thresh[i,j]=nrow(marker.pass); pass.thresh[j,i]=pass.thresh[i,j];
                  
                }
                
              }
              
            }
            
            return(1)
          }
)

setGeneric("merge.clusters.by.de", function(object,genes.use=NULL, min.de.genes = 10, diff.thres=log(2),diff.test='bimod', pval=0.01, tag=NULL, pcs.use=1:10, test.max=NULL,  TPM.mat=NULL, Count.mat=NULL, min.pos.perc=0.1) standardGeneric("merge.clusters.by.de"))
setMethod("merge.clusters.by.de", "scR", 
          function(object,genes.use=NULL, min.de.genes=10, diff.thres=log(2),diff.test='bimod', pval=0.01, tag=NULL, pcs.use=1:10, test.max=NULL,  TPM.mat=NULL, Count.mat=NULL, min.pos.perc=0.1) {
            genes.use = set.ifnull(genes.use, rownames(object@data))
            if (is.null(tag)){
              filename = "CLUSTER_PAIRWISE_MARKERS.txt"
            } else {
              filename = paste0("CLUSTER_PAIRWISE_MARKERS_",tag,".txt")
            }
            zz = file(filename,open="wt")
          

            num.clust=length(levels(object@ident)) 
            print(levels(object@ident))
            print(paste0("Starting with ", num.clust, " clusters"))
            pass.thresh=1e6*data.frame(diag(num.clust)); 
            dist.clust = pass.thresh
            
            for(i in 1:num.clust) {
              print(paste0("Testing Cluster ", i))
              for(j in ((i+1):num.clust)) {
                print(j)
                if (j>num.clust) break
                if (pass.thresh[i,j]==0) {
                  marker=find.markers(object,levels(object@ident)[i],levels(object@ident)[j],genes.use=genes.use,thresh.use=diff.thres,test.use= diff.test, test.max=test.max,  TPM.mat=TPM.mat, Count.mat=Count.mat, min.pos.perc=min.pos.perc)
                  myP2 = p.adjust(marker$myP, method="fdr")
                  marker$myP = myP2
                  marker.pass=subset(marker,myP<pval)
                  sink(zz)
                  print(paste("Clusters ",levels(object@ident)[i],levels(object@ident)[j], "-# DE = ",nrow(marker.pass)))
                  print(head(subset(marker.pass, myDiff > 0),5))
                  print(head(subset(marker.pass, myDiff < 0),5))
                  sink()
                  
                  num.de.genes = min(nrow(subset(marker.pass, myDiff > 0)), nrow(subset(marker.pass, myDiff < 0)))
                  pass.thresh[i,j]=num.de.genes; pass.thresh[j,i]=pass.thresh[i,j];
                  
                }
                
              }
              
            }
            
            colnames(pass.thresh) = levels(object@ident)
            rownames(pass.thresh) = levels(object@ident)
            
            print(pass.thresh)
            
            #iteratively merge clusters
            if (is.null(tag)){
              filename = "CLUSTERS_merge_log.txt"
            } else {
              filename = paste0("CLUSTERS_merge_log_",tag,".txt")
            }
            zz1 = file(filename,open="wt")
            min.val = min(pass.thresh)
            min.val.ind = as.data.frame(which(pass.thresh==min.val, arr.ind = TRUE))
            min.val.ind = min.val.ind[min.val.ind$row < min.val.ind$col,]
            min.val.ind$val = min(pass.thresh)
            rownames(min.val.ind) = 1:nrow(min.val.ind)
            
            
            merge.ind=-1
            while(min.val < min.de.genes) {
              merge.ind=merge.ind+1
              dir.create("Merge_history")
              pdf(paste0("Merge_history/Merge_",merge.ind,".pdf"), w=10,h=10)
              tsplot(object)
              dev.off()
              
              #clust.dists = ComputeClusterDistances(object, reduction.use="pca", pcs.use=pcs.use, dist.type="centroid") 
              clust.dists = ComputeClusterDistances(object, reduction.use="pca", dist.type="centroid", pcs.use=pcs.use) 
              ind.min=which.min(clust.dists[cbind(min.val.ind$row, min.val.ind$col)])
              test.1 = min.val.ind[ind.min,]$row; test.2 = min.val.ind[ind.min,]$col
              
              if (pass.thresh[test.1,test.2]< min.de.genes) {
                object@ident[which(object@ident==test.2)]=test.1
                sink(zz1)
                print(paste("merging",test.1,test.2," with # de.genes=", pass.thresh[test.1, test.2]))
                sink()
#                pass.thresh[test.2,]=rep(1e6,ncol(pass.thresh));
#                pass.thresh[,test.2]=rep(1e6,nrow(pass.thresh))
                pass.thresh = pass.thresh[-test.2,]; pass.thresh = pass.thresh[,-test.2]
                object@ident = drop.levels(object@ident)
                levels(object@ident) = c(1:length(levels(object@ident)))
                object@data.info[,"m"] = object@ident
                num.clust = max(as.numeric(levels(object@ident)))

                #Recompute pairwise markers for merged cluster
                print(paste0("Recomputing pairwise markers for new clust ", test.1))
                for (i in setdiff(c(1:num.clust), test.1)){
                  print(i)
                  marker=find.markers(object,test.1,i,genes.use=genes.use,thresh.use=diff.thres,test.use= diff.test, test.max=test.max, TPM.mat=TPM.mat, Count.mat=Count.mat, min.pos.perc=min.pos.perc)
                  myP2 = p.adjust(marker$myP, method="fdr")
                  marker$myP = myP2
                  marker.pass=subset(marker,myP<pval)
                  pass.thresh[test.1,i]=min(nrow(subset(marker.pass, myDiff>0)),nrow(subset(marker.pass, myDiff<0))); pass.thresh[i,test.1]=pass.thresh[test.1,i];
                  
                }
                
              }
              colnames(pass.thresh) = 1:num.clust
              rownames(pass.thresh) = 1:num.clust

              min.val = min(pass.thresh)
              min.val.ind = as.data.frame(which(pass.thresh==min.val, arr.ind = TRUE))
              min.val.ind = min.val.ind[min.val.ind$row < min.val.ind$col,]
              min.val.ind$val = min(pass.thresh)
              rownames(min.val.ind) = 1:nrow(min.val.ind)

            }
            return(object)
          }
)

setGeneric("merge.clusters.by.de2", function(object,genes.use=NULL, min.de.genes = 10, diff.thres=log(2),diff.test='bimod', pval.cutoff=0.01, tag=NULL, pcs.use=1:10, test.max=NULL, clust.test=NULL, TPM.mat=NULL, Count.mat=NULL, min.pos.perc=0.1,max.neg.perc=1,min.diff=0.1, test.2=FALSE, test.2.use="binom", n.cores=1) standardGeneric("merge.clusters.by.de2"))
setMethod("merge.clusters.by.de2", "scR", 
          function(object,genes.use=NULL, min.de.genes=10, diff.thres=log(2),diff.test='bimod', pval.cutoff=0.01, tag=NULL, pcs.use=1:10, test.max=NULL, clust.test=NULL, TPM.mat=NULL, Count.mat=NULL, min.pos.perc=0.1, max.neg.perc=1,min.diff=0.1, test.2=FALSE, test.2.use="binom", n.cores=1) {
            print(paste0("# PCs = ", length(pcs.use)))
            genes.use = set.ifnull(genes.use, rownames(object@data))
            clust.test = set.ifnull(clust.test, as.numeric(levels(object@ident)))
            
            
            num.clust=length(clust.test) 
            print(clust.test)
            print(paste0("Starting with ", num.clust, " clusters"))
            pass.thresh=1e6*data.frame(diag(length(levels(object@ident)))); 
            
            for (i in setdiff(as.numeric(levels(object@ident)), clust.test)){
              pass.thresh[i,]=1e6; pass.thresh[,i]=1e6;
            } 
            
            dist.clust = pass.thresh
            
            if (n.cores == 1){ 
              if (is.null(tag)){
                filename = "CLUSTER_PAIRWISE_MARKERS.txt"
              } else {
                filename = paste0("CLUSTER_PAIRWISE_MARKERS_",tag,".txt")
              }
              zz = file(filename,open="wt")
                  for(k in 1:num.clust) {
                    i=clust.test[k]
                    print(paste0("Testing Cluster ", i, "against clusters ", i+1,"-",num.clust))
                    for(m in ((k+1):num.clust)) {
                      j=clust.test[m]
                      if (m>num.clust) break
                      if (pass.thresh[i,j]==0) {
                        #print(i)
                        marker=find.markers(object,i,j,genes.use=genes.use,thresh.use=diff.thres,test.use= diff.test, test.max=test.max, TPM.mat=TPM.mat, Count.mat=Count.mat, min.pos.perc=min.pos.perc, min.diff = min.diff, max.neg.perc=max.neg.perc, pval.thresh = 1e-2, verbose=FALSE)
                        pval2 = p.adjust(marker$pval, method="fdr")
                        marker$pval = pval2
                        marker.pass=subset(marker,pval<pval.cutoff)
                        
                        sink(zz)
                        print(paste("For cluster ",i, "vs. cluster", j, " found ",nrow(marker.pass), " DE markers"))
                        #print(head(subset(marker.pass, diff > 0),5))
                        #print(head(subset(marker.pass, diff < 0),5))
                        sink()
                        
                        num.de.genes = min(nrow(subset(marker.pass, diff > 0)), nrow(subset(marker.pass, diff < 0)))
                        print(paste("For cluster ",i, "vs. cluster", j, " found ",nrow(subset(marker.pass, diff > 0)), " +ve DE markers and",
                                    nrow(subset(marker.pass, diff < 0)), " -ve DE markers"))
                        pass.thresh[i,j]=num.de.genes; pass.thresh[j,i]=pass.thresh[i,j];
                        
                      }
                      
                    }
                    
                    
                  }
            } else {
              # not complete
              require(parallel)
              require(foreach)
              require(doParallel)
              for(k in 1:num.clust) {
                clust1=clust.test[k]
                print(paste0("Testing Cluster ", clust1))
                cl = makePSOCKcluster(n.cores, outfile="")
                registerDoParallel(cl, cores=n.cores)
                chunksize = (num.clust - k) / n.cores
                vals = split((k+1):num.clust, ceiling(seq_along(1:((num.clust-k)/chunksize))))
                export.funs = c("find.markers")
                par.out = foreach(run.id = 1:n.cores, .combine=c) %dopar% {
                v = vals[[run.id]]
                do.call(c, lapply(v, function(i) {
                     marker=find.markers(object,i,j,genes.use=genes.use,thresh.use=diff.thres,test.use= diff.test, test.max=test.max, TPM.mat=TPM.mat, Count.mat=Count.mat, min.pos.perc=min.pos.perc, min.diff = min.diff, max.neg.perc=max.neg.perc, pval.thresh = 1e-2, verbose=FALSE)
                    pval2 = p.adjust(marker$pval, method="fdr")
                    marker$pval = pval2
                    marker.pass=subset(marker,pval<pval.cutoff)
                    min(nrow(subset(marker.pass, diff > 0)), nrow(subset(marker.pass, diff < 0)))
                  }))
                }
                
                cat("\nUnregistering parallel backend..")
                stopCluster(cl)
                registerDoSEQ()
                cat(" done\n");
              
              }
              
              
            }
            colnames(pass.thresh) = levels(object@ident)
            rownames(pass.thresh) = levels(object@ident)
            
            write.table(pass.thresh, file=paste0("DE_genes_matrix_",tag, ".txt"), sep="\t", quote=FALSE)
            
            print(pass.thresh)
            
            #iteratively merge clusters
            if (is.null(tag)){
              filename = "CLUSTERS_merge_log.txt"
            } else {
              filename = paste0("CLUSTERS_merge_log_",tag,".txt")
            }
            zz1 = file(filename,open="wt")
            min.val = min(pass.thresh)
            min.val.ind = as.data.frame(which(pass.thresh==min.val, arr.ind = TRUE))
            min.val.ind = min.val.ind[min.val.ind$row < min.val.ind$col,]
            min.val.ind$val = min(pass.thresh)
            rownames(min.val.ind) = 1:nrow(min.val.ind)
            
            
            merge.ind=-1
            while(min.val < min.de.genes) {
              merge.ind=merge.ind+1
              dir.create("Merge_history")
              pdf(paste0("Merge_history/Merge_",merge.ind,".pdf"), w=10,h=10)
              tsplot(object)
              dev.off()
              
              #clust.dists = ComputeClusterDistances(object, reduction.use="pca", pcs.use=pcs.use, dist.type="centroid") 
              clust.dists = ComputeClusterDistances(object, reduction.use="pca", dist.type="centroid", pcs.use=pcs.use) 
              ind.min=which.min(clust.dists[cbind(min.val.ind$row, min.val.ind$col)])
              test.1 = min.val.ind[ind.min,]$row; test.2 = min.val.ind[ind.min,]$col
              
              if (pass.thresh[test.1,test.2]< min.de.genes) {
                object@ident[which(object@ident==test.2)]=test.1
                print(paste("merging",test.1,test.2," with # de.genes=", pass.thresh[test.1, test.2]))
                sink(zz1)
                print(paste("merging",test.1,test.2," with # de.genes=", pass.thresh[test.1, test.2]))
                sink()
                #                pass.thresh[test.2,]=rep(1e6,ncol(pass.thresh));
                #                pass.thresh[,test.2]=rep(1e6,nrow(pass.thresh))
                pass.thresh = pass.thresh[-test.2,]; pass.thresh = pass.thresh[,-test.2]
                old.ident.levels = as.numeric(levels(object@ident))
                old.ident.levels = setdiff(old.ident.levels, test.2)
                clust.test = setdiff(clust.test, test.2)
                
                object@ident = drop.levels(object@ident)
                levels(object@ident) = c(1:length(levels(object@ident)))
                object@data.info[,"m"] = object@ident
                
                new.ident.levels = as.numeric(levels(object@ident))
                names(new.ident.levels) = as.character(old.ident.levels)
                clust.test = new.ident.levels[as.character(clust.test)]
                
                
                
                #Recompute pairwise markers for merged cluster
                print(paste0("Recomputing pairwise markers for new clust ", test.1))
                for (i in setdiff(clust.test, test.1)){
                  print(i)
                  marker=find.markers(object,test.1,i,genes.use=genes.use,thresh.use=diff.thres,test.use= diff.test, test.max=test.max, TPM.mat=TPM.mat, Count.mat=Count.mat, min.pos.perc=min.pos.perc, verbose = FALSE)
                  pval2 = p.adjust(marker$pval, method="fdr")
                  marker$pval = pval2
                  marker.pass=subset(marker,pval<pval.cutoff)
                  pass.thresh[test.1,i]=min(nrow(subset(marker.pass, diff>0)),nrow(subset(marker.pass, diff<0))); pass.thresh[i,test.1]=pass.thresh[test.1,i];
                  
                }
                
              }
              colnames(pass.thresh) = 1:length(levels(object@ident))
              rownames(pass.thresh) = colnames(pass.thresh)
              
              min.val = min(pass.thresh)
              min.val.ind = as.data.frame(which(pass.thresh==min.val, arr.ind = TRUE))
              min.val.ind = min.val.ind[min.val.ind$row < min.val.ind$col,]
              min.val.ind$val = min(pass.thresh)
              rownames(min.val.ind) = 1:nrow(min.val.ind)
              
            }
            return(object)
          }
)

setGeneric("ClusterCentroids", function(object,reduction.use="pca", pcs.use = 1:10, cells.use=NULL, genes.use=NULL) standardGeneric("ClusterCentroids"))
setMethod("ClusterCentroids", "scR", 
          function(object,reduction.use="pca", pcs.use = 1:10, cells.use=NULL, genes.use=NULL) {
            cells.use = set.ifnull(cells.use, colnames(object@data))
            ident.use=object@ident[cells.use]
            if (reduction.use == "tsne"){data.use = object@tsne.rot[cells.use,]}
            if (reduction.use == "pca"){
              pcs.use = set.ifnull(pcs.use, dim(object@pca.rot)[2])
              data.use = object@pca.rot[cells.use,pcs.use]
            }
            if (reduction.use == "orig"){
              genes.use = set.ifnull(genes.use, rownames(object@data))
              data.use = as.data.frame(t(object@data[genes.use,cells.use]))
            }
            
            centroids = c()
            for (i in levels(ident.use)){
                cells.in.cluster = which.cells(object, i)
                cells.in.cluster = cells.in.cluster[cells.in.cluster %in% cells.use]
                centroids = rbind(centroids, colMeans(data.use[cells.in.cluster,]))
            }
            centroids = as.data.frame(centroids)
            colnames(centroids) = colnames(data.use)
            rownames(centroids) = as.numeric(levels(object@ident))
          
            return(centroids)
          }
            

)

setGeneric("ComputeClusterDistances", function(object, reduction.use="pca", dist.type="nn", pcs.use = 1:10, cells.use=NULL, genes.use=NULL) standardGeneric("ComputeClusterDistances"))
setMethod("ComputeClusterDistances", "scR", 
          function(object,reduction.use="pca",dist.type="nn", pcs.use = 1:10, cells.use=NULL, genes.use=NULL) {
            require(RANN)
            cells.use = set.ifnull(cells.use, colnames(object@data))
            ident.use=object@ident[cells.use]
            if (reduction.use == "tsne"){
              data.use = object@tsne.rot[cells.use,]
              centroids = ClusterCentroids(object, reduction.use="tsne", cells.use=cells.use)
            
            }
            if (reduction.use == "pca"){
              pcs.use = set.ifnull(pcs.use, dim(object@pca.rot)[2])
              data.use = object@pca.rot[cells.use,pcs.use]
              centroids = ClusterCentroids(object, reduction.use="pca", pcs.use=pcs.use, cells.use=cells.use)
            }
            if (reduction.use == "orig"){
              genes.use = set.ifnull(genes.use, rownames(object@data))
              data.use = as.data.frame(t(object@data[genes.use,cells.use]))
              centroids = ClusterCentroids(object, reduction.use="orig", genes.use=genes.use, cells.use=cells.use)
            }
            
            if (dist.type=="centroid"){
              clust.dists = as.matrix(dist(centroids, upper=TRUE))
              diag(clust.dists) = 1e6
            }
            
            num.clust = length(levels(ident.use))
            

            if (dist.type == "nn"){
              clust.dists = matrix(0, nrow=num.clust, ncol=num.clust)
              diag(clust.dists) = 1e6
              rownames(clust.dists) = levels(ident.use)
              colnames(clust.dists) = rownames(clust.dists)
                for (i in 1:nrow(clust.dists)){
                   for(j in ((i+1):ncol(clust.dists))){
                      if (j>nrow(clust.dists)) break
                      cells.in.cluster_i = which.cells(object, i)
                      cells.in.cluster_i = cells.in.cluster_i[cells.in.cluster_i %in% cells.use]
                      cells.in.cluster_j = which.cells(object, j)
                      cells.in.cluster_j = cells.in.cluster_j[cells.in.cluster_j %in% cells.use]
                      
                      nnA = nn2(data.use[cells.in.cluster_i,], query = centroids[j,], k=1)
                      nnB = nn2(data.use[cells.in.cluster_j,], query = centroids[i,],k=1)
                      clust.dists[i,j] = min(c(nnA$nn.dists, nnB$nn.dists))
                      
                      clust.dists[j,i] = clust.dists[i,j]
                   }
                }
            }
            
            colnames(clust.dists) = c(1:ncol(clust.dists))
            rownames(clust.dists) = colnames(clust.dists)
            return(clust.dists)
          }
          
          
)


setGeneric("cluster.alpha", function(object,thresh.min=0) standardGeneric("cluster.alpha"))
setMethod("cluster.alpha", "scR", 
          function(object,thresh.min=0) {
            ident.use=object@ident
            data.all=data.frame(row.names = rownames(object@data))
            for(i in sort(unique(ident.use))) {
              temp.cells=which.cells(object,i)
              data.temp=apply(object@data[,temp.cells],1,function(x)return(length(x[x>thresh.min])/length(x)))
              data.all=cbind(data.all,data.temp)
              colnames(data.all)[ncol(data.all)]=i
            }
            colnames(data.all)=sort(unique(ident.use))
            return(data.all)
          }
)

setGeneric("average.pca", function(object) standardGeneric("average.pca"))
setMethod("average.pca", "scR", 
          function(object) {
            ident.use=object@ident
            data.all=data.frame(row.names = colnames(object@pca.rot))
            for(i in levels(ident.use)) {
              temp.cells=which.cells(object,i)
              if (length(temp.cells)==1) data.temp=(object@pca.rot[temp.cells,])
              if (length(temp.cells)>1) data.temp=apply(object@pca.rot[temp.cells,],2,mean)
              data.all=cbind(data.all,data.temp)
              colnames(data.all)[ncol(data.all)]=i
            }
            colnames(data.all)=levels(ident.use)
            return((data.all))
          }
)

setGeneric("average.expression", function(object,genes.use=NULL, fxn=expMean, regress.use=FALSE) standardGeneric("average.expression"))
setMethod("average.expression", "scR", 
          function(object,genes.use=NULL, fxn=expMean, regress.use=FALSE) {
            genes.use=set.ifnull(genes.use,rownames(object@data))
            genes.use=ainb(genes.use,rownames(object@data))
            ident.use=object@ident
            data.all=data.frame(row.names = genes.use)
            if (regress.use){
              data.use=object@regress.data
            } else {
              data.use = object@data
            }
            data.use = data.use[genes.use,]
            for(i in levels(ident.use)) {
              temp.cells=which.cells(object,i)
              if (length(temp.cells)==1) data.temp=(data.use[,temp.cells])
              if (length(temp.cells)>1) data.temp=apply(data.use[,temp.cells],1,fxn)
              data.all=cbind(data.all,data.temp)
              colnames(data.all)[ncol(data.all)]=i
            }
            colnames(data.all)=levels(ident.use)
            return(data.all)
          }
)


setGeneric("print.pca", function(object,pcs.print=1:5,genes.print=30,use.full=FALSE) standardGeneric("print.pca"))
setMethod("print.pca", "scR", 
          function(object,pcs.print=1:5,genes.print=30,use.full=FALSE) {
            pc_scores=object@pca.x
            if (use.full==TRUE) pc_scores = object@pca.x.full
            for(i in pcs.print) {
              code=paste("PC",i,sep="")
              sx=pc_scores[order(pc_scores[,code]),]
              print(code)
              print(rownames(sx[1:genes.print,]))
              print ("")
              
              print(rev(rownames(sx[(nrow(sx)-genes.print):nrow(sx),])))
              print ("")
              print ("")      
            }
          }
)

setGeneric("fetch.data",  function(object, vars.all=NULL,cells.use=NULL,use.imputed=FALSE, use.regress=FALSE, use.raw=FALSE, use.count=FALSE) standardGeneric("fetch.data"))
setMethod("fetch.data","scR",
          function(object, vars.all=NULL,cells.use=NULL,use.imputed=FALSE, use.regress=FALSE,use.raw=FALSE, use.count=FALSE) {
            cells.use=set.ifnull(cells.use,object@cell.names)
            data.return=data.frame(row.names = cells.use)
            data.expression=object@data; if (use.imputed) data.expression=object@imputed[,cells.use]; 
            if (use.regress) data.expression=object@regress.data[,cells.use]
            if (use.raw) data.expression=object@raw.data[,cells.use]
            if (use.count) data.expression=object@count.data[,cells.use]
            var.options0=c("ident","data.ngene","cell.names")
            var.options=c("data.info","pca.rot","ica.rot","tsne.rot",
"dr$cca@cell.embeddings","dr$cca.aligned@cell.embeddings")
            data.expression=t(data.expression)
            vars.in.data = c()
            for (my.var in vars.all) {
              print(my.var)
              data.use=data.frame()
              if (my.var %in% colnames(data.expression)) {
                data.use=data.expression
              } else if (my.var %in% var.options0) {
                eval(parse(text=paste("data.use = object@",my.var,sep=""))) 
                data.use = data.frame(data.use)
                colnames(data.use) = my.var
              } else {
                for(i in var.options) {
                  eval(parse(text=paste("data.use = object@",i,sep="")))
                  if (my.var %in% colnames(data.use)) {
                    break;
                  }
                }
              }
              if (nrow(data.use)==0 | ncol(data.use)==0) {
                print(paste("Error : ", my.var, " not found", sep=""))
                next
                #return(0);
              }
              cells.use=ainb(cells.use,rownames(data.use))
              data.add=data.use[cells.use,my.var]
              data.return=cbind(data.return,data.add) 
              vars.in.data = c(vars.in.data, my.var)
            }
            colnames(data.return)=vars.in.data #unlist(strsplit(vars.in.data,"[|]"))[seq(1,length(vars.in.data),2)]
            rownames(data.return)=cells.use
            return(data.return)
          }
)



setGeneric("viz.pca", function(object,pcs.use=1:5,num.genes=15,use.full=FALSE,font.size=0.5,nCol=NULL, gene.labels="left") standardGeneric("viz.pca"))
setMethod("viz.pca", "scR", 
          function(object,pcs.use=1:5,num.genes=15,use.full=FALSE,font.size=0.5,nCol=NULL,gene.labels="left") {
            pc_scores=object@pca.x
            if (use.full==TRUE) pc_scores = object@pca.x.full
            
            if (is.null(nCol)) {
              nCol=2
              if (length(pcs.use)>6) nCol=3
              if (length(pcs.use)>9) nCol=4
            }         
            num.row=floor(length(pcs.use)/nCol-1e-5)+1
            par(mfrow=c(num.row,nCol))
            
            for(i in pcs.use) {
              code=paste("PC",i,sep="")
              sx=pc_scores[order(pc_scores[,code]),]
              subset.use=sx[c(1:num.genes,(nrow(sx)-num.genes):nrow(sx)),]
              plot(subset.use[,i],1:nrow(subset.use),pch=16,col="blue",xlab=paste("PC",i,sep=""),yaxt="n",ylab="")
              if (gene.labels=="right"){
                  axis(4,at=1:nrow(subset.use),labels = rownames(subset.use),las=1,cex.axis=font.size) 
              } else {
                  axis(2,at=1:nrow(subset.use),labels = rownames(subset.use),las=1,cex.axis=font.size)
              }
            }
            rp()
          }
)

setGeneric("viz.pc.heatmap", function(object,pc.use=1,cells.use=NULL,num.genes=15,font.size=0.5, min.disp=-2.5, max.disp=2.5, use.count=FALSE, cex.row.use=1, cex.col.use=1, suppress.cellnames=TRUE, remove.genes=NULL, do.average=FALSE, thresh.use=0, do.scale=FALSE, ann.id=NULL, do.col.clust=NA) standardGeneric("viz.pc.heatmap"))
setMethod("viz.pc.heatmap", "scR", 
          function(object,pc.use=1,cells.use=NULL,num.genes=15,font.size=0.5,min.disp=-2.5, max.disp=2.5, use.count=FALSE, cex.row.use=1,cex.col.use=1,  suppress.cellnames=TRUE, remove.genes=NULL, do.average=FALSE, thresh.use=0, do.scale=FALSE, ann.id=NULL, do.col.clust=NA) {
            
            require(heatmap3)
            if (length(pc.use) != 1) stop("pc.use must be a single number")
            code=paste("PC",pc.use,sep="")
            cells.use=set.ifnull(cells.use, colnames(object@data))
            
            pc_loads = as.numeric(object@pca.x[,pc.use]); names(pc_loads) = rownames(object@pca.x)
            #Order genes from low to high
            pc_loads=pc_loads[order(pc_loads)]
            genes.use=names(pc_loads)[c(1:num.genes,(length(pc_loads)-num.genes):length(pc_loads))]
            genes.use=setdiff(genes.use, remove.genes)
            pc_scores = as.numeric(object@pca.rot[cells.use,pc.use]); names(pc_scores) = cells.use
            pc_scores=pc_scores[order(pc_scores)]
            cells.use.ord = names(pc_scores)
            
            if (!do.average){
                if (use.count){
                    data.plot = object@count.data[genes.use, cells.use.ord]
                } else {
                  data.plot = object@data[genes.use, cells.use.ord]
                }
            } else {
               data.plot = as.data.frame(t(fetch.data(object, vars.all = genes.use, use.scaled = FALSE, use.count = use.count, use.raw = FALSE)))
            }
            if (do.scale) data.plot = t(scale(t(data.plot)))
            data.plot=minmax(data.plot[rev(genes.use),],min=min.disp,max=max.disp)
            
            if (do.average){
              data.plot.temp = data.plot
              cells.use1 = names(pc_scores)[pc_scores <= thresh.use]
              cells.use2 = names(pc_scores)[pc_scores > thresh.use]
              if (use.count){ 
                stat.fxn = mean
              } else {
                stat.fxn = mean
              }
              data.plot = cbind(apply(data.plot.temp[,cells.use1],1,stat.fxn), apply(data.plot.temp[,cells.use2],1,stat.fxn))
              colnames(data.plot) = c(paste0("PC", pc.use," <= ", thresh.use), paste0("PC", pc.use," >", thresh.use))
              lab.col.use=colnames(data.plot)
            } else {
              if (suppress.cellnames){
                lab.col.use = rep("", ncol(data.plot))
              } else {
                lab.col.use=NULL
              }
            }
            print(paste0("Using min.disp=",min.disp," and max.disp=", max.disp))
            if (!do.average){
              if (is.null(ann.id)){
                heatmap3(as.matrix(data.plot),Colv=NA, Rowv=NA, scale="none", cexCol=cex.col.use, balanceColor = T, cexRow=cex.row.use, labCol=lab.col.use)
              } else {
                object = set.all.ident(object,id=ann.id)
                ColSideAnn=data.frame(object@ident[colnames(data.plot)])
                colnames(ColSideAnn) = "Group"
                heatmap3(as.matrix(data.plot),Colv=do.col.clust, Rowv=NA, scale="none", cexCol=cex.col.use, balanceColor = T, cexRow=cex.row.use, labCol=lab.col.use,
                         ColSideAnn=ColSideAnn,ColSideFun=function(x) showAnn(x),ColSideWidth=0.8)
              }
              } else {
              require(colorRamps)
              if (do.scale){ 
                scale.type = "row"
              } else {
                  scale.type = "none"
                }
              print(aheatmap(data.plot, Rowv = NA, Colv=NA, color=blue2red(20), scale=scale.type))
            }
          }
)

setGeneric("hclust.heatmap", function(object,genes.use=NULL,cells.use=NULL,font.size=0.5, min.disp=-2.5, max.disp=2.5, use.count=FALSE, cex.row.use=1, cex.col.use=1, suppress.cellnames=TRUE,  thresh.use=0, do.scale=FALSE, ann.id=NULL, do.col.clust=TRUE, do.row.clust=TRUE) standardGeneric("hclust.heatmap"))
setMethod("hclust.heatmap", "scR", 
          function(object,genes.use=NULL,cells.use=NULL,font.size=0.5,min.disp=-2.5, max.disp=2.5, use.count=FALSE, cex.row.use=1,cex.col.use=1,  suppress.cellnames=TRUE,  thresh.use=0, do.scale=FALSE, ann.id=NULL, do.col.clust=TRUE, do.row.clust=TRUE) {
            
            require(heatmap3)
            
            cells.use=set.ifnull(cells.use, colnames(object@data))
            genes.use = set.ifnull(genes.use, rownames(object@data))
            genes.use = genes.use[genes.use %in% rownames(object@data)]
            
              if (use.count){
                data.plot = object@count.data[genes.use, cells.use]
              } else {
                data.plot = object@data[genes.use, cells.use]
              }
            
            if (do.scale) data.plot = t(scale(t(data.plot)))
            data.plot=minmax(data.plot[rev(genes.use),],min=min.disp,max=max.disp)
            
            if (suppress.cellnames){
              lab.col.use = rep("", ncol(data.plot))
            } else {
              lab.col.use=NULL
            }
            
            print(paste0("Using min.disp=",min.disp," and max.disp=", max.disp))
          
            if (is.null(ann.id)){
              heatmap3(data.plot,Colv=NA, Rowv=NA, scale="none", cexCol=cex.col.use, balanceColor = T, cexRow=cex.row.use, labCol=lab.col.use, method="average")
            } else {
              object = set.all.ident(object,id=ann.id)
              ColSideAnn=data.frame(object@ident[colnames(data.plot)])
              colnames(ColSideAnn) = "Group"
              heatmap3(data.plot,Colv=do.col.clust, Rowv=do.row.clust, scale="none", cexCol=cex.col.use, balanceColor = T, cexRow=cex.row.use, labCol=lab.col.use,
                         ColSideAnn=ColSideAnn,ColSideFun=function(x) showAnn(x),ColSideWidth=0.8, method="average", distfun = function(x) dist(x, method="euclidean"))
            }

          }
)

set.ifnull=function(x,y) {
  if(is.null(x)) x=y
  return(x)
}

kill.ifnull=function(x,message="Error:Execution Halted") {
  if(is.null(x)) {
    stop(message)
  }
}


expAlpha=function(mu,coefs) {
  logA=coefs$a
  logB=coefs$b
  return(exp(logA+logB*mu)/(1+(exp(logA+logB*mu))))
}

setGeneric("getWeightMatrix", function(object) standardGeneric("getWeightMatrix"))
setMethod("getWeightMatrix", "scR",
          function(object) {
            data=object@data[, rownames(object@drop.coefs)]
            data.humpAvg=apply(data,1,humpMean,min=object@drop.expr)
            wt.matrix=data.frame(t(sapply(data.humpAvg,expAlpha,object@drop.coefs)))
            colnames(wt.matrix)=colnames(data); rownames(wt.matrix)=rownames(data)
            wt.matrix[is.na(wt.matrix)]=0
            object@wt.matrix=wt.matrix
            wt1.matrix=data.frame(sapply(1:ncol(data),function(x)setWt1(data[,x],wt.matrix[,x],min=object@drop.expr)))
            colnames(wt1.matrix)=colnames(data); rownames(wt1.matrix)=rownames(data)
            wt1.matrix[is.na(wt1.matrix)]=0
            object@drop.wt.matrix=wt1.matrix
            return(object)
          }      
)

regression.sig=function(x,score,data,latent,code="rsem") {
  if(var(as.numeric(subc(data,code)[x,]))==0) {
    return(0)
  }
  latent=latent[grep(code,names(data))]
  data=rbind(subc(data,code),vsubc(score,code))
  rownames(data)[nrow(data)]="score"
  data2=data[c(x,"score"),]
  rownames(data2)[1]="fac"
  if (length(unique(latent))>1) {
    mylm=lm(score ~ fac+latent, data = data.frame(t(data2)))
  }
  else {
    mylm=lm(score ~ fac, data = data.frame(t(data2)))
  }
  return(coef(summary(mylm))["fac",3])
}

setGeneric("regulatorScore", function(object, candidate.reg, score.name, cells.use=NULL) standardGeneric("regulatorScore"))
setMethod("regulatorScore", "scR",
          function(object, candidate.reg, score.name, cells.use=NULL) {
            cells.use=set.ifnull(cells.use, colnames(object@data))
            candidate.reg=candidate.reg[candidate.reg%in%rownames(object@data)]
            my.score=retreiveScore(object,score.name)[cells.use]
            my.data=object@data[,cells.use]
            my.ident=object@ident[cells.use]
            reg.score=unlist(lapply(candidate.reg,regressionSig,score = my.score,data = my.data,latent = my.ident,code = "rsem"))
            names(reg.score)=candidate.reg
            return(reg.score)
          }
)





filter.genes.by.effect.size = function(object, cells.1, cells.2, genes.use,thresh.use, fxn.x=expMean){
  
  #print(paste0("Using ", fxn.x, " for calculating the effect size"))
  data.1=apply(object@data[genes.use,cells.1],1,fxn.x)
  data.2=apply(object@data[genes.use,cells.2],1,fxn.x)
  total.diff=abs(data.1-data.2)
  total.diff = total.diff[order(total.diff,decreasing=TRUE)]
  genes.diff = names(which(total.diff>thresh.use))
  if (length(genes.diff)==0){
    print("Warning: thresh.use is too large, no genes found")
    return(genes.diff)
  } 
  return(genes.diff)
  
}

add.genes.by.min.diff = function(object, cells.1, cells.2, genes.use, genes.diff, min.diff){
  
  print(paste0("Considering lowly expressed genes, based on min.diff = ", min.diff, " criterion"))
  print(paste0("Using ", object@is.expr, " as the minimum threshold for expression"))
  data.1 = apply(object@data[genes.use,cells.1],1,function(x) sum(x > object@is.expr) / length(x))
  data.2 = apply(object@data[genes.use,cells.2],1,function(x) sum(x > object@is.expr) / length(x))
  df = data.frame(A=data.1, B=data.2)
  genes.to.add = apply(df, 1, function(x) (x[1] <= 0.01 & x[2] > min.diff) | (x[2] <= 0.01 & x[1] > min.diff) )
  genes.diff2 = names(which(genes.to.add))
  genes.diff = union(genes.diff, genes.diff2)
  return(genes.diff)
  
}

#BIMOD TEST

setGeneric("diffExp.test", function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2), min.diff=NULL, test.max=NULL,...) standardGeneric("diffExp.test"))
setMethod("diffExp.test", "scR",
          function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2), min.diff=NULL, test.max=NULL,...) {
            genes.use=set.ifnull(genes.use,rownames(object@data))
            
            #Genes to test
            genes.diff = filter.genes.by.effect.size(object, cells.1, cells.2, genes.use, thresh.use,...)
            
            if (length(genes.diff)==0) return( data.frame(pval=numeric(0), diff=numeric(0), fracDiff=numeric(0)))
            if (!is.null(test.max)){
              genes.diff = genes.diff[1:min(test.max,length(genes.diff))]
            }
            
            if (!is.null(min.diff)){
              genes.diff = add.genes.by.min.diff(object, cells.1, cells.2, genes.use, genes.diff, min.diff)
            }
            
            to.return=bimod.diffExp.test(object@data[,cells.1],object@data[,cells.2],genes.diff,xmin=object@is.expr,...)
            print(to.return)
            to.return=to.return[order(to.return$pval,-abs(to.return$diff)),]
            return(to.return)
          } 
)

#ROC test

setGeneric("marker.test", function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2), min.diff=NULL, test.max=NULL,...) standardGeneric("marker.test"))
setMethod("marker.test", "scR",
          function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2), min.diff=NULL, test.max=NULL,...) {
            genes.use=set.ifnull(genes.use,rownames(object@data))
            
            #Genes to test
            genes.diff = filter.genes.by.effect.size(object, cells.1, cells.2, genes.use, thresh.use,...)
          
            if (length(genes.diff)==0) return( data.frame(pval=numeric(0), diff=numeric(0), fracDiff=numeric(0)))
            if (!is.null(test.max)){
              genes.diff = genes.diff[1:min(test.max,length(genes.diff))]
            }
            
            if (!is.null(min.diff)){
              genes.diff = add.genes.by.min.diff(object, cells.1, cells.2, genes.use, genes.diff)
            }
            
            to.return=marker.auc.test(object@data[,cells.1],object@data[,cells.2],genes.use,...)
            to.return=to.return[rev(order(abs(to.return$AUC-0.5))),]
            to.return$power=abs(to.return$AUC-0.5)
            return(to.return)
          } 
)

#TOBIT test
setGeneric("tobit.test", function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2), min.diff=NULL, test.max=NULL,...) standardGeneric("tobit.test"))
setMethod("tobit.test", "scR",
          function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2), min.diff=NULL, test.max=NULL,...) {
            genes.use=set.ifnull(genes.use,rownames(object@data))
            
            #Genes to test
            genes.diff = filter.genes.by.effect.size(object, cells.1, cells.2, genes.use, thresh.use,...)
            if (length(genes.diff)==0) return( data.frame(pval=numeric(0), diff=numeric(0), fracDiff=numeric(0)))
            if (!is.null(test.max)){
              genes.diff = genes.diff[1:min(test.max,length(genes.diff))]
            }
            
            #print(genes.diff)
            
            if (!is.null(min.diff)){
              genes.diff = add.genes.by.min.diff(object, cells.1, cells.2, genes.use, genes.diff, min.diff)
            }
            
            to.return=tobit.diffExp.test(object@data[,cells.1],object@data[,cells.2],genes.diff,...)
            to.return=to.return[order(to.return$pval,-abs(to.return$diff)),]
            return(to.return)
          } 
)

setGeneric("diff.t.test", function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2), min.diff=NULL, test.max=NULL,...) standardGeneric("diff.t.test"))
setMethod("diff.t.test", "scR",
          function(object, cells.1,cells.2,genes.use=NULL,thresh.use=log(2), min.diff=NULL, test.max=NULL,...) {
            
            genes.use=set.ifnull(genes.use,rownames(object@data))
            #Genes to test
            genes.diff = filter.genes.by.effect.size(object, cells.1, cells.2, genes.use, thresh.use,...)
            if (length(genes.diff)==0) return( data.frame(pval=numeric(0), diff=numeric(0), fracDiff=numeric(0)))
            if (!is.null(test.max)){
              genes.diff = genes.diff[1:min(test.max,length(genes.diff))]
            }
            
            #print(genes.diff)
            
            if (!is.null(min.diff)){
              genes.diff = add.genes.by.min.diff(object, cells.1, cells.2, genes.use, genes.diff, min.diff)
            }
            
            pval=unlist(lapply(genes.diff,function(x)t.test(object@data[x,cells.1],object@data[x,cells.2])$p.value))
            data.1=apply(object@data[genes.diff,cells.1],1,mean)
            data.2=apply(object@data[genes.diff,cells.2],1,mean)
            diff=(data.1-data.2)[genes.use]
            to.return=data.frame(pval,diff,row.names = genes.use)
            to.return=to.return[with(to.return, order(pval, -abs(diff))), ]
            return(to.return)
          } 
)

setGeneric("binomcount.test", function(object, cells.1,cells.2, thresh.use, TPM.mat, Count.mat,...) standardGeneric("binomcount.test"))
setMethod("binomcount.test", "scR",
          function(object, cells.1,cells.2, thresh.use, TPM.mat, Count.mat,...) {
            
            x=TPM.mat
            y=Count.mat
            
            #Enrichments in cells.1
            m0 = apply(x[, cells.2], 1, mean)
            n0 = apply(x[, cells.1], 1, mean) 
            no = median(apply(y[,cells.1], 2, sum)) / median(apply(y[,cells.2], 2, sum))
            m = m0*no
            n = n0/no
            pv1 = binompval(m/sum(m),sum(n),n)
            d1 <- data.frame(diff=log(n/m),pval=pv1)
            d1 <- subset(d1, diff >= thresh.use)
            d1 <- d1[order(d1$pval,decreasing=FALSE),]
            
            #Enrichments in cells.2
            m = m0*no; n = n0/no;
            pv2 = binompval(n/sum(n),sum(m),m)
            d2 <- data.frame(diff=log(n/m),pval=pv2)
            d2 <- subset(d2, diff <= -thresh.use)
            d2 <- d2[order(d2$pval,decreasing=FALSE),]
            
            d = rbind(d1, d2);
            d = d[order(d$pval, decreasing=FALSE),]
            return(d)
          } 
)

setGeneric("binomcount.digit.test", function(object, cells.1,cells.2, thresh.use, TPM.mat, Count.mat,xmin=0) standardGeneric("binomcount.digit.test"))
setMethod("binomcount.digit.test", "scR",
          function(object, cells.1,cells.2, thresh.use, TPM.mat, Count.mat,xmin=0) {
            
            x=TPM.mat
            y=Count.mat
            
            print(paste0("Using expression TPM threshold, xmin = ", xmin, " for binomial test. Pass xmin=xval argument to find.markers if you want to use a different value"))
            
            #Test for enrichments in cluster #1
            m = apply(x[, cells.2], 1, function(x) sum(x>xmin)) #Number of cells expressing marker in cluster #2
            m1 = m; m1[m==0]=1; # Regularization. Add a pseudo count of 1 to unexpressed markers to avoid false positives
            n = apply(x[, cells.1], 1, function(x) sum(x>xmin)) #Number of cells expressing marker in cluster #1
            #Find the probability of finding n or more cells +ve for marker in cluster 1 given the fraction in cluster 2
            pv1 = pbinom(n, length(cells.1), m1/length(cells.2), lower.tail = FALSE) + dbinom(n, length(cells.1), m1/length(cells.2))
            #m_mean = apply(x[, cells.2], 1, mean)
            #n_mean = apply(x[, cells.1], 1,mean)
            
            log_fold_express = log(n*length(cells.2)/(m*length(cells.1))) #log proportion of expressing cells
            d1 <- data.frame(diff=log_fold_express,pval=pv1)
            d1 <- subset(d1, diff >= thresh.use)
            d1 <- d1[order(d1$pval,decreasing=FALSE),]
            
            #Enrichments in cells.2
            n1 = n; n1[n==0]=1; # Regularization.  Add a pseudo count of 1 to unexpressed markers to avoid false positives
            #Find the probability of finding m or more cells +ve for marker in cluster 2 given the fraction in cluster 1
            pv2 = pbinom(m, length(cells.2), n1/length(cells.1), lower.tail=FALSE) + dbinom(m, length(cells.2), n1/length(cells.1))
            d2 <- data.frame(diff=log_fold_express,pval=pv2)
            d2 <- subset(d2, diff <= -thresh.use)
            d2 <- d2[order(d2$pval,decreasing=FALSE),]
            
            d = rbind(d1, d2);
            d = d[order(d$pval, decreasing=FALSE),]
            return(d)
          } 
)
# setGeneric("clustdiffgenes", function(object,genes.use=NULL,pvalue=0.01, thresh.use=log(2)) standardGeneric("clustdiffgenes"))
# setMethod("clustdiffgenes", "scR",
#           function(object,genes.use=NULL,pvalue=0.01, thresh.use=log(2)) {
#             genes.use=set.ifnull(genes.use,rownames(object@data))
#             if ( ! is.numeric(pvalue) ) stop("pvalue has to be a number between 0 and 1") else if (  pvalue < 0 | pvalue > 1 ) stop("pvalue has to be a number between 0 and 1")
#             to.return <- list()
#             x <- exp(object@data)-1
#             y <- object@count.data[rownames(object@data),colnames(object@data)]
#             ident.use=object@ident
#             
#             med_trans = median(apply(y, 2, sum)); #median number of transcripts
#             
#             for (i in as.numeric(levels(object@ident))){
#                if (sum(ident.use == i) == 0) next
#                m = if (sum(ident.use!=i) > 1) apply(x[, ident.use!=i], 1, mean) else x[, ident.use!=i]
#                n = if (sum(ident.use==i) > 1) apply(x[, ident.use==i], 1, mean) else x[, ident.use==i]
#                no = if (sum(ident.use==i) > 1) median(apply(y[,ident.use==i], 2, sum)) / median(apply(y[,ident.use!=i], 2, sum)) else sum(y[,ident.use == i])/median(apply(y[,ident.use!=i], 2, sum))
#                m = m*no
#                n = n*no
#                pv = binompval(m/sum(m),sum(n),n)
#                d <- data.frame(mean.all=m,mean.cl=n,fc=n/m,pv=pv)[order(pv,decreasing=FALSE),]
#                to.return[[paste("cl",i,sep=".")]] <- d[d$pv < pvalue,]
#             }
#             return(to.return)
#           } 
# )



#Random forest classification for binary categories
setGeneric("find.markers.rf.binary", function(object, ident.1,ident.2=NULL,ident.1.name=NULL,ident.2.name=NULL,genes.use=NULL,thresh.use=log(2), downsample.frac=0.8,min.num.cells=500, n.tree.use=100, genes.print=30) standardGeneric("find.markers.rf.binary"))
setMethod("find.markers.rf.binary", "scR",
          function(object, ident.1,ident.2=NULL,ident.1.name=NULL,ident.2.name=NULL,genes.use=NULL,thresh.use=log(2), downsample.frac=0.8, min.num.cells=500, n.tree.use=100, genes.print=30) {
            genes.use=set.ifnull(genes.use,rownames(object@data))
            ident.use=object@ident
            cells.1=names(ident.use[which(ident.use%in%ident.1)])
            ident.1.name = set.ifnull(ident.1.name, paste0(ident.1, collapse="_"))
            if (is.null(ident.2)) {
              ident.2="rest"
              ident.2.name = set.ifnull(ident.2.name, "rest")
              cells.2=names(ident.use)
              cells.2=cells.2[!(cells.2%in%cells.1)]
            } else {
              ident.2.name = set.ifnull(ident.2.name, paste0(ident.2,collapse="_"))
              cells.2=names(ident.use[which(ident.use%in%ident.2)])
            }
            
            #Establish training and testing sets
            train.1.size = min(min.num.cells, downsample.frac*length(cells.1))
            train.2.size = min(min.num.cells, downsample.frac*length(cells.2))
            train.1 = sample(cells.1, train.1.size); test.1 = setdiff(cells.1, train.1)
            train.2 = sample(cells.2, train.2.size); test.2 = setdiff(cells.2, train.2)
            
            #0 == ident.1.name, 1 == ident.2.name
            labels.train = factor(c(rep(ident.1.name,length(train.1)), rep(ident.2.name, length(train.2))));
            test.labels = factor(c(rep(ident.1.name,length(test.1)), rep(ident.2.name, length(test.2))));
            data.train = as.data.frame(t(object@data[,c(train.1, train.2)]))
            data.test = as.data.frame(t(object@data[,c(test.1, test.2)]))
            
            rF_Model <- randomForest(x=data.train,y=labels.train,ntree=n.tree.use,importance=TRUE, keep.inbag=TRUE,replace=FALSE) 
            print("Training error")
            print(rF_Model$confusion)
            test.prediction = predict(rF_Model, data.test)
            print("Test prediction error")
            print(table(test.labels, test.prediction))
            
            #Class predictors for ident.1
            li<-getLocalIncrements(rF_Model, data.train, binAsReg = FALSE, mcls = ident.1.name)
            fc<-featureContributions(rF_Model, li, data.train, mClass = ident.1.name)
            sort.by.imp1 = sort(colMeans(fc[labels.train==ident.1.name,]), decreasing=TRUE) 
            sort.by.imp1 = sort.by.imp1 * sign(colMeans(data.train[labels.train==ident.1.name,names(sort.by.imp1)]) - colMeans(data.train[labels.train != ident.1.name,names(sort.by.imp1)]))
            sort.by.imp1 = sort.by.imp1 / max(abs(sort.by.imp1))
            
            
            #Class predictors for ident.2
            li<-getLocalIncrements(rF_Model, data.train, binAsReg = FALSE, mcls = ident.2.name)
            fc<-featureContributions(rF_Model, li, data.train, mClass = ident.2.name)
            sort.by.imp2 = sort(colMeans(fc[labels.train==ident.2.name,]), decreasing=TRUE) 
            sort.by.imp2 = sort.by.imp2 * sign(colMeans(data.train[labels.train==ident.2.name,names(sort.by.imp2)]) - colMeans(data.train[labels.train != ident.2.name,names(sort.by.imp2)]))
            sort.by.imp2 = sort.by.imp2/ max(abs(sort.by.imp2))
            
            df1=data.frame(importance=sort.by.imp1[1:30], genes = factor(names(sort.by.imp1[1:30]), levels=names(sort.by.imp1[1:30])) )
            p1=ggplot(df1, aes(y=importance, x= genes)) + geom_bar(stat="identity",fill="red") + ggtitle(paste(ident.1.name, " Discriminatory genes")) + theme_bw() + gg.xax() + gg.yax()
            
            df2=data.frame(importance=sort.by.imp2[1:30], genes = factor(names(sort.by.imp2[1:30]), levels=names(sort.by.imp2[1:30])) )
            p2=ggplot(df2, aes(y=importance, x= genes)) + geom_bar(stat="identity",fill="blue") + ggtitle(paste(ident.2.name, " Discriminatory genes")) + theme_bw() + gg.xax() + gg.yax() +
            
            pdf(paste0(ident.1.name,"_vs_",ident.2.name,"_RF.pdf"), w=12,h=12)
            print(grid.arrange(p1,p2, ncol=1))
            dev.off()
            
            #Assemble top n positive markers for each cell type
            to.return = as.data.frame(matrix(0, nrow=2*genes.print, ncol=4))
            
            
            #Mean number of transcripts per cell
            if (!is.null(attr(object,"count.data"))){
              colnames(to.return) = c("Gene",paste0("nTrans_", ident.1.name), paste0("nTrans_", ident.2.name), "Marker_for")
              genes.1 = names(sort.by.imp1)[sort.by.imp1 > 0][1:genes.print]
              genes.2 = names(sort.by.imp2)[sort.by.imp2 > 0][1:genes.print]
              nTrans.1 = apply(object@count.data[c(genes.1, genes.2), c(cells.1)], 1, function(x) round(mean(x),3))
              nTrans.2 = apply(object@count.data[c(genes.1,genes.2), c(cells.2)], 1, function(x) round(mean(x),3))
              to.return$Gene = c(genes.1,genes.2)
              to.return[, 2] = nTrans.1
              to.return[, 3] = nTrans.2
              to.return[,4] = c(rep(ident.1.name, genes.print), rep(ident.2.name, genes.print))
            }
            
            return(to.return)
          } 
)


#Random forest classification for binary categories
setGeneric("find.all.markers.rf", function(object, clust.include=NULL,genes.use=NULL, downsample.frac=0.8,min.num.cells=200, n.tree.use=100, genes.print=30,...) standardGeneric("find.all.markers.rf"))
setMethod("find.all.markers.rf", "scR",
          function(object, ident.1,ident.2=NULL,clust.include = NULL, downsample.frac=0.8, min.num.cells=200, n.tree.use=100, genes.print=30,...) {
            genes.use=set.ifnull(genes.use,rownames(object@data))
            ident.use=object@ident
            clust.include = set.ifnull(clust.include, unique(levels(object@ident)))
            clust.include = ainb(clust.include,  unique(levels(object@ident)))
            
            cells.train = c(); labels.train = c()
            cells.test = c(); labels.test = c()
            
            for (clust in clust.include){
              cells.in.clust = names(ident.use)[which(ident.use %in% clust)]
              train.size = min(min.num.cells, downsample.frac*length(cells.in.clust))
              train.set = sample(cells.in.clust, train.size); test.set = setdiff(cells.in.clust, train.set)
              cells.train = c(cells.train, train.set); cells.test=c(cells.test, test.set)
              labels.train = c(labels.train, rep(clust, length(train.set))); labels.test = c(labels.test, rep(clust, length(test.set)))
            }
            
            
            labels.train = factor(labels.train); labels.test = factor(labels.test)
            data.train = as.data.frame(t(object@data[,cells.train]))
            data.test = as.data.frame(t(object@data[,cells.test]))
            
            rF_Model <- randomForest(x=data.train,y=labels.train,ntree=n.tree.use,importance=TRUE, keep.inbag=TRUE,replace=FALSE, classwt = length(labels.train)*as.numeric(1/sqrt(table(labels.train)))) 
            print("Training error")
            print(rF_Model$confusion)
            test.prediction = predict(rF_Model, data.test)
            print("Test prediction error")
            print(table(labels.test, test.prediction))
            
            li<-getLocalIncrements(rF_Model, data.train)
            fc<-featureContributions(rF_Model, li, data.train)
            
            
            if (!dir.exists("Figs_RF_signatures")){
              dir.create("Figs_RF_signatures")
            }
            
            #Assemble top n positive markers for each cell type
            to.return = as.data.frame(matrix(0, nrow=length(clust.include)*genes.print, ncol=2+length(clust.include)))
            colnames(to.return) = c("Gene","Marker_for", paste0("nTrans_", clust.include))
            l=1;
            for (i in clust.include){
              print(i)
              sort.by.imp = sort(colMeans(fc[labels.train==i,]), decreasing=TRUE)
              sort.by.imp = sort.by.imp * sign(colMeans(data.train[labels.train==i,names(sort.by.imp)]) - colMeans(data.train[labels.train != i,names(sort.by.imp)]))
              df=data.frame(importance=sort.by.imp[1:genes.print], genes = factor(names(sort.by.imp[1:genes.print]), levels=names(sort.by.imp[1:genes.print])) )
              p=ggplot(df, aes(y=importance, x= genes)) + geom_bar(stat="identity",fill="blue") + theme_bw() + gg.xax() + gg.yax()
              pdf(paste0("Figs_RF_signatures/",i,"_markers.pdf"),w=10,h=6)
              print(p)
              dev.off()
              
              genes = names(sort.by.imp)[sort.by.imp > 0][1:genes.print]
              strtrow = ((l-1)*genes.print+1); endrow = l*genes.print
              to.return[strtrow:endrow,1] = genes
              to.return[strtrow:endrow,2] = i
              
              for (j in clust.include){
                cells.in.cluster = names(ident.use)[which(ident.use == j)]
                vec.exp = apply(object@count.data[genes, cells.in.cluster], 1, function(x) round(mean(x),3)) 
               to.return[strtrow:endrow, paste0("nTrans_",j)] = vec.exp
              }
              
              l=l+1;
            }
            
            
            #Average expression by cluster
            ExpMat = matrix(0, nrow=length(genes.use), ncol = length(ident.levels))
            rownames(ExpMat) = genes.use; colnames(ExpMat) = paste0("nTrans_",ident.levels)
            
            for (i in clust.include){
              cells.in.cluster = cells.use[which(ident.use== i)]
              vec.exp = apply(object@count.data[genes.use, cells.in.cluster], 1, function(x) round(mean(x),3)) 
              ExpMat[, i] = vec.exp
            }
            
            #Mean number of transcripts per cluster
            if (!is.null(attr(object,"count.data"))){
              
              
              genes.2 = names(sort.by.imp2)[sort.by.imp2 > 0][1:genes.print]
              nTrans.1 = apply(object@count.data[c(genes.1, genes.2), c(cells.1)], 1, function(x) round(mean(x),3))
              nTrans.2 = apply(object@count.data[c(genes.1,genes.2), c(cells.2)], 1, function(x) round(mean(x),3))
              to.return$Gene = c(genes.1,genes.2)
              to.return[, 2] = nTrans.1
              to.return[, 3] = nTrans.2
              to.return[,4] = c(rep(ident.1.name, genes.print), rep(ident.2.name, genes.print))
            }
            
            
        
            
            
            return(to.return)
          } 
)



setGeneric("scde.test", function(object, cells.1,cells.2,ident.1, ident.2,genes.use=NULL,thresh.use=log(2), min.diff=NULL, test.max=NULL) standardGeneric("scde.test"))
setMethod("scde.test", "scR",
          function(object, cells.1,cells.2,ident.1,ident.2,genes.use=NULL,thresh.use=log(2), min.diff=NULL, test.max=NULL) {
            genes.use=set.ifnull(genes.use,object@var.genes)
            data.1=apply(object@data[genes.use,cells.1],1,mean)
            data.2=apply(object@data[genes.use,cells.2],1,mean)
            total.diff=abs(data.1-data.2)
            total.diff = total.diff[order(total.diff,decreasing=TRUE)]
            genes.diff = names(which(total.diff>thresh.use))
            
            if (length(genes.diff)==0) genes.diff = names(total.diff)[1]
            
            if (!is.null(test.max)){
              genes.diff = genes.diff[1:min(test.max,length(genes.diff))]
            }
            
            
            if (!is.null(min.diff)){
              data.1 = apply(object@data[genes.use,cells.1],1,function(x) sum(x > 0) / length(x))
              data.2 = apply(object@data[genes.use,cells.2],1,function(x) sum(x > 0) / length(x))
              df = data.frame(A=data.1, B=data.2)
              genes.to.add = apply(df, 1, function(x) (x[1] <= 0.01 & x[2] > min.diff) | (x[2] <= 0.01 & x[1] > min.diff) )
              genes.diff2 = names(which(genes.to.add))
              genes.diff = union(genes.diff, genes.diff2)
            }
            
            print(length(genes.diff))
            
            data.use0 = object@count.data[genes.use, c(cells.1, cells.2)]
            data.use1 = object@count.data[genes.diff, c(cells.1,cells.2)]
            print(dim(data.use0))
            sg = factor(c(rep(ident.1,length(cells.1)), rep(ident.2, length(cells.2))), levels=c(ident.1,ident.2))
            names(sg) = c(cells.1, cells.2);
            print("Cell table: ")
            print(table(sg))
            
            n.cores <- 8
            o.ifm <- scde.error.models(counts=round(data.use0), groups=sg, n.cores=n.cores, threshold.segmentation = T, save.crossfit.plots = F, save.model.plots = F, verbose=1);
            valid.cells <- o.ifm$corr.a > 0;
            print("Valid Cells: ")
            table(valid.cells)
            
            o.ifm <- o.ifm[valid.cells,];
            groups <- sg[valid.cells];
            
            # estimate gene expression prior
            o.prior <- scde.expression.prior(models=o.ifm, counts=round(data.use0[,valid.cells]), length.out = 400, show.plot=T)
            
            #run differential expression tests on all genes
            to.return <- scde.expression.difference(o.ifm, round(data.use1[,valid.cells]), o.prior, groups = groups, n.randomizations = 100, n.cores = n.cores, verbose = 1)
            #to.return <- subset(to.return, (lb > 0 & ub > 0) | (lb < 0 & ub < 0)) #Only significant genes
            return(to.return)
          } 
)



setGeneric("batch.gene", function(object, idents.use,genes.use=NULL,use.imputed=FALSE,auc.cutoff=0.6) standardGeneric("batch.gene"))
setMethod("batch.gene", "scR",
          function(object, idents.use,genes.use=NULL,use.imputed=FALSE,auc.cutoff=0.6) {
            batch.genes=c()
            genes.use=set.ifnull(genes.use,rownames(object@data))
            for(ident in idents.use ) {
              cells.1=names(object@ident)[object@ident==ident]
              cells.2=names(object@ident)[object@ident!=ident]
              if ((length(cells.1)<5)|(length(cells.2)<5)) {
                break;
              }
              markers.ident=marker.test(object,cells.1,cells.2,genes.use)
              batch.genes=unique(c(batch.genes,rownames(subset(markers.ident,myAUC>auc.cutoff))))
            }
            return(batch.genes)
          } 
)




#to remove
setGeneric("retreiveScore", function(object, score.name) standardGeneric("retreiveScore"))
setMethod("retreiveScore", "scR",
          function(object, score.name) {
            my.score=object@gene.scores[,score.name]
            names(my.score)=rownames(object@gene.scores)
            return(my.score)
          } 
)

setGeneric("which.cells", function(object,value=1, id=NULL) standardGeneric("which.cells"))
setMethod("which.cells", "scR",
          function(object, value=1,id=NULL) {
            id=set.ifnull(id,"ident")
            if (id=="ident") {
              data.use=object@ident
            } else {
              if (id %in% colnames(object@data.info)) {
                data.use=object@data.info[,id]; names(data.use)=rownames(object@data.info)
              }
            }
            return(names(data.use[which(data.use%in%value)]))
          } 
)

setGeneric("set.all.ident", function(object,id=NULL) standardGeneric("set.all.ident"))
setMethod("set.all.ident", "scR",
          function(object, id=NULL) {
            id=set.ifnull(id,"orig")
            if (id %in% colnames(object@data.info)) {
              cells.use=rownames(object@data.info)
              ident.use=object@data.info[,id]
              object@ident = factor(ident.use)
              names(object@ident) = rownames(object@data.info)
              #object=set.ident(object,cells.use,ident.use)
            }
            return(object)
          } 
)

setGeneric("set.ident", function(object,cells.use=NULL,ident.use=NULL) standardGeneric("set.ident"))
setMethod("set.ident", "scR",
          function(object, cells.use=NULL,ident.use=NULL) {
            cells.use=set.ifnull(cells.use,object@cell.names)
            if (length(anotinb(cells.use,object@cell.names)>0)) {
              print(paste("ERROR : Cannot find cells ",anotinb(cells.use,object@cell.names)))
            }
          
            ident.new=anotinb(ident.use,levels(object@ident))
            object@ident=factor(object@ident,levels = unique(c(as.character(object@ident),as.character(ident.new))))
            object@ident[cells.use]=ident.use
            object@ident=drop.levels(object@ident)
            return(object)
          } 
)

#likely delete
setGeneric("retreiveCellInfo", function(object, cells.use=NULL,id=NULL) standardGeneric("retreiveCellInfo"))
setMethod("retreiveCellInfo", "scR",
          function(object, cells.use=NULL,id=NULL) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            id=set.ifnull(id,"ident")
            if (id=="ident") {
              data.use=object@ident
            }
            if (id %in% colnames(object@data.info)) {
              data.use=object@data.info[,id]; names(data.use)=rownames(object@data.info)
            }
            return(data.use[cells.use])
          } 
)

#likely delete
setGeneric("retreiveCluster", function(object, names=NULL) standardGeneric("retreiveCluster"))
setMethod("retreiveCluster", "scR",
          function(object, names=NULL) {
            names=set.ifnull(names,colnames(object@data))
            if (names[1] %in% rownames(object@data)) return(object@kmeans.obj[[1]]$cluster[names])
            if (names[1] %in% colnames(object@data)) return(object@kmeans.col[[1]]$cluster[names])
          } 
)

setGeneric("posterior.plot", function(object, name) standardGeneric("posterior.plot"))
setMethod("posterior.plot", "scR",
          function(object, name) {
            post.names=colnames(subc(object@mix.probs,name))
            vlnPlot(object,post.names,inc.first=TRUE,inc.final=TRUE,by.k=TRUE)
            
            
          } 
)

map.cell.score=function(gene,gene.value,insitu.bin,mu,sigma,alpha) {
  code.1=paste(gene,insitu.bin,sep=".")
  mu.use=mu[paste(code.1,"mu",sep="."),1]
  sigma.use=sigma[paste(code.1,"sigma",sep="."),1]
  alpha.use=alpha[paste(code.1,"alpha",sep="."),1]
  bin.prob=unlist(lapply(1:length(insitu.bin),function(x) dnorm(gene.value,mean = mu.use[x],sd = sigma.use[x],log = TRUE) + log(alpha.use[x])))
  return(bin.prob)
}

setGeneric("map.cell",  function(object,cell.name,do.plot=FALSE,safe.use=TRUE,text.val=NULL,do.rev=FALSE) standardGeneric("map.cell"))
setMethod("map.cell", "scR",
          function(object,cell.name,do.plot=FALSE,safe.use=TRUE,text.val=NULL,do.rev=FALSE) {
            insitu.matrix=object@insitu.matrix
            insitu.genes=colnames(insitu.matrix)
            insitu.genes=insitu.genes[insitu.genes%in%rownames(object@imputed)]
            insitu.use=insitu.matrix[,insitu.genes]
            imputed.use=object@imputed[insitu.genes,]
            safe_fxn=sum
            if (safe.use) safe_fxn=log_add
            
            all.needed.cols=unique(unlist(lapply(insitu.genes,function(x) paste(x,insitu.use[,x],"post",sep="."))))
            missing.cols=which(!(all.needed.cols%in%colnames(object@mix.probs)))
            if (length(missing.cols)>0) print(paste("Error : ", all.needed.cols[missing.cols], " is missing from the mixture fits",sep=""))
            all.probs=data.frame(sapply(insitu.genes,function(x) log(as.numeric(object@mix.probs[cell.name,paste(x,insitu.use[,x],"post",sep=".")]))))
            scale.probs=t(t(all.probs)-apply(t(all.probs),1,log_add))
            scale.probs[scale.probs<(-9.2)]=(-9.2)
            #head(scale.probs)
            total.prob=exp(apply(scale.probs,1,safe_fxn))
            total.prob=total.prob/sum(total.prob)
            if (do.plot) {
              #plot(total.prob,main=cell.name)
              par(mfrow=c(1,2))
              txt.matrix=matrix(rep("",64),nrow=8,ncol=8)
              if (!is.null(text.val)) txt.matrix[text.val]="X"
              if (do.rev) scale.probs=scale.probs[unlist(lapply(0:7,function(x)seq(1,57,8)+x)),] 
              aheatmap(matrix(total.prob,nrow=8,ncol=8),Rowv=NA,Colv=NA,txt=txt.matrix,col=bwCols)   
              aheatmap(scale.probs,Rowv=NA,Colv=NA)
              rp()
            }
            return(total.prob)
          } 
)


setGeneric("refined.mapping",  function(object,genes.use) standardGeneric("refined.mapping"))
setMethod("refined.mapping", "scR",
          function(object,genes.use) {
            
            genes.use=ainb(genes.use, rownames(object@imputed))            
            cells.max=t(sapply(colnames(object@data),function(x) exact.cell.centroid(object@final.prob[,x])))
            all.mu=sapply(genes.use,function(gene) sapply(1:64, function(bin) mean(as.numeric(object@imputed[gene,fetch.closest(bin,cells.max,2*length(genes.use))]))))
            all.cov=list(); for(x in 1:64) all.cov[[x]]=cov(t(object@imputed[genes.use,fetch.closest(x,cells.max,2*length(genes.use))]))
            
            mv.probs=sapply(colnames(object@data),function(my.cell) sapply(1:64,function(bin) slimdmvnorm(as.numeric(object@imputed[genes.use,my.cell]),as.numeric(all.mu[bin,genes.use]),all.cov[[bin]])))
            mv.final=exp(sweep(mv.probs,2,apply(mv.probs,2,log_add)))
            object@final.prob=data.frame(mv.final)            
            return(object)
          } 
)


setGeneric("initial.mapping", function(object,cells.use=NULL,do.plot=FALSE,safe=FALSE) standardGeneric("initial.mapping"))
setMethod("initial.mapping", "scR",
          function(object,cells.use=NULL,do.plot=FALSE,safe=FALSE) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            every.prob=sapply(cells.use,function(x)map.cell(object,x,do.plot=FALSE,safe.use=safe))
            object@final.prob=data.frame(every.prob)
            rownames(object@final.prob)=paste("bin.",rownames(object@final.prob),sep="")
            return(object)
          } 
)

setGeneric("calc.insitu", function(object,gene,do.plot=TRUE,do.write=FALSE,write.dir="~/window/insitu/",col.use=bwCols,do.norm=FALSE,cells.use=NULL,do.return=FALSE,probs.min=0,do.log=FALSE, use.imputed=FALSE, bleach.use=0) standardGeneric("calc.insitu"))
setMethod("calc.insitu", "scR",
          function(object,gene,do.plot=TRUE,do.write=FALSE,write.dir="~/window/insitu/",col.use=bwCols,do.norm=FALSE,cells.use=NULL,do.return=FALSE,probs.min=0,do.log=FALSE,use.imputed=FALSE,bleach.use=0) {
            cells.use=set.ifnull(cells.use,colnames(object@final.prob))
            probs.use=object@final.prob
            data.use=exp(object@data)-1  
            if (use.imputed) data.use=exp(object@imputed)-1
            cells.use=cells.use[cells.use%in%colnames(probs.use)]; cells.use=cells.use[cells.use%in%colnames(data.use)]
            #insilico.stain=matrix(unlist(lapply(1:64,function(x) sum(probs.use[x,]*data.use[gene,]))),nrow=8,ncol=8)
            insilico.vector=unlist(lapply(1:64,function(x) sum(as.numeric(probs.use[x,cells.use])*as.numeric(data.use[gene,cells.use]))))
            probs.total=apply(probs.use,1,sum)
            probs.total[probs.total<probs.min]=probs.min
            insilico.stain=(matrix(insilico.vector/probs.total,nrow=8,ncol=8))
            if (do.log) insilico.stain=log(insilico.stain+1)
            if (bleach.use > 0) {
              insilico.stain=insilico.stain-bleach.use
              insilico.stain=minmax(insilico.stain,min=0,max=1e6)
            }
            if (do.norm) insilico.stain=(insilico.stain-min(insilico.stain))/(max(insilico.stain)-min(insilico.stain))
            title.use=gene
            if (gene %in% colnames(object@insitu.matrix)) {
              pred.use=prediction(insilico.vector/probs.total,object@insitu.matrix[,gene],0:1)
              perf.use=performance(pred.use,"auc")
              auc.use=round(perf.use@y.values[[1]],3)
              title.use=paste(gene,sep=" ")
            }
            if (do.write) {
              write.table(insilico.stain,paste(write.dir,gene,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
            }
            if (do.plot) {
              aheatmap(insilico.stain,Rowv=NA,Colv=NA,col=col.use, main=title.use)
            }
            if (do.return) {
              return(as.vector(insilico.stain))
            }
            return(object)
          } 
)

setGeneric("fit.gene.k", function(object, gene, do.k=3,num.iter=2,do.plot=FALSE,genes.use=NULL,start.pct=NULL) standardGeneric("fit.gene.k"))
setMethod("fit.gene.k", "scR",
          function(object, gene, do.k=3,num.iter=2,do.plot=FALSE,genes.use=NULL,start.pct=NULL) {
            data=object@imputed            
            data.use=data[gene,]
            names(data.use)=colnames(data.use)
            scale.data=t(scale(t(object@imputed)))
            genes.use=set.ifnull(genes.use,rownames(scale.data))
            genes.use=genes.use[genes.use%in%rownames(scale.data)]
            scale.data=scale.data[genes.use,]
            #print(genes.use)
            #seed the k-means based on the 0'th component
            data.cut=as.numeric(data.use[gene,])
            cell.ident=as.numeric(cut(data.cut,do.k))
            if (!(is.null(start.pct))) {
              cell.ident=rep(1,length(data.cut))
              cell.ident[data.cut>quantile(data.cut,1-start.pct)]=2
            }
            cell.ident=order(tapply(as.numeric(data.use),cell.ident,mean))[cell.ident]
            #cell.ident=sample(cell.ident)
            ident.table=table(cell.ident)
            #            if (do.plot) {
            #              par(mfrow=c(2,2))
            #              hist(as.numeric(data.use),probability = TRUE,ylim=c(0,1),xlab=gene,main=gene);
            #              for(i in 1:do.k) lines(seq(-10,10,0.01),(ident.table[i]/sum(ident.table)) * dnorm(seq(-10,10,0.01),mean(as.numeric(data.use[cell.ident==i])),sd(as.numeric(data.use[cell.ident==i]))),col=i,lwd=2); 
            #            }
            if (num.iter > 0) {
              for(i2 in 1:num.iter) {
                cell.ident=iter.k.fit(scale.data,cell.ident,data.use)
                ident.table=table(cell.ident)
                #                if (do.plot) {
                #                  hist(as.numeric(data.use),probability = TRUE,ylim=c(0,1),xlab=gene,main=gene);
                #                  for(i in 1:do.k) lines(seq(-10,10,0.01),(ident.table[i]/sum(ident.table)) * dnorm(seq(-10,10,0.01),mean(as.numeric(data.use[cell.ident==i])),sd(as.numeric(data.use[cell.ident==i]))),col=i,lwd=2); 
                #                }
              }
            }
            ident.table=table(cell.ident)
            #for(i in 2:do.k) ident.table[i]=ident.table[1]
            raw.probs=t(sapply(data.use,function(y) unlist(lapply(1:do.k,function(x) ((ident.table[x]/sum(ident.table))*dnorm(y,mean(as.numeric(data.use[cell.ident==x])),sd(as.numeric(data.use[cell.ident==x]))))))))
            norm.probs=raw.probs/apply(raw.probs,1,sum)
            colnames(norm.probs)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"post",sep=".")))
            norm.probs=cbind(norm.probs,cell.ident); colnames(norm.probs)[ncol(norm.probs)]=paste(gene,".ident",sep="")
            new.mix.probs=data.frame(minusc(object@mix.probs,paste(gene,".",sep="")),row.names = rownames(object@mix.probs)); colnames(new.mix.probs)[1]="nGene"
            object@mix.probs=cbind(new.mix.probs,norm.probs)
            
            if (do.plot) {
              nCol=2
              #(object,gene,by.k=TRUE,use.imputed=TRUE)
              num.row=floor((do.k+1)/nCol-1e-5)+1
              #par(mfrow=c(num.row,nCol))         
              #par(mfrow=c(1,1))
              hist(as.numeric(data.use),probability = TRUE,ylim=c(0,1),xlab=gene,main=gene);
              for(i in 1:do.k) lines(seq(-10,10,0.01),(ident.table[i]/sum(ident.table)) * dnorm(seq(-10,10,0.01),mean(as.numeric(data.use[cell.ident==i])),sd(as.numeric(data.use[cell.ident==i]))),col=i,lwd=2); 
              #unlist(lapply(1:do.k,function(x) plot(as.numeric(data.use),norm.probs[,x],ylab=paste("Posterior for Component ",x-1,sep=""),main=gene)))     
              #barplot(ident.table,main=round(ident.table[2]/sum(ident.table),3))
            }
            return(object)
          }
)

iter.k.fit=function(scale.data,cell.ident,data.use) {
  means.all=sapply(sort(unique(cell.ident)),function(x)apply(scale.data[,cell.ident==x],1,mean))
  all.dist=data.frame(t(sapply(1:ncol(scale.data),function(x) unlist(lapply(sort(unique(cell.ident)),function(y)dist(rbind(scale.data[,x],means.all[,y])))))))
  cell.ident=apply(all.dist,1,which.min)
  cell.ident=order(tapply(as.numeric(data.use),cell.ident,mean))[cell.ident]
  return(cell.ident)
}


setGeneric("fit.gene.mix", function(object, gene, do.k=3,use.mixtools=TRUE,do.plot=FALSE,plot.with.imputed=TRUE,min.bin.size=10) standardGeneric("fit.gene.mix"))
setMethod("fit.gene.mix", "scR",
          function(object, gene, do.k=3,use.mixtools=TRUE,do.plot=FALSE,plot.with.imputed=TRUE,min.bin.size=10) {
            require(mixtools)
            data.fit=as.numeric(object@imputed[gene,])
            mixtools.fit=normalmixEM(data.fit,k=do.k)
            comp.order=order(mixtools.fit$mu)
            mixtools.posterior=data.frame(mixtools.fit$posterior[,comp.order])
            colnames(mixtools.posterior)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"post",sep=".")))
            
            #mixtools.mu=data.frame(mixtools.fit$mu[comp.order])
            #mixtools.sigma=data.frame(mixtools.fit$sigma[comp.order])
            #mixtools.alpha=data.frame(mixtools.fit$lambda[comp.order])
            #rownames(mixtools.mu)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"mu",sep=".")))
            #rownames(mixtools.sigma)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"sigma",sep=".")))
            #rownames(mixtools.alpha)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"alpha",sep=".")))
            #object@mix.mu = rbind(minusr(object@mix.mu,gene), mixtools.mu); 
            #object@mix.sigma = rbind(minusr(object@mix.sigma,gene), mixtools.sigma); 
            #o#bject@mu.alpha =rbind(minusr(object@mu.alpha,gene), mixtools.alpha); 
            
            if (do.plot) {
              nCol=2
              num.row=floor((do.k+1)/nCol-1e-5)+1
              par(mfrow=c(num.row,nCol))
              plot.mixEM(mixtools.fit,which=2)
              plot.data=as.numeric(object@imputed[gene,])
              if (!plot.with.imputed) plot.data=as.numeric(object@data[gene,])
              unlist(lapply(1:do.k,function(x) plot(plot.data,mixtools.posterior[,x],ylab=paste("Posterior for Component ",x-1,sep=""),xlab=gene,main=gene)))
            }
            new.mix.probs=data.frame(minusc(object@mix.probs,paste(gene,".",sep="")),row.names = rownames(object@mix.probs)); colnames(new.mix.probs)[1]="nGene"
            object@mix.probs=cbind(new.mix.probs,mixtools.posterior)
            return(object)
          } 
)

lasso.fxn = function(lasso.input,genes.obs,s.use=20,gene.name=NULL,do.print=FALSE,gram=TRUE) {
  lasso.model=lars(lasso.input,as.numeric(genes.obs),type="lasso",max.steps = s.use*2,use.Gram=gram)
  #lasso.fits=predict.lars(lasso.model,lasso.input,type="fit",s=min(s.use,tail(lasso.model$df,2)[1]))$fit
  lasso.fits=predict.lars(lasso.model,lasso.input,type="fit",s=s.use)$fit
  if (do.print) print(gene.name)
  return(lasso.fits)  
}

setGeneric("addImputedScore", function(object, genes.use=NULL,cells.use=NULL,genes.fit=NULL,s.use=20,do.print=FALSE,gram=TRUE, add.imputed=TRUE) standardGeneric("addImputedScore"))
setMethod("addImputedScore", "scR",
          function(object, genes.use=NULL,cells.use=NULL,genes.fit=NULL,s.use=20,do.print=FALSE,gram=TRUE, add.imputed=TRUE) {
            genes.use=set.ifnull(genes.use,object@var.genes)
            genes.fit=set.ifnull(genes.fit,rownames(object@data))
            genes.use=genes.use[genes.use%in%rownames(object@data)]
            genes.fit=genes.fit[genes.fit%in%rownames(object@data)]
            
            cells.use = set.ifnull(cells.use, colnames(object@data))
            
            #genes.use = genes.use[rowSums(object@data[genes.use,cells.use]) > 0]
            #genes.fit = genes.fit[rowSums(object@data[genes.fit,cells.use]) > 0]
            print(paste0("Fitting linear models to ", length(genes.fit), " genes"))
            lasso.fits=data.frame(t(sapply(genes.fit,function(x)lasso.fxn(t(object@data[genes.use[genes.use!=x],]),object@data[x,],s.use=s.use,x,do.print,gram))))
            genes.old=genes.fit[genes.fit%in%rownames(object@imputed)]
            genes.new=genes.fit[!(genes.fit%in%rownames(object@imputed))]
            
            if (length(object@imputed)==0){
              if (length(genes.old)>0) object@imputed[genes.old,]=lasso.fits[genes.old,]
              object@imputed=rbind(object@imputed,lasso.fits[genes.new,])
            } else {
              if (length(genes.old)>0) temp = lasso.fits[genes.old,];
              temp = rbind(temp, lasso.fits[genes.new,])
              object@imputed = cbind(object@imputed, temp)
            }
           
            if (add.imputed){
              return(object)
            } else {
              return(object@imputed)
            }
          }
)    



setGeneric("getNewScore", function(object, score.name,score.genes, cells.use=NULL, score.func=weighted.mean,scramble=FALSE, no.tech.wt=FALSE, biol.wts=NULL,use.scaled=FALSE) standardGeneric("getNewScore"))
setMethod("getNewScore", "scR",
          function(object, score.name,score.genes, cells.use=NULL, score.func=weighted.mean,scramble=FALSE, no.tech.wt=FALSE, biol.wts=NULL,use.scaled=FALSE) {
            data.use=object@data; if (use.scaled) data.use=object@scale.data #minmax(object@scale.data,min = -2,max=2)
            if (!(no.tech.wt)) score.genes=score.genes[score.genes%in%rownames(object@wt.matrix)]
            if (no.tech.wt) score.genes=score.genes[score.genes%in%rownames(data.use)]
            cells.use=set.ifnull(cells.use,colnames(data.use))
            wt.matrix=data.frame(matrix(1,nrow=length(score.genes),ncol = length(cells.use),dimnames = list(score.genes,cells.use))); 
            if (!(no.tech.wt)) wt.matrix=object@drop.wt.matrix[score.genes,cells.use]
            score.data=data.use[score.genes,cell.ids]
            if(scramble) {
              score.data=score.data[,sample(ncol(score.data))]
            }
            #wt.matrix=wt.matrix*(wt.matrix) ???
            if (no.tech.wt) wt.matrix[wt.matrix<1]=1
            biol.wts=set.ifnull(biol.wts,rep(1,nrow(wt.matrix)))
            if (mean(biol.wts)==1) names(biol.wts)=score.genes
            my.scores=unlist(lapply(colnames(score.data),function(x)score.func(score.data[,x],wt.matrix[,x]*biol.wts[score.genes])))
            names(my.scores)=colnames(score.data)
            object@gene.scores[cell.ids,score.name]=my.scores
            return(object)
          }
)

setGeneric("calcNoiseModels", function(object, cells.use=NULL, trusted.genes=NULL,n.bin=20,drop.expr=1,pop.avg=NULL) standardGeneric("calcNoiseModels"))
setMethod("calcNoiseModels","scR",
          function(object, cells.use=NULL, trusted.genes=NULL,n.bin=20,drop.expr=1, pop.avg=NULL) {
            object@drop.expr=drop.expr
            cells.use=set.ifnull(cells.use,colnames(object@data))
            trusted.genes=set.ifnull(trusted.genes,rownames(object@data))
            trusted.genes=trusted.genes[trusted.genes%in%rownames(object@data)]
            object@trusted.genes=trusted.genes
            data=object@data[trusted.genes,cells.use]
            #pop.avg = set.ifnull(pop.avg, apply(data, 1, function(x) expMean(x)))
            
            pop.avg = set.ifnull(pop.avg, apply(data,1,humpMean,min=object@drop.expr))
            trusted.genes = intersect(trusted.genes, names(pop.avg))
            pop.avg = pop.avg[trusted.genes]
            object@pop.avg = pop.avg
            idents=data.frame(data[trusted.genes,1])
            
            code_humpAvg=pop.avg
            code_humpAvg[code_humpAvg>9]=9
            code_humpAvg[is.na(code_humpAvg)]=0
            idents$code_humpAvg=code_humpAvg
            data[data>=object@drop.expr]=1
            data[data<object@drop.expr]=0
            data$bin=cut(code_humpAvg,n.bin)
            data$avg=code_humpAvg
            rownames(idents)=rownames(data)
            my.coefs=data.frame(t(sapply(colnames(data[1:(ncol(data)-2)]),
                                         getAB,data=data,data2=idents,identus="code",code2="humpAvg",hasBin=TRUE,doPlot=FALSE)))
            colnames(my.coefs)=c("a","b")
            object@drop.coefs = my.coefs
            return(object)
          }
)

setGeneric("fitDetectProb", function(object, cells.use=NULL, cell.name=NULL, trusted.genes=NULL,n.bin=20,drop.expr=1,pop.avg=NULL, return.data=FALSE) standardGeneric("fitDetectProb"))
setMethod("fitDetectProb","scR",
          function(object, cells.use=NULL, cell.name=NULL, trusted.genes=NULL,n.bin=20,drop.expr=1, pop.avg=NULL,return.data=FALSE) {
            object@drop.expr=drop.expr
            cells.use=set.ifnull(cells.use,colnames(object@data))
            trusted.genes=set.ifnull(trusted.genes,rownames(object@data))
            trusted.genes=trusted.genes[trusted.genes%in%rownames(object@data)]
            object@trusted.genes=trusted.genes
            #data=object@data[trusted.genes,cells.use]
            data=object@count.data[trusted.genes,cells.use]
            #pop.avg = set.ifnull(pop.avg, apply(data, 1, function(x) expMean(x)))
            cell.name=set.ifnull(cell.name,"Group")
            
          
            pop.avg = set.ifnull(pop.avg, apply(data,1,humpMean,min=drop.expr))
            trusted.genes = intersect(trusted.genes, names(pop.avg))
            pop.avg = pop.avg[trusted.genes]
            object@pop.avg = pop.avg
            data = data[trusted.genes,]
            
            binaryData = data;
            binaryData = 1*(binaryData > drop.expr)
            
            
            probDetect = apply(binaryData,1,function(x) sum(x>0)/length(x))
            cutLocs=cut2(pop.avg,g=n.bin,onlycuts=TRUE)
            bulkBin=cut2(pop.avg,g=n.bin)
            yBin=tapply(probDetect,bulkBin,median)
            yBinSD = tapply(probDetect,bulkBin,mad)
            xBin=tapply(pop.avg,bulkBin,median)
            #wts = sapply(bulkBin,function(x) cutLocs[as.numeric(x)+1]- cutLocs[as.numeric(x)])
            binaryDataVec = c(as.matrix(binaryData))
            bulkDataVec = rep(pop.avg, ncol(binaryData))
            #wts_vec = rep(wts,ncol(binaryData))
            
            options(warn=-1) 
            data.glm=data.frame(binaryData=binaryDataVec,bulkData=bulkDataVec)
            #glm.out = glm(binaryData~bulkData,family=binomial(logit),data=data.glm,weights=wts_vec)
            glm.out = glm(binaryData~1+bulkData+I(bulkData^2),family=binomial(logit),data=data.glm)
            options(warn=0)  
            x_vals=seq(0,10,0.1)
            y_vals=predict(glm.out,data.frame(bulkData=x_vals),type="response")
            plot(xBin,yBin,pch=16,xlab="Average expression",ylab="Probability of detection",xlim=c(0,8),ylim=c(0,1),main=cell.name)
            lines(x_vals,y_vals,lwd=2)
            
            if (return.data){
              to.return=c();
              to.return$glm.fit = glm.out;
              to.return$xdata = xBin;
              to.return$ydata = yBin;
              to.return$ydataSD = yBinSD;
              return(to.return)
            } else {
              return(object)
              
            }
            
          }
)  






setGeneric("signature.vector.pc", function(object,pc.1=1,pc.2=2,sig.genes = NULL,sig.name=NULL,cells.use=NULL,pt.size=3,do.bare=FALSE,cols.use=NULL,reduction.use="pca", min.perc.gene = 0.001, max.perc.gene=1, seg.color="blue", do.return=FALSE) standardGeneric("signature.vector.pc"))
setMethod("signature.vector.pc", "scR", 
          function(object,pc.1=1,pc.2=2,sig.genes = NULL, sig.name=NULL, cells.use=NULL,pt.size=3,do.bare=FALSE,cols.use=NULL,reduction.use="pca", min.perc.gene = 0.001, max.perc.gene=1,seg.color="blue", do.return=FALSE) {
          
            sig.genes = sig.genes[sig.genes %in% rownames(object@data)]
            if (length(sig.genes)==0){
              stop("Signature is either of length zero, or none of the genes are appreciably expressed. Check gene names")
            }
            sig.name = set.ifnull(sig.name,"Signature")
            if (length(sig.genes) == 1){
              score = as.numeric(object@data[sig.genes,])
              names(score) = colnames(object@data)
            } else {
              score = gene.set.score(object, genes.use=sig.genes,score.name = sig.name, min.perc = min.perc.gene, max.perc=max.perc.gene, return.score=TRUE)
            }
            p = pca.plot(object, pc.1=pc.1, pc.2=pc.2, cells.use=cells.use, pt.size=pt.size, do.return=TRUE, do.bare=do.bare, cols.use=cols.use,
                         reduction.use=reduction.use)
            
            if (reduction.use=="pca") {
              dim.code="PC"
            }
            if (reduction.use=="tsne") {
              dim.code="tsne_"
            }
            if (reduction.use=="ica") {
              dim.code="IC"
            }
            
           
            df = fetch.data(object, vars.all=c(paste0(dim.code, c(pc.1, pc.2))))
            r1 = cor(df[,1], score); r2 = cor(df[,2],score)
           
            xend1 = r1*max(sign(r1)*df[,1]); yend1 = r2*max(sign(r2)*df[,2])
            df1 = data.frame(x1=0, y1=0, xend1=xend1, yend1=yend1)
            #p2 = p + geom_segment(aes(x=0,y=0,xend = r1*xmax, yend=r2*ymax), color=seg.color,arrow = arrow(length = unit(0.5,"cm"))) + annotate("text", x=1.1*xmax*r1, y=1.1*ymax*r2, label=sig.name)
            p2 = p + geom_segment(aes(x=x1,y=y1,xend = xend1, yend=yend1),data=df1, color=seg.color,arrow = arrow(length = unit(0.5,"cm"))) + annotate("text", x=1.1*xend1, y=1.1*yend1, label=sig.name)
            
            if (do.return){
              return(p2)
            } else {
              print(p2)
            }
            
            
          }
)

setGeneric("spatial.de", function(object,marker.cells,genes.use=NULL,...) standardGeneric("spatial.de"))
setMethod("spatial.de", "scR", 
          function(object,marker.cells,genes.use=NULL) {
            object=p15
            embed.map=object@tsne.rot
            mult.use=2
            mult.use.far=10
            if ((mult.use.far*length(marker.cells))>nrow(embed.map)) {
              mult.use.far=1
              mult.use=1
            }
            genes.use=set.ifnull(genes.use,object@var.genes)
            marker.pos=apply(embed.map[marker.cells,],2,mean)
            embed.map=rbind(embed.map,marker.pos)
            rownames(embed.map)[nrow(embed.map)]="marker"
            embed.dist=sort(as.matrix(dist((embed.map)))["marker",])
            embed.diff=names(embed.dist[!(names(embed.dist)%in%marker.cells)][1:(mult.use*length(marker.cells))][-1])
            embed.diff.far=names(embed.dist[!(names(embed.dist)%in%marker.cells)][1:(mult.use.far*length(marker.cells))][-1])
            
            diff.genes=rownames(subset(diffExp.test(p15,marker.cells,embed.diff,genes.use=genes.use),myP<(1e-5)))
            diff.genes=subset(diffExp.test(p15,marker.cells,embed.diff,genes.use = diff.genes),myP<(1e-10))
            return(diff.genes)
          }
)


setGeneric("Mclust_dimension", function(object,pc.1=1,pc.2=2,cells.use=NULL,pt.size=4,reduction.use="tsne",G.use=NULL,set.ident=FALSE,seed.use=1,...) standardGeneric("Mclust_dimension"))
setMethod("Mclust_dimension", "scR", 
          function(object,pc.1=1,pc.2=2,cells.use=NULL,pt.size=4,reduction.use="tsne",G.use=NULL,set.ident=FALSE,seed.use=1,...) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            dim.code="PC"
            if (reduction.use=="pca") data.plot=object@pca.rot[cells.use,]
            if (reduction.use=="tsne") {
              data.plot=object@tsne.rot[cells.use,]
              dim.code="tsne_"
            }
            if (reduction.use=="ica") {
              data.plot=object@ica.rot[cells.use,]
              dim.code="IC"
            }
            x1=paste(dim.code,pc.1,sep=""); x2=paste(dim.code,pc.2,sep="")
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            set.seed(seed.use); data.mclust=ds <- dbscan(data.plot[,c("x","y")],eps = G.use,...)
            
            to.set=as.numeric(data.mclust$cluster+1)
            data.names=names(object@ident)
            object@data.info[data.names,"m"]=to.set
            if (set.ident) {
              object@ident=factor(to.set); names(object@ident)=data.names;               
            }
            
            return(object)
          }
)

setGeneric("prune.clust", function(object,min.cells=10,remove.clust=NULL,id=NULL, attributes.yes = NULL, attributes.no=NULL,statfxn=mean,reorder.clust=TRUE,...) standardGeneric("prune.clust"))
setMethod("prune.clust", "scR", 
          function(object,min.cells=10, remove.clust=NULL,id="m",attributes.yes = NULL, attributes.no=NULL,statfxn=mean,reorder.clust=TRUE,...) {
            
            clust.use = names(table(object@ident))[which(table(object@ident) >= min.cells)]
            if (!(is.null(remove.clust))){ 
              print(paste0("Removing clusters ", paste0(remove.clust, collapse=","), " as specified ."))
              cells.use = which.cells(object, setdiff(clust.use,remove.clust))
            } else {
              print(paste0("Removing clusters ", paste0(setdiff(names(table(object@ident)), clust.use), collapse=","), " as these have fewer than min.cells = ", 
                           min.cells, " cells"))
              cells.use = which.cells(object, clust.use)
            }

              
            
            if (!(is.null(attributes.yes)) | !(is.null(attributes.no))){
              print(paste0("Using id = ", id.use, " for pruning"))
              ident.use = object@data.info[,id]
              clusters_to_remove = c()
              df = fetch.data(object, vars.all = c(attributes.yes, attributes.no))
              df_agg = aggregate(df, by=list(ident.use), statfxn)
              attributes.yes = intersect(colnames(df), attributes.yes)
              if (!(is.null(attributes.yes))){
                for (att in attributes.yes){
                clust = rownames(df_agg)[which(df_agg[,att] <= mean(df_agg[,att]) - 1 * sd(df_agg[,att])) ]
                print(paste0("Removing ", paste0(clust, collapse = ","), " based on low ", att))
                clusters_to_remove = union(clusters_to_remove, clust)
                }
              }
              
              if (!(is.null(attributes.no))){
                attributes.no = intersect(colnames(df), attributes.no)
                for (att in attributes.no){
                  clust = rownames(df_agg)[which(df_agg[,att] >= mean(df_agg[,att]) + 1 * sd(df_agg[,att])) ]
                  print(paste0("Removing ", paste0(clust, collapse = ","), " based on high ", att))
                  clusters_to_remove = union(clusters_to_remove, clust)
                }
              }
              
              cells.remove = which.cells(object, clusters_to_remove)
              cells.use = setdiff(cells.use, cells.remove)
            }
            print(paste0("Removing ", length(setdiff(object@cell.names, cells.use)), " cells from the data"))
            
            new.object = subsetData(object, cells.use=cells.use)
            if (reorder.clust){
              levels(new.object@ident) = 1:length(levels(new.object@ident))
            }
            new.object@ident = factor(new.object@ident)
            new.object@data.info[,"m_prune"] = new.object@ident
            
            return(new.object)
          }
)

setGeneric("force.merge", function(object,merge.clust=NULL,new.name=NULL,char=FALSE,set.data.info=NULL,...) standardGeneric("force.merge"))
setMethod("force.merge", "scR", 
          function(object, merge.clust=NULL,new.name=NULL,char=FALSE,set.data.info=NULL,...) {
            
            if (is.null(merge.clust)) return(object)
            if (is.null(new.name)) new.name = merge.clust[1]
            if (!char){
              set.data.info=set.ifnull(set.data.info,"m")
              ident.use=object@ident
              ident.levels = levels(object@ident)
              ident.use[ident.use %in% merge.clust] = new.name
              ident.use = drop.levels(ident.use)
              #ident.levels = ident.levels[-which(ident.levels==merge.clust[2])]
              #ident.use = factor(ident.use,levels = as.numeric(sort(ident.levels)))
              #levels(ident.use) = ident.levels
              #if ((min(as.numeric(ident.use))) == 0){
              #  levels(ident.use) = c(0:(length(levels(ident.use))-1))
              #} else {
              #  levels(ident.use) = c(1:length(levels(ident.use)))
              #}
              object@ident = ident.use
              object@data.info[,set.data.info] = ident.use
            }
            
            if (char){
              set.data.info=set.ifnull(set.data.info,"orig")
              ident.use= as.character(object@ident)
              names(ident.use) = names(object@ident)
              ident.use[ident.use %in% merge.clust] = new.name
              ident.use = factor(ident.use)
              object@data.info[,set.data.info] = ident.use
            }
          
            return(object)
          }
)

setGeneric("pca.sig.genes", function(object,pcs.use,pval.cut=0.1,use.full=TRUE) standardGeneric("pca.sig.genes"))
setMethod("pca.sig.genes", "scR", 
          function(object,pcs.use,pval.cut=0.1,use.full=TRUE) {
            pvals.use=object@jackStraw.empP
            pcx.use=object@pca.x
            if (use.full)  {
              pvals.use=object@jackStraw.empP.full
              pcx.use=object@pca.x.full
            }
            if (length(pcs.use)==1) pvals.min=pvals.use[,pcs.use]
            if (length(pcs.use)>1) pvals.min=apply(pvals.use[,pcs.use],1,min)
            names(pvals.min)=rownames(pvals.use)
            genes.use=names(pvals.min)[pvals.min<pval.cut]
            genes.use=genes.use[genes.use%in%rownames(object@scale.data)]
            return(genes.use)       
          }
)

setGeneric("pca.sig.genes.by.dist", function(object,pcs.use,sd.thres=5,use.full=TRUE) standardGeneric("pca.sig.genes.by.dist"))
setMethod("pca.sig.genes.by.dist", "scR", 
          function(object,pcs.use,sd.thres=2,use.full=TRUE) {
            pvals.use=object@jackStraw.empP
            pcx.use=object@pca.x
            if (use.full)  {
              pvals.use=object@jackStraw.empP.full
              pcx.use=object@pca.x.full
            }
            
            fakePC = object@jackStraw.fakePC
            #Compute mean and sd of fake PCs
            mean.fakePC = apply(fakePC, 2, mean)
            sd.fakePC = apply(fakePC, 2, sd)
            
            zscore.max = apply(pcx.use[,pcs.use], 1, function(x) max (abs(x-mean.fakePC[pcs.use]) / sd.fakePC[pcs.use]))
            genes.use = names(zscore.max)[ zscore.max > sd.thres]
            genes.use=genes.use[genes.use%in%rownames(object@scale.data)]
            return(genes.use)       
          }
)



same=function(x) return(x)




setGeneric("doKMeans", function(object,pcs.use=1,pval.cut=0.1,k.num=NULL,k.seed=1,do.plot=TRUE,clust.cut=2.5,disp.cut=2.5,k.cols=pyCols,do.one=FALSE,do.k.col=FALSE,
                                k.col=NULL,pc.row.order=NULL,pc.col.order=NULL, rev.pc.order=FALSE, cluster.zoom=0, use.full=FALSE,clust.col=TRUE,do.annot=FALSE,
                                only.k.annot=FALSE,do.recalc=TRUE,use.imputed=FALSE,col.annot.show=NULL,genes.use=NULL,print.genes=FALSE) standardGeneric("doKMeans"))

setMethod("doKMeans","scR",
          function(object,pcs.use=1,pval.cut=0.1,k.num=NULL,k.seed=1,do.plot=TRUE,clust.cut=2.5,disp.cut=2.5,k.cols=pyCols,do.one=FALSE,do.k.col=FALSE,
                   k.col=NULL,pc.row.order=NULL,pc.col.order=NULL,rev.pc.order=FALSE, cluster.zoom=0,use.full=FALSE,clust.col=TRUE,do.annot=FALSE,
                   only.k.annot=FALSE,do.recalc=TRUE,use.imputed=FALSE,col.annot.show=NULL,genes.use=NULL,print.genes=FALSE) {
            require(gplots)
            require(NMF)
            
            data.use.orig=object@scale.data
            if (use.imputed) data.use.orig=data.frame(t(scale(t(object@imputed))))
            data.use=minmax(data.use.orig,min=disp.cut*(-1),max=clust.cut)
            pc.row.order=set.ifnull(pc.row.order,pcs.use[1])
            pc.col.order=set.ifnull(pc.col.order,pcs.use[1])
            
            pvals.use=object@jackStraw.empP
            pcx.use=object@pca.x
            if (use.full)  {
              pvals.use=object@jackStraw.empP.full
              pcx.use=object@pca.x.full
            }
            revFxn=same; if (rev.pc.order) revFxn=function(x)max(x)+1-x;
            k.num=set.ifnull(k.num,object@.best)
            #if (length(pcs.use)==1) pvals.min=pvals.use[,pcs.use]
            #if (length(pcs.use)>1) pvals.min=apply(pvals.use[,pcs.use],1,min)
            #names(pvals.min)=rownames(pvals.use)
            #genes.use=set.ifnull(genes.use,names(pvals.min)[pvals.min<pval.cut])
            genes.use=set.ifnull(genes.use,pca.sig.genes(object,pcs.use,pval.cut,use.full))
            
            genes.use=genes.use[genes.use%in%rownames(data.use)]
            cells.use=colnames(data.use)
            kmeans.data=data.use[genes.use,cells.use]      
            if (print.genes) print(genes.use)
            
            k.col=set.ifnull(k.col, k.num)
            if (do.recalc) {
              set.seed(k.seed); kmeans.obj=kmeans(kmeans.data,k.num); kmeans.col=kmeans(t(kmeans.data),k.col)
              kmeans.obj$cluster=as.numeric(revFxn(rank(tapply(pcx.use[genes.use,pc.row.order],as.numeric(kmeans.obj$cluster),mean)))[as.numeric(kmeans.obj$cluster)])
              names(kmeans.obj$cluster)=genes.use
              kmeans.col$cluster=as.numeric(revFxn(rank(tapply(object@pca.rot[cells.use,pc.col.order],kmeans.col$cluster,mean)))[as.numeric(kmeans.col$cluster)])
              names(kmeans.col$cluster)=cells.use
              object@kmeans.obj=list(kmeans.obj)
              object@kmeans.col=list(kmeans.col)
            }    
            kmeans.obj=object@kmeans.obj[[1]]
            kmeans.col=object@kmeans.col[[1]]
            
            if (do.plot) {       
              disp.data=minmax(kmeans.data[order(kmeans.obj$cluster[genes.use]),],min=disp.cut*(-1),max=disp.cut)
              if (do.one)  {
                disp.data=disp.data[order(pcx.use[rownames(disp.data),pcs.use[1]]),order(object@pca.rot[colnames(disp.data),pcs.use[1]])]
              }
              if (do.k.col) {
                disp.data=disp.data[,order(kmeans.col$cluster)]
              }
              if (cluster.zoom[1] != 0) {
                genes.use=names(kmeans.obj$cluster)[which(kmeans.obj$cluster%in%cluster.zoom)]
                disp.data=disp.data[genes.use[order(kmeans.obj$cluster[genes.use])],]
              }  
              row.annot=data.frame(cbind(kmeans.obj$cluster[rownames(disp.data)],pcx.use[rownames(disp.data),pcs.use]))
              colnames(row.annot)=c("K",paste("PC",pcs.use,sep=""))
              k.names=colnames(disp.data)
              pcrot.info=object@pca.rot[k.names,pcs.use]
              pcrot.info=log(pcrot.info+min(pcrot.info)+1)
              col.annot=data.frame(cbind(object@data.ngene[k.names],pcrot.info,object@ident[k.names]))
              #col.annot.data=data.frame(cbind(object@data.ngene[k.names],pcrot.info,object@ident[k.names]),t(object@data.metrics),t(object@gene.scores))
              for(i in c(1:(ncol(col.annot)-1))) {
                col.annot[,i]=as.numeric(as.character(col.annot[,i]))
              }     
              colnames(col.annot)=c("nGene",paste("cPC",pcs.use,sep=""),"ident")     
              do.Colv="TRUE"
              colnames(disp.data)=sub.string(colnames(disp.data),"_rsem","")
              if (!do.annot) {
                col.annot=NULL
                row.annot=NULL
              }
              if (do.k.col) {
                #col.annot=cbind(col.annot,kmeans.col$cluster[colnames(disp.data)])
                #colnames(col.annot)[ncol(col.annot)]="K"
                #print(head(kmeans.col$cluster))
              }
              if (only.k.annot) {
                row.annot=minusc(row.annot,"PC")
                col.annot=minusc(col.annot,"PC")
              }
              if (do.one)  aheatmap(disp.data,annRow = row.annot,annCol=col.annot,Rowv=NA,Colv=NA,col=k.cols)
              if (!(do.one)&& !(do.k.col)) aheatmap(disp.data,annRow = row.annot,annCol=col.annot,Rowv=NA,Colv=clust.col,col=k.cols)
              if (!(do.one)&& (do.k.col)) aheatmap(disp.data,annRow = row.annot,annCol=col.annot,Rowv=NA,Colv=NA,col=k.cols)
              
            }
            return(object)
          }
)



setGeneric("genes.in.cluster", function(object, cluster.num)  standardGeneric("genes.in.cluster"))
setMethod("genes.in.cluster", signature = "scR",
          function(object, cluster.num) {
            print(unlist(lapply(cluster.num,function(x)sort(names(which(object@kmeans.obj[[1]]$cluster==x))))))
          }    
)

setGeneric("cells.in.cluster", function(object, cluster.num)  standardGeneric("cells.in.cluster"))
setMethod("cells.in.cluster", signature = "scR",
          function(object, cluster.num) {
            #            print(sort(names(which(object@kmeans.col[[1]]$cluster==cluster.num))))
            return(sort(unlist(lapply(cluster.num,function(x) names(which(object@kmeans.col[[1]]$cluster==x))))))
          }    
)

setGeneric("cell.cor.matrix", function(object, cor.genes=NULL,cell.inds=NULL, do.k=FALSE,k.seed=1,k.num=4,vis.low=(-1),vis.high=1,vis.one=0.8,pcs.use=1:3,col.use=pyCols)  standardGeneric("cell.cor.matrix"))
setMethod("cell.cor.matrix", signature = "scR",
          function(object, cor.genes=NULL,cell.inds=NULL, do.k=FALSE,k.seed=1,k.num=4,vis.low=(-1),vis.high=1,vis.one=0.8,pcs.use=1:3,col.use=pyCols) {
            cor.genes=set.ifnull(cor.genes,object@var.genes)
            cell.inds=set.ifnull(cell.inds,colnames(object@data))
            cor.genes=cor.genes[cor.genes%in%rownames(object@data)]
            data.cor=object@data[cor.genes,cell.inds]
            cor.matrix=cor((data.cor))
            set.seed(k.seed); kmeans.cor=kmeans(cor.matrix,k.num)
            cor.matrix=cor.matrix[order(kmeans.cor$cluster),order(kmeans.cor$cluster)]
            kmeans.names=rownames(cor.matrix)
            row.annot=data.frame(cbind(kmeans.cor$cluster[kmeans.names],object@pca.rot[kmeans.names,pcs.use]))
            colnames(row.annot)=c("K",paste("PC",pcs.use,sep=""))
            cor.matrix[cor.matrix==1]=vis.one
            cor.matrix=minmax(cor.matrix,min = vis.low,max=vis.high)
            object@kmeans.cell=list(kmeans.cor)
            if (do.k) aheatmap(cor.matrix,col=col.use,Rowv=NA,Colv=NA,annRow=row.annot)
            if (!(do.k)) aheatmap(cor.matrix,col=col.use,annRow=row.annot)
            return(object)
          }    
)

setGeneric("gene.cor.matrix", function(object, cor.genes=NULL,cell.inds=NULL, do.k=FALSE,k.seed=1,k.num=4,vis.low=(-1),vis.high=1,vis.one=0.8,pcs.use=1:3,col.use=pyCols)  standardGeneric("gene.cor.matrix"))
setMethod("gene.cor.matrix", signature = "scR",
          function(object, cor.genes=NULL,cell.inds=NULL, do.k=FALSE,k.seed=1,k.num=4,vis.low=(-1),vis.high=1,vis.one=0.8,pcs.use=1:3,col.use=pyCols) {
            cor.genes=set.ifnull(cor.genes,object@var.genes)
            cell.inds=set.ifnull(cell.inds,colnames(object@data))
            cor.genes=cor.genes[cor.genes%in%rownames(object@data)]
            data.cor=object@data[cor.genes,cell.inds]
            cor.matrix=cor(t(data.cor))
            set.seed(k.seed); kmeans.cor=kmeans(cor.matrix,k.num)
            cor.matrix=cor.matrix[order(kmeans.cor$cluster),order(kmeans.cor$cluster)]
            kmeans.names=rownames(cor.matrix)
            row.annot=data.frame(cbind(kmeans.cor$cluster[kmeans.names],object@pca.x[kmeans.names,pcs.use]))
            colnames(row.annot)=c("K",paste("PC",pcs.use,sep=""))
            cor.matrix[cor.matrix==1]=vis.one
            cor.matrix=minmax(cor.matrix,min = vis.low,max=vis.high)
            object@kmeans.gene=list(kmeans.cor)
            if (do.k) aheatmap(cor.matrix,col=col.use,Rowv=NA,Colv=NA,annRow=row.annot)
            if (!(do.k)) aheatmap(cor.matrix,col=col.use,annRow=row.annot)
            return(object)
          }    
)

setGeneric("calinskiPlot", function(object,pcs.use,pval.cut=0.1,gene.max=15,col.max=25,use.full=TRUE)  standardGeneric("calinskiPlot"))
setMethod("calinskiPlot","scR",
          function(object,pcs.use,pval.cut=0.1,gene.max=15,col.max=25,use.full=TRUE) {
            require(vegan)
            if (length(pcs.use)==1) pvals.min=object@jackStraw.empP.full[,pcs.use]
            if (length(pcs.use)>1) pvals.min=apply(object@jackStraw.empP.full[,pcs.use],1,min)
            names(pvals.min)=rownames(object@jackStraw.empP.full)
            genes.use=names(pvals.min)[pvals.min<pval.cut]
            genes.use=genes.use[genes.use%in%rownames(object@scale.data)]
            
            par(mfrow=c(1,2))
            
            mydata <- object@scale.data[genes.use,]
            wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
            for (i in 1:gene.max) wss[i] <- sum(kmeans(mydata,
                                                       centers=i)$withinss)
            plot(1:gene.max, wss, type="b", xlab="Number of Clusters for Genes",
                 ylab="Within groups sum of squares")
            
            
            mydata <- t(object@scale.data[genes.use,])
            wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
            for (i in 1:col.max) wss[i] <- sum(kmeans(mydata,
                                                      centers=i)$withinss)
            plot(1:col.max, wss, type="b", xlab="Number of Clusters for Cells",
                 ylab="Within groups sum of squares")
            rp()
            return(object)
          }
)

setMethod("show", "scR",
          function(object) {
            cat("An object of class ", class(object), " in project ", object@project.name, "\n", sep = "")
            cat(" ", nrow(object@data), " features by ",
                ncol(object@data), " samples.\n", sep = "")
            invisible(NULL)
          }
)



setGeneric("dot.plot", function(object,genes.plot,cex.use=2,cols.use=NULL,thresh.col=2.5,dot.min=0.05)  standardGeneric("dot.plot"))
setMethod("dot.plot","scR",
          function(object,genes.plot,cex.use=2,cols.use=NULL,thresh.col=2.5,dot.min=0.05) {
            genes.plot=ainb(genes.plot,rownames(object@data))
            object@data=object@data[genes.plot,]
            avg.exp=average.expression(object)
            avg.alpha=cluster.alpha(object)
            cols.use=set.ifnull(cols.use,myPalette(low = "red",high="green"))
            exp.scale=t(scale(t(avg.exp)))
            exp.scale=minmax(exp.scale,max=thresh.col,min=(-1)*thresh.col)
            n.col=length(cols.use)
            data.y=rep(1:ncol(avg.exp),nrow(avg.exp))
            data.x=unlist(lapply(1:nrow(avg.exp),rep,ncol(avg.exp)))
            data.avg=unlist(lapply(1:length(data.y),function(x) exp.scale[data.x[x],data.y[x]]))
            exp.col=cols.use[floor(n.col*(data.avg+thresh.col)/(2*thresh.col)+.5)]
            data.cex=unlist(lapply(1:length(data.y),function(x) avg.alpha[data.x[x],data.y[x]]))*cex.use+dot.min
            plot(data.x,data.y,cex=data.cex,pch=16,col=exp.col,xaxt="n",xlab="",ylab="",yaxt="n")
            axis(1,at = 1:length(genes.plot),genes.plot)
            axis(2,at=1:ncol(avg.alpha),colnames(avg.alpha),las=1)
          }
          
) 

setGeneric("gene.dist.plot", function(object,ident.use=NULL,features.plot,nCol=NULL,ylab.max=12,do.ret=FALSE,do.sort=FALSE,size.x.use=16,size.y.use=16,size.title.use=20, use.imputed=FALSE,adjust.use=1,size.use=1,cols.use=NULL,
                                do.mean=FALSE,use.raw=FALSE,vln.vert=FALSE,...)  standardGeneric("gene.dist.plot"))
setMethod("gene.dist.plot","scR",
          function(object,ident.use=NULL,features.plot,nCol=NULL,ylab.max=12,do.ret=FALSE,do.sort=FALSE,size.x.use=16,size.y.use=16,size.title.use=20,use.imputed=FALSE,adjust.use=1,
                   size.use=1,cols.use=NULL,do.mean=FALSE,use.raw=FALSE,...) {
            if (is.null(nCol)) {
              nCol=2
              if (length(features.plot)>6) nCol=3
              if (length(features.plot)>9) nCol=4
            }
            
            if (is.null(ident.use)){
              cells.use=NULL
            } else {
              cells.use = which.cells(object,ident.use)
            }
            
            data.use=data.frame(t(fetch.data(object,features.plot,use.imputed=use.imputed, use.raw=use.raw, cells.use=cells.use)))
            ident.use=object@ident[cells.use]
            pList=lapply(features.plot,function(x) plot.dist(x,data.use[x,],ident.use,ylab.max,TRUE,do.sort,size.x.use,size.y.use,size.title.use,adjust.use,size.use,cols.use))
            if(do.ret) {
              return(pList)
            }
            else {
              multiplotList(pList,cols = nCol)
              rp()
            }
          }
)

plot.dist=function(gene,data,cell.ident,ylab.max=12,do.ret=FALSE,do.sort=FALSE,size.x.use=16,size.y.use=16,size.title.use=20,
                   adjust.use=1,size.use=1,cols.use=NULL) {
  data$gene=as.character(rownames(data))
  data.use=data.frame(data[gene,])
  if (length(gene)==1) {
    data.melt=data.frame(rep(gene,length(cell.ident))); colnames(data.melt)[1]="gene"
    data.melt$value=as.numeric(data[1,1:length(cell.ident)])
    data.melt$id=names(data)[1:length(cell.ident)]
  }
  #print(head(data.melt))
  
  if (length(gene)>1) data.melt=melt(data.use,id="gene")
  data.melt$ident=cell.ident
  
  data.melt$value=as.numeric(as.character(data.melt$value))
  if(do.sort) {
    data.melt$ident=factor(data.melt$ident,levels=names(rev(sort(tapply(data.melt$value,data.melt$ident,mean)))))
  }
  p=ggplot(data.melt,aes(value))
  p2=p +  xlab("Expression level (log TPM)")
  if(length(unique(cell.ident)) == 1){
    p2 = p2 + geom_histogram(binwidth=0.1,color="black",fill="indianred4")
  } else {
    p2 = p2 + geom_density(aes(color=factor(ident)))
  }
  if (!is.null(cols.use)) {
    p2=p2+scale_fill_manual(values=cols.use)
  }
  p3=p2+theme(legend.position="top")+guides(fill=guide_legend(title=NULL))+ylab("Density")
  p4=p3+ theme(axis.title.x = element_text(face="bold", size=size.x.use), axis.text.x  = element_text(angle=90, vjust=0.5, size=14))+theme_bw()+nogrid
  p5=p4+theme(axis.title.y = element_text(face="bold", size=size.y.use), axis.text.y  = element_text(angle=90, vjust=0.5, size=14))+ggtitle(gene)+theme(plot.title = element_text(size=size.title.use, face="bold"))
  if(do.ret==TRUE) {
    return(p5)
  }
  else {
    print(p5)
  }
}

setGeneric("jackStrawPlot", function(object,plot.lim=0.8,num.pc=5,score.thresh=1e-5,...)  standardGeneric("jackStrawPlot"))
setMethod("jackStrawPlot","scR",
          function(object,plot.lim=0.8,num.pc=5,score.thresh=0.01,...) {
            pAll=object@jackStraw.empP
            nCol=2
            if (num.pc>6) nCol=3
            if (num.pc>9) nCol=4
            num.row=floor(num.pc/nCol-1e-5)+1
            par(mfrow=c(nCol,num.row))
            qq.x=sapply(1:num.pc,function(x) {qqplot(pAll[,x],runif(1000),plot.it=FALSE)$x})
            qq.y=sapply(1:num.pc,function(x) {qqplot(pAll[,x],runif(1000),plot.it=FALSE)$y})
            pc.score=unlist(lapply(1:num.pc, function(x) mean(qq.y[which(qq.x[,x]<=score.thresh),x])))
            pc.score=unlist(lapply(1:num.pc, function(x) prop.test(c(length(which(pAll[,x]<=score.thresh)),floor(nrow(pAll)*score.thresh)),c(nrow(pAll),nrow(pAll)))$p.val))
            
            unlist(lapply(1:num.pc,function(x) {qqplot(pAll[,x],runif(1000),,xlim=c(0,plot.lim),ylim=c(0,plot.lim),xlab=colnames(pAll)[x],pch=16,main=round(pc.score[x],4)); 
                                                lines(seq(0,1,0.01),seq(0,1,0.01),lty=2,lwd=2)}))
            #            unlist(lapply(1:num.pc,function(x) {qqplot(pAll[,x],runif(1000),,xlim=c(0,plot.lim),ylim=c(0,plot.lim),xlab=colnames(pAll)[x],pch=16); 
            #                                                lines(seq(0,1,0.01),seq(0,1,0.01),lty=2,lwd=2)}))
            rp()
          }
)

setGeneric("addMetaData", function(object,metadata)  standardGeneric("addMetaData"))
setMethod("addMetaData","scR",
          function(object,metadata) {
            cols.add=colnames(metadata)
            object@data.info[,cols.add]=metadata[rownames(object@data.info),]
            return(object)
          }
)

setGeneric("jackStrawPlot2", function(object,plot.lim=0.4,num.pc=5,...)  standardGeneric("jackStrawPlot2"))
setMethod("jackStrawPlot2","scR",
          function(object,plot.lim=0.4,num.pc=5,...) {
            pAll=object@emp.pval
            num.pc=ncol(pAll)
            num.row=floor(num.pc/2-1e-5)+1
            par(mfrow=c(2,num.row))
            unlist(lapply(1:num.pc,function(x) {qqplot(pAll[,x],runif(1000),,xlim=c(0,plot.lim),ylim=c(0,plot.lim),xlab=colnames(pAll)[x],pch=16); 
                                                lines(seq(0,1,0.01),seq(0,1,0.01),lty=2,lwd=2)}))
            rp()
          }
)


setGeneric("geneScorePlot", function(object, gene1, score.name, cell.ids=NULL,col.use=NULL,nrpoints.use=Inf,pch.use=16,cex.use=2,...)  standardGeneric("geneScorePlot"))
setMethod("geneScorePlot","scR",
          function(object, gene1, score.name, cell.ids=NULL,col.use=NULL,nrpoints.use=Inf,pch.use=16,cex.use=2,...) {
            cell.ids=set.ifnull(cell.ids,colnames(object@data))
            g1=as.numeric(object@data[gene1,cell.ids])
            my.score=retreiveScore(object,score.name)
            s1=as.numeric(my.score[cell.ids])
            col.use=set.ifnull(as.numeric(as.factor(object@ident[cell.ids])))
            gene.cor=round(cor(g1,s1),2)
            smoothScatter(g1,s1,xlab=gene1,ylab=score.name,col=col.use,nrpoints=nrpoints.use,cex=cex.use,main=gene.cor,pch=pch.use)
          }
)



setGeneric("jackStraw", function(object,num.pc=8,num.replicate=100,prop.freq=0.01,do.print=FALSE)  standardGeneric("jackStraw"))
setMethod("jackStraw","scR",
          function(object,num.pc=5,num.replicate=100,prop.freq=0.01,do.print=FALSE) {
            pc.genes=rownames(object@pca.x)
            if (length(pc.genes)<200) prop.freq=max(prop.freq,0.015)
            md.x=as.matrix(object@pca.x)
            md.rot=as.matrix(object@pca.rot)
            if (!(do.print)) fake.pcVals.raw=sapply(1:num.replicate,function(x)jackRandom(scaled.data=object@scale.data[pc.genes,],prop=prop.freq,r1.use = 1,r2.use = num.pc,seed.use=x),simplify = FALSE)
            if ((do.print)) fake.pcVals.raw=sapply(1:num.replicate,function(x){ print(x); jackRandom(scaled.data=object@scale.data[pc.genes,],prop=prop.freq,r1.use = 1,r2.use = num.pc,seed.use=x)},simplify = FALSE)
            
            fake.pcVals=sapply(1:num.pc,function(x)as.numeric(unlist(lapply(1:num.replicate,function(y)fake.pcVals.raw[[y]][,x]))))
            object@jackStraw.fakePC = data.frame(fake.pcVals)
            object@jackStraw.empP=data.frame(sapply(1:num.pc,function(x)unlist(lapply(abs(md.x[,x]),empP,abs(fake.pcVals[,x])))))
            colnames(object@jackStraw.empP)=paste("PC",1:ncol(object@jackStraw.empP),sep="")
            return(object)
          }
)

jackRandom=function(scaled.data,prop.use=0.01,r1.use=1,r2.use=5, seed.use=1) {
  set.seed(seed.use); rand.genes=sample(rownames(scaled.data),nrow(scaled.data)*prop.use)
  data.mod=scaled.data
  data.mod[rand.genes,]=shuffleMatRow(scaled.data[rand.genes,])
  fake.pca=prcomp(data.mod)
  fake.x=fake.pca$x
  fake.rot=fake.pca$rotation
  return(fake.x[rand.genes,r1.use:r2.use])
}

setGeneric("jackStrawFull", function(object,num.pc=5,num.replicate=100,prop.freq=0.01)  standardGeneric("jackStrawFull"))
setMethod("jackStrawFull","scR",
          function(object,num.pc=5,num.replicate=100,prop.freq=0.01) {
            pc.genes=rownames(object@pca.x)
            if (length(pc.genes)<200) prop.freq=max(prop.freq,0.015)
            md.x=as.matrix(object@pca.x)
            md.rot=as.matrix(object@pca.rot)
            real.fval=sapply(1:num.pc,function(x)unlist(lapply(pc.genes,jackF,r1=x,r2=x,md.x,md.rot)))
            rownames(real.fval)=pc.genes
            object@real.fval=data.frame(real.fval)
            
            fake.fval=sapply(1:num.pc,function(x)unlist(replicate(num.replicate,
                                                                  jackStrawF(prop=prop.freq,data=object@scale.data[pc.genes,],myR1 = x,myR2 = x),simplify=FALSE)))
            rownames(fake.fval)=1:nrow(fake.fval)
            object@fake.fval=data.frame(fake.fval)
            
            object@emp.pval=data.frame(sapply(1:num.pc,function(x)unlist(lapply(object@real.fval[,x],empP,object@fake.fval[,x]))))
            
            rownames(object@emp.pval)=pc.genes
            colnames(object@emp.pval)=paste("PC",1:ncol(object@emp.pval),sep="")
            return(object)
          }
)

setGeneric("mean.var.plot", function(object, cells.use=NULL, genes.use = NULL, fxn.x=humpMean, fxn.y=sd,do.plot=TRUE,set.var.genes=TRUE,do.text=TRUE,
                                     x.low.cutoff=4,x.high.cutoff=8,y.cutoff=2,y.high.cutoff=12,cex.use=0.5,cex.text.use=0.5,do.spike=FALSE, 
                                     pch.use=16, col.use="black", spike.col.use="red",use.imputed=FALSE,plot.both=FALSE,do.contour=TRUE,
                                     contour.lwd=3, contour.col="white", contour.lty=2,num.bin=20, do.ident=FALSE) standardGeneric("mean.var.plot"))
setMethod("mean.var.plot", signature = "scR",
          function(object,cells.use=NULL, genes.use=NULL, fxn.x=humpMean, fxn.y=sd,do.plot=TRUE,set.var.genes=TRUE,do.text=TRUE,
                   x.low.cutoff=4,x.high.cutoff=8,y.cutoff=1,y.high.cutoff=12,cex.use=0.5,cex.text.use=0.5,do.spike=FALSE, 
                   pch.use=16, col.use="black", spike.col.use="red",use.imputed=FALSE,plot.both=FALSE,do.contour=TRUE,
                   contour.lwd=3, contour.col="white", contour.lty=2,num.bin=20, do.ident=FALSE) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            genes.use=set.ifnull(genes.use, rownames(object@data))
            data=object@data[genes.use,cells.use]
            print(dim(data))
            data.x=apply(data,1,fxn.x); data.y=apply(data,1,fxn.y)
            data.norm.y=meanNormFunction(data,fxn.x,fxn.y,num.bin)
            data.norm.y[is.na(data.norm.y)]=0
            names(data.norm.y)=names(data.x)
            pass.cutoff=names(data.x)[which(((data.x>x.low.cutoff) & (data.x<x.high.cutoff)) & (data.norm.y>y.cutoff) & (data.norm.y < y.high.cutoff))]
            mv.df=data.frame(data.x,data.y,data.norm.y)
            rownames(mv.df)=rownames(data)
            object@mean.var=mv.df
            if (do.spike) spike.genes=rownames(subr(data,"^ERCC"))
            if (do.plot) {
              if (plot.both) {
                par(mfrow=c(1,2)) 
                smoothScatter(data.x,sqrt(data.y),pch=pch.use,cex=cex.use,col=col.use,xlab="Average expression",ylab="SD (expression)",nrpoints=Inf)
                
                if (do.contour) {
                  data.kde=kde2d(data.x,data.y)
                  contour(data.kde,add=TRUE,lwd=contour.lwd,col=contour.col,lty=contour.lty)
                }
                if (do.spike) points(data.x[spike.genes],data.y[spike.genes],pch=16,cex=cex.use,col=spike.col.use)
                if(do.text) text(data.x[pass.cutoff],sqrt(data.y)[pass.cutoff],pass.cutoff,cex=cex.text.use)
              }
              smoothScatter(data.x,data.norm.y,pch=pch.use,cex=cex.use,col=col.use,xlab="Average expression",ylab="Dispersion",nrpoints=Inf)
              if (do.contour) {
                data.kde=kde2d(data.x,data.norm.y)
                contour(data.kde,add=TRUE,lwd=contour.lwd,col=contour.col,lty=contour.lty)
              }
              if (do.spike) points(data.x[spike.genes],data.norm.y[spike.genes],pch=16,cex=cex.use,col=spike.col.use,nrpoints=Inf)
              if(do.text) text(data.x[pass.cutoff],data.norm.y[pass.cutoff],pass.cutoff,cex=cex.text.use)
              
              if (do.ident) {
                identify(data.x,data.y,labels = names(data.x))
              }
            }
            if (set.var.genes) { 
              object@var.genes=pass.cutoff
              return(object)
            }
            if (!set.var.genes) return(pass.cutoff)
          }
          
)



#' Identify variable genes
#'
#' Identifies genes that are outliers on a 'mean variability plot'. First, uses
#' a function to calculate average expression (fxn.x) and dispersion (fxn.y)
#' for each gene. Next, divides genes into num.bin (deafult 20) bins based on
#' their average expression, and calculates z-scores for dispersion within each
#' bin. The purpose of this is to identify variable genes while controlling for
#' the strong relationship between variability and average expression.
#'
#' Exact parameter settings may vary empirically from dataset to dataset, and
#' based on visual inspection of the plot.
#' Setting the y.cutoff parameter to 2 identifies genes that are more than two standard
#' deviations away from the average dispersion within a bin. The default X-axis function
#' is the mean expression level, and for Y-axis it is the log(Variance/mean). All mean/variance
#' calculations are not performed in log-space, but the results are reported in log-space -
#' see relevant functions for exact details.
#'
#' @param object Seurat object
#' @param fxn.x Function to compute x-axis value (average expression). Default
#' is to take the mean of the detected (i.e. non-zero) values
#' @param fxn.y Function to compute y-axis value (dispersion). Default is to
#' take the standard deviation of all values/
#' @param do.plot Plot the average/dispersion relationship
#' @param set.var.genes Set object@@var.genes to the identified variable genes
#' (default is TRUE)
#' @param do.text Add text names of variable genes to plot (default is TRUE)
#' @param x.low.cutoff Bottom cutoff on x-axis for identifying variable genes
#' @param x.high.cutoff Top cutoff on x-axis for identifying variable genes
#' @param y.cutoff Bottom cutoff on y-axis for identifying variable genes
#' @param y.high.cutoff Top cutoff on y-axis for identifying variable genes
#' @param cex.use Point size
#' @param cex.text.use Text size
#' @param do.spike FALSE by default. If TRUE, color all genes starting with ^ERCC a different color
#' @param pch.use Pch value for points
#' @param col.use Color to use
#' @param spike.col.use if do.spike, color for spike-in genes
#' @param plot.both Plot both the scaled and non-scaled graphs.
#' @param do.contour Draw contour lines calculated based on all genes
#' @param contour.lwd Contour line width
#' @param contour.col Contour line color
#' @param contour.lty Contour line type
#' @param num.bin Total number of bins to use in the scaled analysis (default
#' is 20)
#' @param do.recalc TRUE by default. If FALSE, plots and selects variable genes without recalculating statistics for each gene.
#' @importFrom MASS kde2d
#' @return Returns a Seurat object, placing variable genes in object@@var.genes.
#' The result of all analysis is stored in object@@mean.var
#' @export
setGeneric("MeanVarPlot", function(object, genes.use = NULL, fxn.x=expMean, fxn.y=logVarDivMean,do.plot=TRUE,set.var.genes=TRUE,do.text=TRUE,
                                   x.low.cutoff=0.1,x.high.cutoff=8,y.cutoff=2,y.high.cutoff=Inf,cex.use=0.5,cex.text.use=0.5,do.spike=FALSE,
                                   pch.use=16, col.use="black", spike.col.use="red",plot.both=FALSE,do.contour=TRUE,
                                   contour.lwd=3, contour.col="white", contour.lty=2,num.bin=20,do.recalc=TRUE) standardGeneric("MeanVarPlot"))
#' @export
setMethod("MeanVarPlot", signature = "scR",
          function(object, genes.use = NULL, fxn.x=expMean, fxn.y=logVarDivMean, do.plot=TRUE,set.var.genes=TRUE,do.text=TRUE,
                   x.low.cutoff=0.1,x.high.cutoff=8,y.cutoff=1,y.high.cutoff=Inf,cex.use=0.5,cex.text.use=0.5,do.spike=FALSE,
                   pch.use=16, col.use="black", spike.col.use="red",plot.both=FALSE,do.contour=TRUE,
                   contour.lwd=3, contour.col="white", contour.lty=2,num.bin=20,do.recalc=TRUE) {
            
            genes.use = set.ifnull(genes.use, rownames(object@data))
            data=object@data[genes.use,]
            
            if (do.recalc) {    
              data.x=rep(0,length(genes.use)); names(data.x)=genes.use; data.y=data.x; data.norm.y=data.x;
              
              bin.size <- 1000
              max.bin <- floor(length(genes.use)/bin.size) + 1
              cat("Calculating gene dispersion", file = stderr())
              pb <- txtProgressBar(min = 0, max = max.bin, style = 3, file = stderr())
              for(i in 1:max.bin) {
                my.inds <- ((bin.size * (i - 1)):(bin.size * i - 1))+1
                my.inds <- my.inds[my.inds <= length(genes.use)]
                genes.iter=genes.use[my.inds]; data.iter=data[genes.iter,]
                data.x[genes.iter]=apply(data.iter,1,fxn.x); data.y[genes.iter]=apply(data.iter,1,fxn.y)
                setTxtProgressBar(pb, i)  
              }
              close(pb)
              data.y[is.na(data.y)]=0
              data.x[is.na(data.x)]=0
              
              data_x_bin=cut(data.x,num.bin)
              names(data_x_bin)=names(data.x)
              mean_y=tapply(data.y,data_x_bin,mean)
              sd_y=tapply(data.y,data_x_bin,sd)
              data.norm.y=(data.y-mean_y[as.numeric(data_x_bin)])/sd_y[as.numeric(data_x_bin)]
              #data.x=apply(data,1,fxn.x); data.y=apply(data,1,fxn.y); data.x[is.na(data.x)]=0
              #data.norm.y=meanNormFunction(data,fxn.x,fxn.y,num.bin)
              
              data.norm.y[is.na(data.norm.y)]=0
              names(data.norm.y)=names(data.x)
              
              mv.df=data.frame(data.x,data.y,data.norm.y)
              rownames(mv.df)=rownames(data)
              #object@mean.var=mv.df
            }
            data.x=mv.df[,1]; data.y=mv.df[,2]; data.norm.y=mv.df[,3]; 
            names(data.x)=names(data.y)=names(data.norm.y)=rownames(object@data[genes.use,])
            
            pass.cutoff=names(data.x)[which(((data.x>x.low.cutoff) & (data.x<x.high.cutoff)) & (data.norm.y>y.cutoff) & (data.norm.y < y.high.cutoff))]
            if (do.spike) spike.genes=rownames(subr(data,"^ERCC"))
            if (do.plot) {
              if (plot.both) {
                par(mfrow=c(1,2))
                smoothScatter(data.x,data.y,pch=pch.use,cex=cex.use,col=col.use,xlab="Average expression",ylab="Dispersion",nrpoints=Inf)
                
                if (do.contour) {
                  data.kde=kde2d(data.x,data.y)
                  contour(data.kde,add=TRUE,lwd=contour.lwd,col=contour.col,lty=contour.lty)
                }
                if (do.spike) points(data.x[spike.genes],data.y[spike.genes],pch=16,cex=cex.use,col=spike.col.use)
                if(do.text) text(data.x[pass.cutoff],data.y[pass.cutoff],pass.cutoff,cex=cex.text.use)
              }
              smoothScatter(data.x,data.norm.y,pch=pch.use,cex=cex.use,col=col.use,xlab="Average expression",ylab="Dispersion",nrpoints=Inf)
              if (do.contour) {
                data.kde=kde2d(data.x,data.norm.y)
                contour(data.kde,add=TRUE,lwd=contour.lwd,col=contour.col,lty=contour.lty)
              }
              if (do.spike) points(data.x[spike.genes],data.norm.y[spike.genes],pch=16,cex=cex.use,col=spike.col.use,nrpoints=Inf)
              if(do.text) text(data.x[pass.cutoff],data.norm.y[pass.cutoff],pass.cutoff,cex=cex.text.use)
            }
            if (set.var.genes) {
              object@var.genes=pass.cutoff
              return(object)
            } else {
              return(pass.cutoff)
            }
          }
)

#' Find Variable genes using a negative binomial null model
#' cut.quantile - quantile ceiling on max counts
NB.var.genes <- function(
  object = object, 
  cells.use=NULL, 
  min.cells = 200,
  do.idents=NULL,
  genes.use = NULL, 
  do.plot=TRUE,
  set.var.genes=TRUE,
  x.low.cutoff=0.005, 
  x.high.cutoff=3, 
  diffCV.cutoff=NULL,
  num.sd=NULL, 
  cex.use=0.5,
  cex.text.use=0.5,
  do.spike=FALSE,
  pch.use=16, 
  col.use="black", 
  spike.col.use="red",
  do.ident=FALSE, 
  do.text=TRUE, 
  cut.quantile=0.99,
  rem.mt.rp = FALSE,
  max.cor = 0.5) {
  
  require(MASS)
  print("Identifying variable genes based on UMI Counts. Warning - use this only for UMI based data")
  
  genes.use=set.ifnull(genes.use, rownames(object@data))
  
  if (is.null(do.idents)){
      ident.label = "all"
      cells.use=set.ifnull(cells.use,colnames(object@data))  
      multi.idents = FALSE
  } else {
    if (is.logical(do.idents) & do.idents == T){
      do.idents = levels(object@ident)
    }
    ident.label = do.idents[do.idents %in% levels(object@ident)]
    cells.use = NULL
    if (length(ident.label) >0){ multi.idents = TRUE }
  }
  
  genes.return = c()
  cells.final = c()
  for (ident in ident.label){
    print(1)
    if (multi.idents){ 
      print(paste0("Extracting variable genes for sample : ", ident))
      }
    cells.use = set.ifnull(cells.use, which.cells(object, ident))
    cells.final = c(cells.final, cells.use)
    if (length(cells.use) < min.cells){
      print(paste0("Skipping sample ", ident, " as there are fewer than ", min.cells, " cells"))
      cells.use=NULL
      next
    }
    
    count.data=object@count.data[genes.use,cells.use]
    
    # Empirical mean, var and CV
    mean_emp = Matrix::rowMeans(count.data)
    var_emp = Matrix::rowMeans(count.data^2) - mean_emp^2
    genes.use=names(mean_emp)[mean_emp > 0]
    
    mean_emp = mean_emp[genes.use]
    var_emp = var_emp[genes.use]
    cv_emp = sqrt(var_emp) / mean_emp
    
    # NB sampling
    a=Matrix::colSums(count.data)
    a = a[a <= quantile(a,cut.quantile)]
    size_factor =  a/ mean(a)
    fit=fitdistr(size_factor, "Gamma")
    hist(size_factor, 50, probability=TRUE, xlab="N_UMI/<N_UMI>")
    
    curve(dgamma(x, shape=fit$estimate[1], rate=fit$estimate[2]),from=0, to=quantile(size_factor, 0.99), add=TRUE, col="red",
          main="Gamma dist fit for size factor")
    text(0.8*max(size_factor),0.6, paste("shape = ", round(fit$estimate[1],2)))
    text(0.8*max(size_factor),0.5, paste("rate = ", round(fit$estimate[2],2)))
    
    # Gamma distributions of individual genes are just scaled versions. If X ~ Gamma(a,b)
    # then cX ~ Gamma(a, b/c)
    a_i = rep(fit$estimate[1], length(mean_emp)); names(a_i) = names(mean_emp)
    b_i = fit$estimate[2] / mean_emp; names(b_i) = names(mean_emp)
    mean_NB = a_i / b_i; var_NB = a_i*(1+b_i) / (b_i^2)
    cv_NB = sqrt(var_NB)/mean_NB
    diffCV = log(cv_emp) - log(cv_NB)
    
    hist(diffCV,500, main="Select a delta-logCV cutoff for variable gene: ", xlab="delta-logCV", ylim = c(0,1500), xlim = c(0, quantile(diffCV,0.99)))
    
    if (!is.null(num.sd)){
      diffCV.cutoff = mean(diffCV) + num.sd*sd(diffCV)
    }
    
    if (is.null(diffCV.cutoff)){
      diffCV.cutoff = readline("Select a delta-logCV cutoff (genes with a higher value will be considered):")
      diffCV.cutoff = as.numeric(diffCV.cutoff)
    }
    
    
    print(paste0("Using diffCV = ", round(diffCV.cutoff,2), " as the cutoff"))
    abline(v=diffCV.cutoff, col="red")
    Sys.sleep(3)
    
    print(paste0("Considering only genes with mean counts less than ", x.high.cutoff, " and more than ", x.low.cutoff))
    pass.cutoff=names(diffCV)[which(diffCV > diffCV.cutoff & (mean_emp > x.low.cutoff & mean_emp < x.high.cutoff))]
    print(paste0("Found ", length(pass.cutoff), " variable genes"))
    if (multi.idents){
      genes.return = union(genes.return, pass.cutoff)
      print(paste0("Found ", length(genes.return), " unique genes"))
    } else {
      genes.return = union(genes.return, pass.cutoff)
    }
    if (multi.idents){
      mv.df=data.frame(mean_emp,cv_emp)
      rownames(mv.df)=names(mean_emp)
      mv.df$gene = rownames(mv.df)
      mv.df$ident = ident
      
      if (nrow(object@hvg.info)==0){
        object@hvg.info=mv.df
      } else {
        object@hvg.info = rbind(object@hvg.info, mv.df)
      }
      
    } else {
      
      mv.df=data.frame(mean_emp,cv_emp)
      mv.df$ident = "all"
      rownames(mv.df)=names(mean_emp)
      object@hvg.info=mv.df
    }
    
    if (do.spike) spike.genes=grep("^ERCC", rownames(count.data), value=TRUE)
    if (do.plot) {
      
      plot(mean_emp,cv_emp,pch=pch.use,cex=cex.use,col="black",xlab="Mean Counts",ylab="CV (counts)", log="xy")
      curve(sqrt(1/x), add=TRUE, col="red", log="xy", lty=2, lwd=2)
      or = order(mean_NB)
      lines(mean_NB[or], cv_NB[or], col="magenta", lwd=2)
      points(mean_emp[pass.cutoff], cv_emp[pass.cutoff], col="blue", pch=16, cex=cex.use)
      
      if (do.spike) points(mean_emp[spike.genes],cv_emp[spike.genes],pch=16,cex=cex.use,col=spike.col.use)
      if(do.text) text(mean_emp[pass.cutoff],cv_emp[pass.cutoff],pass.cutoff,cex=cex.text.use)
      
    }  
    cells.use = NULL
    
  }
  
  if (is.null(do.idents) ){
    if (do.ident) {
      identify(mean_emp,cv_emp,labels = names(mean_emp))
    }
  }
  
  if (rem.mt.rp){
    print(paste0("Removing genes that have a correlation greater than ", max.cor," with mt.score or rp.score"))
    cor.vec1 = rep(0, length(genes.return))
    cor.vec2 = rep(0, length(genes.return))
    if ("mt.score" %in% colnames(object@data.info)){
      cor.vec1 = apply(object@data[genes.return,cells.final],1,function(x) abs(cor(x, object@data.info[cells.final,"mt.score"])))
    }
    
    if ("rp.score" %in% colnames(object@data.info)){
    cor.vec2 = apply(object@data[genes.return,cells.final],1,function(x) abs(cor(x, object@data.info[cells.final,"rp.score"])))
    }
    var.genes = genes.return[cor.vec1 < max.cor & cor.vec2 < max.cor]
    print(paste0("Removed ", length(genes.return) - length(var.genes), " variable genes"))
    genes.return = var.genes
  }
  
  if (set.var.genes) { 
    object@var.genes=genes.return
    return(object)
  }
  if (!set.var.genes) return(genes.return)
          
}
          

setGeneric("prune.var.genes", function(object, cells.use=NULL, genes.use = NULL, pos.corr.thresh=0.3, neg.corr.thresh=-0.4, n.thresh=4, n.pos.thresh=1, set.var.genes=TRUE, do.remove=TRUE, do.print=FALSE) standardGeneric("prune.var.genes"))
setMethod("prune.var.genes", signature = "scR",
          function(object,cells.use=NULL, genes.use=NULL, pos.corr.thresh=0.3, neg.corr.thresh=-0.4, n.thresh=3, n.pos.thresh=1, set.var.genes=TRUE, do.remove=TRUE, do.print=FALSE) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            genes.use=set.ifnull(genes.use, rownames(object@data))
            data.use=object@data[genes.use,cells.use]
            
            var.genes0 = object@var.genes
            print(length(var.genes0))
            if (length(var.genes0)==0){
              print("ERROR: First set variable genes")
              return(object)
            }
            
            data.var.genes0 = data.use[var.genes0, ]
            
            #First expand variable gene list by including correlated genes 
            var.genes1 = var.genes0
            test.genes = setdiff(genes.use, var.genes0)
            nblocks = ceiling(length(test.genes) /  1000)
            rem = length(test.genes) %% 1000
            
            count=0
            for (i in 1:nblocks){
              
              if (count+1000 < length(test.genes)){
                count_plus = count+1000
              } else {
                count_plus = length(test.genes)
              }
              if (do.print){print(paste0("Testing genes ", count+1, "-", count_plus))}
              data1 = data.use[test.genes[(count+1):(count_plus)],]
              genes.in.block = rownames(data1)
              CorrMat = cor(t(data1), t(data.var.genes0))
              
              PosCorr = rowSums(1 * (CorrMat > pos.corr.thresh))
              NegCorr = rowSums(1 * (CorrMat < neg.corr.thresh))
              
              pass.ind = which((PosCorr >= n.pos.thresh) & (PosCorr + NegCorr >= n.thresh))
              for (k in pass.ind){
                if (do.print){
                print(paste0("Adding ", genes.in.block[k], " with ", PosCorr[k], " +ve correlations and ", NegCorr[k], " -ve correlations within old set"))
                }
                var.genes1 = c(var.genes1, genes.in.block[k])
              }    
              count=count_plus
            }
            print(length(var.genes1))
            
            #Prune var.genes
            if (do.remove){
              data.var.genes1 = data.use[var.genes1, ]
              var.genes2 = c()
              
              nblocks = ceiling(length(var.genes1) /  100)
              rem = length(var.genes1) %% 100
              
              count=0
              for (i in 1:nblocks){
                
                if (count+100 < length(var.genes1)){
                  count_plus = count+100
                } else {
                  count_plus = length(var.genes1)
                }
                if (do.print){
                print(paste0("Testing genes ", count+1, "-", count_plus))
                }
                data1 = data.var.genes1[var.genes1[(count+1):(count_plus)],]
                genes.in.block = rownames(data1)
                CorrMat = cor(t(data1), t(data.var.genes1))
                
                PosCorr = rowSums(1 * (CorrMat > pos.corr.thresh))
                NegCorr = rowSums(1 * (CorrMat < neg.corr.thresh))
                
                pass.ind = which(PosCorr + NegCorr > 1)
                print(length(pass.ind))
                nopass.ind = which(PosCorr + NegCorr <= 1)
                print(length(nopass.ind))
                var.genes2 = c(var.genes2, genes.in.block[pass.ind])
                
                for (k in nopass.ind){
                  if (do.print){
                  print(paste0("Removing ", genes.in.block[k], " with no correlations within expanded set"))
                  }
                }
                count=count_plus
              }
              } else {
                var.genes2=var.genes1
            }  
            
            if (set.var.genes) { 
              object@var.genes=var.genes2
              return(object)
            } else { return(var.genes2)}
            
          
            }
          
)



setGeneric("jackstraw.permutation.test", function(object, genes.use=NULL, use.full=FALSE, n.resamp=100, p.thres=0.05, do.verbose=TRUE, seed.val=NULL, max.pc=1000, n.cores=1) standardGeneric("jackstraw.permutation.test"))
setMethod("jackstraw.permutation.test", "scR",
          function(object, genes.use=NULL, use.full=FALSE, n.resamp=100, p.thres=0.05, do.verbose=TRUE, seed.val=NULL, max.pc=1000, n.cores=1) {
            genes.use=set.ifnull(genes.use,object@var.genes)
            data = as.matrix(object@scale.data[genes.use, ])
            if (dim(data)[1] < dim(data)[2]) {
              data = t(data)
            }
            perm.obj = permutationPA(data, B=n.resamp, threshold=p.thres, verbose=do.verbose, seed=seed.val, max.pc=max.pc, n.cores=n.cores)
            return(perm.obj)
          })
          


setGeneric("doGOseq", function(object,bg.genes=NULL,fg.genes=NULL, genome=NULL, id = NULL, geneLengths=NULL, p.adj.method="BH", pval.thres=1e-4,... ) standardGeneric("doGOseq"))
setMethod("doGOseq", "scR", 
          function(object,bg.genes=NULL,fg.genes=NULL, genome=NULL, id = NULL, geneLengths=NULL,p.adj.method="BH", pval.thres=1e-4,... ) {
                bg.genes = set.ifnull(bg.genes, rownames(object@data))
                fg.genes = set.ifnull(fg.genes, object@var.genes)
                kill.ifnull(genome, "Error: Need to specify Genome")
                kill.ifnull(id, "Error: Need to specify gene identifier")
                
                gene.vector <- as.integer(bg.genes%in%fg.genes)
                names(gene.vector) <- bg.genes
                
                pwf <- nullp(gene.vector, genome,id,bias.data=geneLengths) #Empirical null distribution with length correction
                GO.wall=goseq(pwf,genome,id,...)
                enriched.GO = GO.wall[p.adjust(GO.wall$over_represented_pvalue, method=p.adj.method) < pval.thres, 
                                      c("category","over_represented_pvalue", "term","ontology", "numDEInCat", "numInCat")]
                rownames(enriched.GO) = enriched.GO$category
                enriched.GO = enriched.GO[,-1]
                
                return(enriched.GO)
            
          }
)



setGeneric("PlotMetaData", function(object,features=NULL,adjust.use=2,...) standardGeneric("PlotMetaData"))
setMethod("PlotMetaData", "scR", 
          function(object,features=NULL,adjust.use=2,...) {
            
            features = set.ifnull(features, colnames(object@data.info))
            features = union(features,"orig"); features = setdiff(features,"m")
            
            meta.data = object@data.info[,features]
            
            pList = list(); l=1;
            for (feature in colnames(meta.data)){
              if (!(is.numeric(object@data.info[,feature]))) next
              print(l)
              x = as.numeric(meta.data[,feature])
              pList[[l]] = ggplot(meta.data, aes_string(x=feature)) + theme_bw() + xlab(feature) + ylab("Freq")
              
              if (max(x) > 5e4){
                pList[[l]] = pList[[l]] + geom_histogram(aes(fill=orig), binwidth=0.2, alpha=0.6) + scale_x_log10()
                  
              } else {
                pList[[l]] = pList[[l]] + geom_histogram(aes(fill=orig), binwidth=max(x)/20, alpha=0.6)
              }
              
              l=l+1;
              
            }
            multiplotList(pList, cols=2) 
            rp()
          }
)


AddMetaData <- function(object, metadata, col.name = NULL) {
  if (typeof(x = metadata) != "list") {
    metadata <- as.data.frame(x = metadata)
    if (is.null(x = col.name)) {
      stop("Please provide a name for provided metadata")
    }
    colnames(x = metadata) <- col.name
  }
  cols.add <- colnames(x = metadata)
  meta.add <- metadata[rownames(x = object@data.info), cols.add]
  if (all(is.null(x = meta.add))) {
    stop("Metadata provided doesn't match the cells in this object")
  }
  object@data.info[, cols.add] <- meta.add
  return(object)
}

setGeneric("CorrelateMetaData", function(object,features=NULL,pcs.use = c(1,2,3),...) standardGeneric("CorrelateMetaData"))
setMethod("CorrelateMetaData", "scR", 
          function(object,features=NULL,pcs.use = c(1,2,3),...) {
            
            features = set.ifnull(features, c("nReads","nGene","nTranscripts") )
            if (ncol(object@data.info) > 3){
              features = c(features, colnames(object@data.info)[4:ncol(object@data.info)])
            }
            features = features[features %in% colnames(object@data.info)]
            features = setdiff(features,"orig")
            features = unique(features)
            
            MetaData = object@data.info[,features]
            
            PC = object@pca.rot[,pcs.use]
            
            pc_vec = sapply(pcs.use, function(x) paste0("PC",x))
            df = data.frame(PC=pc_vec)
            print(features)
            print(pc_vec)
            for (f in features){
               v = c()
               for (pc in pc_vec){
                 v = c(v, -log10(cor.test(PC[,pc], MetaData[,f])$p.value))
               }
               df = cbind(df, v)
            }
            colnames(df)[2:ncol(df)] = features
            df$PC = factor(as.character(df$PC), levels=as.character(df$PC))
            
            #melt
            df.m = melt(df, id.vars = "PC", variable_name="feature")
            p = ggplot(df.m, aes(PC, value)) + geom_bar(stat="identity", aes(fill=feature),position="dodge", color="black") + theme_bw()
            p = p + xlab("PC") + ylab("-log10(P-val)") + ylim(0,5) + ggtitle("Correlation between PCs and QC metrics")
            p = p + theme(axis.title.y = element_text(face="bold", colour="#990000", size=20), axis.text.y  = element_text(size=15),
                          axis.title.x = element_text(face="bold", colour="#990000", size=20), axis.text.x  = element_text(size=15),
                          legend.text = element_text(size = 15), legend.title= element_text(size = 20,colour="#990000"))
            
            print(p)
            
          }
)


setGeneric("QCnorm", function(object,qc.features=NULL,var.thresh=80,return.scaled=FALSE,n.bin.use=10,...) standardGeneric("QCnorm"))
setMethod("QCnorm", "scR", 
          function(object,qc.features=NULL,var.thresh=80,return.scaled=FALSE,n.bin.use=10,...) {
            
            
            qc.features = set.ifnull(qc.features, c("nReads","nGene","nTranscripts") )
            if (ncol(object@data.info) > 3){
              qc.features = c(qc.features, colnames(object@data.info)[4:ncol(object@data.info)])
            }
            qc.features = qc.features[qc.features %in% colnames(object@data.info)]
            qc.features = setdiff(qc.features,c("orig","m"))
            qc.features = unique(qc.features)
            
            print("Correcting for features")
            print(qc.features)
            
            MetaData = object@data.info[,qc.features]
            MetaData.scaled = scale(MetaData, center=TRUE, scale=TRUE)
            
            pca.obj = prcomp(t(MetaData.scaled),center=FALSE, scale=FALSE)
            pca.var = 100*pca.obj$sdev^2 / sum(pca.obj$sdev^2)
            var.explained=0;
            for (i in 1:dim(MetaData)[2]){
              var.explained = var.explained + pca.var[i]
              if (var.explained > var.thresh){
                break
              }
            }
            
            pcs.use = c(1:i)
            print(paste0(i, " significant PCs for quality metrics"))
            
            n.bins = n.bin.use
            pc.bins = data.frame(PC1=cut(pca.obj$rotation[,1], n.bins))
            for (pc in setdiff(pcs.use,1)){
              pc.bins = cbind(pc.bins,cut(pca.obj$rotation[,pc],n.bins))
              colnames(pc.bins)[pc] = paste0("PC",pc)
            }
            
            #Normalization
            ExpMat = object@data
            #ExpMat0 = object@data
            #MedMat0 = apply(object@data,1,median)
           #MedMat0 = as.data.frame(replicate(ncol(object@data), MedMat0)); rownames(MedMat0) = rownames(object@data); colnames(MedMat0) = colnames(object@data)
            
            for (pc in colnames(pc.bins)){
              for (i in levels(pc.bins[,pc])){
                  cells.in.bin = which(pc.bins[,pc]==i)
                  if(length(cells.in.bin) == 0) next
                  if (length(cells.in.bin) > 1){
                     MedMat0 = apply(ExpMat,1,median)
                     MedMat0 = as.data.frame(replicate(ncol(ExpMat), MedMat0)); rownames(MedMat0) = rownames(ExpMat); colnames(MedMat0) = colnames(ExpMat)
                     MedMat = apply(ExpMat[,cells.in.bin],1,median);
                     MedMat = as.data.frame(replicate(length(cells.in.bin), MedMat)); rownames(MedMat) = rownames(ExpMat); colnames(MedMat) = colnames(ExpMat)[cells.in.bin]    
                     ExpMat[,cells.in.bin] = ExpMat[,cells.in.bin] - MedMat + MedMat0[,cells.in.bin]
                  } else {
                    MedMat = ExpMat[,cells.in.bin]
                    MedMat0 = apply(ExpMat,1,median)
                    ExpMat[,cells.in.bin] = ExpMat[,cells.in.bin] - MedMat + MedMat0
                  }
              }
            }
            
            if (return.scaled){
              ExpMat.scaled=t(scale(t(ExpMat),center=TRUE,scale=TRUE))
              ExpMat.scaled=ExpMat.scaled[complete.cases(ExpMat.scaled),]
              return(ExpMat.scaled)
            } else {
              return(ExpMat)
            }
            
          }
)

setGeneric("Avg.by.ident", function(object,features.use=NULL, use.count=FALSE, ident.use=NULL, return.scaled=TRUE,fxn.x=expMean,...) standardGeneric("Avg.by.ident"))
setMethod("Avg.by.ident", "scR", 
          function(object,features.use=NULL, use.count=FALSE, ident.use=NULL, return.scaled=TRUE,fxn.x=expMean,...) {

    features.use=set.ifnull(features.use, object@var.genes)
    features.use = ainb(features.use, rownames(object@data))
    #features.use = ainb(features.use, colnames(object@data.info))
    ident.use=set.ifnull(ident.use, levels(object@ident))
    
    ExpMat = matrix(0, nrow=length(features.use), ncol = length(ident.use))
    rownames(ExpMat) = features.use; colnames(ExpMat) = ident.use
    
    if (use.count){
      data.use = object@count.data[features.use,]
      #vec.exp = apply(object@count.data[features.use, cells.in.cluster], 1, function(x) fxn.x(x)) 
    } else {
      data.use = object@data[features.use,]
      #data.use = t(object@data.info[,features.use])
      #vec.exp = apply(object@data[features.use, cells.in.cluster], 1, function(x) fxn.x(x)) 
    }
    data.use = as.matrix(data.use)
    for (i in ident.use){
      cells.in.cluster = names(object@ident)[which(object@ident== i)]
      vec.exp = apply(data.use[features.use, cells.in.cluster], 1, function(x) fxn.x(x)) 
      ExpMat[, i] = vec.exp
    }
    
    if (return.scaled){
      ExpMat.scaled = t(scale(t(ExpMat), scale=TRUE, center=TRUE))
      return(ExpMat.scaled)    
    } else {
      return(ExpMat)
    }
}
)



setGeneric("PercExpDf", function(object,features.use=NULL, ident.use=NULL, thresh.use=NULL, perc.floor = 0.3,exp.floor = 1,...) standardGeneric("PercExpDf"))
setMethod("PercExpDf", "scR", 
          function(object,features.use=NULL,ident.use=NULL, thresh.use=NULL, perc.floor = 0.3,exp.floor = 1,...) {
            
            features.use=set.ifnull(features.use, object@var.genes[1:20])
            features.use=features.use[features.use %in% rownames(object@data)]
            ident.use=set.ifnull(ident.use, levels(object@ident))
            thresh.use=set.ifnull(thresh.use,object@is.expr)
            
            #Matrix of percent expressing cells
            PercMat = matrix(0, nrow=length(features.use), ncol = 0)
            rownames(PercMat) = features.use; 
            
            #Matrix of average transcript levels
            ExpMat = PercMat;
            
            #Count mat
            Count.mat = object@count.data[features.use, colnames(object@data)]
            
            
            for (i in ident.use){
              cells.in.cluster = names(object@ident)[which(object@ident== i)]
              vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) sum(x>thresh.use)/length(x)) 
              PercMat = cbind(PercMat,vec.exp)
              
              vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) if (sum(x>0) > 1){ mean(x[x>0]) } else {sum(x)})
              ExpMat = cbind(ExpMat, vec.exp)
            }
            colnames(ExpMat) = ident.use
            colnames(PercMat) = ident.use
            
            if (!is.null(perc.floor)){
              print(paste0("Zeroing all percentage values less than ", perc.floor, " percent." ))
             PercMat[PercMat < perc.floor] = 0
            }
            
            if (!is.null(exp.floor)){ 
                print(paste0("Zeroing all percentage values less than ", exp.floor, "." ))
                 ExpMat[ExpMat < exp.floor] = 0
            }
            
            df = data.frame(row.names=ident.use)
            df[,"clust"] = ident.use
            for (f in features.use){
              df[,paste0("Exp_",f)] = ExpMat[f,]
              df[,paste0("Perc_",f)] = PercMat[f,]
            }

            
          return(df)  
          
          }
)

setGeneric("Perc.pos.ident.matrix", function(object,feature, ident1.use="m",ident2.use="orig", thresh.use=0,do.plot=TRUE,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,cols.use=gray.colors(10),max.size=10,min.perc=0,...) standardGeneric("Perc.pos.ident.matrix"))
setMethod("Perc.pos.ident.matrix", "scR", 
          function(object,feature,ident1.use="m",ident2.use="orig", thresh.use=0,do.plot=TRUE,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,cols.use=gray.colors(10),max.size=10,min.perc=0,...) {
            
            
            feature=feature[feature %in% rownames(object@data)]
            ident1=object@data.info[,ident1.use]; names(ident1) = rownames(object@data.info)
            ident2=object@data.info[,ident2.use]; names(ident2) = rownames(object@data.info)
            
            #Matrix of percent expressing cells
            PercMat = matrix(0, nrow=length(levels(ident1)), ncol = length(levels(ident2)))
            colnames(PercMat) = levels(ident2)
            rownames(PercMat) = levels(ident1) 
            
            #Matrix of average transcript levels
            ExpMat = PercMat;
            
            #Count mat
            Count.mat = object@count.data[feature, colnames(object@data)]
            
            
            for (i in 1:length(levels(ident1))){
              for (j in 1:length(levels(ident2))){
                
                cells.in.cluster = intersect(rownames(object@data.info)[which(ident1== levels(ident1)[i])],rownames(object@data.info)[which(ident2== levels(ident2)[j])])
                x=Count.mat[cells.in.cluster]
                PercMat[i,j] = sum(x > thresh.use)/length(x)
                if (length(cells.in.cluster) > 0){
                  ExpMat[i,j] = mean(x[x>0])
                } else {
                  ExpMat[i,j] = x
                }
              }
              
            }
            
            
            if (!is.null(max.val.perc)) PercMat[PercMat > max.val.perc] = max.val.perc
            if (!is.null(max.val.exp)) ExpMat[ExpMat > max.val.exp] = max.val.exp
            
            if (do.plot){
              
              ExpVal = melt(ExpMat)
              PercVal = melt(PercMat)
              colnames(ExpVal) = c("m","orig","nTrans")
              ExpVal$percExp = PercVal$value*100
              
              if (!do.transpose){
                ExpVal$m = factor(ExpVal$m, levels=rev(levels(ident1)))
                ExpVal$orig = factor(ExpVal$orig, levels= levels(ident2))
                p=ggplot(ExpVal, aes(y = factor(m),  x = factor(orig))) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                  scale_color_gradient(low ="blue",   high = "red", limits=c(1, max(ExpVal$nTrans) ))+scale_size(range = c(1, max.size))+   theme_bw() +nogrid
                p = p + xlab("m") + ylab("orig") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
                  theme(axis.text.y=element_text(size=12, face="italic"))  
                print(p)
              } else {
                ExpVal$m = factor(ExpVal$m, levels=ident1)
                ExpVal$orig = factor(ExpVal$orig, levels= rev(levels(ident2)))
                p=ggplot(ExpVal, aes(y = factor(orig),  x = factor(m))) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                  scale_color_gradient(low ="blue",   high = "red", limits=c( 1, max(ExpVal$nTrans) ))+scale_size(range = c(1, max.size))+   theme_bw() +nogrid
                p = p + ylab("orig") + xlab("m") + theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1)) + 
                  theme(axis.text.y=element_text(size=12, face="italic"))
                print(p)
                
              }
              
            }else {
              to.return=list()
              to.return$ExpMat = ExpMat;
              to.return$PercMat = PercMat;
              return(to.return)
            }
          }
)


setGeneric("clust.composition", function(object, sample.use=NULL, col.ident="m", row.ident="orig", scale.by.col=TRUE, do.scale=TRUE,max.size=20,...) standardGeneric("clust.composition"))
setMethod("clust.composition", "scR", 
          function(object,sample.use=NULL,col.ident="m",row.ident="orig", scale.by.col=TRUE, do.scale=TRUE,max.size=20,...) {
            
            sample.use=set.ifnull(sample.use, levels(object@data.info$orig))
            counts = table(object@data.info[,row.ident], object@data.info[,col.ident])
            
            
            
            if (do.scale){
              if (scale.by.col){counts.scale = scale(counts, center=FALSE, scale=colSums(counts))*100
              print(1)
              PercVal = melt(counts.scale)
              colnames(PercVal) = c("sample","cluster","perc_cluster")
              PercVal$cluster = factor(PercVal$cluster)
              PercVal$sample = factor(PercVal$sample, levels= sample.use)
              p=ggplot(PercVal, aes(y = factor(sample),  x = factor(cluster))) + geom_point(aes(size=perc_cluster))
              }
              if (!scale.by.col){
                counts.scale = t(scale(t(counts), center=FALSE, scale=rowSums(counts)))*100
                PercVal = melt(counts.scale)
                colnames(PercVal) = c("sample","cluster","perc_sample")
                PercVal$cluster = factor(PercVal$cluster)
                PercVal$sample = factor(PercVal$sample, levels= sample.use)
                p=ggplot(PercVal, aes(y = factor(sample),  x = factor(cluster))) + geom_point(aes(size=perc_sample))
              }
             
            } else{
              PercVal = melt(counts)
              colnames(PercVal) = c("sample","cluster","counts")
              PercVal$cluster = factor(PercVal$cluster)
              PercVal$sample = factor(PercVal$sample, levels= sample.use)
              p=ggplot(PercVal, aes(y = factor(sample),  x = factor(cluster))) + geom_point(aes(size=counts))
            }
              

            
            p= p+scale_size(range = c(0, max.size))+   theme_bw() +
              xlab("Cluster") + ylab("Sample") + 
              theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
              theme(axis.text.y=element_text(size=12, face="italic"))  
            
            clust.size = as.data.frame(colSums(counts))
            colnames(clust.size) = "clust_size"
            clust.size$cluster = factor(rownames(clust.size), levels=rownames(clust.size))
            p2 = ggplot(clust.size, aes(x = factor(cluster), y="Size")) + geom_point(aes(size=clust_size), color="blue") + theme_bw() + nogrid + 
              scale_size(range = c(0, 0.7*max.size)) + 
              theme(axis.text.x=element_blank(), axis.ticks = element_blank()) + 
              theme(axis.text.y=element_blank()) + theme(legend.position = "top", legend.direction="horizontal") +
              xlab("") + ylab("# Cells")
            p <- ggplot_gtable(ggplot_build(p))
            p2 <- ggplot_gtable(ggplot_build(p2))
            p2$widths <- p$widths
            grid.arrange(p2,p, ncol=1, heights=c(1.5,3))
              
          
          }
)

setGeneric("gene.set.score", function(object,genes.use=NULL, use.scaled=FALSE, score.name=NULL,use.weights=FALSE, return.score=FALSE, verbose=TRUE,min.perc=0.05,max.perc=0.2,scale.by.perc=FALSE,...) standardGeneric("gene.set.score"))
setMethod("gene.set.score", "scR", 
          function(object,genes.use=NULL,use.scaled=FALSE, score.name=NULL,use.weights=FALSE, return.score=FALSE, verbose=TRUE,min.perc=0.05, max.perc=0.2,scale.by.perc=FALSE,...) {
            
            if (is.null(score.name)){
              score.name = readline(prompt="Provide a name for this score (no spaces): ")
            }
            
            genes.use0 = genes.use
            genes.use = ainb(genes.use, rownames(object@data))
            kill.ifnull(genes.use, "Error : Must provide a list of genes that are present in the data")
            
            if (verbose){
              if (length(genes.use0) > length(genes.use)){
                n = length(genes.use0) - length(genes.use)
                print(paste0("Note that ", n, " genes (", round(n*100 / length(genes.use0)), " %) in the gene-set  are not expressed in this dataset"))
              }
            }
            

            
            if (use.scaled){
              data.use = object@scale.data
            } else {
              data.use = object@data
            }
            
            percpos=Perc.pos.by.ident(object, features.use=genes.use, do.plot=FALSE)
            if (!is.null(dim(percpos$PercMat))){
                #genes.use = genes.use[apply(percpos$PercMat,1, function(x) mean(x) >= min.perc & mean(x) <= max.perc)]
                genes.use = genes.use[apply(percpos$PercMat,1, function(x) max(x) >= min.perc)]
            } else {
                #genes.use = genes.use[percpos$PercMat >= min.perc & percpos$PercMat <= max.perc]
                genes.use = genes.use[percpos$PercMat >= min.perc]
            }
            print(paste0("Removing genes that are expressed in less than ", 100*min.perc, " % of cells per cluster"))
            print(paste0("Final signature consists of ", length(genes.use), " genes"))
            print(genes.use)
            
            scores=NULL
            if (length(genes.use) > 0 & !scale.by.perc){
            data.use1 = data.use[genes.use,];
            scores = Matrix::colMeans(data.use1) - Matrix::colMeans(data.use)
            }
            
            if (length(genes.use) > 0 & scale.by.perc){
              weights = apply(percpos$PercMat[genes.use,],1,function(x){y=sort(x); max(y) / mean(y[1:10])})
              data.use1 = data.use[genes.use,];
              scores = apply(data.use1,2, function(x) sum(x*weights)/sum(weights)) - Matrix::colMeans(data.use)
            }
            
            if (return.score){
              scores = scores[rownames(object@data.info)]
              return(scores)
            } else {
              object@data.info[, score.name] = scores
              return(object)
            }
              
            }
  )


setGeneric("signature.spatial", function(object,genesets.list=NULL, use.weights=FALSE, reduction.use="tsne",verbose=FALSE, alpha=0.33, min.perc=0.3, min.genes=5, random.genesets.size=c(5,10,20,50,100,200),n.random = 1000,...) standardGeneric("signature.spatial"))
setMethod("signature.spatial", "scR", 
          function(object,genesets.list=NULL,use.weights=FALSE, reduction.use="tsne",  verbose=FALSE, alpha=0.33, min.perc=0.3, min.genes=5, random.genesets.size=c(5,10,20,50,100,200), n.random = 1000,...) {
            
            
            #Spatial coordinates scaled
            coord.use=reduction.space(object, reduction.use=reduction.use, pc.1=1,pc.2=2)
            # kernel matrix
            K = exp(-as.matrix(dist(coord.use, method="euclidean", upper=TRUE))^2 / alpha^2)
            
            #Random genesets
            rand.dist = data.frame(size=random.genesets.size, mean=0, sigma=0)
            index=1
            for (size in random.genesets.size){
              consis.temp = c()
               print(paste0("Evaluating random gene sets of size ", size))
                for (iter in c(1:n.random)){
                    genes.rand = sample(rownames(object@scale.data), size)
                    
                    score = gene.set.score(object, genes.use=genes.rand, use.scaled=TRUE, score.name="Rand", use.weights=use.weights, return.score=TRUE, verbose=verbose)
                    
                    #Expected score 
                    exp_score = (K %*% score) / rowSums(K)
                    
                    #consistency
                    consis.temp = c(consis.temp,median(abs(score-exp_score)))
                }
              rand.dist[index, "mean"] = mean(consis.temp)
              rand.dist[index,"sigma"] = sd(consis.temp)
              index=index+1;

            }
            
            l=0
            ind.tested = c()
            pval.tested = c();
            sig = 0;
            for (i in 1:length(genesets.list$genesets)){
              
              if ( i %% 1000 == 0) {
                print(paste0("Genesets considered:", i, "; Tested:", l, "; Significant:", sig))
              }
              genes = genesets.list$genesets[[i]]
              name = genesets.list$geneset.names[[i]]
              name = gsub("%","_",name)
              name = name.shorten(name)
              
              genes = genes[genes != ""]
              genes0=genes
              genes = ainb(genes, rownames(object@scale.data))
              
              #Check if many genes have been removed
              
              if ( (length(genes)/length(genes0) < min.perc ) | length(genes) < min.genes ){
                next;
              } else {
                ind.tested = c(ind.tested, i)
                l=l+1
                score = gene.set.score(object, genes.use=genes, use.scaled=TRUE, score.name=name, use.weights=use.weights, return.score=TRUE, verbose=verbose)
                
                #Expected score 
                exp_score = K %*% score / rowSums(K)
              
                #consistency
                consis = median(abs(score-exp_score))
                
                #Statistical significance
                index = which.min(abs(random.genesets.size-length(genes)))
                pval = pnorm(consis, mean=rand.dist[index, "mean"], sd = rand.dist[index, "sigma"], lower.tail=FALSE)
                pval.tested = c(pval.tested, pval)
                
                if (pval < 1e-3){
                  sig=sig+1;
                }
              }
              
            }
            
            
          }
)

name.shorten = function(name){
    name = gsub("DIFFERENTIATION","DIFF",name)  
    name = gsub("NEGATIVE","NEG",name)
    name = gsub("POSITIVE","POS",name)
    name = gsub("PROCESS","PROC",name)
    name = gsub("REGULATION","REG",name)
    name = gsub("STRUCTURAL","STRUC",name)
    name = gsub("STRUCTURE","STRUC",name)
    name = gsub("ACTIVITY","ACT",name)
    name = gsub("ACTIVATION","ACT",name)
    name = gsub("DEVELOPMENT","DEV",name)
    name = gsub("TRANSMEMBRANE","TRANSMEM",name)
    name =  str_replace_all(name,"[^[:alnum:]]", " ")
}

setGeneric("reduction.space", function(object,reduction.use="tsne",pc.1=1,pc.2=2,do.scale=TRUE,...) standardGeneric("reduction.space"))
setMethod("reduction.space", "scR", 
          function(object,reduction.use="tsne",pc.1=1,pc.2=2,do.scale=TRUE...) {
            
            if (reduction.use=="pca") {
              data.plot=object@pca.rot[,c(pc.1, pc.2)]
            }
            if (reduction.use=="tsne") {
              data.plot=object@tsne.rot
            }
            if (reduction.use=="ica") {
              data.plot=object@ica.rot[,c(pc.1, pc.2)]
            }
            
            if (do.scale) data.plot = scale(data.plot)
            return(as.data.frame(data.plot))
            
          }
)

setGeneric("doBatchCorrection", function(object, genes.use=NULL,  batch.cov=NULL, max.val=6, set.data=FALSE, set.scale.data=FALSE) standardGeneric("doBatchCorrection"))
# Batch Correction using ComBat
# Memory heavy - consider running on a cluster
setMethod("doBatchCorrection","scR",
          function(object,  genes.use=NULL,  batch.cov=NULL, max.val=6, set.data=FALSE, set.scale.data=FALSE) {
            require(sva)
            print("Using ComBat for batch correction")
            genes.use=set.ifnull(genes.use, rownames(object@data))
            pheno=data.frame(batch=batch.cov); rownames(pheno) = colnames(object@data)
            modcombat = model.matrix(~1, data=pheno)
            correct.data = ComBat(dat=object@data[genes.use,],batch=batch.cov,mod=modcombat, prior.plots=TRUE, par.prior=TRUE)
            correct.data[correct.data > max.val] = max.val
            
            if (set.data) object@data[genes.use,] = as.data.frame(correct.data[genes.use,])
            if (set.scale.data) object@scale.data[genes.use,] = t(scale(t(correct.data), center=TRUE, scale=TRUE))
            if (set.data | set.scale.data){
               return(object)
            } else {
              print("Returning unscaled batch corrected matrix")
              return(correct.data)
            }
          }
)

setGeneric("doLMBatchCorrection", function(object, genes.use=NULL,  batch.cov=NULL, max.val=6, set.data=FALSE, set.scale.data=FALSE) standardGeneric("doLMBatchCorrection"))
# Batch Correction using a linear model
# Memory heavy - consider running on a cluster
setMethod("doLMBatchCorrection","scR",
          function(object,  genes.use=NULL,  batch.cov=NULL, max.val=6, set.data=FALSE, set.scale.data=FALSE) {
            require(sva)
            print("Using ComBat for batch correction")
            genes.use=set.ifnull(genes.use, rownames(object@data))
            pheno=data.frame(batch=batch.cov); rownames(pheno) = colnames(object@data)
            modcombat = model.matrix(~1, data=pheno)
            correct.data = ComBat(dat=object@data[genes.use,],batch=batch.cov,mod=modcombat, prior.plots=TRUE, par.prior=TRUE)
            correct.data[correct.data > max.val] = max.val
            
            if (set.data) object@data[genes.use,] = as.data.frame(correct.data[genes.use,])
            if (set.scale.data) object@scale.data[genes.use,] = t(scale(t(correct.data), center=TRUE, scale=TRUE))
            if (set.data | set.scale.data){
              return(object)
            } else {
              print("Returning unscaled batch corrected matrix")
              return(correct.data)
            }
          }
)



setGeneric("summarize.clusters", function(object, markers.table=NULL, tier.pos.min = c(0.3, 0.8, 0.6, 0.5), tier.neg.max = c(0.05, 0.3, 0.3, 0.1), min.cells.clust=10, max.other.clust=c(0, 5, 8, 60)) standardGeneric("summarize.clusters"))
setMethod("summarize.clusters","scR",
          function(object, markers.table=NULL, tier.pos.min = c(0.3, 0.8, 0.6, 0.5), tier.neg.max = c(0.05, 0.3,0.3, 0.1), min.cells.clust=10, max.other.clust=c(0, 5, 8,60)) {
            ident.use=object@ident
            print(paste0("Considering clusters that have greater than ", min.cells.clust, " cells"))
            clust.use = names(table(ident.use))[table(ident.use) > min.cells.clust]
            
            print(paste0("Using threshold for expression as ", exp(object@is.expr) - 1))
            expmat = Perc.pos.by.ident(object, features.use=rownames(object@data), do.plot = FALSE, ident.use=clust.use, thresh.use = exp(object@is.expr) - 1)
            
            clust.summary = data.frame(clust=as.character())
            if (length(tier.pos.min) != length(tier.neg.max)) stop("Error : tier.pos.min must be the same length as tier.pos.max")
            if (length(tier.pos.min) != length(max.other.clust)) stop("Error : tier.pos.min must be the same length as max.other.clust")
            l=1
            for (tier in 1:length(tier.pos.min)){
              clust.summary[,paste0("tier",l)] = as.character()
              l=l+1
            }
            
            if (is.null(markers.table)) stop("Error : run find_all_markers and provide output of markers table")
            for (clust in clust.use){
              subtable = subset(markers.table, cluster==clust)
              l=1 
              genes.use=subtable$gene
              temp = data.frame(clust=clust)
              for (tier in 1:length(tier.pos.min)){
                print(paste0("For tier", l, " Require gene to be expressed in at >= ", tier.pos.min[tier]*100, "% of cluster and > ", 
                             tier.neg.max[tier]*100, "% in atmost ", max.other.clust[tier], " clusters"))
                genes1 = genes.use[(expmat$PercMat[genes.use, clust] >= tier.pos.min[tier]) & apply(expmat$PercMat[genes.use,],1, function(x) sum(x>tier.neg.max[tier]) < max.other.clust[tier] + 1 )]
                temp[,paste0("tier",l)] = paste(genes1, collapse=",")
                genes.use=setdiff(genes.use,genes1)
                if (length(genes.use)==1) break
                l=l+1
              }
              clust.summary = rbind(clust.summary, temp)
            }
            
            return(clust.summary)
            
          }
)



#' Regress out technical effects and cell cycle
#'
#' Remove unwanted effects from scale.data
#'
#'
#' @param object Seurat object
#' @param latent.vars effects to regress out
#' @param latent.vars.neg effects to model but not regress out
#' @param genes.regress gene to run regression for (default is all genes)
#' @param model.use Use a linear model or generalized linear model (poisson, negative binomial) for the regression. Options are 'linear' (default), 'poisson', and 'negbinom'
#' @param use.umi Regress on UMI count data. Default is FALSE for linear modeling, but automatically set to TRUE if model.use is 'negbinom' or 'poisson' 
#' @param do.scale Whether to scale the data. 
#' @param do.center Whether to center the data.
#' @param scale.max Max value to accept for scaled data. The default is 10. Setting this can help 
#' reduce the effects of genes that are only expressed in a very small number of cells. 
#' @return Returns Seurat object with the scale.data (object@scale.data) genes returning the residuals from the regression model
#' @import Matrix
#' @export
setGeneric("RegressOut", function(object,latent.vars=NULL, latent.vars.neg=NULL,genes.regress=NULL, model.use="linear", use.umi=F, scale.max = 10, set.regress.data=TRUE, exclude.genes.prior=NULL, scale.factor=NULL) standardGeneric("RegressOut"))
#' @export
setMethod("RegressOut", "scR",
          function(object,latent.vars=NULL,latent.vars.neg=NULL,genes.regress=NULL, model.use="linear", use.umi=F, scale.max = 10, set.regress.data=TRUE,exclude.genes.prior=NULL, scale.factor=NULL) {
            genes.regress=set.ifnull(genes.regress,rownames(object@data))
            genes.regress=ainb(genes.regress,rownames(object@data))
            latent.data=fetch.data(object,c(latent.vars,latent.vars.neg))
            
            if (!is.null(latent.vars.neg) & model.use != "linear"){
              print("WARNING : latent.vars.neg is not supported for non-linear models")
            }
            bin.size <- 100;
            if (model.use=="negbinom") bin.size=5;


            if (!is.null(exclude.genes.prior)){
              exclude.genes.prior = ainb(exclude.genes.prior, rownames(object@data))
              print(paste0("Re-doing the normalization by excluding ", length(exclude.genes.prior), " genes"))
              genes.use = setdiff(rownames(object@data), exclude.genes.prior)
              genes.regress = setdiff(genes.regress, exclude.genes.prior)
              count.data = object@count.data[genes.use,]
              if (is.null(scale.factor)){
                print(paste0("Performing median normalization"))
                num.trans = Matrix::colSums(count.data)
                scale.factor = median(num.trans)
              }
              normalized.data <-  scale.factor * t(t(count.data) / num.trans)
              normalized.data <- log(normalized.data + 1)
              data.use = normalized.data[genes.regress, , drop=FALSE]
              rm(normalized.data)
              
            } else {
              data.use=object@data[genes.regress, , drop=FALSE];
            }
            print(1)
            if (is.null(latent.vars)){
              # No regression, just storing renormalized matrix in regress.data
              object@regress.data = as.matrix(object@data)
              object@regress.data[rownames(data.use),] = as.matrix(data.use[,object@cell.names])
              object@regress.data = Matrix(object@regress.data, sparse=TRUE)
              return(object)
            }
            
            bin.ind <- ceiling(1:length(genes.regress)/bin.size)
            max.bin <- max(bin.ind)
            
            print(paste("Regressing out",c(latent.vars, latent.vars.neg)))
            pb <- txtProgressBar(min = 0, max = max.bin, style = 3, file = stderr())
            data.resid=c()
            
            if (model.use != "linear") {
              use.umi=T
            }
            
            if (!is.null(latent.vars.neg)){
              print("Adding back effects of latent.vars.neg")
            }
            
            if (use.umi) data.use=object@count.data[genes.regress,object@cell.names, drop=FALSE]
            for(i in 1:max.bin) {
              genes.bin.regress <- rownames(data.use)[bin.ind == i]
              gene.expr <- as.matrix(data.use[genes.bin.regress, , drop=FALSE])
              new.data <- do.call(rbind, lapply(genes.bin.regress, function(x) {
                regression.mat = cbind(latent.data, gene.expr[x,])
                colnames(regression.mat) <- c(colnames(latent.data), "GENE")
                fmla=as.formula(paste("GENE ", " ~ ", paste(latent.vars,collapse="+"),sep=""));
                if (model.use=="linear"){
                  lmodel = lm(fmla,data = regression.mat)
                  if (is.null(latent.vars.neg)){
                    return(lmodel$residuals)
                  } else {
                    mm = model.matrix(lmodel)
                    
                    coefs = c()
                    for (var in latent.vars.neg){
                      coefs = c(coefs, grep(paste0("^",var), colnames(mm),value=TRUE))
                    }
                    
                    to.return = apply(mm[,coefs],1,function(x) sum(x * lmodel$coefficients[coefs]))
                    return(to.return + lmodel$residuals)
                  }
                  
                  # get coefficients for latent.vars.neg
                  
                  
                }
                if (model.use=="poisson") return(residuals(glm(fmla,data = regression.mat,family = "poisson"), type='pearson'))
                if (model.use=="negbinom") return(nb.residuals(fmla, regression.mat, x))
              }))
              if (i==1) data.resid=new.data
              if (i>1) data.resid=rbind(data.resid,new.data)
              setTxtProgressBar(pb, i)
            }
            print(2)
            close(pb)
            rownames(data.resid) <- genes.regress
            if (use.umi) {
              data.resid=sweep(data.resid,MARGIN = 1,apply(data.resid,1,min),"-")
            }
            
            if (!use.umi){
              print("Recalculating object@data")
              
              if (!set.regress.data){
                if (is.matrix(object@data[genes.regress,object@cell.names])){ 
                  object@data[genes.regress,object@cell.names]=data.resid
                  object@data[is.na(object@data)]=0
                } else{
                  
                  if (ncol(object@data) < 35000){
                    object@data = as.matrix(object@data)
                    object@data[genes.regress,object@cell.names]=Matrix(data.resid,sparse=TRUE)
                    object@data[is.na(object@data)]=0
                    object@data = Matrix(object@data, sparse=TRUE)
                  } else {
                    object@regress.data = Matrix(data.resid, sparse=TRUE)
                  }
                  
                }
              } else {
                print("Storing new matrix in object@regress.data")
                
                if (ncol(object@data) < 35000){
                  object@regress.data = as.matrix(object@data)
                  object@regress.data[rownames(data.resid),] = as.matrix(data.resid[,object@cell.names])
                  object@regress.data = Matrix(object@regress.data, sparse=TRUE)
                } else {
                  object@regress.data = Matrix(data.resid, sparse=TRUE)
                }
              }
            } else {
              print("Recalculating object@data. Median re-normalizing (CAUTION).")
              C = as.matrix(object@count.data)
              print(paste0("old median counts was ", median(Matrix::colSums(C)), " UMIs per cell"))
              C[genes.regress, object@cell.names] = data.resid
              C[is.na(C)]=0
              med_trans = median(Matrix::colSums(C))
              print(paste0("New median counts is ", med_trans, " UMIs per cell"))
              TPM.mat = med_trans*t(t(C) / Matrix::colSums(C))
              
              if (is.matrix(object@data[genes.regress,object@cell.names])){ 
                object@data[genes.regress,object@cell.names]=log(TPM.mat+1)
              } else{
                object@data = Matrix(log(TPM.mat+1), sparse=TRUE)
              }
            }
           
            return(object)
          }
)


MNNCorrect <- function(object,genes.regress=NULL,ident.use="orig", 
                       ref.batch=NULL, set.regress.data=TRUE,
                       k.use=40, sigma.use=1.5) {
  
  genes.regress=set.ifnull(genes.regress,rownames(object@data))
  genes.regress=ainb(genes.regress,rownames(object@data))
  if (ident.use %in% colnames(object@data.info)){
    batch_ids = levels(object@data.info[,ident.use])
  } else {
    stop("Error : ident.use is invalid as specified")
  }
  object = set.all.ident(object, id = ident.use)
  
  ref.batch = set.ifnull(ref.batch,batch_ids[1])
  if (!(ref.batch %in% batch_ids)) stop("Error : ref.batch is invalid")
  other.batches = setdiff(batch_ids, ref.batch)
  
  
  Bref = as.matrix(object@data[genes.regress, which.cells(object, ref.batch)])
  l=0
  for (batch in other.batches){
    l=l+1
    B = as.matrix(object@data[genes.regress, which.cells(object, batch)])
    #Xmnn <- mnnCorrect(Bref, B, k=k.use, sigma=sigma.use,  cos.norm=TRUE, svd.dim=2, hvg.genes = rgc@var.genes,inquiry.genes = rgc@var.genes)
    Xmnn <- mnnCorrect(Bref, B, k=k.use, sigma=sigma.use, svd.dim=2, subset.row = genes.regress)
    if (l==1){
      mnnData = cbind(Xmnn$corrected[[1]],Xmnn$corrected[[2]])
      colnames(mnnData) = c(which.cells(object, ref.batch), which.cells(object, batch))
      rownames(mnnData) = genes.regress
    } else {
      mnnData1 = Xmnn$corrected[[2]]
      colnames(mnnData1) = which.cells(object, batch)
      rownames(mnnData1) = genes.regress
      mnnData = cbind(mnnData, mnnData1)
    }
    
    mnnData_mnnzscore = mnnData[,colnames(object@data)]
  }
  # MNN correction
  A = Matrix(mnnData_mnnzscore, sparse=TRUE)
  
  # Renormalize A
  Ae=A
  min_rowvals = apply(Ae,1,min)
  Ae = sweep(Ae,1,min_rowvals,"-")
  norm_factors0 = Matrix::rowSums(Ae)
  norm_factors = Matrix::rowSums(object@data[rownames(A),])
  Ae = t(scale(t(Ae), center=F, scale= norm_factors0 / norm_factors))
  if (set.regress.data){
    object@regress.data = as.matrix(object@data)
    object@regress.data[rownames(Ae),] = as.matrix(Ae)
    object@regress.data = Matrix(object@regress.data, sparse=TRUE)
    return(object)
  } else {
    return(Ae)
  }
  
  
}

#' Finding clusters to subcluster
#' @param min.cells.no.cluster = default (40), minimum number of cells below which subclustering is not done
#' @param cells.tsne to avoid oversensitivity  we downsample the number of points on tSNE space. 
#' # Dip statistic

FindClustersToSubCluster = function(object,
                                    min.cells.no.cluster = 40,
                                    min.cells.for.var.genes = 200,
                                    cells.tsne = 400,
                                    cluster.test = NULL){
      require(diptest)
      # Find clusters that possibly have subclusters
      clust_to_subcluster = c()
      if (is.null(cluster.test)) cluster.test = levels(object@ident)
      for (clust in cluster.test){
        
        print(paste0("Testing cluster ", clust))
        cells.all.clust = which.cells(object,clust)
        if (length(cells.all.clust) < min.cells.no.cluster) next
        cells.clust = sample(cells.all.clust, min(cells.tsne, length(cells.all.clust))) # use a small number of cells to test on tSNE map
        
        if (length(cells.clust) > min.cells.for.var.genes){ 
          # PCA on the tSNE plot
          ydata = as.matrix(object@tsne.rot[cells.clust,])
          pca_obj = prcomp(ydata)
          a_tsne=dip.test(pca_obj$x[,1])
          tsne_pvalue = a_tsne$p.value
          # PCA on its own (use all cells)
          clust.obj = subsetData(object, cells.use = cells.all.clust)
          sink("/dev/null")
          clust.obj@var.genes = NB.var.genes(clust.obj, set.var.genes = FALSE, num.sd = 1, x.high.cutoff = 12, do.idents = TRUE, x.low.cutoff = 50/length(clust.obj@cell.names), rem.mt.rp = TRUE)
          clust.obj <- RegressOut(clust.obj,latent.vars = c( "batch", "nTranscripts","mt.score","rp.score"), genes.regress = clust.obj@var.genes, set.regress.data = TRUE) 
          clust.obj = pca2(clust.obj,pcs.store = 10,pcs.print = 5,genes.print = 5, maxit=300)
          sink()
          pca_pvalue = c()
          for (pc in c(1:10)){
            pca_pvalue = c(pca_pvalue, dip.test(clust.obj@pca.rot[,pc])$p.value)
          }
          
          pca_pvalue = min(pca_pvalue)
        } else {
          pca_pvalue=1
        }
        if (tsne_pvalue < 1e-2 | pca_pvalue < 1e-2){
          print(paste0("Cluster ", clust, " may possibly have substructure. tSNE pval < ", tsne_pvalue, " and PCA pval < ", pca_pvalue))
          clust_to_subcluster = c(clust_to_subcluster, clust)
          
        }
      }
      
      return(clust_to_subcluster)

}


FindSubClusters = function(object,
                           clust_to_subcluster=NULL,
                           batch.corr.method = "RegressOut",
                           cluster.method = "Infomap",
                           thresh.use=0.5,
                           min.de.genes = 3,
                           min.cells.per.cluster=100,
                           do.tsne=FALSE,
                           de.test = "MAST",
                           gene.hi.cutoff.x = 1,
                           num.sd=0.8,
                           latent.vars = c( "batch", "nTranscripts","mt.score","rp.score"),
                           analogizer.batch = "batch",
                           do.sink=FALSE,
                           min.on.off=2,
                           rem.mt.rp=TRUE,
                           genes.exclude=NULL,
                           exclude.genes.prior = NULL,
                           return.obj = FALSE,
                           ...){
  
  if (is.null(clust_to_subcluster)){
    clust_to_subcluster = levels(object@ident)
    cells.use = object@cell.names
  } else {
    cells.use = which.cells(object, clust_to_subcluster)
    if (do.tsne){ plot.cluster(object, cells.plot = cells.use)}
    if (length(cells.use) == 0){
      stop("Error : Provide valid clust_to_subcluster")
    }
  }
  orig.ident = object@ident[cells.use]
  new.ident = as.numeric(orig.ident); names(new.ident) = names(orig.ident)
  new.clust = max(new.ident) + 1
  for (clust in clust_to_subcluster){
    print(paste0("Testing cluster ", clust))
    cells.in.clust = which.cells(object, clust)
    if (length(cells.in.clust) < min.cells.per.cluster){
      print("Fewer than 50 cells in clust. skipping ...")
      next
    }
    clust.obj = subsetData(object, cells.use = cells.in.clust)
    if (do.sink) sink("/dev/null")
    if (length(clust.obj@cell.names) >= 100){
      clust.obj@var.genes = NB.var.genes(clust.obj, set.var.genes = FALSE, num.sd = num.sd,min.cells = 100, x.high.cutoff = gene.hi.cutoff.x, do.idents = TRUE, x.low.cutoff = 25/length(clust.obj@cell.names), rem.mt.rp = rem.mt.rp)
      var.genes = clust.obj@var.genes
    } else {
      clust.obj@var.genes = object@var.genes
    }
    
    if (is.null(clust.obj@var.genes)) stop("Error : Must provide valid variable genes as object@var.genes")
    if (batch.corr.method == "RegressOut"){
      latent.vars = latent.vars[latent.vars %in% colnames(clust.obj@data.info)]
      if (length(levels(drop.levels(clust.obj@data.info$batch))) < 2 ){
        
        clust.obj <- RegressOut(clust.obj,latent.vars = setdiff(latent.vars,"batch"), genes.regress = clust.obj@var.genes, set.regress.data = TRUE,exclude.genes.prior = exclude.genes.prior) 
        
      } else {
        print(1)
        clust.obj <- RegressOut(clust.obj,latent.vars = latent.vars, genes.regress = clust.obj@var.genes, set.regress.data = TRUE, exclude.genes.prior = exclude.genes.prior) 
        print(2)
      }
      clust.obj =pca2(clust.obj,pcs.store = 50,pcs.print = 5,genes.print = 5, maxit=300, use.regress = TRUE)
      n_K = KEst_TracyWidom(clust.obj)
      if (do.tsne){ clust.obj  = run_tsne(clust.obj, pcs.use = c(1:n_K), do.fast = TRUE); tsplot(clust.obj) }
      do.regress.for.dendro=FALSE
      do.regress = FALSE
    } else if (batch.corr.method == "Analogizer"){
      clust.obj = set.all.ident(clust.obj, id=analogizer.batch)
      dge_list = list()
      l=1
      for (i in levels(clust.obj@ident)){
        dge_list[[paste0("name",l)]] = clust.obj@count.data[rownames(clust.obj@data),which.cells(clust.obj,i)]
        l=l+1
      }
      analogy = Analogizer(dge_list) #Can also pass in more than 2 datasets
      analogy = Analogizer::normalize(analogy)
      analogy@var.genes = clust.obj@var.genes
      clust.obj =pca2(clust.obj,pcs.store = 50,pcs.print = 5,genes.print = 5, maxit=300, use.regress = FALSE)
      n_K = KEst_TracyWidom(clust.obj)
      analogy = Analogizer::scaleNotCenter(analogy)
      analogy = optimizeALS(analogy,k=n_K) 
      analogy = quantile_align_SNF(analogy) #SNF clustering and quantile alignment
      analogy = run_tSNE(analogy)
      clust.obj@tsne.rot[,1] = analogy@tsne.coords[clust.obj@cell.names,1]
      clust.obj@tsne.rot[,2] = analogy@tsne.coords[clust.obj@cell.names,2]
      clust.obj@dr$tsne@cell.embeddings[,1] = analogy@tsne.coords[clust.obj@cell.names,1]
      clust.obj@dr$tsne@cell.embeddings[,2] = analogy@tsne.coords[clust.obj@cell.names,2]
      clust.obj@dr$tsne@key = "tsne_"
      clust.obj@pca.rot = as.data.frame(analogy@H.norm)
      colnames(clust.obj@pca.rot) = paste0("PC",c(1:ncol(analogy@H.norm)))
      clust.obj@dr$pca@cell.embeddings = as.matrix(clust.obj@pca.rot)
      do.regress.for.dendro=FALSE
      do.regress = FALSE
      
    } else {
      do.regress.for.dendro=FALSE
      do.regress = FALSE
      clust.obj =pca2(clust.obj,pcs.store = 50,pcs.print = 5,genes.print = 5, maxit=300, use.regress = FALSE)
      n_K = KEst_TracyWidom(clust.obj)
      if (do.tsne){ clust.obj  = run_tsne(clust.obj, pcs.use = c(1:n_K), do.fast = TRUE); tsplot(clust.obj) }
    }
    
    
    if (n_K <= 2) next
   
    do.jaccard = FALSE
    if (cluster.method == "Louvain") do.jaccard = TRUE
    clust.obj = doGraph_clustering(clust.obj, pcs.use = 1:max(n_K,2), num.nn=round(log2(length(clust.obj@cell.names)*2)), method=cluster.method, do.jaccard = do.jaccard)
    
    if (length(clust.obj@ident) > 1){
      clust.obj = buildClusterTree(clust.obj, regress.use = do.regress, genes.use = clust.obj@var.genes, ...)
      clust.obj = MergeClustersByDendro(clust.obj, var.genes = clust.obj@var.genes, thresh.use = thresh.use,
                                        min.de.genes = min.de.genes, test.use = de.test, 
                                        remove.mito_ribo = TRUE,min_nTrans_ratio = 2, regress.use.for.dendro = do.regress.for.dendro,
                                        regress.use.for.de = FALSE,min.on.off = min.on.off, do.tsne=do.tsne, genes.exclude=genes.exclude, ...)
      if (do.tsne){ tsplot(clust.obj) }
    }
    if (do.sink) sink()
    if (length(levels(clust.obj@ident)) > 1){
      print(paste0("Old Cluster ", clust, " has ", length(levels(clust.obj@ident)), " subclusters. Changing identity ... "))
      l = 0
      for (subclust in levels(clust.obj@ident)){
        if (l==0){
          l=1;
          next
        }
        cells.subclust = which.cells(clust.obj, subclust)
        new.ident[cells.subclust] = new.clust;
        new.clust = new.clust + 1
      }
    }
    
  }
  
  if (return.obj) return(clust.obj)
  
  return(new.ident)
}

Ligerize = function(object,batch_id="batch", var.genes = NULL,  do.clustering=TRUE, max.pcs=90, do.umap=FALSE, ...){
  
  library(liger)
  object = set.all.ident(object, id=batch_id)
  if (is.null(var.genes)){
      var.genes_each_no_mt = NB.var.genes(object, set.var.genes = FALSE, num.sd = 1, x.high.cutoff = 3, do.idents = TRUE, x.low.cutoff = 0.005, rem.mt.rp = TRUE, do.text = FALSE, ...)
      object@var.genes = var.genes_each_no_mt
    
  } else {
    object@var.genes = var.genes
  }
  
  object=pca2(object,pcs.store = max.pcs,pcs.print = 5,genes.print = 5, maxit=300, use.regress = FALSE) # regressed
  n_K = KEst_TracyWidom(object)
  if (n_K == max.pcs){print("Warning : Re-run by increasing max PCs")}
  object  = run_tsne(object, pcs.use = c(1:n_K), do.fast = TRUE, fast.method="FFT")
  object@dr$tsne@key = "tsne_"
  tsplot(object)
  orig_pca = object@pca.rot
  print(1)
  dge_list = list()
  l=1
  for (i in levels(object@data.info[,batch_id])){
    dge_list[[paste0("name",l)]] = object@count.data[rownames(object@data),which.cells(object,i)]
    l=l+1
  }
  print(2)
  analogy = createLiger(dge_list) #Can also pass in more than 2 datasets
  analogy = liger::normalize(analogy)
  analogy <- liger::selectGenes(analogy)
  if(!is.null(var.genes)){
    analogy@var.genes = intersect(analogy@var.genes, var.genes)
   }
  analogy = liger::scaleNotCenter(analogy)
  analogy = optimizeALS(analogy,k=n_K+3) 
  analogy = quantileAlignSNF(analogy) #SNF clustering and quantile alignment
  object@var.genes = analogy@var.genes
 
  analogy = runTSNE(analogy, dims.use = 1:n_K)
  orig_tsne = object@tsne.rot
  object@tsne.rot[,1] = analogy@tsne.coords[object@cell.names,1]
  object@tsne.rot[,2] = analogy@tsne.coords[object@cell.names,2]
    object@dr$tsne@cell.embeddings[,1] = analogy@tsne.coords[object@cell.names,1]
    object@dr$tsne@cell.embeddings[,2] = analogy@tsne.coords[object@cell.names,2]
    object@dr$tsne@key = "tsne"
  
  if (do.umap){  
    analogy = runUMAP(analogy, dims.use = 1:n_K)
    colnames(analogy@tsne.coords) = c("UMAP1","UMAP2")
    diff.obj <- new(Class = "dim.reduction",
                    gene.loadings = matrix(rand(10),2,5),
                    cell.embeddings = as.matrix(analogy@tsne.coords),
                    dim.var = c(1,2,3),
                    key = "UMAP")
    object@dr$UMAP = diff.obj
  }
  
  
  object@pca.rot = as.data.frame(analogy@H.norm)
  colnames(object@pca.rot) = paste0("PC",c(1:ncol(analogy@H.norm)))
  object@dr$pca@cell.embeddings = as.matrix(object@pca.rot)
  if (do.clustering){
    object = doGraph_clustering(object, pcs.use = 1:n_K, num.nn=25, method="Louvain", do.jaccard=TRUE, use.reduction = "pca")
  }
  return(object)
}
