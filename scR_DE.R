# Suite of functions for evaluating differential expression
library(pbapply)
require(spdep)

#' Differential Expression analysis between two clusters, uses many tests
#' @param object scR object
#' @param ident.1 = Required, identifier for first cluster
#' @param ident.2 = Optional, identifier for second cluster (if NULL, ident.1 is compared with the rest)
#' @param cells.1 = Optional, cells in first group, if provided, ident.1 is overriden
#' @param cells.2 = Optional, cells in second group, if provided, ident.2 is overriden
#' @param genes.use = Optional, set of genes to be considerered (if NULL, all genes are used)
#' @param thresh.use = Optional, minimum average log-fold expression difference between the two clusters for a gene to be considered
#' @param test.use = Optional, which test to use (options are)
#' wilcox : Wilcoxon rank sum test (default)
#' bimod : Likelihood ratio test for single cell gene expression (a la McDavid et al., Bioinformatics, 2013)
#' roc : Standard ROC classifier
#' t : Student's t-test
#' tobit : Tobit-test for differential expression (Trapnell et al., Nat. Biotech, 2014)
#' poisson : Likelihood ratio test assuming an underlying poisson distribution. Use only for UMI-based data
#' negbinom : Likelihood ratio test assuming an underlying negative binomial distribution. Use only for UMI-based data
#' MAST : GLM-framework that treats cellular detection rate as a covariate (Finak et al., Genome Biology, 2015)
#' Tweede-seq
#' binom.digit : Binomial test - good for very lowly expressed genes. Will miss genes that change in an analog fashion
#' @param min.perc = Optional, only genes that are expressed in at least 100*min.pos.perc % cells in either cluster
#' @param min.diff = fraction<1, Optional, minimum difference in the proportion of expressing cells for the gene to be considered between ident.1 and ident.1
#' set to -1 by default
#' @param only.pos = Only return positive markers (FALSE by default)
#' @param test.max = integer, Optional, test only the first test.max genes sorted by effect size
#' @param TPM.mat = matrix, optional, matrix/data.frame of normalized counts. Only relevant for binom.digit test. Speeds up calculations
#' @param verbose = Print a progress bar (default = TRUE)  
#' @param max.cells.per.ident = downsample each identity class to a max number. default is no downsampling
#' @param random.seed = random seed for downsampling (default 42)
#' @param latent.vars = Variables to test
#' @param min.cells - Minimum number of cells expressing the gene in at least one of the two groups
#' @param pseudocount.use Pseudocount to add to averaged expression values when
#' @param pval.thresh - p-value threshold
#' @param min.exp.thresh - value above which a gene is considered expressed
#' @param append.all.pct - Boolean, if multiple clusters are included in any group then the percentage of expressing cells in each cluster are appended
#' @return.genes - if TRUE, DE is skipped and genes that pass the pre-filtering are returned (useful for GO analysis)
#'@param \dots Additional parameters to pass to specific DE functions
#' calculating logFC. 1 by default.
#' TO DO : Downsampling

find.markers =  function(
  object, 
  ident.1=NULL,
  ident.2=NULL,
  cells.1=NULL,
  cells.2=NULL,
  ident.use=NULL,
  genes.use=NULL,
  thresh.use=0.25,
  test.use="wilcox", 
  min.perc=0.1, 
  min.diff=-1, 
  only.pos=FALSE, 
  test.max=NULL, 
  TPM.mat=NULL, 
  verbose=TRUE,
  max.cells.per.ident = Inf,
  random.seed = 42,
  latent.vars = "nTranscripts",
  min.cells = 10,
  pseudocount.use = 0.3,
  pval.thresh=1e-2,
  min.exp.thresh = NULL,
  append.all.pct = TRUE,
  return.genes = FALSE,
  min.genes = 2,
  use.regress.data = FALSE,
  ...
  ){

            
            genes.use=set.ifnull(genes.use,rownames(object@data))
            if (use.regress.data){
              data.use = object@regress.data[genes.use,]
            } else {
              data.use = object@data[genes.use,]
            }
            
            if (!is.null(ident.use)){
              if (ident.use %in% colnames(object@data.info)){
                 orig.idents = object@ident
                 object = set.all.ident(object, id = ident.use)
                 ident.use = object@ident
                 object@ident = orig.idents
              } else {
                print("Error : ident.use is not valid")
                return(NULL)
              } 
              } else {
                ident.use=object@ident
              }
            
            # Find the cells to do DE
            if (is.null(ident.1)){
              ident.1 = "clust1"
            }
            cells.1=set.ifnull(cells.1,names(ident.use[which(ident.use%in%ident.1)]))
            cells.1 = ainb(cells.1, object@cell.names)
            
            if (is.null(ident.2)) {
              if (is.null(cells.2)){ 
                ident.2="rest"
                cells.2=set.ifnull(cells.2,names(ident.use))
                cells.2=cells.2[!(cells.2%in%cells.1)]
              } else {
                ident.2 = "clust2"
                cells.2 = ainb(cells.2, object@cell.names)
              }
            } else {
              cells.2=set.ifnull(cells.2,names(ident.use[which(ident.use%in%ident.2)]))
            }

            # Error checking
            if (length(cells.1)==0) {
              cat(paste("Error:Clusters ", ident.1, "vs", ident.2, "Error: Cell group 1 is empty - no cells with identity class", ident.1,"\n"),file=stderr())
              return(NULL)
            }
            if (length(cells.2)==0) {
              cat(paste("Error:Clusters ", ident.1, "vs", ident.2, "Error:Cell group 2 is empty - no cells with identity class", ident.2,"\n"),file=stderr())
              return(NULL)
            }
            
            if (length(cells.1) < min.cells) {
              cat(paste("Error:Clusters ", ident.1, "vs", ident.2, ":Cell group 1 has fewer than", as.character(min.cells), "cells in identity class", ident.1,"\n"),file=stderr())
              return(NULL)
            }
            if (length(cells.2) < min.cells) {
              cat(paste("Error:Clusters ", ident.1, "vs", ident.2, ": Cell group 2 has fewer than", as.character(min.cells), " cells in identity class", ident.2,"\n"),file=stderr())
              return(NULL)
            }
            
            if (verbose){
              print(paste0("DE : Cluster ", ident.1, " vs Cluster ", ident.2," : Group 1 has ", length(cells.1), " cells"))
              print(paste0("DE : Cluster ", ident.1, " vs Cluster ", ident.2," : Group 2 has ", length(cells.2), " cells"))
              if (length(intersect(cells.1, cells.2)) > 0){
                print(paste0("WARNING: There are ",length(intersect(cells.1, cells.2)), " common cells between the two groups." ))
              }
            }
            
            print(paste0("number of cells in group 2 is ", length(cells.2)))
            # Downsample?
            if (max.cells.per.ident < Inf) {
              if (verbose) print(paste0("Downsampling each cluster to ", max.cells.per.ident, " cells"))
              set.seed(seed = random.seed)
              if (length(cells.1) > max.cells.per.ident) cells.1 = sample(x = cells.1, size = max.cells.per.ident)
              if (length(cells.2) > max.cells.per.ident) cells.2 = sample(x = cells.2, size = max.cells.per.ident)
            }
            
        
            # Gene selection (based on percent expressed)
            
            min.exp.thresh=set.ifnull(min.exp.thresh,object@is.expr)
            if (verbose) print(paste0("Considering genes expressed in at least ", 100*min.perc, " % of cluster, and at min.diff = ", min.diff))
            if (verbose) print(paste0("Genes with expression value greater than ", min.exp.thresh, " are considered expressed"))
            
            data.use = data.use[,c(cells.1,cells.2)]
            data.temp1 = round(Matrix::rowSums(data.use[genes.use, cells.1] > min.exp.thresh) / length(cells.1),3)
            data.temp2 = round(Matrix::rowSums(data.use[genes.use, cells.2] > min.exp.thresh) / length(cells.2),3)
            #data.temp1=round(apply(data.use[genes.use,cells.1],1,function(x)return(length(x[x>min.exp.thresh])/length(x))),3)
            #data.temp2=round(apply(data.use[genes.use,cells.2],1,function(x)return(length(x[x>min.exp.thresh])/length(x))),3)
            
            
            data.pct=cbind(data.temp1,data.temp2); colnames(data.pct)=c("pct.1","pct.2")
            pct.min = apply(data.pct,1,max); names(pct.min)=rownames(data.pct);
            genes.use <- names(which(pct.min > min.perc))
            
            if (length(genes.use) < min.genes){ 
              cat(paste0("Error:Clusters ", ident.1, " vs ", ident.2,": No genes passed min.perc threshold\n"),file=stderr()); return(NULL)
            }
            
            pct.diff = pct.min - apply(data.pct,1,min); 
            genes.use <- names(which(pct.min > min.perc & pct.diff > min.diff))
            
            if (length(genes.use) < min.genes){ 
              cat("No genes passed min.diff threshold",file=stderr()); return(NULL)
            }
            
            # Gene selection (based on average difference) 
            data.1=apply(data.use[genes.use,cells.1],1,function(x) log(mean(expm1(x)) + pseudocount.use))
            data.2=apply(data.use[genes.use,cells.2],1,function(x) log(mean(expm1(x)) + pseudocount.use))
            total.diff=(data.1-data.2)
            
            if (!only.pos){ 
              genes.use <- names(which(abs(total.diff) >= thresh.use))
              
              if (length(genes.use) < min.genes) {
                cat(paste0("Error:No genes pass expression threshold, thresh.use = ", thresh.use,"\n"),file=stderr())
                return(NULL)
              }
              
              if (verbose){
                print(paste0("Considering ", length(genes.use), " for DE"))
              }
              
            }
            
            if (only.pos){
              genes.use <- names(which(total.diff >= thresh.use))
              
              if (length(genes.use) < min.genes) {
                cat(paste0("Error:No genes pass expression threshold, thresh.use = ", thresh.use, " in Cluster ", ident.1,"\n"),file=stderr())
                return(NULL)
              }
              if (verbose){
                print(paste0("Considering ", length(genes.use), " genes that are enriched in Cluster ", ident.1, " for DE"))
              }
             
            }
            
            if (return.genes){
              if (verbose) print("Returning genes that pass the prefilters without performing DE")
              return(genes.use)
            }
            

          
            # Perform DE
            tested = F
            # Wilcoxon Rank Sum Test (also Mann-Whitney Test)
            if (test.use == "wilcox") {
              if (verbose) print("Using Wilcoxon Rank Sum test")
              to.return <- WilcoxDETest(
                object = object,
                cells.1 = cells.1,
                cells.2 = cells.2,
                genes.use = genes.use,
                verbose = verbose,
                ...
              )
              tested = T
            }
            
            if (test.use == "bimod"){
              if (verbose) print("Using Bimod test")
              to.return=DiffExpTest(
              object = object,
              cells.1 = cells.1,
              cells.2 = cells.2,
              genes.use = genes.use,
              verbose = verbose,
              xmin = min.exp.thresh
              ) 
              tested = T
            }
            
            if (test.use == "tobit") {
              if (verbose) print("Using Tobit test")
              to.return <- TobitTest(
                object = object,
                cells.1 = cells.1,
                cells.2 = cells.2,
                genes.use = genes.use,
                verbose = verbose
              )
              tested = T
            }
            
            if (test.use == "binom.digit"){ 
              if (verbose) print("Using binomial count test (ON/OFF)")
              
              if (is.null(TPM.mat)){
                print("Consider providing a precomputed TPM.mat for speeding up")
                data1 = as.matrix(data.use)
                TPM.mat = exp(data1)-1
                TPM.mat = TPM.mat[genes.use, c(cells.1, cells.2)]
              } else {
                TPM.mat = TPM.mat[genes.use, c(cells.1, cells.2)]
              }
              
              to.return <- BinomTest(
                object = object, 
                cells.1 = cells.1, 
                cells.2 = cells.2, 
                genes.use = genes.use, 
                TPM.mat = TPM.mat, 
                xmin = min.exp.thresh,
                ...)
              tested = T
            }
            
            if (test.use == "t") {
              if (verbose) print("Using Student's t-test")
              to.return <- DiffTTest(
                object = object,
                cells.1 = cells.1,
                cells.2 = cells.2,
                genes.use = genes.use,
                verbose = verbose
              )
              tested = T
            }
            
            if (test.use == "roc") {
              if (verbose) print("Using ROC test")
              to.return <- ROCTest(
                object = object,
                cells.1 = cells.1,
                cells.2 = cells.2,
                genes.use = genes.use,
                verbose = verbose
              )
              tested = T
            }
            
            if (test.use == "MAST") {
              require(MAST)
              if (verbose) print("Using MAST test")
              to.return <- MASTDETest(
                object = object,
                cells.1 = cells.1,
                cells.2 = cells.2,
                genes.use = genes.use,
                latent.vars = latent.vars
              )
              tested = T
            }
            
            if (test.use == "edgeRQLF") {
              require(edgeR)
              if (verbose) print("Using edgeR QLF test")
              to.return <- edgeRQLFTest(
                object = object,
                cells.1 = cells.1,
                cells.2 = cells.2,
                genes.use = genes.use
              )
              tested = T
            }
            
            if (test.use == "tweeDEseq") {
              require(tweeDEseq)
              if (verbose) print("Using tweeDE seq test")
              to.return <- tweeDEseqTest(
                object = object,
                cells.1 = cells.1,
                cells.2 = cells.2,
                genes.use = genes.use
              )
              tested = T
            }
            
            if (!tested){
              cat("Error : DE test not available. Please select a valid option\n", file=stderr())
              return(NULL)
              
            }
            
            # Deprecated
            #if (test.use=="scde") to.return=scde.test(object,cells.1,cells.2,paste(ident.1,collapse="_"), paste(ident.2,collapse="_"),genes.use,thresh.use, min.diff,test.max,...) 
            
            # if (test.use=="negbinom"){ 
            #   if (verbose) print("Using NB regression test (use only with UMI data)")
            #   if (is.null(Count.mat)){
            #     if (verbose) print("Consider providing a precomputed Count.mat for speeding up")
            #     Count.mat = set.ifnull(Count.mat, object@count.data[genes.use,c(cells.1, cells.2)])
            #   }
            #   to.return=NegBinomDETest(object,cells.1,cells.2,genes.use,Count.mat[genes.use,])
            #  if (verbose) print(head(to.return))
            # }
            # 
            
            
            to.return[, "diff"] <- total.diff[rownames(to.return)]
            df1 = data.pct[rownames(to.return), ]
            colnames(df1) = c(paste0("pct_clust", 1), paste0("pct_clust", 2))
            to.return <- cbind(to.return, df1)
            to.return[, "test"] <- test.use
            
            if (test.use == "roc"){
              to.return <- to.return[order(-to.return$power,-to.return$diff), ]
            } else {
              to.return$pval_adj = p.adjust(p = to.return$pval,method = "bonferroni")
              to.return <- to.return[order(-to.return$diff,to.return$pval_adj), ]
            }
            
            
            if ((length(ident.1) > 1 | length(ident.2) > 1) & append.all.pct){
              mat=Perc.pos.by.ident(object, features.use=genes.use,ident.use = c(ident.1,ident.2), do.plot = FALSE, thresh.use = min.exp.thresh)
              colnames(mat$PercMat) = paste0("pct_", colnames(mat$PercMat))
              to.return <- cbind(to.return, round(mat$PercMat[rownames(to.return),], digits = 3))
            }
            
            #Mean number of transcripts per cell
            if (!is.null(attr(object,"count.data"))){
              if (length(ident.1) == 1){
              nTrans.1 = apply(object@count.data[rownames(to.return), cells.1], 1, function(x) round(mean(x),3))
              to.return[,paste0("nTrans_", ident.1)] = nTrans.1
              } else {
                for (clust in ident.1){
                  nTrans.1 = apply(object@count.data[rownames(to.return), which.cells(object,clust)], 1, function(x) round(mean(x),3))
                  to.return[,paste0("nTrans_", clust)] = nTrans.1
                }
              }
              
              if (length(ident.2) == 1){
                nTrans.2 = apply(object@count.data[rownames(to.return), cells.2], 1, function(x) round(mean(x),3))
                to.return[,paste0("nTrans_", ident.2)] = nTrans.2
              } else {
                for (clust in ident.2){
                  nTrans.2 = apply(object@count.data[rownames(to.return), which.cells(object,clust)], 1, function(x) round(mean(x),3))
                  to.return[,paste0("nTrans_", clust)] = nTrans.2
                }
              }
            }
            
            # if (test.use == "negbinom"){
            #    to.return$diff = total.diff[rownames(to.return)]
            # }
            
            if (only.pos) {
              to.return <- subset(to.return, subset = diff > 0)
            }
            
            if (test.use == "roc"){
              if (pval.thresh < 0.1) pval.thresh = 0.7
              to.return=subset(to.return, myAUC >= pval.thresh | myAUC <= 1-pval.thresh)     
            } else {
              to.return=subset(to.return,pval_adj <= pval.thresh)
            }
            
            if (verbose) print(paste0("Found ", nrow(to.return), " genes with p-value < ", pval.thresh))
            if (verbose) print(head(to.return[,c(1:5)]))
            
            return(to.return)
} 



#' Gene expression markers for all identity classes
#'
#' Finds markers (differentially expressed genes) for each of the identity classes in a dataset
#'
#' @inheritParams FindMarkers
#' @param print.bar Print a progress bar once expression testing begins (uses pbapply to do this)
#' @param max.cells.per.ident Down sample each identity class to a max number. Default is no downsampling.
#' @param random.seed Random seed for downsampling
#' @param return.thresh Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)
#' @param do.print FALSE by default. If TRUE, outputs updates on progress.
#' @param min.cells Minimum number of cells expressing the gene in at least one of the two groups
#' @param latent.vars remove the effects of these variables
#' @param assay.type Type of assay to perform DE for (default is RNA)
#' @param \dots Additional parameters to pass to specific DE functions
#'
#' @return Matrix containing a ranked list of putative markers, and associated
#' statistics (p-values, ROC score, etc.)
#'
#' @export
#' @examples
#' all_markers <- FindAllMarkers(object = pbmc_small)
#' head(x = all_markers)
#'
find.all.markers <- function(
    object, 
    thresh.use=0.25,
    genes.use = NULL,
    idents.use = NULL,
    test.use="wilcox",
    min.perc = 0.1,
    min.diff = -1,
    TPM.mat=NULL, 
    Count.mat=NULL,
    verbose = TRUE,
    do.print = TRUE,
    only.pos = TRUE,
    max.cells.per.ident = Inf,
    return.thresh = 1e-2,
    random.seed = 42,
    pseudocount.use = 1,
    min.cells = 10,
    latent.vars = "nTranscripts",
    min.exp.thresh = NULL,
    add.to.slot=TRUE,
    ...
){
  
  genes.use=set.ifnull(genes.use,rownames(object@data))
  data.use = object@data[genes.use,]
  pval.thresh = set.ifnull(return.thresh, 1e-2)
  
  if ((test.use == "roc") && (pval.thresh < 0.1)) {
    pval.thresh = 0.7
  }
  
  idents.use = set.ifnull(idents.use, levels(object@ident))
  if (sum(!(idents.use %in% levels(object@ident))) != 0){
    cat(paste0("Warning : Group(s) ", idents.use[!(idents.use %in% levels(object@ident))], 
               " is/are not present. Check!! \n" ), file=stderr())
  }
  
  genes.de <- list()

  
  if (test.use == "binom.digit" & is.null(TPM.mat)){
    TPM.mat = exp(object@data) - 1
  }
  
  
  for (i in 1:length(x = idents.use)) {
    if (verbose) {
      print(paste("Calculating DE genes for cluster", idents.use[i]))
    }
    genes.de[[i]] <- tryCatch(
      {
        find.markers(
          object = object,
          ident.1 = idents.use[i],
          ident.2 = NULL,
          genes.use = genes.use,
          thresh.use = thresh.use,
          test.use = test.use,
          min.perc = min.perc,
          min.diff = min.diff,
          only.pos = only.pos,
          verbose = verbose,
          TPM.mat = TPM.mat,
          max.cells.per.ident = max.cells.per.ident,
          random.seed = 42,
          latent.vars = latent.vars,
          min.cells = min.cells,
          pseudocount.use = pseudocount.use,
          pval.thresh = pval.thresh,
          min.exp.thresh = min.exp.thresh,
          ...
        )
      },
      error = function(cond){
        return(NULL)
      }
    )

  }
  
  gde.all <- data.frame()
  #for (i in 1:length(idents.use)) {
  for (i in 1:length(genes.de)) {
    if (is.null(genes.de[[i]])) {
      next
    }
    gde <- genes.de[[i]]
    if (nrow(x = gde) > 0) {
      if (test.use == "roc") {
        gde <- subset(
          x = gde,
          subset = (power > return.thresh | power < (1 - return.thresh))
        )
      } else {
        gde <- gde[order(gde$pval, -gde$diff), ]
        gde <- subset(x = gde, subset = pval_adj < return.thresh)
      }
      if (nrow(gde) > 0) {
        gde$cluster <- idents.use[i]
        gde$gene <- rownames(gde)
      }
      if (nrow(x = gde) > 0) {
        colnames(gde)[7] = "nTrans_clust"
        gde.all <- rbind(gde.all, gde)
      }
    }
  }
  
  if (only.pos) {
   gde.all = subset(x = gde.all, subset = diff > 0)
  }
  rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
  
  if (add.to.slot){
    testid = paste0(test.use,"_","thresh",thresh.use,"minperc",min.perc,"_DE")
    if (testid %in% names(object@de.list)){
      new.names = make.unique(c(names(object@de.list), testid),sep = "_")
      print(paste0("Adding results to the slot @de with the name ", testid))
      object@de.list[[new.names[length(new.names)]]] = gde.all
    } else {
      object@de.list[[testid]] = gde.all
    }
    return(object)
  } else {
    
    return(gde.all)
  }
 
  

}

#' Find markers corresponding to a specific node in the dendrogram
#' @param object scR object
#' @param node one of the nodes in the scR dendrogram
#' @param genes.use list of genes to perform DE on
#' @param thresh.use  
#' @param test.use
#' @param min.perc
#' @param verbose

find.markers.node = function(
  object,
  node,
  genes.use=NULL,
  thresh.use=0.5,
  test.use="wilcox", 
  min.perc=0.3, 
  verbose=TRUE, 
  only.pos=FALSE,
  pval.thresh=1e-2,
  max.cells.per.ident=Inf,
  append.all.pct = TRUE,
  ...) {
  
  
  genes.use=set.ifnull(genes.use,rownames(object@data))
  tree=object@cluster.tree[[1]]
  ident.order=tree$tip.label
  
  nodes.1=ident.order[getLeftDecendants(tree,node)]
  nodes.2=ident.order[getRightDecendants(tree,node)]
  
  if (verbose){
    print(paste0("Node ", node, " splits cluster(s) ", paste0(nodes.1, collapse=","), " vs cluster(s) ", paste0(nodes.2, collapse=",")))
    print(paste0("Performing ", test.use, " differential expression test"))
  }
  #print(nodes.1)
  #print(nodes.2)
  to.return=find.markers(object,
                         ident.1 = nodes.1,
                         ident.2=nodes.2,
                         genes.use=genes.use,
                         thresh.use=thresh.use,
                         test.use=test.use, 
                         min.perc = min.perc,
                         only.pos = only.pos,
                         verbose = FALSE,
                         pval.thresh = pval.thresh,
                         max.cells.per.ident = max.cells.per.ident,
                         append.all.pct = append.all.pct,
                         ...)
  
  if (verbose){
    print(head(to.return[,c(1:5)]))
  }
  return(to.return)
  
} 

#' Find markers for all nodes in the dendrogram below a particular node
#' Requires the @cluster.tree slot filled. Else use BuildClusterTree

find.all.markers.node = function(
  object,
  node,
  thresh.use=0.5,
  test.use="wilcox",
  genes.use=NULL,
  min.perc=0.3,
  min.diff=0.1,
  verbose=TRUE,
  only.pos=FALSE,
  max.cells.per.ident=Inf,
  return.thresh=1e-2,
  random.seed=42,
  min.cells=10,
  append.all.pct=TRUE,
  add.to.slot=TRUE,
  ...) {

  if (length(object@cluster.tree) == 0) {
    stop("Tree hasn't been built yet. Run BuildClusterTree to build.")
  }
  genes.use = set.ifnull(genes.use, rownames(object@data))
  tree.use=object@cluster.tree[[1]]
  node <- set.ifnull(node, tree.use$edge[1,1])
  
  print("Check to make sure that object@ident is consistent with the nodes. Else use set.all.ident to fix")
  ident.use=object@ident
  
  node.descendants <- getDescendants(tree.use,node)
  node.descendants <- node.descendants[node.descendants > tree.use$Nnode+1]
  node.descendants <- c(node, node.descendants)
  
  
  if ((test.use=="roc") && (return.thresh < 1e-2)) return.thresh=0.7
  
  #if (verbose) print(paste0("Finding markers that define the nodes ", new.nodes))
  
  genes.de=list()

  for(i in node.descendants) {
    if (verbose) print(paste0("Computing markers for Node ", i))
    genes.de[[i]]=find.markers.node(
      object,
      i,
      genes.use=genes.use,
      thresh.use = thresh.use,
      test.use = test.use,
      min.perc=min.perc,
      verbose=verbose,
      only.pos=only.pos,
      pval.thresh=return.thresh,
      max.cells.per.ident = max.cells.per.ident,
      append.all.pct = TRUE,
      random.seed = random.seed,
      ...)
    }

  gde.all=data.frame()

  for(i in ((tree.use$Nnode+2):max(tree.use$edge))) {
    if (is.null(x = unlist(genes.de[i]))) {
      next
    }
    gde=genes.de[[i]]

    if (nrow(gde)>0) {
      if (test.use=="roc"){
        gde=subset(gde,(myAUC>return.thresh | myAUC<(1-return.thresh)))
      } else {
        gde=gde[order(gde$pval,-gde$diff),]
        gde=subset(gde,pval<return.thresh)
      }

      colnames(gde)[c(3,4)] = c("pct_left", "pct_right")
      gde=cbind(rownames(gde), gde)
      gde = cbind(i, gde)
      colnames(gde)[c(1,2)] = c("node", "gene")
      
      left_desc = tree.use$tip.label[getLeftDecendants(tree.use,i)]
      left_desc = paste0(left_desc,collapse=",")
      right_desc = tree.use$tip.label[getRightDecendants(tree.use,i)]
      right_desc = paste0(right_desc,collapse=",")

      gde$left_desc = left_desc
      gde$right_desc = right_desc
      gde = gde[,c("node","gene", "diff","pval","pval_adj","pct_left","pct_right","test","left_desc","right_desc")]
      gde.all=rbind(gde.all,gde)
    }
    
    
    
  }
  
  rownames(gde.all) <- make.unique(names = as.character(gde.all$gene))
  
  if (add.to.slot){
    testid = paste0("NodeTest", node,"_",test.use,"_","thresh",thresh.use,"minperc",min.perc,"_DE")
    if (testid %in% names(object@de.list)){
      new.names = make.unique(c(names(object@de.list), testid),sep = "_")
      print(paste0("Adding results to the slot @de with the name ", testid))
      object@de.list[[new.names[length(new.names)]]] = gde.all
    } else {
      object@de.list[[testid]] = gde.all
    }
    return(object)
  } else {
    
    return(gde.all)
  }
}

#' Merge clusters based on the dendrogram
#' Iteratively checks the closest clusters and merges them if they don't show sufficient differential expression
#' @param num.nodes.test Number of internal nodes to test at a time
#' @param remove.mito_ribo TRUE/FALSE, whether mitochondrial and ribosomal genes should be removed from the list of valid DE genes

MergeClustersByDendro = function(
  object,
  var.genes=NULL, 
  min.de.genes=5, 
  num.nodes.test=3,
  thresh.use=0.3,
  test.use='wilcox', 
  pval.thresh=1e-2, 
  min.perc=0.2, 
  ident.use="m",
  stat.fxn = expMean,
  remove.mito_ribo=TRUE,
  linkage.method = "complete", 
  dist.fun = "euclidean",
  min_nTrans_ratio = 2,
  min_perc_ratio=1,
  verbose=TRUE,
  regress.use.for.dendro=FALSE,
  regress.use.for.de=FALSE,
  do.scale.dendro=FALSE,
  merge_id = "m_merge",
  min.on.off = 2,
  do.tsne=FALSE,
  genes.exclude = NULL,
  ...) {
  
  
  if (!(ident.use %in% colnames(object@data.info))) stop("Error : ident.use is not valid")
  object = set.all.ident(object, id = ident.use)
  # check if cluster tree is present
  if (length(object@cluster.tree) == 0){
    cat("Building dendrogram since it is not provided!!!\n", file=stderr())
    
    if (is.null(var.genes)){
      cat("Using all the genes to build the dendrogram since var.genes is not specified. This might be slow!!!\n", file=stderr())
      var.genes = object@var.genes
    }
    if (verbose) print(paste0("Building dendrogram using N=", length(var.genes), " genes, linkage = ", linkage.method, " and ", dist.fun, " distance"))
    object = buildClusterTree(object, genes.use=var.genes, stat.fxn = stat.fxn, linkage.method = linkage.method, dist.fun = dist.fun, regress.use = regress.use.for.dendro, do.scale=do.scale.dendro)
  }
  
  min.de.genes = set.ifnull(min.de.genes, 5)
  if (verbose) print(paste0("Using DE test ", test.use))
  if (verbose) print(paste0("Requiring at least ", min.de.genes, " genes to be enriched on each side, at log-fold change = ", 
                thresh.use,
                " and pval < ",
                pval.thresh,
                " to avoid merging"))

  tree.use =  object@cluster.tree[[1]]
  
  orig_ident = object@ident
  l = -1 # Indicator variable that will switch to +1 if no nodes are merged
  iter = 1
  while(l == -1){
    
    if (verbose) print(paste0("Round ", iter, " of merging. N=", length(levels(object@ident))," clusters are present."))
    if (do.tsne) tsplot(object)
    iter = iter+1
    if (verbose){ cat("=======================================================\n"); cat("\n")}
    # Find closest leaves in the dendrogram, and a set of clusters to test for DE
    clusters_to_test = ClosestClustersInDendro(tree.use, num.nodes.test = num.nodes.test);
   
    num_merge = 0
    # Merging (for each row in clusters_to_test)
    # use find.markers to find DE genes (based on choices provided)
    # Decide to merge or not
    # use force.merge to merge clusters
    # if > 0 pairs of clusters have been merged then the indicator variable l remains -1, else change it to +1
    # Deleting a cluster?
    for (idx in 1:dim(clusters_to_test)[1]){
      cluster_1 = clusters_to_test[idx,1] 
      cluster_2 = clusters_to_test[idx,2] 
      if (verbose) print(paste0("Testing clusters ", cluster_1, " vs ", cluster_2))
      DE_gene_12 = find.markers(object, 
                                ident.1=cluster_1, 
                                ident.2=cluster_2, genes.use = NULL, 
                                thresh.use=thresh.use, 
                                test.use=test.use,
                                pval.thresh = 1e-2,
                                min.perc=min.perc,
                                only.pos = FALSE, 
                                append.all.pct = FALSE,
                                min.genes = 2*min.de.genes,
                                use.regress.data = regress.use.for.de)

      if (length(DE_gene_12$pval_adj) >= 2*min.de.genes){
        if (verbose) print(paste0("Found ", dim(subset(DE_gene_12, diff>0))[1], 
                                  " genes enriched in cluster ",
                                  cluster_1, " and ",
                                  dim(subset(DE_gene_12, diff<0))[1],
                                  " genes enriched in cluster ",
                                  cluster_2))
        
        if (!is.null(genes.exclude)){
          print("Removing censored genes")
          DE_gene_12 = DE_gene_12[!(rownames(DE_gene_12) %in% genes.exclude), ]
        }
        
        if (remove.mito_ribo){
          # Remove mitochondrial genes and ribosomal genes
          DE_gene_12 = DE_gene_12[!grepl("mt-",rownames(DE_gene_12)),]
          DE_gene_12 = DE_gene_12[!grepl("MT-",rownames(DE_gene_12)),]
          DE_gene_12 = DE_gene_12[!grepl("^Rps",rownames(DE_gene_12)),]
          DE_gene_12 = DE_gene_12[!grepl("^Rpl",rownames(DE_gene_12)),]
          DE_gene_12 = DE_gene_12[!grepl("^RPS",rownames(DE_gene_12)),]
          DE_gene_12 = DE_gene_12[!grepl("^RPL",rownames(DE_gene_12)),]
        
          if (verbose) print(paste0("Removing mitochondrial and ribosomal genes"))
          if (verbose) print(paste0("Remaining ", dim(subset(DE_gene_12, diff>0))[1], 
                                    " genes enriched in cluster ",
                                    cluster_1, " and ",
                                    dim(subset(DE_gene_12, diff<0))[1],
                                    " genes enriched in cluster ",
                                    cluster_2))
        }
        
        if (length(DE_gene_12$pval) >= 2*min.de.genes){
          if (verbose) print(paste0("Ensuring that enrichments satisfy min_perc_ratio ", min_perc_ratio))
          gene_Inc1 = DE_gene_12[DE_gene_12$diff>0 & DE_gene_12[,3]/DE_gene_12[,4] >= min_perc_ratio,]
          gene_Dec1 = DE_gene_12[DE_gene_12$diff<0 & DE_gene_12[,4]/DE_gene_12[,3] >= min_perc_ratio,]
          
          if (verbose) print(paste0("Ensuring that enrichments satisfy min_nTrans_ratio ", min_nTrans_ratio))
          gene_Inc2 = DE_gene_12[DE_gene_12$diff>0 & DE_gene_12[,7]/DE_gene_12[,8] >= min_nTrans_ratio,]
          gene_Dec2 = DE_gene_12[DE_gene_12$diff<0 & DE_gene_12[,8]/DE_gene_12[,7] >= min_nTrans_ratio,]
          
          gene_Inc = DE_gene_12[intersect(rownames(gene_Inc1), rownames(gene_Inc2)),]
          gene_Dec = DE_gene_12[intersect(rownames(gene_Dec1), rownames(gene_Dec2)),]
          if (verbose) print(paste0("Remaining ", nrow(gene_Inc), 
                                    " genes enriched in cluster ",
                                    cluster_1, " and ",
                                    nrow(gene_Dec),
                                    " genes enriched in cluster ",
                                    cluster_2))
          # Check how many ON-OFF genes are there
          nDE_on = nrow(subset(DE_gene_12, pct_clust1 > 0.5 & pct_clust2 < 0.1))
          nDE_off = nrow(subset(DE_gene_12, pct_clust2 > 0.5 & pct_clust1 < 0.1))
          
          if (nDE_on < min.on.off | nDE_off < min.on.off){ 
            if (((dim(gene_Inc)[1] + dim(gene_Dec)[1]) < 2*min.de.genes)) {
              num_merge = num_merge + 1
              object = force.merge(object, merge.clust = c(cluster_1, cluster_2))
              print(cat(paste0("Cluster ", cluster_1, " and cluster ", cluster_2, " are merged!\n")), file=stderr())
            } 
          }
          
        } else {
          num_merge = num_merge + 1
          object = force.merge(object, merge.clust = c(cluster_1, cluster_2))
          print(cat(paste0("Cluster ", cluster_1, " and cluster ", cluster_2, " are merged!\n")), file=stderr())
        }
        
      } else {
        if (verbose) print(paste0("Did not find sufficiently DE genes between Cluster ", cluster_1, " and cluster ", cluster_2))
        num_merge = num_merge + 1
        object = force.merge(object, merge.clust = c(cluster_1, cluster_2))
        print(cat(paste0("Cluster ", cluster_1, " and cluster ", cluster_2, " are merged!\n")), file=stderr())
      }
      
    }
    
    # If merging has happened
    # Renumber clusters 
    # Spit out a diagnostic statement
    # build dendrogram
    # repeat until convergence
    if (num_merge>0){
      if(verbose) print("Before re-build cluster tree")
      if(verbose) print(table(object@ident))
      object@ident=mapvalues(object@ident, from = levels(object@ident), to = c(1:length(levels(object@ident))))
      object@data.info[, merge_id] = object@ident[object@cell.names]
      if (length(levels(object@ident)) > 1){
        print(length(var.genes))
        object = buildClusterTree(object, genes.use=var.genes, stat.fxn = stat.fxn, linkage.method = linkage.method, dist.fun = dist.fun, regress.use = regress.use.for.dendro, do.scale=do.scale.dendro) 
      } else {
        l=1
      }
      tree.use = object@cluster.tree[[1]]
      if(verbose) print("after re-build cluster tree")
      if(verbose) print(table(object@ident))
    } else{
      l=1
    }
    
  }
  
  print(paste0("Recording the merged identities as ", merge_id))
  new_ident = object@ident
  object@data.info[names(new_ident), merge_id] = new_ident
  
  
  # Compare the current clusters with original and summarize merged clusters
  if (verbose) CompareClusterings(orig_ident, new_ident)
  
  return(object)
}

########################## UTILITIES FOR DE #####################

#' Returns an Nx2 matrix of clusters cluster pairs in a dendrogram
ClosestClustersInDendro = function(tree.use, num.nodes.test = num.nodes.test){
  
  internal_nodes = unique(tree.use$edge[,1])
  edge_sums = c()
  for (node in internal_nodes){
    # Calculate edge sums
    # First find the internal descendents of this node
    node_descendants = getDescendants(tree.use, node)
    node_descendants = node_descendants[node_descendants > tree.use$Nnode + 1]
    #Find the edge ids that have the node or its internal descendants
    edge_ids = which(tree.use$edge[,1] %in% c(node, node_descendants))
    # sum the edge lengths
    edge_sums = c(edge_sums, sum(tree.use$edge.length[edge_ids]))
  }
  
  # Find 3 nodes with the closest descendants and get the corresponding cluster pars
  min_node = internal_nodes[order(edge_sums)[1:3]]
  clusters_to_test = matrix(0,nrow=0, ncol=2)
  for (node in min_node){
    clusters = getDescendants(tree.use, node)
    if (length(clusters) == 2){
      clusters_to_test = rbind(clusters_to_test, clusters)
    } else {
      next
    }
  }
  return(clusters_to_test)
  
}

#' Differential expression using Wilcoxon Rank Sum Test
#'
#' Identifies differentially expressed genes between two groups of cells using
#' a Wilcoxon Rank Sum test
#' @param object scR object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param min.cells Minimum number of cells expressing the gene in at least one
#' of the two groups
#' @param genes.use Genes to use for test
#' @param print.bar Print a progress bar
#' @param ... Extra parameters passed to wilcox.test
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @importFrom pbapply pbsapply
#' @importFrom stats wilcox.test
#'
#' @export
#'
#' @examples
#' pbmc_small
#' WilcoxDETest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, ident = 1),
#'             cells.2 = WhichCells(object = pbmc_small, ident = 2))
#'
WilcoxDETest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  verbose = TRUE,
  ...
) {
  data.test <- object@data
  genes.use <- set.ifnull(genes.use, rownames(data.test))
  
  # check that the gene and cells made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(data.test)]
  cells.1 <- cells.1[cells.1 %in% colnames(data.test)]
  cells.2 <- cells.2[cells.2 %in% colnames(data.test)]
  
  coldata <- object@data.info[c(cells.1, cells.2), ]
  coldata[cells.1, "group"] <- "Group1"
  coldata[cells.2, "group"] <- "Group2"
  coldata$group <- factor(coldata$group)
  coldata$wellKey <- rownames(coldata)
  
  countdata.test <- data.test[genes.use, rownames(coldata)]
  if (nrow(countdata.test) == 0) return(data.frame())
  funapply <- if (verbose) {pbsapply} else {sapply}
  pval <- funapply(
    X = 1:nrow(countdata.test),
    FUN = function(x) {
      return(wilcox.test(countdata.test[x, ] ~ coldata$group, ...)$p.value)
    }
  )
  genes.return <- rownames(countdata.test)
  to.return <- data.frame(pval, row.names = genes.return)
  return(to.return)
}

#' ROC-based marker discovery
#'
#' Identifies 'markers' of gene expression using ROC analysis. For each gene,
#' evaluates (using AUC) a classifier built on that gene alone, to classify
#' between two groups of cells.
#'
#' An AUC value of 1 means that expression values for this gene alone can
#' perfectly classify the two groupings (i.e. Each of the cells in cells.1
#' exhibit a higher level than each of the cells in cells.2). An AUC value of 0
#' also means there is perfect classification, but in the other direction. A
#' value of 0.5 implies that the gene has no predictive power to classify the
#' two groups.
#'
#' @inheritParams FindMarkers
#' @inheritParams DiffExpTest
#' @param object Seurat object
#'
#' @return Returns a 'predictive power' (abs(AUC-0.5)) ranked matrix of
#' putative differentially expressed genes.
#'
#' @export
#'
#' @examples
#' pbmc_small
#' MarkerTest(pbmc_small, cells.1 = which.cells(object = pbmc_small,  1),
#'             cells.2 = which.cells(object = pbmc_small, 2))
#'
ROCTest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  verbose = TRUE
) {
  
  data.test <- object@data
  genes.use <- set.ifnull(genes.use, rownames(data.test))
  
  # check that the gene and cells made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(data.test)]
  cells.1 <- cells.1[cells.1 %in% colnames(data.test)]
  cells.2 <- cells.2[cells.2 %in% colnames(data.test)]
  
  to.return <- AUCMarkerTest(
    data1 = data.test[, cells.1],
    data2 = data.test[, cells.2],
    genes.use = genes.use,
    verbose = verbose
  )
  to.return$power <- abs(x = to.return$myAUC - 0.5) * 2
  #print(head(to.return))
  return(to.return)
}

# internal function to calculate AUC values
AUCMarkerTest <- function(data1, data2, genes.use, verbose = TRUE) {
  myAUC <- unlist(x = lapply(
    X = genes.use,
    FUN = function(x) {
      return(DifferentialAUC(
        x = as.numeric(x = data1[x, ]),
        y = as.numeric(x = data2[x, ])
      ))
    }
  ))
  
  myAUC[is.na(x = myAUC)] <- 0
  if (verbose) {
    funapply <- pblapply
  } else {
    funapply <- lapply
  }
  diff <- unlist(x = funapply(
    X = genes.use,
    FUN = function(x) {
      return(
        expMean(
          x = as.numeric(x = data1[x, ])
        ) - expMean(
          x = as.numeric(x = data2[x, ])
        )
      )
    }
  ))
  toRet <- data.frame(cbind(myAUC, diff), row.names = genes.use)
  toRet <- toRet[rev(x = order(toRet$myAUC)), ]
  return(toRet)
}

# internal function to calculate AUC values
#' @importFrom ROCR prediction performance
DifferentialAUC <- function(x, y) {
  prediction.use <- prediction(
    predictions = c(x, y),
    labels = c(rep(x = 1, length(x = x)), rep(x = 0, length(x = y))),
    label.ordering = 0:1
  )
  perf.use <- performance(prediction.obj = prediction.use, measure = "auc")
  auc.use <- round(x = perf.use@y.values[[1]], digits = 3)
  return(auc.use)
}


#' Likelihood ratio test for zero-inflated data
#'
#' Identifies differentially expressed genes between two groups of cells using
#' the LRT model proposed in McDavid et al, Bioinformatics, 2013
#'
#' @inheritParams FindMarkers
#' @param object scR object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @export
#' @examples
#' pbmc_small
#' DiffExpTest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, ident = 1),
#'             cells.2 = WhichCells(object = pbmc_small, ident = 2))
#'
DiffExpTest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  verbose = TRUE,
  xmin = 0
) {
  data.test <- object@data
  # check that the gene and cells made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(data.test)]
  cells.1 <- cells.1[cells.1 %in% colnames(data.test)]
  cells.2 <- cells.2[cells.2 %in% colnames(data.test)]
  data.test <- data.test[genes.use, c(cells.1, cells.2)]
  if (verbose) {
    funapply <- pblapply
  } else {
    funapply <- lapply
  }
  
  
  #internal function
  #
  DifferentialLRT <- function(x, y, xmin = 0) {
    lrtX <- bimodLikData(x = x, xmin=xmin)
    lrtY <- bimodLikData(x = y, xmin = xmin)
    lrtZ <- bimodLikData(x = c(x, y), xmin = xmin)
    lrt_diff <- 2 * (lrtX + lrtY - lrtZ)
    return(pchisq(q = lrt_diff, df = 3, lower.tail = F))
  }
  
  #internal function to run mcdavid et al. DE test
  ##
  bimodLikData <- function(x, xmin = 0) {
    x1 <- x[x <= xmin]
    x2 <- x[x > xmin]
    xal <- min(
      length(x2) / length(x),
      min = 1e-5,
      max = (1 - 1e-5)
    )
    likA <- length(x1) * log(1 - xal)
    if (length(x2) < 2) {
      mysd <- 1
    } else {
      mysd <- sd(x2)
    }
    likB <- length(x2) *
      log(xal) +
      sum(dnorm(x2, mean = mean(x = x2), sd = mysd, log = TRUE))
    return(likA + likB)
  }
  
  
  pval <- unlist(
    x = funapply(
      X = genes.use,
      FUN = function(x) {
        return(
          DifferentialLRT(
            x = as.numeric(data.test[x, cells.1]),
            y = as.numeric(data.test[x, cells.2]),
            xmin = xmin
          )
        )
      }
    )
  )
  to.return <- data.frame(pval, row.names = genes.use)
  return(to.return)
}

#' Differential expression using MAST
#'
#' Identifies differentially expressed genes between two groups of cells using
#' a hurdle model tailored to scRNA-seq data. Utilizes the MAST package to run
#' the DE testing.
#'
#' @references Andrew McDavid, Greg Finak and Masanao Yajima (2017). MAST: Model-based
#' Analysis of Single Cell Transcriptomics. R package version 1.2.1.
#' https://github.com/RGLab/MAST/
#' @param object scR object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param min.cells Minimum number of cells expressing the gene in at least one
#' of the two groups
#' @param genes.use Genes to use for test
#' @param latent.vars Confounding variables to adjust for in DE test. Default is
#' "nTranscripts", which adjusts for cellular depth (i.e. cellular detection rate). For
#' non-UMI based data, set to nGene instead.
#' @param \dots Additional parameters to zero-inflated regression (zlm) function
#' in MAST
#' @details
#' To use this method, please install MAST, using instructions at https://github.com/RGLab/MAST/
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @importFrom stats relevel
#' @importFrom utils installed.packages
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   pbmc_small
#'   MASTDETest(pbmc_small, cells.1 = which.cells(object = pbmc_small, ident = 1),
#'               cells.2 = which.cells(object = pbmc_small, ident = 2))
#' }
#'
MASTDETest <- function(
  object,
  cells.1,
  cells.2,
  min.cells = 10,
  genes.use = NULL,
  latent.vars = NULL,
  ...
) {
  # Check for MAST
  if (!'MAST' %in% rownames(x = installed.packages())) {
    stop("Please install MAST - learn more at https://github.com/RGLab/MAST")
  }
  data.test <- object@data
  # check that the gene and cells made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(data.test)]
  cells.1 <- cells.1[cells.1 %in% colnames(data.test)]
  cells.2 <- cells.2[cells.2 %in% colnames(data.test)]
  data.test <- data.test[genes.use, c(cells.1, cells.2)]
  print(1)
  my.latent <- fetch.data(
    object = object,
    vars.all = latent.vars,
    cells.use = c(cells.1, cells.2),
    use.count  = TRUE
  )
  print(2)
  if (length(latent.vars) > 0) {
    my.latent <- scale(my.latent)
  }
  coldata <- data.frame(my.latent, row.names = c(cells.1, cells.2))
  coldata[cells.1, "group"] <- "Group1"
  coldata[cells.2, "group"] <- "Group2"
  coldata$group <- factor(coldata$group)
  coldata$wellKey <- rownames(coldata)
  latent.vars <- c("condition", latent.vars)
  countdata.test <- data.test[genes.use, rownames(coldata)]
  fdat <- data.frame(rownames(countdata.test))
  colnames(fdat)[1] <- "primerid"
  rownames(fdat) <- fdat[, 1]
  sca <- MAST::FromMatrix(
    exprsArray = as.matrix(countdata.test),
    cData = coldata,
    fData = fdat
  )
  cond <- factor(SummarizedExperiment::colData(sca)$group)
  cond <- relevel(cond, ref = "Group1")
  SummarizedExperiment::colData(sca)$condition <- cond
  fmla <- as.formula(
    object = paste0(" ~ ", paste(latent.vars, collapse = "+"))
  )
  zlmCond <- MAST::zlm(formula = fmla, sca = sca, ...)
  summaryCond <- summary(object = zlmCond, doLRT = 'conditionGroup2')
  summaryDt <- summaryCond$datatable

  # fcHurdle <- merge(
  #   summaryDt[contrast=='conditionGroup2' & component=='H', .(primerid, `Pr(>Chisq)`)], #hurdle P values
  #   summaryDt[contrast=='conditionGroup2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid'
  # ) #logFC coefficients
  # fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  pval <- subset(summaryDt, component=="H")[,4]$`Pr(>Chisq)`
  #p_val <- summaryDt[summaryDt[, "component"] == "H", 4]
  genes.return <- as.character(subset(summaryDt, component=="H")[,1]$primerid)
  #genes.return <- summaryDt[summaryDt[, "component"] == "H", 1]
  # p_val <- subset(summaryDt, component == "H")[, 4]
  # genes.return <- subset(summaryDt, component == "H")[, 1]
  to.return <- data.frame(pval, row.names = genes.return)
  return(to.return)
}


#' Differential expression using edgeR QLF Test
#'
#' @param object scR object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param genes.use Genes to use for test.
#' @param \dots Additional parameters to edgeR
#' @details
#' To use this method, please install edgeR, using instructions at https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @importFrom stats relevel
#' @importFrom utils installed.packages
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   pbmc_small
#'   edgeRQLFTest(pbmc_small, cells.1 = which.cells(object = pbmc_small, ident = 1),
#'               cells.2 = which.cells(object = pbmc_small, ident = 2))
#' }
#'
edgeRQLFTest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL
) {
  # Check for MAST
  if (!'edgeR' %in% rownames(x = installed.packages())) {
    stop("Please install edgeR - learn more at http://bioconductor.org/packages/release/bioc/html/edgeR.html")
  }
  data.test <- object@data
  # check that the gene and cells made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(data.test)]
  cells.1 <- cells.1[cells.1 %in% colnames(data.test)]
  cells.2 <- cells.2[cells.2 %in% colnames(data.test)]
  data.test <- data.test[genes.use, c(cells.1, cells.2)]
  data.test <- exp(data.test) - 1
  
  coldata = c(rep("Group1", length(cells.1)), rep("Group2", length(cells.2)))
  coldata <- data.frame(coldata, row.names = c(cells.1, cells.2))
  colnames(coldata) = "group"
  coldata$group <- factor(coldata$group)
  y <- edgeR::DGEList(counts = data.test, group=coldata$group)
  y <- calcNormFactors(y)
  design <- model.matrix(~coldata$group)
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  pval = qlf$table$PValue
  to.return <- data.frame(pval, row.names = rownames(qlf$table))
  return(to.return)
}


#' Differential expression using tweeDEseq Test
#'
#' @param object scR object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @param genes.use Genes to use for test.
#' @param \dots Additional parameters to edgeR
#' @details
#' To use this method, please install tweeDEseq, using instructions at https://bioconductor.org/packages/release/bioc/vignettes/tweeDEseq/inst/doc/tweeDEseq.pdf
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @importFrom stats relevel
#' @importFrom utils installed.packages
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   pbmc_small
#'   tweeDEseqTest(pbmc_small, cells.1 = which.cells(object = pbmc_small, ident = 1),
#'               cells.2 = which.cells(object = pbmc_small, ident = 2))
#' }
#'
tweeDEseqTest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL
) {
  # Check for tweeDEseq
  if (!'tweeDEseq' %in% rownames(x = installed.packages())) {
    stop("Please install tweeDEseq - learn more at https://bioconductor.org/packages/release/bioc/html/tweeDEseq.html")
  }
  data.test <- object@count.data
  # check that the gene and cells made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(data.test)]
  cells.1 <- cells.1[cells.1 %in% colnames(data.test)]
  cells.2 <- cells.2[cells.2 %in% colnames(data.test)]
  data.test <- as.matrix(data.test[genes.use, c(cells.1, cells.2)])
  
  # Define groups
  coldata = c(rep("Group1", length(cells.1)), rep("Group2", length(cells.2)))
  coldata <- data.frame(coldata, row.names = c(cells.1, cells.2))
  colnames(coldata) = "group"
  coldata$group <- factor(coldata$group)
  
  normCounts = tweeDEseq::normalizeCounts(data.test, group=coldata$group)
  
  # Fitting distribution (Diagnosis)
  #xf = normCounts["Tmeff2",cells.1]
  #compareCountDist(xf)
  
  resPT <- tweeDE(normCounts, group = coldata$group)
  pval <- resPT$pval
  
  to.return <- data.frame(pval, row.names = rownames(resPT))
  return(to.return)
}

#' Differential expression testing using Tobit models
#'
#' Identifies differentially expressed genes between two groups of cells using
#' Tobit models, as proposed in Trapnell et al., Nature Biotechnology, 2014
#'
#' @inheritParams FindMarkers
#' @inheritParams DiffExpTest
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @export
#'
#'@examples
#' pbmc_small
#' \dontrun{
#' TobitTest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, ident = 1),
#'             cells.2 = WhichCells(object = pbmc_small, ident = 2))
#' }
#'
TobitTest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  verbose = TRUE
) {
  
  
  data.test <- object@data
  # check that the gene and cells made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(data.test)]
  cells.1 <- cells.1[cells.1 %in% colnames(data.test)]
  cells.2 <- cells.2[cells.2 %in% colnames(data.test)]
  data.test <- data.test[genes.use, c(cells.1, cells.2)]
  
  to.return <- TobitDiffExpTest(
    data1 = data.test[, cells.1],
    data2 = data.test[, cells.2],
    genes.use = genes.use,
    verbose = verbose
  )
  return(to.return)
}


#internal function to run Tobit DE test
TobitDiffExpTest <- function(data1, data2, genes.use, verbose) {
  
  if (verbose) {
    funapply <- pblapply
  } else {
    funapply <- lapply
  }
  
  
  pval <- unlist(funapply(
    X = genes.use,
    FUN = function(x) {
      return(DifferentialTobit(
        x1 = as.numeric(x = data1[x, ]),
        x2 = as.numeric(x = data2[x, ])
      ))}
  ))
  
  pval[is.na(pval)] <- 1
 
  to.return <- data.frame(pval, row.names = genes.use)
  return(to.return)
}

#internal function to run Tobit DE test
#
#' @importFrom stats pchisq logLik
#
DifferentialTobit <- function(x1, x2, lower = 1, upper = Inf) {
  df <- data.frame(
    c(x1, x2),
    c(rep(x = 0, length(x = x1)), rep(x = 1, length(x = x2)))
  )
  colnames(x = df) <- c("Expression", "Stat")
  #model.v1=vgam(Expression~1,family = tobit(Lower = lower,Upper = upper),data = my.df)
  model.v1 <- TobitFitter(
    x = df,
    modelFormulaStr = "Expression~1",
    lower = lower,
    upper = upper
  )
  #model.v2=vgam(Expression~Stat+1,family = tobit(Lower = lower,Upper = upper),data = my.df)
  model.v2 <- TobitFitter(
    x = df,
    modelFormulaStr = "Expression~Stat+1",
    lower = lower,
    upper = upper
  )
  # if (is.null(x = model.v1) == FALSE && is.null(x = model.v2) == FALSE) {
  if (! is.null(x = model.v1) && ! is.null(x = model.v2)) {
    p <- pchisq(
      q = 2 * (logLik(object = model.v2) - logLik(object = model.v1)),
      df = 1,
      lower.tail = FALSE
    )
  } else {
    p <- 1
  }
  return(p)
}

#internal function to run Tobit DE test
#credit to Cole Trapnell for this
#
#' @importFrom VGAM vgam tobit
#' @importFrom stats as.formula
#
TobitFitter <- function(x, modelFormulaStr, lower = 1, upper = Inf){
  tryCatch(
    expr = return(suppressWarnings(expr = vgam(
      formula = as.formula(object = modelFormulaStr),
      family = tobit(Lower = lower, Upper = upper),
      data = x
    ))),
    #warning = function(w) { FM_fit },
    error = function(e) { NULL }
  )
}

#' Binomial Test (suitable for detecting lowly expressed genes in UMI based data)
#'
#'
#' @inheritParams find.markers
#' @param object scR object
#' @param cells.1 Group 1 cells
#' @param cells.2 Group 2 cells
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @export
#' @examples
#' pbmc_small
#' DiffExpTest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, ident = 1),
#'             cells.2 = WhichCells(object = pbmc_small, ident = 2))
#'


BinomTest = function(
          object, 
          cells.1,
          cells.2, 
          genes.use = NULL, 
          TPM.mat, 
          xmin=0,
          ...) {
            
            
            #Test for enrichments in cluster #1
            m = apply(TPM.mat[, cells.2], 1, function(x) sum(x>xmin)) #Number of cells expressing marker in cluster #2
            m1 = m; m1[m==xmin]=1; # Regularization. Add a pseudo count of 1 to unexpressed markers to avoid false positives
            n = apply(TPM.mat[, cells.1], 1, function(x) sum(x>xmin)) #Number of cells expressing marker in cluster #1
            
            #Find the probability of finding n or more cells +ve for marker in cluster 1 given the fraction in cluster 2
            pv1 = pbinom(n, length(cells.1), m1/length(cells.2), lower.tail = FALSE) + dbinom(n, length(cells.1), m1/length(cells.2))
            
            log_fold_express = log(n*length(cells.2)/(m*length(cells.1))) #log proportion of expressing cells
            d1 <- data.frame(logFoldPerc=log_fold_express,pval=pv1)
            d1 <- subset(d1, logFoldPerc >= 0)
            d1 <- d1[order(d1$pval,decreasing=FALSE),]
            
            #Enrichments in cells.2
            n1 = n; n1[n==xmin]=1; # Regularization.  Add a pseudo count of 1 to unexpressed markers to avoid false positives
            #Find the probability of finding m or more cells +ve for marker in cluster 2 given the fraction in cluster 1
            pv2 = pbinom(m, length(cells.2), n1/length(cells.1), lower.tail=FALSE) + dbinom(m, length(cells.2), n1/length(cells.1))
            d2 <- data.frame(logFoldPerc=log_fold_express,pval=pv2)
            d2 <- subset(d2, logFoldPerc < 0)
            d2 <- d2[order(d2$pval,decreasing=FALSE),]
            
            d = rbind(d1, d2);
            d = d[order(d$pval, decreasing=FALSE),]
            to.return = data.frame(pval = d$pval, row.names = rownames(d))
            return(d)
          
}


#' Differential expression testing using Student's t-test
#'
#' Identify differentially expressed genes between two groups of cells using
#' the Student's t-test
#'
#' @inheritParams FindMarkers
#' @inheritParams DiffExpTest
#'
#' @return Returns a p-value ranked matrix of putative differentially expressed
#' genes.
#'
#' @importFrom stats t.test
#' @importFrom pbapply pblapply
#'
#' @export
#'
#' @examples
#' pbmc_small
#' DiffTTest(pbmc_small, cells.1 = WhichCells(object = pbmc_small, ident = 1),
#'             cells.2 = WhichCells(object = pbmc_small, ident = 2))
DiffTTest <- function(
  object,
  cells.1,
  cells.2,
  genes.use = NULL,
  verbose = TRUE
) {
  
  data.test <- object@data
  # check that the gene and cells made it through the any filtering that was done
  genes.use <- genes.use[genes.use %in% rownames(data.test)]
  cells.1 <- cells.1[cells.1 %in% colnames(data.test)]
  cells.2 <- cells.2[cells.2 %in% colnames(data.test)]
  data.test <- data.test[genes.use, c(cells.1, cells.2)]
  
  if (verbose) {
    funapply <- pblapply
  } else {
    funapply <- lapply
  }
  
  pval <- unlist(
    x = funapply(
      X = genes.use,
      FUN = function(x) {
        t.test(x = data.test[x, cells.1], y = data.test[x, cells.2])$p.value
      }
    )
  )
  to.return <- data.frame(pval,row.names = genes.use)
  return(to.return)
}


# Depth first traversal path of a given tree
#
# @param tree              Tree object (from ape package)
# @param node              Internal node in the tree
# @param path              Path through the tree (for recursion)
# @param include.children  Include children in the output path
# @param only.children     Only include children in the output path
# @return                  Returns a vector representing the depth first
#                          traversal path
#
DFT <- function(
  tree,
  node,
  path = NULL,
  include.children = FALSE,
  only.children = FALSE
) {
  if (only.children) {
    include.children = TRUE
  }
  children <- which(x = tree$edge[, 1] == node)
  child1 <- tree$edge[children[1], 2]
  child2 <- tree$edge[children[2], 2]
  if (child1 %in% tree$edge[, 1]) {
    if(! only.children){
      path <- c(path, child1)
    }
    path <- DFT(
      tree = tree,
      node = child1,
      path = path,
      include.children = include.children,
      only.children = only.children
    )
  } else {
    if (include.children) {
      path <-c(path, child1)
    }
  }
  if (child2 %in% tree$edge[, 1]) {
    if (! only.children) {
      path <- c(path, child2)
    }
    path <- DFT(
      tree = tree,
      node = child2,
      path = path,
      include.children = include.children,
      only.children = only.children
    )
  } else {
    if (include.children) {
      path <- c(path, child2)
    }
  }
  return(path)
}



#################################### OLD STUFF - need to revamp ########################################################

# t-test
setGeneric("node.dot.plot", function(object, node=NULL, test.use="t", thresh.use=log(2), min.pos.perc=0.4, num.markers=5, force.genes=NULL,max.val.exp = 5, max.val.perc = 1, norm.exp=NULL) standardGeneric("node.dot.plot"))
setMethod("node.dot.plot", "scR",
          function(object, node=NULL, test.use="t",thresh.use=log(2), min.pos.perc=0.4, num.markers=5, force.genes=NULL,max.val.exp = 5, max.val.perc = 1,norm.exp=NULL) {
            
            if (is.null(force.genes)){
              markers = find.markers.node(object,node,test.use=test.use,thresh.use=thresh.use, min.pos.perc=min.pos.perc)
              ident = c(object@cluster.tree[[1]]$tip.label[getLeftDecendants(object@cluster.tree[[1]],node)],rgcs_only@cluster.tree[[1]]$tip.label[getRightDecendants(object@cluster.tree[[1]],node)])
              markers.1 = subset(markers,diff>0)
              if (nrow(markers.1) > 0){
                markers.1 = markers.1[order(-markers.1$diff),]
                genes = c(rownames(markers.1)[1:min(nrow(markers.1),num.markers)])
              }
              
              markers.2 = subset(markers,diff<0)
              if (nrow(markers.2) > 0){
                markers.2 = markers.2[order(-markers.2$diff),]
                genes = c(genes,c(rownames(markers.2)[1:min(nrow(markers.2),num.markers)]))
              }
            } else{
              genes = force.genes
            }
            
            
            Perc.pos.by.ident(object,genes, ident.use = ident, max.val.exp = max.val.exp, max.val.perc = max.val.perc, norm.exp = norm.exp)
            return(1)
          } 
)

# t-test
setGeneric("marker.code", function(object, clust=NULL,marker.lists=NULL, genes.use=NULL, length.code=2, thresh.use=0.2, min.pos.perc=0.3, tests.use=c("bimod","binom","t"), TPM.mat=NULL, Count.mat=NULL, max.markers=NULL) standardGeneric("marker.code"))
setMethod("marker.code", "scR",
          function(object, clust=NULL,marker.lists=NULL, genes.use=NULL, length.code=2,thresh.use=0.2, min.pos.perc=0.3,tests.use=c("bimod","binom","t"), TPM.mat=NULL, Count.mat=NULL, max.markers=NULL) {
            
            genes = c()
          if (is.null(genes.use)){
            if (is.null(marker.lists)){ 
              markers = list()
              l=1
              for (test in tests.use){
                markers[[l]] = find.markers(object, cluster, thresh.use = thresh.use, test.use=test, min.pos.perc = min.pos.perc, max.neg.perc = 1, only.pos=TRUE, TPM.mat=TPM.mat, Count.mat=Count.mat)
                markers[[l]] = subset(markers[[l]],diff>0)
                genes = union(genes, rownames(markers[[l]]))
              }
           } else {
             markers_big = marker.lists[[1]]
             for (l in 2:length(marker.lists)){
               if (l >= length(marker.lists)) break
               markers_big = rbind(markers_big, marker.lists[[l]])
             }
             markers_big = subset(markers_big, cluster==clust & diff > 0)
             markers_big = markers_big[order(-markers_big$diff, -markers_big$nTrans_cluster),]
             genes = unique(markers_big$gene)
           }
            
            if (is.null(max.markers)){
               if (length(genes) > 100){
                 print(paste0("Considering top ", 100, " markers. Specifiy max.markers to increase..."))
                 genes=genes[1:100]
               }
            } else {
              genes = genes[1:min(length(genes),max.markers)]
            }
            
          } else {
            genes = genes.use
          }
            df = PercExpDf(object,genes, perc.floor = 0,exp.floor = 0)
            # filter genes
            col.select = "clust"
            genes.filt = c()
            for (g in genes){
              if (df[df$clust == clust,paste0("Perc_",g)] >= min.pos.perc){
                col.select = c(col.select,paste0(c("Exp_","Perc_"), g))
                genes.filt = c(genes.filt, g)
              }
            }
            genes = genes.filt
            df = df[,col.select]
            
            df_norm = df
            for (c in 2:ncol(df)){
              df_norm[,c] = df[,c] / max(df[,c])
            }
            
            marker.comb = combn(genes,length.code)
            is.ok.perc = c()
            
            for (comb in 1:ncol(marker.comb)){
              col.select = c("clust")
              for (j in 1:length.code){
                col.select = c(col.select, grep(marker.comb[j,comb], colnames(df_norm), value=TRUE))
              }
              df1 = df_norm[,col.select]
              df2 = dplyr::select(df1, contains("clust"),contains("Perc"),contains("Exp"))
              for (j in marker.comb[,comb]){
                if (median(df2[,paste0("Perc_",j)]) < 0.3){
                      df2 = df2[df2[,paste0("Perc_",j)] > 0.3, ]
                } else {
                      df2 = df2[df2[,paste0("Exp_",j)] > 0.3, ]
                }
              }
              if (nrow(df2) == 1){
                if (rownames(df2) == clust) is.ok.perc = c(is.ok.perc, 1)
              } else {
                is.ok.perc = c(is.ok.perc,0)
              }
            }
            
           marker.comb.success = marker.comb[,is.ok.perc==1]
           if (sum(is.ok.perc) == 0){
             print("No combinations found. Returning ...")
             return(NULL)
           }
           df.return = as.data.frame(matrix(nrow=0, ncol=length.code+6))
           
           
           for (c in 1:ncol(marker.comb.success)){
             markers.name = marker.comb.success[,c]
             perc.markers = df[as.character(clust), paste0("Perc_", markers.name)]
             exp.markers = df[as.character(clust), paste0("Exp_", markers.name)]
             #perc.ent = 0
             perc.vec = df[,paste0("Perc_", markers.name)]
             perc.vec = scale(perc.vec, center=FALSE, scale=colSums(perc.vec))
             perc.ent = -perc.vec * log(perc.vec)
             perc.ent[is.nan(perc.ent)] = 0
             perc.ent = max(colSums(perc.ent))
             #for (m in markers.name){
               
            #   perc.vec = perc.vec[perc.vec > 0] / sum(perc.vec)
            #   perc.ent = perc.ent - sum(perc.vec * log(perc.vec))
            # }
             
            # exp.ent = 0
             #for (m in markers.name){
            #   exp.vec = df[,paste0("Exp_", m)]
            #   exp.vec = exp.vec[exp.vec > 0] / sum(exp.vec)
            #   exp.ent = exp.ent - sum(exp.vec * log(exp.vec))
            # }
             exp.vec = df[,paste0("Exp_", markers.name)]
             exp.vec = scale(exp.vec, center=FALSE, scale=colSums(exp.vec))
             exp.ent = -exp.vec * log(exp.vec)
             exp.ent[is.nan(exp.ent)] = 0
             exp.ent = max(colSums(exp.ent))
             
             df.return = rbind(df.return, rep(1,ncol(df.return)))
             df.return[c,c(1:length.code)] = markers.name
             df.return[c,c((length.code+1):(2*length.code))] = perc.markers
             df.return[c,c((2*length.code+1):(3*length.code))] = exp.markers
             df.return[c,3*length.code+1] = perc.ent / (log(nrow(df)))
             df.return[c,3*length.code+2] = exp.ent / (log(nrow(df)))
             
           }
           colnames(df.return) = c(paste0("marker_", c(1:length.code)), paste0("PercClust_",c(1:length.code)),
                                   paste0("ExpClust_",c(1:length.code)), "Perc_entropy","Exp_entropy")
            df.return = plyr::arrange(df.return, Perc_entropy, Exp_entropy)
           
            return(df.return)
          } 
)

# t-test
setGeneric("single.marker", function(object, clust=NULL,marker.lists=NULL, thresh.use=0.2, min.pos.perc=0.3, tests.use=c("bimod","binom","t"), TPM.mat=NULL, Count.mat=NULL) standardGeneric("single.marker"))
setMethod("single.marker", "scR",
          function(object, clust=NULL,marker.lists=NULL, thresh.use=0.2, min.pos.perc=0.3,tests.use=c("bimod","binom","t"), TPM.mat=NULL, Count.mat=NULL) {
            
            genes = c()
            if (is.null(marker.lists)){ 
              markers = list()
              l=1
              for (test in tests.use){
                markers[[l]] = find.markers(object, cluster, thresh.use = thresh.use, test.use=test, min.pos.perc = min.pos.perc, max.neg.perc = 1, only.pos=TRUE, TPM.mat=TPM.mat, Count.mat=Count.mat)
                markers[[l]] = subset(markers[[l]],diff>0)
                genes = union(genes, rownames(markers[[l]]))
              }
            } else {
              markers_big = marker.lists[[1]]
              for (l in 2:length(marker.lists)){
                if (l >= length(marker.lists)) break
                markers_big = rbind(markers_big, marker.lists[[l]])
              }
              markers_big = subset(markers_big, cluster==clust & diff > 0)
              markers_big = markers_big[order(-markers_big$diff, -markers_big$nTrans_cluster),]
              genes = unique(markers_big$gene)
            }
            
            df = PercExpDf(object,genes, perc.floor = 0,exp.floor = 0)
            df_norm = df
            for (c in 2:ncol(df)){
              df_norm[,c] = df[,c] / max(df[,c])
            }
            
            
            df.return = as.data.frame(matrix(nrow=0, ncol=5))
            l=1
            for (g in genes){
                df.return = rbind(df.return, rep(1,ncol(df.return)))
                df.return[l,1] = g
                
                perc.vec = df[,paste0("Perc_", g)]; names(perc.vec) = df$clust
                if (names(perc.vec)[which.max(perc.vec)] == clust){
                  df.return[l,2] = perc.vec[which.max(perc.vec)] / sort(perc.vec,decreasing = TRUE)[2]
                } else {
                  df.return[l,2] = perc.vec[as.character(clust)] / perc.vec[which.max(perc.vec)]
                }
                
                exp.vec = df[,paste0("Exp_", g)]; names(exp.vec) = df$clust
                if (names(exp.vec)[which.max(exp.vec)] == clust){
                  df.return[l,3] = exp.vec[which.max(exp.vec)] / sort(exp.vec,decreasing = TRUE)[2]
                } else {
                  df.return[l,3] = exp.vec[as.character(clust)] / exp.vec[which.max(exp.vec)]
                }
                
        
                perc.ent = 0
                perc.vec = perc.vec[perc.vec > 0] / sum(perc.vec)
                perc.ent = perc.ent - sum(perc.vec * log(perc.vec))
                
                df.return[l,4] = perc.ent / log(nrow(df))
                
                exp.ent = 0
                exp.vec = exp.vec[exp.vec > 0] / sum(exp.vec)
                exp.ent = exp.ent - sum(exp.vec * log(exp.vec))
                df.return[l,5] = exp.ent / log(nrow(df))
                l=l+1
            }
            
            colnames(df.return) = c("marker", "Perc_diff","Exp_diff", "Perc_entropy","Exp_entropy")
            df.return = plyr::arrange(df.return, Perc_entropy, Exp_entropy)
            
            return(df.return)
          } 
)



# Neg Binom DE test

setGeneric("NegBinomDETest", function(object, cells.1,cells.2,genes.use=NULL,Count.mat=NULL, progress.bar=TRUE) standardGeneric("NegBinomDETest"))
setMethod("NegBinomDETest", "scR",
          function(object, cells.1,cells.2,genes.use=NULL, Count.mat=NULL, progress.bar=TRUE) {
            genes.use=set.ifnull(genes.use, rownames(object@data))
            if (is.null(Count.mat)){ 
              print("Need UMI data for NB test")
              return(NULL)
              }
            to.test.data=(Count.mat[genes.use,c(cells.1,cells.2)]); 
            to.test=data.frame(object@data.info[c(cells.1, cells.2),"nTranscripts"]); rownames(to.test) = c(cells.1,cells.2)
            to.test[cells.1,"group"]="A";
            to.test[cells.2,"group"]="B";
            to.test$group=factor(to.test$group)
            colnames(to.test) = c("nTranscripts","group")
            latent.vars=c("group","nTranscripts")
            iterate.fxn=lapply; if (progress.bar) iterate.fxn=pblapply
            pval=unlist(iterate.fxn(genes.use,function(x) {
              to.test[,"GENE"]=as.numeric(to.test.data[x,]);
              fmla=as.formula(paste("GENE ", " ~ ", paste(latent.vars,collapse="+"),sep=""));
              return(tryCatch(summary(glm.nb(fmla,data = to.test))$coef[2,4], error = function(e) NA))
            }))
          
            to.return=data.frame(pval,row.names = genes.use)
            remove.genes = rownames(subset(to.return, is.na(pval)))
            print(paste0("Neg binom regression did not converge for genes: ", paste(remove.genes, sep = ",")))
            to.return = subset(to.return, !is.na(pval))
            genes.use=rownames(to.return)
            to.return = data.frame(to.return[order(to.return$pval, decreasing=FALSE),], row.names=genes.use)
            colnames(to.return) = "pval"
            return(to.return)
            
          }
)

setGeneric("PoissonDETest", function(object, cells.1,cells.2,genes.use=NULL,progress.bar=TRUE) standardGeneric("PoissonDETest"))
#' @export
setMethod("PoissonDETest", "scR",
          function(object, cells.1,cells.2,genes.use=NULL,progress.bar=TRUE) {
            genes.use=set.ifnull(genes.use, rownames(object@data))
            if (is.null(Count.mat)){ 
              print("Need UMI data for NB test")
              return(NULL)
            }
            to.test.data=(Count.mat[genes.use,c(cells.1,cells.2)]); 
            to.test=data.frame(object@data.info[c(cells.1, cells.2),"nTranscripts"]); rownames(to.test) = c(cells.1,cells.2)
            to.test[cells.1,"group"]="A";
            to.test[cells.2,"group"]="B";
            to.test$group=factor(to.test$group)
            colnames(to.test) = c("nTranscripts","group")
            latent.vars=c("group","nTranscripts")
            
            iterate.fxn=lapply; if (progress.report) iterate.fxn=pblapply
            pval=unlist(iterate.fxn(genes.use,function(x) {
              to.test[,"GENE"]=as.numeric(to.test.data[x,])
              fmla=as.formula(paste("GENE ", " ~ ", paste(latent.vars,collapse="+"),sep=""));
              return(summary(glm(fmla,data = to.test,family = "poisson"))$coef[2,4])
            }))
            to.return=data.frame(pval,row.names = genes.use)
            to.return = to.return[order(to.return$pval, decreasing=FALSE),]
            return(to.return)
          }
)



# Computes the spatial autocorrelation according to the Moran I and Geary C statistic
# Used for ranking genes that show spatial structure

setGeneric("spatial.autocorr", function(object, clust=NULL,genes.use=NULL, num.nn=10, verbose=FALSE, min.perc=0.1, max.perc=0.8) standardGeneric("spatial.autocorr"))
setMethod("spatial.autocorr", "scR",
          function(object, clust=NULL, genes.use=NULL, num.nn=10, verbose=FALSE, min.perc=0.1, max.perc=0.8) {
            
            genes.use = set.ifnull(genes.use, object@var.genes)
            cells.use = which.cells(object, clust)
            cells.use = set.ifnull(cells.use, object@cell.names)
            coord.use = object@tsne.rot[cells.use,]
            
            kk = knn2nb(knearneigh(as.matrix(coord.use),num.nn))
            df.perc = PercExpDf(object,genes.use, perc.floor = 0,exp.floor = 0) %>% dplyr::select(contains("clust"),contains("Perc"))
            df.perc = df.perc[,c(1, 1+which(df.perc[as.character(clust), 2:ncol(df.perc)] > min.perc))]
            df.perc = df.perc[,c(1, 1+which(df.perc[as.character(clust), 2:ncol(df.perc)] < max.perc))]
            df.return = as.data.frame(matrix(nrow=0, ncol=4))
            l=1
            for (g in sapply(colnames(df.perc)[2:ncol(df.perc)], getStat2)){
              if (verbose) print(paste0("Testing gene: ", g))
              df.return = rbind(df.return, rep(1,ncol(df.return)))
              temp = moran.test(as.numeric(object@data[g,cells.use]), nb2listw(kk, style="C"))
              temp2 = geary.test(as.numeric(object@data[g,cells.use]), nb2listw(kk, style="C"))
              df.return[l,1] = g
              df.return[l,2] = df.perc[as.character(clust),paste0("Perc_",g)]
              df.return[l,3] = temp$estimate[1]
              df.return[l,4] = temp2$estimate[1]
              l=l+1
            }
            
            colnames(df.return) = c("marker","PercClust","MoranI","GearyC")
            df.return = plyr::arrange(df.return, GearyC, -MoranI)
            
            return(df.return)
          } 
)
