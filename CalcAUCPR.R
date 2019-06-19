#' For a list of markers (singlets, pairs, triplets ...), computes the AUC value under the precision recall curve
#' If genes.list is a character vector, then the function treats each element as a singlet
#' If a single gene combination is given, then the function also plots a PR curve
CalcAUCPR = function (object, ident.1, ident.2 = NULL, cells.1 = NULL, cells.2 = NULL, genes.list = NULL, 
          print.bar = TRUE, thresh.use=NULL, do.plot=FALSE, max.cells.per.ident=NULL) 
{
  library(pbapply)
  thresh.use = set.ifnull(thresh.use,0)
  if (is.null(genes.list)) stop("Error: genes.list is empty")

  if (is.character(genes.list)){
    genes.list=genes.list[genes.list %in% rownames(object@data)]
    all_features = genes.list
    feature.rownames = all_features
  } else if(is.list(genes.list)){
    all_features = c(); feature.rownames = c();
    for (l in 1:length(genes.list)){
      genes.list[[l]]=genes.list[[l]][genes.list[[l]] %in% rownames(object@data)]
      feature.rownames = c(feature.rownames, paste0(genes.list[[l]], collapse="__"))
      all_features = union(all_features, genes.list[[l]])
    }
  }
 
  ident.use = object@ident
  
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
  
  if (length(x = cells.1) == 0) {
    print(paste("Cell group 1 is empty - no cells with identity class", 
                ident.1))
    return(NULL)
  }
  if (length(x = cells.2) == 0) {
    print(paste("Cell group 2 is empty - no cells with identity class", 
                ident.2))
    return(NULL)
  }

  if (length(genes.list)==1 && !is.null(do.plot)){
    do.plot=TRUE;
    main.use = paste0(genes.list, " (Group ", ident.1, " vs Group ", ident.2, ")" )
  } else {
    do.plot=FALSE
    main.use=NULL
  }
  
  if (!is.null(max.cells.per.ident)){
    set.seed(42)
    if (length(cells.1) > max.cells.per.ident) cells.1 = sample(cells.1, max.cells.per.ident)
    
    if (length(cells.2) > max.cells.per.ident)  cells.2 = sample(cells.2, max.cells.per.ident)
  }

 
  to.return <- MarkerTestPR(object = object, cells.1 = cells.1, 
                            cells.2 = cells.2, genes.list = genes.list, print.bar = print.bar, 
                            thresh.use=thresh.use, feature.rownames=feature.rownames, do.plot=do.plot, main.use=main.use)


  return(to.return)
}



MarkerTestPR = function (object, cells.1, cells.2, genes.list = genes.list, print.bar = TRUE, 
                         thresh.use=0, feature.rownames=NULL,  do.plot=FALSE, main.use=NULL) 
{
  feature.rownames = set.ifnull(feature.rownames, genes.list)
  to.return <- AUCPRMarkerTest(data1 = object@data[, cells.1], 
                             data2 = object@data[, cells.2], genes.list = genes.list, 
                             print.bar = print.bar, thresh.use=0, feature.rownames=feature.rownames,
                             do.plot=do.plot, main.use=main.use)
  return(to.return)
}

AUCPRMarkerTest = function (data1, data2, genes.list, print.bar = TRUE, thresh.use=0, feature.rownames = NULL,  do.plot=FALSE, main.use=NULL) 
{
  feature.rownames = set.ifnull(feature.rownames, genes.list)
    AUCPR <- unlist(x = lapply(X = genes.list, FUN = function(x) {
    
    if (length(x) == 1){
      in_data1 = data1[x,]
      in_data2 = data2[x,]
      return(DifferentialAUCPR(x = as.numeric(x = in_data1), 
                           y = as.numeric(x = in_data2), do.plot = do.plot, main.use=main.use))
    } else {
      x1 = data1[x,]; x2 = data2[x,];
      temp1 = apply(x1 > thresh.use, 2, prod); temp2 = apply(x2 > thresh.use, 2, prod); 
      in_data1 = apply(x1, 2, mean) * temp1; in_data2 = apply(x2, 2, mean) * temp2
      
      return(DifferentialAUCPR(x = as.numeric(x=in_data1),
                               y = as.numeric(x=in_data2),do.plot = do.plot, main.use=main.use))
    }
  }))
  
  AUCPR[is.na(x = AUCPR)] <- 0
  if (print.bar) {
    iterate.fxn <- pblapply
  } else {
    iterate.fxn <- lapply
  }
  avg_diff <- unlist(x = iterate.fxn(X = genes.list, FUN = function(x) {
    
    if (length(x) == 1){
      in_data1 = data1[x,]
      in_data2 = data2[x,]
      return(expMean(x = as.numeric(x = in_data1)) - expMean(x = as.numeric(x = in_data2)))
    } else {
      x1 = data1[x,]; x2 = data2[x,];
      temp1 = apply(x1 > thresh.use, 2, prod); temp2 = apply(x2 > thresh.use, 2, prod); 
      in_data1 = apply(x1, 2, mean) * temp1; in_data2 = apply(x2, 2, mean) * temp2
      
      return(expMean(x = as.numeric(x = in_data1)) - expMean(x = as.numeric(x = in_data2)))
    }
    
  }))
  
  pct_clust1 <- unlist(x = iterate.fxn(X = genes.list, FUN = function(x) {
    
    if (length(x) == 1){
      in_data1 = data1[x,]
      return(sum(in_data1 > thresh.use) / length(in_data1))
    } else {
      x1 = data1[x,];
      temp1 = apply(x1 > thresh.use, 2, prod); 
      return(sum(temp1)/length(temp1))
    }
    
  }))
  
  pct_clust2 <- unlist(x = iterate.fxn(X = genes.list, FUN = function(x) {
    
    if (length(x) == 1){
      in_data2 = data2[x,]
      return(sum(in_data2 > thresh.use) / length(in_data2))
    } else {
      x2 = data2[x,]; 
      temp2 = apply(x2 > thresh.use, 2, prod); 
      return(sum(temp2)/length(temp2))
    }
    
  }))
  
  toRet <- data.frame(cbind(cbind(AUCPR, avg_diff), cbind(pct_clust1, pct_clust2)), row.names = feature.rownames)
  toRet <- toRet[rev(x = order(toRet$AUCPR)), ]
  return(toRet)
}

DifferentialAUCPR = function (x, y, do.plot=FALSE, main.use=NULL) 
{
  prediction.use <- ROCR::prediction(predictions = c(x, y), labels = c(rep(x = 1, 
                                                                     length(x = x)), rep(x = 0, length(x = y))), label.ordering = 0:1)
  
  perf1 <- ROCR::performance(prediction.use, "prec", "rec")

  is.nanvals = is.nan(perf1@x.values[[1]]) | is.nan(perf1@y.values[[1]])
  require(caTools)
  auc = trapz(perf1@x.values[[1]][!is.nanvals], 
        perf1@y.values[[1]][!is.nanvals])
  auc.use <- round(x = auc, digits = 3)
  if (do.plot){
    main.use = paste0(main.use, ", AUC = ", auc.use)
    print(plot(perf1, main=main.use))
  }
  return(auc.use)
}