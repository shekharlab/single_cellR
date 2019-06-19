set.ifnull=function(x,y) {
  if(is.null(x)) x=y
  return(x)
}
nogrid=theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

Dotplot_multiple = function(object,features.use=NULL, Count.mat = NULL, ident.use=NULL, thresh.use=0,do.plot=TRUE,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,cols.use=gray.colors(10),max.size=10,min.perc=0,norm.exp=NULL,...) {
  
  features.use=set.ifnull(features.use, object@var.genes[1:20])
  if (is.character(features.use)){
    features.use=features.use[features.use %in% rownames(object@data)]
    all_features = features.use
    feature.rownames = all_features
  } else if(is.list(features.use)){
    all_features = c(); feature.rownames = c();
    for (l in 1:length(features.use)){
      #if (length(features.use[[l]]) > 2){
      #  stop("Error: Must specify gene sets of length 2")
      #}
      features.use[[l]]=features.use[[l]][features.use[[l]] %in% rownames(object@data)]
      feature.rownames = c(feature.rownames, paste0(features.use[[l]], collapse="_"))
      all_features = union(all_features, features.use[[l]])
    }
  }
  
  ident.use=set.ifnull(ident.use, levels(object@ident))

  #Matrix of percent expressing cells
  PercMat = matrix(0, nrow=length(features.use), ncol = 0)
  rownames(PercMat) = feature.rownames; 
  
  #Matrix of average transcript levels
  ExpMat = PercMat;
  
  #Count mat
  Count.mat = set.ifnull(Count.mat, exp(object@data[all_features, colnames(object@data)]) - 1)
  
  for (i in ident.use){
    cells.in.cluster = names(object@ident)[which(object@ident== i)]
    if (is.character(features.use)){
      vec.perc = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) sum(x>thresh.use)/length(x)) 
      vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) if (sum(x>0) > 1){ mean(x[x>0]) } else {sum(x)})
    } else {
      vec.perc = c(); vec.exp=c()
      for (l in 1:length(features.use)){
        genes = features.use[[l]]
        if (length(genes) > 1){
            x1 = Count.mat[genes, cells.in.cluster]
            temp = apply(x1 > thresh.use, 2, prod)
            vec.perc = c(vec.perc, sum(temp)/length(temp))
            cells.all = cells.in.cluster[temp==1]
            
             if (length(cells.all)<2){
               vec.exp = c(vec.exp, 0)
             } else {
               vec.exp = c(vec.exp,mean(apply(x1[,cells.all],2, mean)))
             }
        } else {
          
          x1 = Count.mat[genes, cells.in.cluster]
          vec.perc = c(vec.perc, sum(x1 > thresh.use)/length(x1))
          vec.exp = c(vec.exp, if (sum(x1>0) > 1) { mean(x1[x1>0]) } else {sum(x1)} ) 
          
        }
      }
      names(vec.perc) = feature.rownames
      names(vec.exp) = feature.rownames
    }
    
    PercMat = cbind(PercMat,vec.perc)
    ExpMat = cbind(ExpMat, vec.exp)
  }
  colnames(ExpMat) = ident.use
  colnames(PercMat) = ident.use

  if (!is.null(norm.exp)){
    if (norm.exp < 0 | norm.exp > 1){
      print("Warning: norm.exp should be a value between (0,1). Skipping normalization")
      next
    } else{
      quant.vals = apply(ExpMat,1, function(x) quantile(x, norm.exp))
      ExpMat = t(scale(t(ExpMat), center=FALSE, scale=quant.vals))
    }
  }
  
  rows.use = rownames(PercMat)[apply(PercMat, 1, function(x) max(x) >= min.perc)];
  PercMat = PercMat[rows.use,]
  ExpMat = ExpMat[rows.use,]
  feature.rownames = rows.use
  if (!is.null(max.val.perc)) PercMat[PercMat > max.val.perc] = max.val.perc
  if (!is.null(max.val.exp)) ExpMat[ExpMat > max.val.exp] = max.val.exp
  
  if (do.plot){
    
    ExpVal = melt(ExpMat)
    PercVal = melt(PercMat)
    colnames(ExpVal) = c("gene","cluster","logExp")
    ExpVal$percExp = PercVal$value*100
    
    if (!do.transpose){
      ExpVal$gene = factor(ExpVal$gene, levels=rev(feature.rownames))
      ExpVal$cluster = factor(ExpVal$cluster, levels= ident.use)
      p=ggplot(ExpVal, aes(y = factor(gene),  x = factor(cluster))) + geom_point(aes(colour = logExp,  size =percExp)) + 
        scale_color_gradient(low ="blue",   high = "red", limits=c(0, max(ExpVal$logExp) ))+scale_size(range = c(1, max.size))+   theme_bw() +nogrid
      p = p + xlab("Cluster") + ylab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
        theme(axis.text.y=element_text(size=12, face="italic"))  
      print(p)
    } else {
      ExpVal$gene = factor(ExpVal$gene, levels=feature.rownames)
      ExpVal$cluster = factor(ExpVal$cluster, levels= rev(ident.use))
      p=ggplot(ExpVal, aes(y = factor(cluster),  x = factor(gene))) + geom_point(aes(colour = logExp,  size =percExp)) + 
        scale_color_gradient(low ="black",   high = "olivedrab2", limits=c( 0, max(ExpVal$logExp) ))+scale_size(range = c(1, max.size))+   theme_bw() +nogrid
      p = p + ylab("Cluster") + xlab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1)) + 
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
