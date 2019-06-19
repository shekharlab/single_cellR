source("~/Dropbox/CompRepos/sc_R_sparse/Oct2017/PermutationPA.R")
#source("./PermutationPA.R")

#nmf.options(grid.patch=TRUE)
nogrid=theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
sort.column=function(x, col) {
  return(x[order(x[,col]),])
}

# Function to map values in a vector `v` as defined in `from`` to the values
# defined in `to`.
#
# @param v     vector of values to map
# @param from  vector of original values
# @param to    vector of values to map original values to (should be of equal
#              length as from)
# @return      returns vector of mapped values
#
MapVals <- function(v, from, to){
  if (length(from) != length(to)) {
    stop("from and to vectors are not the equal length.")
  }
  vals.to.match <- match(v, from)
  vals.to.match.idx  <- !is.na(vals.to.match)
  v[vals.to.match.idx] <- to[vals.to.match[vals.to.match.idx]]
  return(v)
}

tsplot=function(object,x=1,cex.use=0.6, use.reduction="tsne", pc.1=1, pc.2=2, cols.use=NULL, id.use=NULL, 
                text.cex=1.25, xlim=NULL, ylim=NULL,do.legend=FALSE, do.segments=FALSE, segment_coord = NULL, RFlabels = NULL, cells.use=NULL) {
   
  cells.use = set.ifnull(cells.use, colnames(object@data))          
  
  if (is.null(cols.use)){
      cols.use=rainbow(length(levels(object@ident[cells.use]))); cols.use[x]="lightgrey"
  }
  order = sample(c(1:length(object@ident[cells.use])), replace=FALSE)
  ident.use = object@ident[cells.use][order]

  if (!is.null(id.use)){
    object = set.all.ident(object, id=id.use)
  }
  
  #ident.use= as.numeric(object@ident[cells.use])
 
  if (use.reduction == "tsne"){
    data.use0 = object@tsne.rot[cells.use, c(pc.1, pc.2)]
    data.use = data.use0[order,]
    xlab.use = paste0("tSNE",pc.1)
    ylab.use =  paste0("tSNE",pc.2)
  }
  if (use.reduction == "pca"){
    data.use0 = object@pca.rot[cells.use, c(pc.1, pc.2)]
    data.use = data.use0[order,]
    xlab.use = paste0("PC",pc.1)
    ylab.use =  paste0("PC",pc.2)
  } 

  #Plot
  plot(data.use[,1],data.use[,2],col=cols.use[as.integer(ident.use)],pch=16,xlab=xlab.use,ylab=ylab.use,cex=cex.use, xlim=xlim, ylim=ylim)

  #plot(object@tsne.rot[order,1],object@tsne.rot[order,2],col="lightgrey",pch=20,xlab="TSNE_1",ylab="TSNE_2",cex=cex.use)
  
  if (!is.null(RFlabels)){
    plot(data.use[,1],data.use[,2],col="gray",pch=16,xlab=xlab.use,ylab=ylab.use,cex=0.01, xlim=xlim, ylim=ylim)  
    l=1
    cols.use=rainbow(length(levels(RFlabels)))
    for(i in levels(RFlabels)){
      cells.clust = names(RFlabels)[RFlabels==i]
      text(data.use[cells.clust,1],
           data.use[cells.clust,2],
           as.character(RFlabels[cells.clust]),col=cols.use[l],cex=1,font=4)
      l=l+1
    }
  }
  k.centers=t(sapply(levels(object@ident),function(x) apply(object@tsne.rot[intersect(cells.use,which.cells(object,x)),],2,median)))
  points(k.centers[,pc.1],k.centers[,pc.2],cex=1.3,col="white",pch=16); 
  text(k.centers[,pc.1],k.centers[,pc.2],levels(object@ident[cells.use]),text.cex=1.25)
  
  
  if (do.legend){
    legend("topright", 
           legend = levels(object@ident), 
           col = cols.use, 
           pch = 16, 
           pt.cex = 2, 
           cex = 1.2, 
           text.col = "black", 
           horiz = F )
    
  }
  
  # Line segments for MNN diagnosis
  if (do.segments){
    # segment_coord is an Nx2 or Nx4 matrix
    if ("in.clust" %in% colnames(segment_coord)){
      col.segments = rep("black", nrow(segment_coord))
      col.segments[segment_coord[,"in.clust"] == 0] = "red"
    } else {
      col.segments = rep("black", nrow(segment_coord))
    }
    segments(data.use0[segment_coord[,1],1],data.use0[segment_coord[,1],2],
             data.use0[segment_coord[,2],1],data.use0[segment_coord[,2],2], col = col.segments)
  }
}

# diffplot=function(object,x=1,cex.use=3, use.reduction="tsne", pc.1=1, pc.2=2,pc.3=3, cols.use=NULL) {
#   if (is.null(cols.use)){
#     cols.use=rainbow(length(levels(object@ident))); cols.use[x]="lightgrey"
#   }
#   
#   
#   ident.use=as.numeric(object@ident)
#   if (use.reduction=="tsne") {
#     order1 = sample(c(1:dim(object@tsne.rot)[1]), replace=FALSE)
#     if (dim(object@tsne.rot)[2] < 3){
#       stop("Error : include at least 3 dimensions")
#     }
#     open3d()
#     plot3d(object@tsne.rot[order1,1],object@tsne.rot[order1,2],object@tsne.rot[order1,3],col=cols.use[as.integer(object@ident[order1])],xlab="phi_1",ylab="phi_2", zlab="phi_3",size=cex.use)
#     #plot(object@tsne.rot[order,1],object@tsne.rot[order,2],col="lightgrey",pch=20,xlab="TSNE_1",ylab="TSNE_2",cex=cex.use)
#     k.centers=t(sapply(levels(object@ident),function(x) apply(object@tsne.rot[which.cells(object,x),],2,mean)))
#     points3d(k.centers[,1],k.centers[,2],k.centers[,3],col="white", size=10); text3d(k.centers[,1],k.centers[,2],k.centers[,3],levels(object@ident),adj=c(1,1))
#     #for ( i in 1:max(ident.use) ){
#     #   if ( sum(ident.use == i) > 0 ) text(object@tsne.rot[ident.use == i,1],object@tsne.rot[ident.use == i,2],i,col=cols.use[i],cex=.75,font=4)
#     # }
#   }
#   
#   if (use.reduction=="pca") {
#     order = sample(c(1:dim(object@pca.rot)[1]), replace=FALSE)
#     if (dim(object@pca.rot)[2] < 3){
#       stop("Error : include at least 3 dimensions")
#     }
#     open3d()
#     plot3d(object@pca.rot[order,pc.1],object@pca.rot[order,pc.2],object@pca.rot[order,pc.3],col=cols.use[as.integer(object@ident[order])],xlab="phi_1",ylab="phi_2", zlab="phi_3",size=cex.use)
#     #plot(object@tsne.rot[order,1],object@tsne.rot[order,2],col="lightgrey",pch=20,xlab="TSNE_1",ylab="TSNE_2",cex=cex.use)
#     k.centers=t(sapply(levels(object@ident),function(x) apply(object@pca.rot[which.cells(object,x),c(pc.1,pc.2,pc.3)],2,mean)))
#     points3d(k.centers[,1],k.centers[,2],k.centers[,3],col="white", size=10); text3d(k.centers[,1],k.centers[,2],k.centers[,3],levels(object@ident),adj=c(1,1))
#     #for ( i in 1:max(ident.use) ){
#     #   if ( sum(ident.use == i) > 0 ) text(object@tsne.rot[ident.use == i,1],object@tsne.rot[ident.use == i,2],i,col=cols.use[i],cex=.75,font=4)
#     # }
#   }
# }

diffplot3d <- function(object, dims.use = c(1,2,3), reduction.use="dmap", cells.use=NULL, color.by = "ident", ...){
  cells.use = set.ifnull(cells.use, colnames(object@data))
  if (length(dims.use) != 2 & length(dims.use) !=3){
    stop("Error : dims.use must be of length 2 or 3")
  }
  data.plot <- as.data.frame(GetCellEmbeddings(object, 
                                  reduction.use = reduction.use, 
                                  dims.use = dims.use, 
                                  cells.use = cells.use))
  
  if (!(is.null(color.by))){
    if (color.by != "ident"){
       color.vals = fetch.data(object, vars.all = color.by, cells.use = cells.use, ...)
       color.vals = color.vals[,1]
       
  } else {
       color.vals = object@ident
       require(randomcoloR)
       nColors = length(unique(object@ident[cells.use]))
       col.use = unname(distinctColorPalette(nColors))
  }
  } else {
    color.vals=NULL
  }
  
  if (length(dims.use) == 3){
    
    if (color.by != "ident"){
      plot_ly(as.data.frame(data.plot),x = data.plot[,1], y = data.plot[,2], z = data.plot[,3], color = color.vals) %>% 
        layout(title = color.by, scene = list(xaxis = list(title = colnames(data.plot)[1]), 
                            yaxis = list(title = colnames(data.plot)[2]),
                            zaxis = list(title = colnames(data.plot)[3])))
    } else {
      plot_ly(as.data.frame(data.plot),x = data.plot[,1], y = data.plot[,2], z = data.plot[,3], color = color.vals, colors = col.use) %>% 
        layout(title = color.by, scene = list(xaxis = list(title = colnames(data.plot)[1]), 
                                              yaxis = list(title = colnames(data.plot)[2]),
                                              zaxis = list(title = colnames(data.plot)[3])))
    }
  } else {
    if (color.by != "ident"){
      plot_ly(as.data.frame(data.plot),x = data.plot[,1], y = data.plot[,2], color = color.vals[,1]) %>% 
        layout(title = color.by, scene = list(xaxis = list(title = colnames(data.plot)[1]), 
                                              yaxis = list(title = colnames(data.plot)[2])))
    } else {
      plot_ly(as.data.frame(data.plot),x = data.plot[,1], y = data.plot[,2],color = color.vals, colors = col.use) %>% 
        layout(title = color.by, scene = list(xaxis = list(title = colnames(data.plot)[1]), 
                                              yaxis = list(title = colnames(data.plot)[2])))
    }
  }


}

getLeftDecendants=function(tree,node) {
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  if (daughters[1] <= (tree$Nnode+1)) return(daughters[1])
  daughter.use=getDescendants(tree,daughters[1])
  daughter.use=daughter.use[daughter.use<=(tree$Nnode+1)]
  return(daughter.use)
}


getRightDecendants=function(tree,node) {
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  if (daughters[2] <= (tree$Nnode+1)) return(daughters[2])
  daughter.use=getDescendants(tree,daughters[2])
  daughter.use=daughter.use[daughter.use<=(tree$Nnode+1)]
  return(daughter.use)
}

getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w)) 
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}

situ3d=function(data, label=NULL, ...) {
  # Call Seurat function to get the in situ values out.
  exp.1=data
  exp.1=(exp.1-min(exp.1))/(max(exp.1)-min(exp.1))
  # Reformat them into an expression matrix as expected by the plotting function
  expression.matrix <- data.frame(matrix(exp.1, nrow=8, ncol=8))
  rownames(expression.matrix) <- c("24-30", "17-23", "13-16", "9-12", "7-8", "5-6", "3-4", "1-2")
  names(expression.matrix) <- c("1-4", "5-8", "9-12", "13-16", "17-20", "21-24", "25-28", "29-32")
  
  # Call the plotting function.
  zf.insitu.side(expression.matrix)
  par3d(windowRect=c(0, 0, 800, 800))
  
  # Label or not and then set the view.
  if (!(is.null(label))) {
    text3d(x=0, y=0, z=1.5, text=label, cex=3)
  }
  view3d(zoom=.75, theta=0, phi=-90, fov=0)
}

aucFxn=function(preds,truth,do.plot=FALSE,lab.main="",...) {
  pred.use=prediction(preds,truth,0:1)
  perf.use1=performance(pred.use,"tpr","fpr")
  perf.use=performance(pred.use,"auc")
  auc.use=round(perf.use@y.values[[1]],3)
  if (do.plot) plot(perf.use1@x.values[[1]],perf.use1@y.values[[1]],type="l",xlab="FPR",ylab="TPR",main=paste(lab.main, auc.use))
  return(auc.use)
}

humpMean=function(x, min=0) {
  return(mean(x[x>min]))
}

dispMeasure = function(x) return(var(x) / mean(x))
CoeffVar = function(x) return(sd(x) / mean(x))

humpdispMeasure = function(x, min=0) return(var(x[x>min]) / mean(x[x>min]))

WeightMatrix=function(n, sigma=10){
  a=c(1:n);b=a;
  W = exp(-(outer(b,a,FUN="-"))^2 /(n*sigma))
  return(W);
}

BackSPIN = function(X=X, max_splits=max_splits,SPINiter=SPINiter, min.cells=10, fxn.x=expMean, B.iter=10){
  
  splits=0;
  X0=X;
  N=dim(X0)[2]; # cells
  M=dim(X0)[1]; # genes
  isFinal = rep(0,N); #Indicator vector for: Is the cluster identity frozen?
  names(isFinal) = colnames(X0);
  
  col.ident = rep(0,N); names(col.ident)=colnames(X0); # cell cluster 
  row.ident = rep(0,M); names(row.ident)=rownames(X0); # gene cluster
  
  sort.order.col = colnames(X0);
  sort.order.row = rownames(X0);
  
  clust.ord = c(0); #Final cluster ordering to make things look pretty
  
  while ((splits < max_splits) && (sum(isFinal) < N)) {
    splits = splits+1
    col.ident0=col.ident;
    print(paste0("Split # ", splits))
    print(paste0("Currently ", length(unique(col.ident0)), " clusters"))
    #Iterate through all the current clusters
    for (clust in sort(unique(col.ident0))){
      #Order the sub-expression matrix (Genes by Cells) corresponding to this cluster
      cells.in.clust = names(which(col.ident==clust)) 
      cells.in.clust=cells.in.clust[order(match(cells.in.clust,sort.order.col))]
      genes.in.clust = names(which(row.ident==clust))
      genes.in.clust = genes.in.clust[order(match(genes.in.clust, sort.order.row))];
      
      #Freeze the cluster if the number of cells is already too low and move on to next cluster
      if (length(cells.in.clust) < 2*min.cells){
       # sort.order.row = c(sort.order.row, genes.in.clust)
        isFinal[cells.in.clust] = 1;
        next;
      }
      
      #Extract sub-expression matrix and compute the correlation matrices for cells and genes
      X = X0[genes.in.clust,cells.in.clust];
      Cc = cor(X)
      Cr = cor(t(X))
      
      #Check if the cluster is already frozen
      if (sum(isFinal[cells.in.clust]) == 0){
        #If not frozen order the cells by SPIN, order the matrix by SPIN
        sort.order.temp = SPIN_ordering(1-Cc, iter=SPINiter)
        sort.order.col[match(cells.in.clust,sort.order.col)] = cells.in.clust[sort.order.temp]
        
        #Find the optimal splits
        print(paste0("Finding split of cluster ", clust))
        colSplit = FindSplittingPoint(Cc[sort.order.temp,sort.order.temp], min.cells=min.cells, B = B.iter)
        if (colSplit == -1){ #No more splitting
          isFinal[cells.in.clust]=1;
        } else {
          print("Found a new cluster .. ")
          new.clust = max(col.ident) + 1;
          idx.curr = which(clust.ord==clust)
          clust.ord = append(clust.ord,new.clust, after=idx.curr)
          
          col.ident[cells.in.clust[sort.order.temp[c((colSplit+1):length(cells.in.clust))]]] = new.clust;
          
          #Get gene order
          genes.clust = BS.gene.split(Xc=X0[genes.in.clust,cells.in.clust[sort.order.temp]],col.ident.split = col.ident[cells.in.clust[sort.order.temp]], fxn.x = fxn.x);
          sort.order.row[match(genes.in.clust,sort.order.row)] = names(genes.clust)
          row.ident[names(genes.clust)] = genes.clust
          
        }
      }
      
      
    }
  }
  
  final.rows = c(); final.cols = c();
  cell.ident=col.ident; gene.ident = row.ident;
  l=1
  for (clust in clust.ord){
    
    cells.in.clust = names(col.ident)[col.ident==clust]
    genes.in.clust = names(row.ident)[row.ident==clust]
    cell.ident[cells.in.clust]=l; gene.ident[genes.in.clust]=l;
    
    final.cols = c(final.cols, cells.in.clust[order(match(cells.in.clust, sort.order.col))])
    final.rows = c(final.rows, genes.in.clust[order(match(genes.in.clust, sort.order.row))])
    l=l+1;                            
    
  }
  
  to.return=list()
  to.return$cell.order = final.cols
  to.return$gene.order = final.rows
  to.return$cell.clust = cell.ident
  to.return$gene.clust = gene.ident
  
  return(to.return)
  
}

BS.gene.split = function(Xc=NULL,col.ident.split=NULL, fxn.x=expMean){
   #Assign genes
   genes.use=rownames(Xc);
   ExpMat=matrix(0,nrow=length(genes.use),ncol=0)
   for (i in unique(col.ident.split)){
     cells.in.clust=names(col.ident.split)[col.ident.split==i];
     vec = apply(Xc[,cells.in.clust],1, fxn.x)
     ExpMat = cbind(ExpMat,vec)
   }
   rownames(ExpMat) = genes.use; colnames(ExpMat) = paste0("clust_",unique(col.ident.split));
   
   max.clust =  apply(ExpMat,1, function(x) unique(col.ident.split)[which.max(x)])
   max.clust = sort(max.clust)
   
   return(max.clust)
   
 }

FindSplittingPoint = function(R, min.cells=10, B=10){
  Fullsum=sum(R);
  N=dim(R)[1];
  mat.split = calcSplit(R, Fullsum, min.cells)
    
    Smax_vals = c()
    for (i in c(1:B)){
      reord = sample(N)
      temp.split = calcSplit(R[reord, reord], Fullsum, min.cells)
      Smax_vals = c(Smax_vals, temp.split$Smax)
    }
    
    pval = sum(Smax_vals > 0.9*mat.split$Smax) / B;
    #print(paste0("pval = ", pval))
  
    if ( pval < 0.05){
      return(mat.split$Isplit)
    } else {
      return(-1)
    }
  
  
}

calcSplit = function(R, Fullsum, min.cells){
  
  N=dim(R)[1];
  
  S= c();
  
  for (i in c(1:(N-1))){
    leftSum = sum(R[1:i,1:i]);
    rightSum = sum(R[(i+1):N,(i+1):N]);
    S = c(S, (leftSum+rightSum)/(i^2 + (N-i)^2));
  }
  Isplit = which.max(S);
  if (min(Isplit, N-Isplit) < min.cells){ Isplit = -1; S = -1}
  to.return = list();
  to.return$Isplit = Isplit;
  to.return$Smax = max(S);
  return(to.return)
  
}

#Sleft = c(); Sright = c();
#Sleft = c(Sleft, (leftSum/i^2) / (Fullsum/N^2));
#Sright = c(Sright,(rightSum / (N-i)^2) / (Fullsum/N^2) );

SPIN_ordering = function(D,iter=10){
  
  n = nrow(D);
  sigma=0.4*n;
  flag=1;
  sort.order = c(1:n);
  
  while (flag > 0){
    W = WeightMatrix(n, sigma=sigma)
    row.argsort = SPIN_neighborhood_sort(D,W,itermax=iter);
    sort.order = sort.order[row.argsort]
    Dnew = D[row.argsort,row.argsort]
    rm(D)
    D=Dnew;
    sigma=max(1,sigma/2);
    
    if (sigma==1){
      flag=flag-0.5;
    }
  }
  return(sort.order)
}

SPIN_neighborhood_sort = function(D,W,itermax=50){
  n=nrow(D);
  #Init
  t=0;
  flag=0;
  iter=0;
  M0 = D %*% W;
  px0 = c(1:n);
  
  M=M0;
  
  while (flag==0 & iter < itermax){
    iter=iter+1;
    #Solve linear assignment problem
    px = LinearAssignment(M)
    
    a1 = sum(diag(M[px,]))
    a0 = sum(diag(M[px0,]))
    
    #Compare
    if (a0==a1){
      flag=1;
    } else {
      px0=px;
      M0=M;
      M = D %*% W[order(px),]
    }
    
  }
  return(px)
}

humpVar=function(x, min=0) {
  return(var(x[x>min]))
} 

debugdmvnorm=function (x, mean = rep(0, p), sigma = diag(p), log = FALSE) 
{
  if (is.vector(x)) 
    x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  dec <- tryCatch(chol(sigma), error = function(e) e)
  tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
  rss <- colSums(tmp^2)
  logretval <- -(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
  names(logretval) <- names(x)
  logretval
}

slimdmvnorm=function (x, mean = rep(0, p), sigma = diag(p), log = FALSE) 
{
  x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  dec <- tryCatch(chol(sigma), error = function(e) e)
  tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
  rss <- colSums(tmp^2)
  logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
  names(logretval) <- rownames(x)
  logretval
}

slimdmvnorm_nosum=function (x, mean = rep(0, p), sigma = diag(p), log = FALSE) 
{
  x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  dec <- tryCatch(chol(sigma), error = function(e) e)
  tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
  rss <- colSums(tmp^2)
  logretval <- -(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
  names(logretval) <- rownames(x)
  logretval
}

compareBins=function(object,cell.use,bin.1,bin.2,bins.mu,bins.cov) {
  num.genes=ncol(bins.mu)
  genes.use=colnames(bins.mu)
  to.par=floor(sqrt(num.genes))+1
  par(mfrow=c(to.par,to.par))
  data.use=object@imputed[genes.use,cell.use]; names(data.use)=genes.use
  lik.diff=sort(unlist(lapply(genes.use,function(g)dnorm(data.use[g],bins.mu[bin.1,g],sqrt(bins.cov[[bin.1]][g,g]))/dnorm(data.use[g],bins.mu[bin.2,g],sqrt(bins.cov[[bin.2]][g,g])))))
  for(g in names(lik.diff)) {
    plot(0,0,type="n",xlim=c(0,8),ylim=c(0,2),main=paste(g, round(lik.diff[g],2)))
    lines(density(rnorm(10000,bins.mu[bin.1,g],sqrt(bins.cov[[bin.1]][g,g]))),lwd=2,col="black")
    lines(density(rnorm(10000,bins.mu[bin.2,g],sqrt(bins.cov[[bin.2]][g,g]))),lwd=2,col="red")
    points(data.use[g],0,cex=1.6,col="darkgreen",pch=16)
  }
  #rp()
}

subr=function(data,code) {
  return(data[grep(code,rownames(data)),])
}

cv=function(x)sd(x)/mean(x)

humpCt=function(x, min=0) {
  return(length(x[x>min]))
}


crushRNANewNoLog=function(mynames,mycodes) {
  r=c()
  for(i in 1:length(mynames)) {
    name=mynames[i]
    file = paste("~/big/",name,"/",name,".rsem.all.6.res",sep="")
    input=read.table(file,header=TRUE)
    #input[,8:ncol(input)]=log(input[,8:ncol(input)]+1)
    code=mycodes[i]
    colnames(input)[8:ncol(input)]=paste(code,colnames(input)[8:ncol(input)],sep="")
    if (i==1) {
      r=input
    }
    if (i>1) {
      r=cbind(r,input[,8:ncol(input)])
    }
    print(name)
  }
  return(r)
}

log_add=function(x) {
  mpi=max(x)
  return(mpi+log(sum(exp(x-mpi))))
}

minusr=function(data,code) {
  matchCode=rownames(data)[grep(code,rownames(data))]
  toIgnore=which(rownames(data) %in% matchCode)
  if (length(toIgnore)==0) return(data)
  toRet=data.frame(data[-toIgnore,])
  rownames(toRet)=rownames(data)[-toIgnore]
  colnames(toRet)=colnames(data)
  return(toRet)                   
}

minusc=function(data,code) {
  matchCode=colnames(data)[grep(code,colnames(data))]
  toIgnore=which(colnames(data) %in% matchCode)
  if (length(toIgnore)==0) return(data)
  return(data[,-toIgnore])
}

getUMI=function(proj,code,cols=4,fullproj="none",nameCol=7,nameCode="rsem",subCode="umi",rmLast=FALSE, verbose=TRUE) {
  myCmd = paste("find ~/big/",proj,"/bam -name ", code, sep="")  
  if (fullproj != "none") {
    myCmd = paste(fullproj," -name ", code, sep="")  
  }
  myUMI=system(myCmd,intern=TRUE)
  if (verbose) print(myUMI)
  myNames=unlist(lapply(myUMI,function(x)gsub(nameCode,subCode,strsplit(x,"/")[[1]][nameCol])))
  data=read.table(myUMI[1],sep="\t",fill=TRUE)[,c(1,cols)]
  if (rmLast) {
    data=data[-nrow(data),]
  }
  colnames(data)[ncol(data)]=myNames[1]
  rownames(data)=as.character(data$V1)
  for(i in 2:length(myUMI)) {
    if (verbose) print(paste(myUMI[i],i))
    #print(data)
    if ((file.info(myUMI[i])$size) > 0) {
      data2=read.table(myUMI[i],sep="\t",fill=TRUE)[,c(1,cols)]
      if (rmLast) {
        data2=data2[-nrow(data2),]
      }
      data=merge(data,data2,by="V1",all=TRUE)
      colnames(data)[ncol(data)]=myNames[i]
    }
  }
  #data2=data
  data[is.na(data)]=0
  data2=aggregate(data[,-1],list(Gene = toupper(data$V1)),sum)
  rownames(data2)=toupper(data2$Gene)
  colnames(data2)=unlist(lapply(colnames(data2),function(x)gsub("-","_",x)))
  data2=subc(data2,subCode)
  return(data2)
}

ainb=function(a,b) {
  a2=a[a%in%b]
  return(a2)
}
meanNormFunction=function(data,myfuncX,myfuncY,nBin=20) {
  data_x=apply(data,1,myfuncX)
  data_y=apply(data,1,myfuncY)
  data_x_bin=cut(data_x,nBin)
  names(data_x_bin)=names(data_x)
  mean_y=tapply(data_y,data_x_bin,mean)
  sd_y=tapply(data_y,data_x_bin,sd)
  return((data_y-mean_y[as.numeric(data_x_bin)])/sd_y[as.numeric(data_x_bin)])
}

shift.cell =function(bin,x,y) {
  bin.y=(bin-1)%/%8+1
  bin.x=(bin-1)%%8+1
  new.x=minmax(bin.x+x,min = 1,max=8)
  new.y=minmax(bin.y+y,min = 1,max=8)
  new.bin=8*(new.y-1)+new.x
  return(new.bin)
}

empP=function(x,nullval) {
  return(length(which(nullval>x))/length(nullval))
}

neighbor.cells=function(bin) {
  return(unique(c(bin,shift.cell(bin,0,1),shift.cell(bin,1,0),shift.cell(bin,-1,0),shift.cell(bin,0,-1))))
}

all.neighbor.cells=function(bin,dist=1) {
  all.comb=expand.grid(rep(list(-dist:dist), 2)) 
  return(unique(unlist(lapply(1:nrow(all.comb),function(x)shift.cell(bin,all.comb[x,1],all.comb[x,2])))))
}

no.legend.title=theme(legend.title=element_blank())
ggplot.legend.text=function(x=12,y="bold") return(theme(legend.text = element_text(size = x, face = y)))
gg.legend.pts=function(x=6) guides(colour = guide_legend(override.aes = list(size=x)))
gg.xax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.x = element_text(face=z, colour=y, size=x), 
                                                              axis.text.x  = element_text(angle=90, vjust=0.5, size=x2)))

gg.yax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.y = element_text(face=z, colour=y, size=x), 
                                                              axis.text.y  = element_text(angle=0, vjust=0.5, size=x2)))

sub.string=function(x,s1,s2) return(unlist(lapply(x,function(y)gsub(s1,s2,y))))

jackStrawF=function(prop=0.1,myR1,myR2=3,data=smD) {
  randGenes=sample(rownames(data),nrow(data)*prop)
  smD.mod=data
  smD.mod[randGenes,]=shuffleMatRow(data[randGenes,])
  fmd.pca=prcomp(smD.mod)
  fmd.x=fmd.pca$x
  fmd.rot=fmd.pca$rotation
  fakeF=unlist(lapply(randGenes,jackF,r1=myR1,r2=myR2,x=fmd.x,rot=fmd.rot))
  return(fakeF)
}

jackF=function(gene,r1=1,r2=2,x=md.x,rot=md.rot) {
  if (r2==1) { #assuming r1, r2=1
    mod.x=x[,r1]
    mod.x[gene]=0
    return(var.test((x[,r1]%*%t(rot[,r1])),(mod.x%*%t(rot[,r1])))$statistic)
  }
  mod.x=x[,1:r2]
  mod.x[gene,r1:r2]=rep(0,r2-r1+1)
  return(var.test((x[,1:r2]%*%t(rot[,1:r2])),(mod.x[,1:r2]%*%t(rot[,1:r2])))$statistic)
}

shuffleMatRow=function(x) {
  x2=x
  x2 <- t(x)
  ind <- order(c(col(x2)), runif(length(x2)))
  x2 <- matrix(x2[ind], nrow=nrow(x), ncol=ncol(x), byrow=TRUE)
  return(x2)
}

GenesTranscriptsHist=function(ExpMat, UMI.data=TRUE, is.expr=0, do.bar=FALSE) {
  num.transcripts = colSums(ExpMat)
  num.genes = colSums(ExpMat > is.expr)
  num.cells = rowSums(ExpMat > is.expr)
  
  data = data.frame(nGene = num.genes, nTrans = num.transcripts, cell.names = colnames(ExpMat), stringsAsFactors = FALSE)
  data.gene = data.frame(nCells = num.cells)
  
  if (!(do.bar)){
    if (UMI.data){
      #Data consists of UMIs
      pList = list()
      pList[[1]]=ggplot(data, aes(x=nTrans)) + geom_histogram(color="black",fill="red") + theme_bw() + #scale_x_log10() +
        xlab("Num. Transcripts per cell") + ylab("Frequency")
      pList[[2]]=ggplot(data, aes(x=nGene)) + geom_histogram(color="black",fill="blue") +  theme_bw() + #scale_x_log10() +
        xlab("Num. Genes per cell") + ylab("Frequency")
      pList[[3]] = ggplot(data, aes(x=nTrans, y=nGene)) + geom_point(fill="black") +  theme_bw() + #scale_x_log10() + scale_y_log10() +
        xlab("Num. Transcripts per cell") + ylab("Num. Genes per cell")
      pList[[4]] = ggplot(data.gene, aes(x=nCells)) + geom_histogram(color="black",fill="blue") + theme_bw() +  #scale_x_log10() +
        xlab("Num. Cells per gene") + ylab("Frequency")
    } else {
      #Data consists of read counts
      pList = list()
      pList[[1]]=ggplot(data, aes(x=nTrans)) + geom_histogram(color="black",fill="red")  + theme_bw() + #scale_x_log10() +
        xlab("Num. Reads per cell") + ylab("Frequency")
      print(1)
      pList[[2]]=ggplot(data, aes(x=nGene)) + geom_histogram(color="black",fill="blue") +  theme_bw() + #scale_x_log10() +
        xlab("Num. Genes per cell") + ylab("Frequency")
      pList[[3]] = ggplot(data, aes(x=nTrans, y=nGene)) + geom_point(fill="black")  + theme_bw() + #scale_y_log10() +  scale_x_log10() +
        xlab("Num. Reads per cell") + ylab("Num. Genes per cell")
      pList[[4]] = ggplot(data.gene, aes(x=nCells)) + geom_histogram(color="black",fill="blue") +  theme_bw() + #scale_x_log10() +
        xlab("Num. Cells per gene") + ylab("Frequency")
      
    }
    
    multiplotList(pList, cols=2)
    rp()
  } else {
    if (UMI.data){
      print("#Data consists of UMIs")
      pList = list()
      
      data1 = transform(data, cell.names = reorder(cell.names, nTrans))
      pList[[1]]=ggplot(data1, aes(x=cell.names,y=nTrans)) + geom_bar(stat="identity", fill="red",color="black",width=0.9) + theme_bw() + #scale_y_log10() + 
        xlab("Cells (Low to high transcripts)") + ylab("Num. Transcripts per cell")
      if (dim(data1)[1] > 10){
        pList[[1]] = pList[[1]] + scale_x_discrete(breaks=NULL)
      }
      
      pList[[2]]=ggplot(data1, aes(x=cell.names, y=nGene)) + geom_bar(stat="identity", fill="blue",color="black",width=0.9)  + theme_bw() + 
        xlab("Cells (Low to high transcripts)") + ylab("Num. Genes per cell")
      if (dim(data1)[1] > 10){
        pList[[2]] = pList[[2]] + scale_x_discrete(breaks=NULL)
      }
      
      pList[[3]] = ggplot(data, aes(x=nTrans, y=nGene)) + geom_point(fill="black") + theme_bw() + # scale_x_log10() + scale_y_log10() +
        xlab("Num. Transcripts per cell") + ylab("Num. Genes per cell")
      
      pList[[4]] = ggplot(data.gene, aes(x=nCells)) + geom_histogram(color="black", fill="blue") +  theme_bw() + #scale_x_log10() +
        xlab("Num. Cells per gene") + ylab("Frequency")
    } else {
      #Data consists of read counts
      pList = list()
      data1 = transform(data, cell.names = reorder(cell.names, nTrans))
      pList[[1]]=ggplot(data1, aes(x=cell.names,y=nTrans)) + geom_bar(stat="identity", fill="red",color="black",width=0.9) + theme_bw() + #scale_y_log10() + 
        xlab("Cells (Low to high reads)") + ylab("Num. Reads per cell")
      if (dim(data1)[1] > 10){
        pList[[1]] = pList[[1]] + scale_x_discrete(breaks=NULL)
      }
      
      pList[[2]]=ggplot(data1, aes(x=cell.names, y=nGene)) + geom_bar(stat="identity", fill="blue",color="black",width=0.9) + theme_bw() + 
        xlab("Cells (Low to high reads)") + ylab("Num. Genes per cell")
      if (dim(data1)[1] > 10){
        pList[[2]] = pList[[2]] + scale_x_discrete(breaks=NULL)
      }
      
      pList[[3]] = ggplot(data, aes(x=nTrans, y=nGene)) + geom_point(fill="black") +  theme_bw() + #scale_x_log10() + scale_y_log10() +
        xlab("Num. Reads per cell") + ylab("Num. Genes per cell")
      pList[[4]] = ggplot(data.gene, aes(x=nCells)) + geom_histogram(color="black", fill="blue")  + theme_bw() + #scale_x_log10() +
        xlab("Num. Cells per gene") + ylab("Frequency")
      
    }
    
    multiplotList(pList, cols=2)
    rp()
     
  }
  
}




logMeanMinus= function(x)log(mean(exp(as.numeric(x))-1)+1)
logVarMinus= function(x)(mean(var(as.numeric(x))-1)+1)

logVarMinus2= function(x)(var(exp(as.numeric(x))-1)+1)

quickRNAHuman=function(x) {
  dataFile=paste("~/big/",x,"/summary/",x,".expMatrix.txt",sep="")
  data=log(read.table(dataFile,header=TRUE,sep="\t")[,-1]+1)
  data=subc(data,"rsem")
  return(data)
}

expVar=function(x) {
  return(log(var(exp(x)-1)+1))
}

expSD=function(x) {
  return(log(sd(exp(x)-1)+1))
}

expMean=function(x) {
  return(log(mean(exp(x)-1)+1))
}

perc.pos=function(x, thres=0) {
  return(sum(x>thresh)/length(x))
}

expMedian=function(x) {
  return(log(median(exp(x)-1)+1))
}

quickRNAZfish=function(x) {
  zdict=read.table("~/window/annotate/danRer7ensgene110512.txt.short.txt.tbl.txt",sep="\t",header=TRUE)
  rownames(zdict)=as.character(zdict$name)
  #zcan=subset(zdict,iscanonical==TRUE)
  #zdt=zcan[!duplicated(toupper(zcan$geneSymbol)),]
  #rownames(zdt)=toupper(zdt$geneSymbol)
  #zdt=zdt[,1:8]
  dataFile=paste("~/big/",x,"/",x,".rsem.iso.all.3.res",sep="")
  rawFin=read.table(dataFile,sep="\t",header=TRUE)
  rawFin$gene=toupper(zdict[as.character(rawFin$V1),"geneSymbol"])
  rF=rawFin
  rF=subc(rawFin,"rsem|gene")
  rF=aggregate(subc(rawFin,"rsem"),list(Gene = rawFin$gene),sum)
  rownames(rF)=rF$Gene
  fin=log(subc(rF,"rsem")*1e6+1)
  return(fin)
}

normal.sample=function(x) {
  return(rnorm(10000,mean(x),sd(x)))
  
}

fetch.closest=function(bin,all.centroids,num.cell) {
  bin.y=(bin-1)%/%8+1
  bin.x=(bin-1)%%8+1
  all.centroids=rbind(all.centroids,c(bin.x,bin.y))
  all.dist=as.matrix(dist(all.centroids))
  return(names(sort(all.dist[nrow(all.dist),]))[2:(num.cell+2)])
}

fetch.mincells=function(bin,cells.max,min.cells) {
  for(i in 1:5) {
    my.names=names(ainb(cells.max,all.neighbor.cells(bin,i)))
    if (length(my.names) > min.cells) break;
  }
  return(my.names)
}

cell.centroid=function(cell.probs) {
  centroid.x=round(sum(sapply(1:64,function(x)(x-1)%%8+1)*cell.probs))
  centroid.y=round(sum(sapply(1:64,function(x)(x-1)%/%8+1)*cell.probs))
  centroid.bin=8*(centroid.y-1)+centroid.x
  return(centroid.bin)
}

cell.centroid.x=function(cell.probs) {
  return(centroid.x=round(sum(sapply(1:64,function(x)(x-1)%%8+1)*cell.probs)))
}

cell.centroid.y=function(cell.probs) {
  return(centroid.y=round(sum(sapply(1:64,function(x)(x-1)%/%8+1)*cell.probs)))
}

exact.cell.centroid=function(cell.probs) {
  centroid.x=(sum(sapply(1:64,function(x)(x-1)%%8+1)*cell.probs))
  centroid.y=(sum(sapply(1:64,function(x)(x-1)%/%8+1)*cell.probs))
  return(c(centroid.x,centroid.y))
}


marker.auc.test=function(data1,data2,mygenes) {
  myAUC=unlist(lapply(mygenes,function(x)diffAUC(as.numeric(data1[x,]),as.numeric(data2[x,]))))
  myAUC[is.na(myAUC)]=0
  myDiff=unlist(lapply(mygenes,function(x)(expMean(as.numeric(data1[x,]))-expMean(as.numeric(data2[x,])))))
  toRet=data.frame(cbind(myAUC,myDiff),row.names=mygenes)
  toRet=toRet[rev(order(toRet$myAUC)),]
  return(toRet)
}  

#credit to Cole Trapnell for this
tobit_fitter <- function(x, modelFormulaStr, lower=1, upper=Inf){
  tryCatch({
    FM_fit <-  suppressWarnings(vgam(as.formula(modelFormulaStr), family=tobit(Lower=lower, Upper=upper),data = x))
    FM_fit
  }, 
  #warning = function(w) { FM_fit },
  error = function(e) { NULL }
  )
}

anotinb=function(x,y) {
  x2=x[!x%in%y]
  return(x2)
}
diffTobit=function(x1,x2,lower=1,upper=Inf) {
  my.df=data.frame(c(x1,x2),c(rep(0,length(x1)),rep(1,length(x2))))
  colnames(my.df)=c("Expression","Stat")
  #model.v1=vgam(Expression~1,family = tobit(Lower = lower,Upper = upper),data = my.df)
  model.v1=tobit_fitter(my.df,"Expression~1",lower,upper)
  #model.v2=vgam(Expression~Stat+1,family = tobit(Lower = lower,Upper = upper),data = my.df)
  model.v2=tobit_fitter(my.df,"Expression~Stat+1",lower,upper)
  p=1
  if (is.null(model.v1) == FALSE && is.null(model.v2) == FALSE) {
      (p <- pchisq(2 * (logLik(model.v2) - logLik(model.v1)), df = 1, lower.tail = FALSE))
  }
  return(p)
}

tobit.diffExp.test=function(data1,data2,mygenes) {
  myP=unlist(lapply(mygenes,function(x)diffTobit(as.numeric(data1[x,]),as.numeric(data2[x,]))))
  myP[is.na(myP)]=1
  myDiff=unlist(lapply(mygenes,function(x)(expMean(as.numeric(data1[x,]))-expMean(as.numeric(data2[x,])))))
  toRet=data.frame(cbind(myP,myDiff),row.names=mygenes)
  toRet=toRet[order(toRet$myP),]
  return(toRet)
}

bimod.diffExp.test=function(data1,data2,mygenes, xmin=1) {
  pval=unlist(lapply(mygenes,function(x) diffLRT(as.numeric(data1[x,]),as.numeric(data2[x,]), xmin=xmin)))
  pval[is.na(pval)]=1
  diff=unlist(lapply(mygenes,function(x)(expMean(as.numeric(data1[x,]))-expMean(as.numeric(data2[x,])))))
  toRet=data.frame(cbind(pval,diff),row.names=mygenes)
  toRet=toRet[order(toRet$pval),]
  return(toRet)
}

diffAUC = function(x,y) {
  prediction.use=prediction(c(x,y),c(rep(1,length(x)),rep(0,length(y))),0:1)
  perf.use=performance(prediction.use,"auc")
  auc.use=round(perf.use@y.values[[1]],3)
  return(auc.use)
}

diffLRT = function(x,y,xmin=1) {
  lrtX=bimodLikData(x, xmin=xmin)
  lrtY=bimodLikData(y, xmin=xmin)
  lrtZ=bimodLikData(c(x,y), xmin=xmin)
  lrt_diff=2*(lrtX+lrtY-lrtZ)
  return(1-pchisq(lrt_diff,3))
}

bimodLikData=function(x,xmin=0) {
  x1=x[x<=xmin]
  x2=x[x>xmin]
  xal=minmax(length(x2)/length(x),min=1e-5,max=(1-1e-5))
  likA=length(x1)*log(1-xal)
  mysd=sd(x2)
  if(length(x2)<2) {
    mysd=1
  }
  likB=length(x2)*log(xal)+sum(dnorm(x2,mean(x2),mysd,log=TRUE))
  return(likA+likB)
}

makeAlnPlot=function(proj) {
  alnFile=paste("~/big/",proj,"/summary/",proj,".all.aln.metrics.txt",sep="")
  alnData=read.table(alnFile)
  par(mar=c(10,5,4,1))
  mymax=max(100,max(apply(alnData[,4:7],1,sum)))
  x=barplot(as.matrix(t(alnData[,4:7])),names.arg=alnData$V1,las=2,col=1:4,ylim=c(0,mymax))
  text(x,mymax-10,alnData$V2)
  rp()
}

getCoefs=function(data) {
  my_stats=data.frame(data[,1])
  my_stats$code_humpAvg=apply(data,1,humpMean,min=1)
  rownames(my_stats)=rownames(data)
  my_coefs=data.frame(t(sapply(colnames(data),getAB,data=data,data2=my_stats,status="code",code2="humpAvg",doPlot=FALSE)))
  colnames(my_coefs)=c("a","b")
  return(my_coefs)
}

meanVarPlot=function(x,...) {
  myMean=apply(x,1,logMeanMinus)
  myVar=apply(x,1,logVarMinus)
  plot(myMean,myVar)
}

vsubc=function(data,code) {
  return(data[grep(code,names(data))])
}

vminusc=function(data,code) {
  matchCode=names(data)[grep(code,names(data))]
  toIgnore=which(names(data) %in% matchCode)
  if (length(toIgnore)==0) return(data)
  return(data[-toIgnore])
}

plotVln=function(gene,data=dc2,code="rsem",mmax=12,getStat=getStat1,doRet=FALSE,doSort=FALSE) {
  data$GENE=as.character(rownames(data))
  a1=data[gene,c(colnames(data)[grep(code,colnames(data))],"GENE")]
  a2=melt(a1,id="GENE")
  a2$stat=unlist(lapply(as.character(a2$variable),getStat))
  noise <- rnorm(length(a2$value))/100000
  a2$value=a2$value+noise
  if(doSort) {
    a2$stat=factor(a2$stat,levels=names(rev(sort(tapply(a2$value,a2$stat,mean)))))
  }
  p=ggplot(a2,aes(factor(stat),value))
  p2=p + geom_violin(scale="width",adjust=0.75,aes(fill=factor(stat))) + ylab("Expression level (log TPM)")
  p3=p2+theme(legend.position="top")+guides(fill=guide_legend(title=NULL))+geom_jitter(height=0)+blackbg+xlab("Cell Type")
  p4=p3+ theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=90, vjust=0.5, size=12))
  p5=(p4+theme(axis.title.y = element_text(face="bold", colour="#990000", size=16), axis.text.y  = element_text(angle=90, vjust=0.5, size=12))+ggtitle(gene)+theme(plot.title = element_text(size=20, face="bold")))
  if(doRet==TRUE) {
    return(p5)
  }
  else {
    print(p5)
  }
}

getStat1=function(x)return(strsplit(x,"_")[[1]][1])

getStat=function(x,y=1) return(strsplit(x,"_")[[1]][y])
getStat2=function(x,y=2) return(strsplit(x,"_")[[1]][y])
getStat3=function(x,y=3) return(strsplit(x,"_")[[1]][y])


multiplotList <- function(plots, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  #plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

subc=function(data,code) {
  return(data[,grep(code,colnames(data))])
}

minmax=function(data,min,max) {
  data2=data
  data2[data2>max]=max
  data2[data2<min]=min
  return(data2)
}

arrange <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
  dots <- list(...)
  n <- length(dots)
  if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
  if(is.null(nrow)) { nrow = ceiling(n/ncol)}
  if(is.null(ncol)) { ncol = ceiling(n/nrow)}
  ## NOTE see n2mfrow in grDevices for possible alternative
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
  ii.p <- 1
  for(ii.row in seq(1, nrow)){
    ii.table.row <- ii.row
    if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
    for(ii.col in seq(1, ncol)){
      ii.table <- ii.p
      if(ii.p > n) break
      print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
      ii.p <- ii.p + 1
    }
  }
}

vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
  dots <- list(...)
  n <- length(dots)
  if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
  if(is.null(nrow)) { nrow = ceiling(n/ncol)}
  if(is.null(ncol)) { ncol = ceiling(n/nrow)}
  ## NOTE see n2mfrow in grDevices for possible alternative
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
  ii.p <- 1
  for(ii.row in seq(1, nrow)){
    ii.table.row <- ii.row
    if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
    for(ii.col in seq(1, ncol)){
      ii.table <- ii.p
      if(ii.p > n) break
      print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
      ii.p <- ii.p + 1
    }
  }
}


subSort = function(vdat, my, fNum, sortBy) {
  vNames=colnames(vdat)
  mycol=which(vNames==sortBy)
  v2 = vdat[order(vdat[,mycol]),]
  vsort=v2
  return(v2)
}

calcMedians = function(p, start, end) {
  medians=c()
  for(i in start:end) {
    scores=p[which(p$factorNum==i),i+6]
    medians=c(medians,mean(scores))
  }
  return(medians)
}

myPalette=
  function (low = "white", high = c("green", "red"), mid = NULL,
            k = 50)
  {
    low <- col2rgb(low)/255
    high <- col2rgb(high)/255
    if (is.null(mid)) {
      r <- seq(low[1], high[1], len = k)
      g <- seq(low[2], high[2], len = k)
      b <- seq(low[3], high[3], len = k)
    }
    if (!is.null(mid)) {
      k2 <- round(k/2)
      mid <- col2rgb(mid)/255
      r <- c(seq(low[1], mid[1], len = k2), seq(mid[1], high[1],
                                                len = k2))
      g <- c(seq(low[2], mid[2], len = k2), seq(mid[2], high[2],
                                                len = k2))
      b <- c(seq(low[3], mid[3], len = k2), seq(mid[3], high[3],
                                                len = k2))
    }
    rgb(r, g, b)
  }

bwCols=myPalette(low = "white",high="black",k = 50)


comparePCA = function(a, b) {
  inds=c(6:26)
  pa=prcomp(a[,inds],scale=TRUE,center=TRUE)
  pb=prcomp(b[,inds],scale=TRUE,center=TRUE)
  print(summary(pa))
  print(summary(pb))
}

getSmooth = function(vsort, myBin, n=0, smooth=0, overlap=0.9, type=0) {
  if (smooth==0) {
    delta=1/(1-overlap)
    nbin=round((length(vsort[,1])-delta)/(n-1+delta))
    smooth=(nbin+1)*delta
    
  }
  return(smooth)
}

condInt = function(vsort, myCtrl, myX, myY, n=0, smooth=0, over=0.9, type=0) {
  overlap=over
  vNames=colnames(vsort)
  myCol=which(vNames==myCtrl)
  xCol=which(vNames==myX)
  yCol=which(vNames==myY)
  binTotal = c()
  nbin=round((length(vsort[,myCol])-smooth)/(n-1)-0.5)
  if (smooth==0) {
    delta=1/(1-overlap)
    nbin=round((length(vsort[,myCol])-delta)/(n-1+delta))
    smooth=(nbin+1)*delta
    
  }
  if (n==0) {
    n=round((length(vsort[,myCol])-smooth)/(smooth*(1-overlap)))
    nbin=round((length(vsort[,myCol])-smooth)/(n-1)-0.5)
  }
  vsort=vsort[order(vsort[,myCtrl]),]
  for (i in 1:n) {
    q1=nbin*(i-1)+1
    q2=q1+smooth
    if (q2 > length(vsort[,1])) {
      q2= length(vsort[,1])
    }
    if (type==0) {
      w=vsort[q1:q2,]
      w=w[order(w[,myX]),]
      x=calcBin(w,myX,n=100,overlap=over)
      y=calcBin(w,myY,n=100,overlap=over)
      val=cor(x,y)
      binTotal=c(binTotal,val)
    }
    if (type==7) {
      w=vsort[q1:q2,]
      w=w[order(w[,myX]),]
      x=calcBin(w,myX,n=100,overlap=over)
      y=calcBin(w,myY,n=100,overlap=over)
      ct=cor.test(x,y)
      width=(ct$conf.int[2]-ct$conf.int[1])/4
      binTotal=c(binTotal,width)
    }
    if (type==5) {
      w=vsort[q1:q2,]
      w=w[order(w[,myX]),]
      x=calcBin(w,myX,n=100)
      y=calcBin(w,myY,n=100)
      y2=calcBin(w,vNames[yCol+1],n=100)
      val=cor(x,y/(y+y2))
      binTotal=c(binTotal,val)
    }
    if (type==6) {
      w=vsort[q1:q2,]
      w=w[order(w[,myX]),]
      x=calcBin(w,myX,n=100)
      y=calcBin(w,myY,n=100,type=2)
      val=cor(x,y)
      binTotal=c(binTotal,val)
    }
  }
  return(binTotal)
}


calcBin = function(vsort, myBin, n=0, smooth=0, overlap=0.9, type=0,cut=0) {
  vNames=colnames(vsort)
  myCol=which(vNames==myBin)
  binTotal = c()
  nbin=round((length(vsort[,myCol])-smooth)/(n-1)-0.5)
  if (smooth==0) {
    delta=1/(1-overlap)
    nbin=round((length(vsort[,myCol])-delta)/(n-1+delta))
    smooth=(nbin+1)*delta
    
  }
  if (n==0) {
    n=round((length(vsort[,myCol])-smooth)/(smooth*(1-overlap)))
    nbin=round((length(vsort[,myCol])-smooth)/(n-1)-0.5)
  }
  
  for (i in 1:n) {
    q1=nbin*(i-1)+1
    q2=q1+smooth
    if (q2 > length(vsort[,1])) {
      q2= length(vsort[,1])
    }
    if (type==0) {
      binTotal=c(binTotal,mean(vsort[q1:q2,myCol]))
    }
    if (type==1) {
      binTotal=c(binTotal,median(vsort[q1:q2,myCol]))
    }
    if (type==2) {
      binTotal=c(binTotal,length(which(vsort[q1:q2,myCol]>cut))/(q2-q1+1))
    }
    if (type==3) {
      binTotal=c(binTotal,sd(vsort[q1:q2,myCol])/sqrt(q2-q1))
    }
    if (type==4) {
      myBin=vsort[q1:q2,]
      myBin=myBin[order(myBin[,myCol]),]
      yv=calcBin(myBin,"CagPkCons",n=100,overlap=0.9)
      xv=calcBin(myBin,vNames[myCol],n=100,overlap=0.9)
      lmod=coef(lm(yv ~ xv))
      slope=cor(xv,yv)
      binTotal=c(binTotal,slope)
    }
    if (type==5) {
      myBin=vsort[q1:q2,]
      myBin=myBin[order(myBin[,myCol]),]
      yv=calcBin(myBin,"CagPkCons",n=100,overlap=0.7)
      xv=calcBin(myBin,vNames[myCol],n=100,overlap=0.7)
      lmod=coef(lm(yv ~ xv))
      ct=cor.test(xv,yv)
      width=(ct$conf.int[2]-ct$conf.int[1])/4
      binTotal=c(binTotal,width)
    }
    
  }
  return(binTotal)
}

genCols=function(al=50) {
  cols=c("darkblue","darkred","darkgreen", "black","orange","purple","khaki","grey","gold4","seagreen3","chocolate")
  tcols=c()
  for (i in 1:6) {
    tcols=c(tcols,rgb (t (col2rgb (cols[i])), alpha = al, maxColorValue = 255))
  }
  return(tcols)
}

plos = function(xvals,yvals,dev,colNum,lwid=4,al=50,ltty=1) {
  #if (length(xvals)>100) {
  #  inds=seq(1,length(xvals),round(length(xvals)/60))+1
  #	xvals=xvals[inds]
  #	yvals=yvals[inds]
  #	dev=dev[inds]
  #}
  cols=c("blue","red","green", "black","orange","purple","khaki","grey","gold4","seagreen3","chocolate")
  tcols=genCols(al)
  lines(xvals,yvals,col=cols[colNum],lwd=lwid,lty=ltty)
  x=c(xvals,rev(xvals))
  y=c(yvals-dev,rev(yvals+dev))
  if (ltty==1) {
    polygon(x,y,col=tcols[colNum],border=NA)
  }
}



getBin = function(vsort, i, n=0, smooth=0,overlap) {
  vNames=colnames(vsort)
  binTotal = c()
  nbin=round((length(vsort[,myCol])-smooth)/(n-1)-0.5)
  if (smooth==0) {
    delta=1/(1-overlap)
    nbin=round((length(vsort[,myCol])-delta)/(n-1+delta))
    smooth=(nbin+1)*delta
    
  }
  if (n==0) {
    n=round((length(vsort[,myCol])-smooth)/(smooth*(1-overlap)))
    nbin=round((length(vsort[,myCol])-smooth)/(n-1)-0.5)
  }
  q1=nbin*(i-1)+1
  q2=q1+smooth
  #print(q1)
  #??print(q2)
  return(vsort[q1:q2,])
}


rHeatMap=function(vdat,ctrlDat,xax,yax,cond,nr=50,nc=50,xwid=0.25,xstep=0.075,xstart=-2,minNum=50,ov=0.9,maxVal=1,ret=1,minVal=0) {
  
}

plosHeatMap=function(vdat,ctrlDat,xax,yax,cond,nr=50,nc=50,xwid=0.25,xstep=0.075,xstart=-2,minNum=50,ov=0.9,maxVal=1,ret=1,minVal=0) {
  
  cNames=colnames(vdat)
  condCol=which(cNames==cond)
  xCol=which(cNames==xax)
  yCol=which(cNames==yax)
  
  ctrl=1:nc
  for(j in 1:nc) {
    start=xstart+xstep*(j-1)
    end=start+xwid
    mySpot= ctrlDat[which(ctrlDat[,xCol]>start),]
    mySpot=mySpot[which(mySpot[,xCol]<=end),]
    ctrl[j]=mean(mySpot[,yCol])
  }
  sum=0
  x=matrix(nrow=nr,ncol=nc,0)
  xDiff=x
  vsort=vdat[order(vdat[,condCol]),]
  numSpots=0
  for(i in 1:nr) {
    myBin=getBin(vsort,i,n=nr,overlap=ov)
    print(paste(i, sum, numSpots, sum/numSpots))
    for (j in 1:nc) {
      start=xstart+xstep*(j-1)
      end=start+xwid
      mySpot=myBin[which(myBin[,xCol]>start),]
      mySpot=mySpot[which(mySpot[,xCol]<=end),]
      myPct=mean(mySpot[,yCol])
      myLen=length(mySpot[,1])
      if (myLen > minNum) {
        myVar=sqrt(myPct*(1-myPct)/myLen)
        x[i,j]=myPct
        obs=myPct*myLen
        
        
        exp=myLen*ctrl[j]
        sum = sum + ((obs-exp)*(obs-exp)/exp)
        xDiff[i,j]=((obs-exp)*(obs-exp)/(exp*myLen))
        #sum = sum + abs(myLen*log(myPct/ctrl[j]))
        numSpots=numSpots+1
        
      }
    }
  }
  x[1,nc]=maxVal
  x[which(x==0)]=1000
  x[which(x<minVal)]=minVal
  x[which(x==1000)]=minVal-0.001
  print(sum)
  print(numSpots)
  if(ret==2) {
    rlist=list(1,2,sum,numSpots)
    rlist[[1]]=x
    rlist[[2]]=xDiff
    return(rlist)
  }
  else {
    return(x)
  }
}


processCutoffs = function(vsort, strongCut, weakCut, colP1, colP2) {
  len=length(vsort[,1])
  hasStrong=rep(0,len)
  hasWeak=hasStrong
  consWeak=hasStrong-1
  fNum = vsort$factorNum
  vNames=colnames(vsort)
  myCol=which(vNames==colP1)
  myCol2=which(vNames==colP2)
  pVals = vsort[,myCol]
  p2Vals=vsort[,myCol2]
  for(i in 1:len) {
    if (((fNum[i] != 5) && (fNum[i] > 0)) && !(is.na(pVals[i]))){
      if (pVals[i]<weakCut[fNum[i]]) {
        hasWeak[i]=1
      }
      if (pVals[i]<strongCut[fNum[i]]) {
        hasWeak[i]=1
        hasStrong[i]=1
      }
      if (hasWeak[i]==1) {
        if (p2Vals[i]<=weakCut[fNum[i]]) {
          consWeak[i]=1
        }
        else {
          consWeak[i]=0
        }
      }
    }
  }
  cutData=cbind(hasStrong,hasWeak, consWeak)
  return(cutData)
}

setupAll = function(strongCut=c(-11.2,-13.5,-11.2,-11,-11,-10),
                    weakCut=c(-9,-11.65,-9.4,-9,-11,-8)) {
  vplus=setup()
  
  vret=cbind(vplus,processCutoffs(vplus,strongCut,weakCut,"PWMP1","PWMP2"),processCutoffs(vplus,strongCut,weakCut,"DistP1","DistP2"))
  return(vret)
}


myCor=function(vdat, toSort, toCalc, n=30,type=1,over=0.9,off=0,ymin=-1,ymax=-1,cut=0) {
  vNames=colnames(vdat)
  mycol=which(vNames==toSort)
  col2=which(vNames==toCalc)
  
  adat=vdat[order(vdat[,mycol]),]
  mdat=calcBin(adat,toCalc,n=n,overlap=over,type=type,cut=cut)
  if (type==2) {
    m1=calcBin(adat,toCalc,n=n,overlap=over)
    m2=calcBin(adat,vNames[col2+1],n=n,overlap=over)
    mdat=(m1/(m1+m2))
  }
  if (type==3) {
    mdat=calcBin(adat,toCalc,n=n,type=2,cut=cut)
  }
  xvals=calcBin(adat,toSort,n=n,overlap=over)
  if(off==0) {
    plot(xvals,mdat,xlab=toSort,ylab=toCalc)
  }
  if(off==1) {
    if (ymax != (-1)) {
      plot(xvals,mdat,xlab="",ylab="",type="n",ylim=c(ymin,ymax))
    }
    else {
      plot(xvals,mdat,xlab="",ylab="",type="n")
    }
  }
  #print(cor.test(xvals,mdat))
  #plot(mdat)
  return(cor(xvals,mdat))
}


plosCor=function(vdat, toSort, toCalc, n=30,col=1,type=0,over=0.9,flip=0,pct=1,xov=0.9,plty=1) {
  vNames=colnames(vdat)
  mycol=which(vNames==toSort)
  col2=which(vNames==toCalc)
  
  adat=vdat[order(vdat[,mycol]),]
  mdat=calcBin(adat,toCalc,n=n,overlap=over,type=type)
  if (flip==1) {
    mdat=1-mdat
  }
  nsmth=getSmooth(adat,toCalc,n=n,overlap=over,type=type)
  if (pct==1) {
    adev=sqrt(mdat*(1-mdat)/nsmth)
  }
  if(pct==0) {
    adev=calcBin(adat,toCalc,n=n,overlap=over,type=3)
  }
  xvals=calcBin(adat,toSort,n=n,overlap=xov)
  plos(xvals,mdat,adev,col,lwid=lwid,ltty=plty)
  #print(cor.test(xvals,mdat))
  #plot(mdat)
  return(cor(xvals,mdat))
}

rd= function(dat,n=2) {
  return(round(dat, n))
}
src= function() {
  source("rahulFxns.R")
  source("RobFig1.txt")
  source("RobFig2.txt")
  #source("RobFig3.txt")
  source("RobFig4.txt")
  #source("RobFig5.txt")
  #source("RobFig6.txt")
  #source("gelShift.R.txt")
}

setupP = function() {
  p1=read.table("pocketOut/summits.1FDR.50",header=TRUE)
  p2=read.table("mappedMotifs/summits.1FDR.9.2.100",header=TRUE)
  p2b=read.table("mappedMotifs/summits.1FDR.11.5.100",header=TRUE)
  p3=read.table("mappedTTM/TAGTM01.summits.1FDR.9.8.500",header=TRUE)
  pca=prcomp(p1[,6:29],scale=TRUE,center=TRUE)
  pcs=-pca$x[,1:5]
  pcmel=-pca$x[,1]
  p=cbind(p1,p2,p2b,p3,pcmel,pcs)
}

setupQ = function() {
  q1=read.table("pocketOut/summits.25ORC.200",header=TRUE)
  q2=read.table("mappedMotifs/summits.25FDR.9.2.100",header=TRUE)
  q2b=read.table("mappedMotifs/summits.25FDR.11.5.100",header=TRUE)
  q2c=read.table("mappedMotifs/summits.25FDR.9.25.100",header=TRUE)
  q2d=read.table("mappedMotifs/summits.25FDR.11.55.100",header=TRUE)
  
  q3=read.table("mappedTTM/TAGTM01.summits.25FDR.9.8.500",header=TRUE)
  q4=read.table("mappedTTM/TAGTM27.summits.25FDR.9.500",header=TRUE)
  
  
  cn=colnames(q1)
  
  tfs=cn[6:29]
  tfs[24]="Z_Stage14"
  tfMeans=c()
  tfVars=c()
  bindingData=matrix(nrow=nrow(q1),ncol=24)
  melRatio=c()
  for(i in 1:24) {
    tmp=subset(q1,factorName==tfs[i])
    tfMeans[i]=mean(tmp[,i+5])
    tfVars[i]=sqrt(var(tmp[,i+5]))
    newm=(tmp[,i+5]-tfMeans[i])
    newsd=tfVars[i]
    melRatio=c(melRatio,newm/newsd);
    bindingData[,i]=(q[,i+5]-tfMeans[i])/tfVars[i]
  }
  pca=prcomp(q1[,6:29],scale=TRUE,center=TRUE)
  pcs=-pca$x[,1:5]
  pcmel=-pca$x[,1]
  
  ap=subset(q1,factorNum<6)
  pca=prcomp(ap[,6:11],scale=TRUE,center=TRUE)
  pcsAP=-pca$x[,1:5]
  pcAP=rep(0,nrow(q1))
  pcAP[1:nrow(ap)]=-pca$x[,1]
  pcAP[(nrow(ap)+1):nrow(q1)]=0
  
  q=cbind(q1,q2,q2b,q2c,q2d,q3,q4,pcmel,pcs,pcAP)
  p=cbind(q,melRatio)
  
}

setupV = function() {
  v1=read.table("summitData/summits.chipSeq",header=TRUE)
  v2=read.table("mappedMotifs/summits.chipSeq.9.2.100",header=TRUE)
  v3=read.table("mappedTTM/TAGTM01.summits.chipSeq.9.8.250",header=TRUE)
  v=cbind(v1,v2,v3)
}

gF1 = function() {
  pdf("1.pdf",height=6,width=8)
  genFig1()
  dev.off()
}

gF2 = function() {
  pdf("2.pdf",height=6,width=8)
  genFig2()
  dev.off()
}

gF3 = function() {
  pdf("33.pdf",height=9,width=5)
  genFig3.3()
  dev.off()
}
gF4 = function() {
  pdf("4.pdf",height=9,width=11)
  genFig4.1()
  dev.off()
}	

gF6 = function() {
  pdf("6.pdf",height=6,width=12)
  genFig6()
  dev.off()
}	


combineIntersection = function (prefix, peakList, toCombine, nDelim=3, returnBinary=0, returnColumn=11, returnReadCount=0,addSuffix=".bed.intersect.bed", promGeneName=0) {
  intersect = read.table(paste(prefix,"/",toCombine,addSuffix,sep=""))
  returnVals=rep(0,nrow(peakList))
  
  uniqueInt = intersect[!duplicated(intersect$V4),]
  spec1IDs = (uniqueInt$V6)
  spec2IDs = as.character(uniqueInt$V10)
  
  if (returnBinary==1) {
    returnVals[spec1IDs]=1
    return(returnVals)
  }
  
  returnVals[spec1IDs]=uniqueInt[,returnColumn]
  return(returnVals)
}

meanmin=function(x,min=1,val=0) {
  if (length(x)>min) {
    return(mean(x))
  } 
  return(val)
}


medmin=function(x,min=1,val=0) {
  if (length(x)>min) {
    return(median(x))
  } 
  return(val)
}

combineIntersection = function (prefix, peakList, toCombine, nDelim=3, returnBinary=0, returnColumn=10, returnReadCount=0,addSuffix=".bed.intersect.bed", promGeneName=0, uniqueCol=4,id1Col=6) {
  intersect = read.table(paste(prefix,"/",toCombine,addSuffix,sep=""))
  returnVals=rep(0,nrow(peakList))
  
  uniqueInt = intersect[!duplicated(intersect[,uniqueCol]),]
  spec1IDs = (uniqueInt[,id1Col])
  spec2IDs = as.character(uniqueInt$V10)
  
  if (returnBinary==1) {
    returnVals[spec1IDs]=1
    return(returnVals)
  }
  
  spec2IDList = unlist(strsplit(spec2IDs,"_"))
  spec2PeakNums=as.numeric(spec2IDList[seq(nDelim,length(spec2IDList),nDelim)])
  
  
  if (promGeneName == 1) {
    geneNames=substr(as.character(spec2IDs),1,10)
    returnVals[spec1PeakNums]=geneNames
    return(returnVals)
  }
  if (returnReadCount > 0) {
    summitFile = read.table(paste("../",toCombine,"_summits.bed",sep=""))
    returnVals[spec1PeakNums]=summitFile$V5[spec2PeakNums]
    return(returnVals)
  }
  
  
  returnVals[spec1IDs]=uniqueInt[,returnColumn]
  return(returnVals)
}

init2 = function() {
  library(ggplot2)
  
  
  opt <-  opts(legend.title = theme_blank(), # switch off the legend title
               legend.text = theme_text(size=12,face="bold"),        
               legend.key.size = unit(2.5, "lines"),
               legend.key = theme_blank(),
               axis.title.x = theme_text(size = 14, vjust = -0.5),
               axis.title.y = theme_text(size = 14, angle = 90),
               axis.text.x = theme_text(size = 12),
               axis.text.y = theme_text(size = 12),
               plot.title = theme_text(size = 18,vjust=2.5,face="bold"),
               plot.margin = unit(c(2,2,0.75,0.75), "lines"))
  
  lwid=2
}

writ.table=function(a, b) {
  write.table(a,b,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
}


sp=function(data,x,y,z,n=50,min=10,minval=0,maxval=0,qx=0.975,qy=0.975,func="mean") {
  xs=seq(quantile(data[,x],1-qx),quantile(data[,x],qx),length.out=n)
  ys=seq(quantile(data[,y],1-qy),quantile(data[,y],qy),length.out=n)
  t=data
  t$x=as.numeric(factor(cut(t[,x],xs)))
  t$y=as.numeric(factor(cut(t[,y],ys)))
  t=subset(t,!is.na(x)&!is.na(y))
  t$code=paste(t$x,t$y,sep="_")
  res=matrix(0,nrow=n,ncol=n)
  mylev=levels(factor(t$code))
  myres=tapply(t[,z],t$code,meanmin,min)
  if (func=="median") {
    myres=tapply(t[,z],t$code,medmin,min)
  }
  for(i in 1:length(mylev)) {
    inds = as.numeric(unlist(strsplit(mylev[i],"_")))
    res[inds[1],inds[2]]=myres[i]
  }
  mapCols=myPalette(low="white",high="black")
  mapCols[1]="white"
  res[which(res<minval)]=minval
  
  res[which(res==0)]=min(res)-0.01
  if (maxval ==0) {
    maxval=max(res)
  }
  res[which(res>maxval)]=maxval
  n=n-1
  yqs=quantile(t[,y],seq(0,1,1/n))
  xqs=quantile(t[,x],seq(0,1,1/n))
  labRow=as.character(round(xs,2))
  labCol=as.character(round(ys,2))
  #heatmap.2((res),Rowv=NA,Colv=NA,scale="none",col=mapCols,ylab=x,xlab=y,main=z,,,density.info="none",trace="none",keysize=1)
  #,scales=list(x=list(at=1:n, labels=labRow)), x=list(at=1:n, labels=labCol)
  #print(levelplot(res[1:(n-1),1:(n-1)],col.regions=mapCols,xlab=list(label=y,cex=1.5),ylab=list(label=x,cex=1.5),main=list(label=z,cex=2)),scales=list(x=list(at=1:(n-1), labels=labRow[1:(n-1)]), y=list(at=1:(n-1), labels=labCol[1:(n-1)])))
  print(levelplot(res[1:n,1:n],col.regions=mapCols,,xlab=list(label=x,cex=1.3),ylab=list(label=y,cex=1.3),main=list(label=z,cex=2),scales=list(x=list(at=1:n, labels=labRow),y=list(at=1:n, labels=labCol))))
  #levelplot(res)
  return(list(t,res))
}

calcTP = function(cutoff,data,score,real,nTP) {
  return(length(which((data[,score]>cutoff)&(data[,real]>0)))/nTP)
}

calcFP = function(cutoff,data,score,real,nFP) {
  return(length(which((data[,score]>cutoff)&(data[,real]==0)))/nFP)
}

auc=function(data,score,real,n=20) {
  totalPos=length(which(data[,real]==1))
  totalNeg=length(which(data[,real]==0))
  scores=data[,score]
  data$myScore=(scores+min(scores))/(max(scores)+min(scores))
  tp=unlist(lapply(seq(-0.0001,0.9999,1/n),calcTP,data,"myScore",real,totalPos))
  fp=unlist(lapply(seq(-0.0001,0.9999,1/n),calcFP,data,"myScore",real,totalNeg))
  plot(c(fp,1),c(tp,1),xlim=c(0,1),ylim=c(0,1))
  x1=c(1,fp)
  x2=c(1,tp)
  print(sum(diff(rev(x2))*diff(rev(x1)))/2+sum(diff(rev(x1))*rev(x2[-1])))
  return(list(c(1,fp),c(1,tp)))
}

bplot=function(s,CR,q=0.98,save=0) {
  s$curdata=s[,CR]
  org = s[1,"Org"]
  x=qplot(curdata,data=s,geom="density",col=stat,xlim=c(0,quantile(s$curdata,q,na.rm=TRUE)),main=paste(CR,org),xlab="",size=I(1))+opt
  print(x)
  return(x)
}

bplotSave=function(CR) {
  a=bplot(s,CR)
  ggsave(file=paste("pic/",CR,".jpg",sep=""))
  print(paste("pic/",CR,".jpg",sep=""))
}

bcplot=function(s,CR,q=0.98,...) {
  s$curdata=s[,CR]
  qplot(curdata,data=subset(s,cg!="NA"),geom="density",col=stat,xlim=c(0.01,quantile(s$curdata,q),xlab="",size=I(1)),main=CR,...)+facet_grid(. ~ cg )+ opt
}

gea = function(data,file) {
  data[,1]=paste("chr",data[,1],sep="")
  writ.table(data[,1:3],file) 
}


bbplot = function(data,name,start=0.2,doSort=TRUE,...) {
  if (doSort) {
    x=boxplot(split(data[,name],data$stat)[order(tapply(data[,name],data$stat,median,na.rm=TRUE))],...)
    text(1:length(summary(data$stat)),start,summary(data$stat)[order(tapply(data[,name],data$stat,median,na.rm=TRUE))])
  }
  else {
    x=boxplot(split(data[,name],data$stat),...)
    text(1:length(summary(data$stat)),start,summary(data$stat))
  }
}

pyCols=myPalette(low = "magenta",high = "yellow",mid = "black")

rp=function() {par(mfrow=c(1,1))}

calcResidLog=function(x1,y1,mcut=30,toAdd=1) {
  touse=which((x1>mcut)&(y1>mcut))
  x=log(x1+toAdd,2)
  y=log(y1+toAdd,2)
  a=x[touse]
  b=y[touse]
  myLM=lm(b~a)
  myResid=y-predict(myLM,data.frame(a=x))
  return(myResid)
}

logVarDivMean=function(x) return(log(var(exp(x)-1)/mean(exp(x)-1)))

calcResid=function(x1,y1,mcut=30,toAdd=1) {
  touse=which((x1>mcut)&(y1>mcut))
  x=x1
  y=y1
  a=x[touse]
  b=y[touse]
  myLM=lm(b~a)
  myResid=y-predict(myLM,data.frame(a=x))
  return(myResid)
}

gtCut=function(x, cutoff=1) {
  return(length(which(x>cutoff)))
}

findNGene=function(data,is.expr=1) {
  toRet=unlist(lapply(1:ncol(data),function(x)gtCut(data[,x],cutoff=is.expr)))  
  names(toRet)=colnames(data)
  return(toRet)
}

getCoefs=function(data,nbin=20,mycut=1) {
  my_stats=data.frame(data[,1])
  code_humpAvg=apply(data,1,humpMean,min=mycut)
  code_humpAvg[code_humpAvg>9]=9
  code_humpAvg[is.na(code_humpAvg)]=0
  my_stats$code_humpAvg=code_humpAvg
  data[data>mycut]=1
  data[data<mycut]=0
  data$bin=cut(code_humpAvg,nbin)
  data$avg=code_humpAvg
  rownames(my_stats)=rownames(data)
  my_coefs=data.frame(t(sapply(colnames(data[1:(ncol(data)-2)]),getAB,data=data,data2=my_stats,status="code",code2="humpAvg",hasBin=TRUE,doPlot=FALSE)))
  colnames(my_coefs)=c("a","b")
  return(my_coefs)
}

makeScorePlot2=function(allscores,getStatFxn=getStat2,mytitle="Title") {
  alls=data.frame(allscores)
  alls$stat=unlist(lapply(names(allscores),getStatFxn))
  p=ggplot(alls,aes(factor(stat),allscores))
  p2=p + geom_violin(scale="width",adjust=0.75,aes(fill=factor(stat))) + ylab("Expression level (log TPM)")
  p3=p2+theme(legend.position="top")+guides(fill=guide_legend(title=NULL))+geom_jitter(height=0)+blackbg+xlab("Cell Type")
  p4=p3+ theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=90, vjust=0.5, size=12))
  print(p4+theme(axis.title.y = element_text(face="bold", colour="#990000", size=16), axis.text.y  = element_text(angle=90, vjust=0.5, size=12))+ggtitle(mytitle)+theme(plot.title = element_text(size=20, face="bold")))
}

getNewScore=function(mygenes,data=mcell,wt_matrix,myfunc=weighted.mean,code="rsem",scramble=FALSE) {
  my_data=subc(data,code)
  if(scramble) {
    my_data=my_data[,sample(ncol(my_data))]
  }
  myscores=unlist(lapply(colnames(my_data),function(x)myfunc(my_data[mygenes,x],wt_matrix[mygenes,x])))
  names(myscores)=colnames(my_data)
  return(myscores)
}

corCellWeightFast=function(cell1,cell2,wt_matrix,data=subc(cell,"lps_t2"),spear=FALSE) {
  cell1Drop=wt_matrix[rownames(data),cell1]
  cell2Drop=wt_matrix[rownames(data),cell2]
  my_weights=cell1Drop*cell2Drop
  cell1Data=an(data[,cell1])
  cell2Data=an(data[,cell2])
  # print(my_weights)
  my_weights[is.na(my_weights)]=0
  ret=corr(d=as.matrix(cbind(cell1Data,cell2Data)),w=my_weights)
  if (spear==TRUE) {
    ret=corr(d=as.matrix(cbind(rank(cell1Data),rank(cell2Data))),w=my_weights)
  }
  return(ret)
}

covCellWeightFast=function(cell1,cell2,wt_matrix,data=subc(cell,"lps_t2"),spear=FALSE) {
  cell1Drop=wt_matrix[rownames(data),cell1]
  cell2Drop=wt_matrix[rownames(data),cell2]
  my_weights=cell1Drop*cell2Drop
  cell1Data=an(data[,cell1])
  cell2Data=an(data[,cell2])
  
  # print(my_weights)
  ret=wtCov(cell1Data,cell2Data,my_weights)
  if (spear==TRUE) {
    ret=wtCov(rank(cell1Data),rank(cell2Data),my_weights)
  }
  return(ret)
}

scaleSCMatrix2=function(data,wt_matrix,code="rsem") {
  wtX=unlist(lapply(rownames(data),function(x)(sum(data[x,]*wt_matrix[x,])/sum(wt_matrix[x,]))))
  sdX=unlist(lapply(rownames(data),function(x)(wtCov(data[x,],data[x,],wt_matrix[x,]))))
  sData=(data-wtX)/sqrt(sdX)
  return(sData)
}

wtCov=function(x,y,w) {
  w=w/sum(w)
  wtX=sum(x*w)
  wtY=sum(y*w)
  wt_cov=sum(w*(x-wtX)*(y-wtY))
  return(wt_cov)
}
expAlpha=function(mu,coefs) {
  logA=coefs$a
  logB=coefs$b
  return(exp(logA+logB*mu)/(1+(exp(logA+logB*mu))))
}

setWt1=function(x,wts,min=1) {
  wts[x>min]=1
  return(wts)
}

gtCut=function(x, cutoff=1) {
  return(length(which(x>cutoff)))
}

findNGene=function(data,mycut=1) {
  toRet=unlist(lapply(1:ncol(data),function(x)gtCut(data[,x],cutoff=mycut)))  
  names(toRet)=colnames(data)
  return(toRet)
}

setWtMatrix1=function(data,wt_matrix,mycut=1) {
  wt1_matrix=sapply(1:ncol(data),function(x)setWt1(data[,x],wt_matrix[,x],min=mycut))
  colnames(wt1_matrix)=colnames(data)
  return(wt1_matrix)
}


getAB=function(cn="lps_t1_S1_rsem",code="lps_t1",data=cell,data2=cs,code2="avg",status="",ncut=25,hasBin=FALSE,doPlot=FALSE,myfunc=plot,func2=lines,...) {
  #Check this
  #if (status == "") {
  # status=getStatus(cn)
  #}
  #Check this
  if (!(hasBin)) {
    data[data>1]=1
    data[data<1]=0
    data$avg=data2[rownames(data),paste(status,"_",code2,sep="")]
    data[data>9]=9
    data$bin=cut(data$avg,20)
  }
  data$val=data[,cn]
  x1=(tapply(data[,cn],data$bin,mean))
  x2=(tapply(data[,"avg"],data$bin,mean))
  #glm.out = glm(val~avg,data=data,family=binomial)
  glm.out = glm(x1~x2,family=binomial)
  
  if (doPlot==TRUE) {
    pred=(predict(glm.out,data.frame(avg=as.numeric(x2)),type="response"))
    myfunc(x2,x1,pch=16,xlab="Average Expression",ylab="P(detection) or 1-FNR",main=cn,ylim=c(0,1))
    func2(x2,pred,...)
  } 
  return(glm.out$coefficients)
}




sensitivityCurve=function(cellName,scData,bulkData,mycut=1,mycex=1,mynew=TRUE,...) {
  cutLocs=cut2(bulkData,g=100,onlycuts=TRUE)
  bulkBin=cut2(bulkData,g=100)
  binaryData=scData[,cellName]
  binaryData[binaryData >= mycut]=1
  binaryData[binaryData<mycut]=0
  yBin=tapply(binaryData,bulkBin,mean)
  xBin=tapply(bulkData,bulkBin,mean)
  #glm.out = glm(val~avg,data=data,family=binomial)
  options(warn=-1) #otherwise glm throws an unnecessary error
  glm.out = glm(binaryData~bulkData,family=binomial)
  options(warn=0)  
  x_vals=seq(0,10,0.1)
  y_vals=predict(glm.out,data.frame(bulkData=x_vals),type="response")
  if (mynew) {
    plot(xBin,yBin,pch=16,xlab="Average expression",ylab="Probability of detection",...)
  }
  lines(x_vals,y_vals,lwd=2,...)
}

getAdjMatrix=function(X,nn=10,edge.weights=FALSE,do.jaccard=TRUE,do.sparse=TRUE,full.eval=FALSE,...) {
  require(RANN)
  if (dim(X)[1] < 500){
    searchtype.use="standard"
    treetype.use = "kd"
  } else {
    searchtype.use="priority"
    treetype.use = "bd"
  }
  print(dim(X))
  nearest=nn2(X,X,k=nn+1, treetype = treetype.use, searchtype=searchtype.use)
  print("Found nearest neighbors")
  nearest$nn.idx = nearest$nn.idx[,-1]
  nearest$nn.dists = nearest$nn.dists[,-1] #Convert to a similarity score
  
  if (edge.weights){
    nearest$nn.dists = nearest$nn.dists / median(nearest$nn.dists)
    nearest$nn.sim = 1 / nearest$nn.dists  #Convert to a similarity score
  } else { 
    nearest$nn.sim = 1*(nearest$nn.dists >= 0 )
  }
  
  edges = melt(t(nearest$nn.idx)); colnames(edges) = c("B", "A", "C"); edges = edges[,c("A","B","C")]
  edges$B = edges$C; edges$C=1
  #Remove repetitions
  edges = unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))
  
  if (!full.eval){
    if (do.jaccard){
      
      NN = nearest$nn.idx
      jaccard_dist = apply(edges, 1, function(x) length(intersect(NN[x[1], ],NN[x[2], ]))/length(union(NN[x[1], ], NN[x[2], ])) )
      
      edges$C = jaccard_dist
      #edges$B = rownames(X)[edges$B]
      #edges$A = rownames(X)[edges$A]
      edges = subset(edges, C != 0)
      edges$C = edges$C/max(edges$C)
    } else {
      
      #     edges = data.frame()
      #     edge_temp = data.frame()
      #     for (i in 1:dim(X)[1]){
      #       
      #       temp = data.frame(A=rownames(X)[i], B=rownames(X)[nearest$nn.idx[i,]],C=nearest$nn.sim[i,])
      #       edge_temp = rbind(edge_temp,temp)
      #       
      #       if (i %% 100 == 0) {
      #         if (i %% 1000 == 0) print(i)
      #         edges = rbind(edges,edge_temp)
      #         edge_temp = data.frame()
      #       }
      #       
      #     }
      #     edges = rbind(edges,edge_temp)
    }
    
    Adj = Matrix(0, nrow=nrow(X), ncol=nrow(X), sparse=do.sparse)
    Adj[cbind(edges$A,edges$B)] = edges$C
    Adj[cbind(edges$B,edges$A)] = edges$C
    rownames(Adj) = rownames(X); colnames(Adj) = rownames(X)
    return(Adj)
    
  }
  
  if (full.eval){
    if (do.jaccard){
      Adj = Matrix(0, nrow=nrow(X), ncol=nrow(X), sparse=do.sparse)
      rownames(Adj) = rownames(X)
      colnames(Adj) = rownames(X)
      NN = nearest$nn.idx
      
      for (i in 1:nrow(Adj)){
        Adj[i,] = sapply(c(1:nrow(Adj)), function(x) length(intersect(NN[i, ],NN[x, ]))/length(union(NN[i, ], NN[x, ])))
        Adj[i,i] = 0
      }
      
      Adj = Adj + t(Adj)
      Adj = Adj/max(Adj)
      
      
      
      
    } else {
      
      Adj = Matrix(0, nrow=nrow(X), ncol=nrow(X), sparse=do.sparse)
      rownames(Adj) = rownames(X)
      colnames(Adj) = rownames(X)
      Adj[cbind(edges$A,edges$B)] = edges$C
      Adj[cbind(edges$B,edges$A)] = edges$C
      
      return(Adj)
    }
    
    
    
 }
}  
  
plotAdjMatrix=function(object,clusters.use=NULL, title=NULL, nn.use=100, pcs.use=1:50) {
    
    cells.use = c()
    
    for (i in clusters.use){
      cells = which.cells(object, i)
      cells.use = c(cells.use, cells[sample(length(cells))])
    }
    Adj = getAdjMatrix(object@pca.rot[cells.use,pcs.use],nn=nn.use,edge.weights=FALSE, do.jaccard=FALSE)
    
    #   k=1
    #   for (i in clusters.use){
    #     l = length(which.cells(object, i))
    #     Adj[k,c(k:(k+l-1))] = 5
    #     Adj[k+l-1,c(k:(k+l-1))] = 5
    #     Adj[c(k:(k+l-1)),k] = 5
    #     Adj[c(k:(k+l-1)),k+l-1] = 5
    #     
    #     k=k+l
    #   }
    
    #Adj = Adj[length(cells.use):1,]
    
    
    
    image(Adj, useRaster=TRUE, xlab="Cells (ordered by cluster)", ylab="Cells (ordered by cluster)")
    
    
  }
  
  graphAdjMatrix=function(Adj,cluster.ident=NULL,clusters.use=NULL,order.by.cluster=TRUE, cluster.annot=TRUE, annot.ident=NULL, title=NULL) {
    
    cells.use = rownames(Adj)
    if (!is.null(cluster.ident)){
      if (order.by.cluster){
        cluster.ident=sample(cluster.ident, length(cluster.ident))
        cluster.ident = cluster.ident[order(cluster.ident)]
      }
      cells.use = names(cluster.ident)
      if (!is.null(clusters.use)){
        cluster.ident = cluster.ident[cluster.ident %in% clusters.use]
        cells.use = names(cluster.ident)
      }
      cluster.sizes = as.numeric(table(cluster.ident))
    }
    
    if (!is.null(annot.ident)){
      annot.ident = annot.ident[cells.use]
    }
    
    X = Adj[cells.use,cells.use]
    
    if (order.by.cluster & !is.null(cluster.ident)){
      start_ind = 1; end_ind = cumsum(cluster.sizes)
      colsep.use =cumsum(cluster.sizes)
      rowsep.use = cumsum(cluster.sizes)
      title = set.ifnull(title,paste0("N = ", sum(cluster.sizes), ", max_size = ", max(cluster.sizes), ", min_size = ", min(cluster.sizes)))
    } else {
      colsep.use = NULL
      rowsep.use = NULL
      title= set.ifnull(title, paste0("N = ", length(cells.use)))
    }
    
    if (!is.null(cluster.ident) & cluster.annot){
      annot.ident = set.ifnull(annot.ident,cluster.ident)
      annot.cols.list = sample(rainbow(length(unique(annot.ident))))
      annot.ident.new=annot.ident
      l=1
      for (i in unique(annot.ident)){
        annot.ident.new[annot.ident==i] = l
        l=l+1
      }
      annot.cols = annot.cols.list[as.numeric(annot.ident.new)]
      heatmap.2(as.matrix(X), Rowv=FALSE,Colv=FALSE,dendrogram="none", trace="none",labRow=NA,labCol=NA,col=grey(seq(0,1,length=10)), 
                colsep=colsep.use, rowsep=colsep.use, RowSideColors = annot.cols, ColSideColors = annot.cols, main=title)
    } else {
      heatmap.2(as.matrix(X), Rowv=FALSE,Colv=FALSE,dendrogram="none", trace="none",labRow=NA,labCol=NA,col=grey(seq(0,1,length=10)), 
                colsep=colsep.use, rowsep=colsep.use, main=title)
    }
    
    
  }
  
  getInfomapNetwork=function(network,X,nn=10,edge.weights=FALSE,do.jaccard=TRUE,do.sparse=TRUE,is.symm=TRUE, do.noise=FALSE,...) {
    require(RANN)
    nearest=nn2(X,X,k=nn+1, treetype = "bd", searchtype="priority")
    print("Found nearest neighbors")
    nearest$nn.idx = nearest$nn.idx[,-1]
    nearest$nn.dists = nearest$nn.dists[,-1] #Convert to a similarity score
    
    if (edge.weights){
      nearest$nn.dists = nearest$nn.dists / median(nearest$nn.dists)
      nearest$nn.sim = 1 / nearest$nn.dists  #Convert to a similarity score
    } else { 
      nearest$nn.sim = 1*(nearest$nn.dists >= 0 )
    }
    
    edges = melt(t(nearest$nn.idx)); colnames(edges) = c("B", "A", "C"); edges = edges[,c("A","B","C")]
    edges$B = edges$C; edges$C=1
    #Remove repetitions
    edges = unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))
    
    #Add 10% spurious edges
    
    if (do.noise){
      print("Add Noise")
      l=0
      n_edge = dim(edges)[1]
      to_remove = sample(1:n_edge, round(0.025*n_edge))
      nodes1 = sample(1:dim(X)[1], round(0.025*n_edge)*100, replace = TRUE);
      nodes2 = sample(1:dim(X)[1], round(0.025*n_edge)*100, replace = TRUE);
      new_edges = data.frame(A=nodes1, B=nodes2, C=1)
      new_edges = unique(transform(new_edges, A = pmin(A,B), B=pmax(A,B)))
      distinct_edges = anti_join(new_edges, edges)
      edges = rbind(edges,distinct_edges[1:round(0.025*n_edge),])
      edges <- edges[-to_remove,]
    }
    rownames(edges) = 1:nrow(edges)
    if (do.jaccard){
      print("do Jaccard Correction")
      
      NN = nearest$nn.idx
      jaccard_dist = apply(edges, 1, function(x) length(intersect(NN[x[1], ],NN[x[2], ]))/length(union(NN[x[1], ], NN[x[2], ])) )
      
      edges$C = jaccard_dist
      #edges$B = rownames(X)[edges$B]
      #edges$A = rownames(X)[edges$A]
      edges = subset(edges, C != 0)
      edges$C = edges$C/max(edges$C)
    }
    
    if (do.noise){
      print("Multiplicative Noise")
      mult.noise = exp(runif(length(edges$C), min = -0.2, max = 0.2))
      edges$C = edges$C * mult.noise 
    }
    
    
    
    
    for (i in 1:nrow(edges)){
      network$addLink(edges$A[i] - 1, edges$B[i] - 1, edges$C[i])
      if (is.symm){
        network$addLink(edges$B[i] - 1, edges$A[i] - 1, edges$C[i])
      }
    }
    print("Done")
    return(network)
    
  }
  
  
  get_nn_graph=function(X=NULL, k = 100, write.file = NULL,downsample.frac=1, do.jaccard=TRUE, do.edge.weights=FALSE,...) {
    require(RANN)
    if (is.null(X)){
      stop("Error : need expression matrix to compute NN graph")
    }
    
    if (downsample.frac != 1){
      cells.use = sample.int(dim(X)[1],round( downsample.frac*dim(X)[1] ))
      data.use = X[cells.use,]
    } else {
      data.use = X
    }
    nearest=nn2(data.use,data.use,k+1)
    nearest$nn.idx = nearest$nn.idx[,-1]
    nearest$nn.dists = nearest$nn.dists[,-1]
    
    if (!do.jaccard){
      #nearest$nn.dists = nearest$nn.dists/quantile(nearest$nn.dists,0.95) #normalize distances
      #nearest$nn.sim = 1- nearest$nn.dists #Convert to a similarity score
      #nearest$nn.sim[nearest$nn.sim < 0] = 0
      #nearest$nn.sim = nearest$nn.sim / max(nearest$nn.sim)
      if (!do.edge.weights){
        edges = melt(t(nearest$nn.idx)); colnames(edges) = c("B", "A", "C"); edges = edges[,c("A","B","C")]
        edges$B = edges$C; edges$C=1
        #Remove repetitions
        edges = unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))
        edges$B = rownames(data.use)[edges$B]
        edges$A = rownames(data.use)[edges$A]
        edges = subset(edges, C!=0)
      } else {
        nearest$nn.dists = nearest$nn.dists/quantile(nearest$nn.dists,0.95) #normalize distances
        nearest$nn.sim = 1- nearest$nn.dists #Convert to a similarity score
        nearest$nn.sim[nearest$nn.sim < 0] = 0
        nearest$nn.sim = nearest$nn.sim / max(nearest$nn.sim)
        
        edges = melt(t(nearest$nn.sim)); colnames(edges) = c("B", "A", "C"); edges = edges[,c("A","B","C")]
        edges$B = melt(t(nearest$nn.idx))[,3];
        #Remove repetitions
        edges = unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))
        edges$B = rownames(data.use)[edges$B]
        edges$A = rownames(data.use)[edges$A]
        edges = subset(edges, C !=0)
      }
    }
    
    if (do.jaccard){
      edges = melt(t(nearest$nn.idx)); colnames(edges) = c("B", "A", "C"); edges = edges[,c("A","B","C")]
      edges$B = edges$C; edges$C=1
      #Remove repetitions
      edges = unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))
      
      NN = nearest$nn.idx
      jaccard_dist = apply(edges, 1, function(x) length(intersect(NN[x[1], ],NN[x[2], ]))/length(union(NN[x[1], ], NN[x[2], ])) )
      
      edges$C = jaccard_dist
      edges$B = rownames(data.use)[edges$B]
      edges$A = rownames(data.use)[edges$A]
      edges = subset(edges, C != 0)
    }
    
    if (is.null(write.file)){
      return(edges)
    } else { 
      write.table(edges, file=write.file, quote=F, sep="\t", row.names=F, col.names=F)
    }
    
  }
  
  run_MCL=function(nn_graph_file=NULL, out_file = "MCL_clusters.txt", inflation.use=1.5,...) {
    if (is.null(nn_graph_file)){
      stop("Error : need nearest neighbor graph to proceed")
    }
    CMD = paste0("mcl ", nn_graph_file, " --abc -I ", inflation.use, " -te 16 -o ", out_file)
    system(CMD)
    return(1)
  } 
  
  read_MCL_clusters=function(mcl_out_file=NULL, cells.return=NULL,min.points = 1,...) {
    if (is.null(mcl_out_file)){
      stop("Error : need MCL cluster assignments to proceed")
    }
    
    mcl_clusters = read.table(file=mcl_out_file,sep="\t",fill=TRUE)
    mcl_clusters <- data.frame(lapply(mcl_clusters, as.character), stringsAsFactors=FALSE)
    clust.ident = c(); cell.names = c()
    for (i in 1:dim(mcl_clusters)[1]){
      a = mcl_clusters[i,]
      a = a[a!=""]
      if (length(a) >= min.points){
        clust.ident = c(clust.ident, rep(i, length(a)))
        cell.names = c(cell.names, a)
      } else {break}
    }
    
    names(clust.ident) = cell.names
    if (!is.null(cells.return)){
      clust.ident[setdiff(cells.return, names(clust.ident))] = 0
      return(clust.ident[cells.return])
    } else {return(clust.ident)}
  }
  
  ############### FUNCTIONS FOR SATURATION ANALYSIS #################################
  
  nGene_downsample = function(reads.matrix, downsample.level=0.5, nboot=10){
    
    gene.names = rownames(reads.matrix)
    num.genes = c()
    for (cell in colnames(reads.matrix)){
      all.reads = rep(gene.names, reads.matrix[,cell])
      num.vec = c()
      for (i in 1:nboot){
        num.vec = c(num.vec, length(unique(sample(all.reads, round(downsample.level*length(all.reads)), replace=FALSE))))
      }
      
      num.genes = c(num.genes, median(num.vec))       
    }
    return(num.genes)
    
  }
  
  nTrans_downsample = function(reads.matrix, counts.matrix, downsample.level=0.5, nboot=10){
    
    
    gene.names = rownames(counts.matrix)
    num.trans = c()
    l=0;
    for (cell in colnames(counts.matrix)){
      
      num.reads.cell = as.numeric(reads.matrix[,cell]); num.trans.cell = as.numeric(counts.matrix[,cell])
      num.reads.cell = num.reads.cell[num.reads.cell > 0] ; num.trans.cell = num.trans.cell[num.trans.cell > 0]
      num.reads.per.trans = round(rep(num.reads.cell, num.trans.cell) / rep(num.trans.cell, num.trans.cell))
      
      all.trans = paste0(rep(gene.names, counts.matrix[,cell]), '_', c(1:(sum(counts.matrix[,cell]))))
      all.reads = rep(all.trans, num.reads.per.trans)
      num.vec = c()
      for (i in 1:nboot){
        num.vec = c(num.vec, length(unique(sample(all.reads, round(downsample.level*length(all.reads)), replace=FALSE))))
      }
      
      num.trans = c(num.trans, median(num.vec))
    }
    return(num.trans)
    
  }
  
  
  max_modularity = function(X=NULL,NN_graph=NULL, clust.ident.df = NULL, inflation.vals=NULL, do.plot=TRUE){
    NN_graph = set.ifnull(NN_graph, get_nn_graph(X=X,k=100))
    
    A = clust.ident.df
    
    #Adjacency matrix
    Adj = matrix(0, nrow=nrow(A), ncol=nrow(A))
    rownames(Adj) = rownames(A)
    colnames(Adj) = rownames(A)
    Adj[cbind(as.character(NN_graph$A),as.character(NN_graph$B))] = 1
    m = sum(Adj) #Number of edges
    modularity_vec = c()
    
    for (k in 1:dim(A)[2]){
      print(k)
      ident = A[,k]
      clust.names = sort(unique(ident))
      modularity = 0
      for (i in clust.names){
        
        cells.in.clust = names(ident)[ident == i]
        Adj_clust = Adj[cells.in.clust, cells.in.clust]
        
        lc = sum(Adj_clust)
        dc = sum(Adj[cells.in.clust,]) + sum(Adj[,cells.in.clust])
        modularity = modularity + lc/m - (dc/m)^2
      }
      
      modularity_vec = c(modularity_vec, modularity)
      
    }
    
    if (do.plot){
      plot(inflation_vals,modularity_vec, xlab="Inflation",ylab="Modularity", pch=20)
    } else{
      return(min(inflation_vals[m_vec > 0.95*max(m_vec)])) #Pick the minimum value of inflation that achives the best modularity
    }
    
  }
  
  
  ### Graph clustering
  cluster_jl <- function(M, s=1, self=1, debug=0, ddebug=0, verbose=TRUE){
    #INPUTS : 
    #M : adjacency matrix (the matrix is symmetrized) - should be sparse if big. Note that M need not necessarily be the original
    #adjacency matrix. 
    #s : 1 = recursive computation, 0 = one level computation
    #self : 1 = use self weights
    #       0 = Do not use self weights
    #debug : 1 = outputs some debug messages
    #verbose : 1 = outputs some messages
    
    # Output :
    # COMTY, structure with the following information
    # for each level i :
    #   COMTY$COM[[i]] : vector of community IDs (sorted by community sizes)
    #   COMTY$SIZE[[i]] : vector of community sizes
    #   COMTY$MOD[i] : modularity of clustering
    #   COMTY$Niter[i] : Number of iteration before convergence
    
    
    
    S = dim(M);
    N = S[1]; #Number of nodes
    
    to.return = list();
    to.return[[2]] = 0;
    
    #Symmetrize matrix
    M = M + t(M);
    
    if (!self){
      diag(M) = 0
    }
    M2=M; 
    diag(M2) = 0;
    
    m = sum(M); #total sum of weights or # of edges
    Niter = 1;
    
    if (m == 0 | N == 1){
      stop("No more decomposition possible")
      to.return[[1]] = 0
      to.return[[2]] = 1
      return(to.return);
    }
    
    #Main loop
    K = colSums(M); #Sum of weight incident to node i
    sumTot = K; #SigmaTot in the paper
    sumIn = diag(M); #Sum of weight inside community i
    COM = c(1:S[1]); #Each node is its own community 
    print("Finding Nearest Neighbors")
    Neighbor = lapply(1:N, function(x) {if (x %% 500 == 0) {print(x)};which(M2[x,] > 0)})
    
    sCost=10
    gain=1;
    
    while (gain == 1){
      Cost = rep(0, N);
      gain = 0;
      #Now loop through each node
      for (i in 1:N){
        Ci = COM[i]
        NB = Neighbor[[i]]
        G = rep(0, N); #Gain vector
        best_increase = -1;
        Cnew = Ci;
        COM[i] = -1;
        #Decrease in the weights of links incident to Ci from outside
        sumTot[Ci] = sumTot[Ci] - K[i];
        CNi1 = which(COM==Ci); #Other members of community Ci 
        #Decrease in the weights of links inside C
        sumIn[Ci] = sumIn[Ci] - 2 * sum(M[i, CNi1]) - M[i, i];
        
        for (j in 1:length(NB)){
          Cj = COM[NB[j]]; #Community of neighbor j of i
          if (length(NB) == 0){
            break
          }
          if (G[Cj] == 0){
            CNj = which(COM == Cj); # Members of community Cj
            Ki_in = 2 * sum(M[i, CNj]); #Sum of weights of the links from node i to nodes in Cj
            G[Cj] = Ki_in / m - 2 * K[i] * sumTot[Cj] / (m*m) #Change in Cost function
            if (ddebug){
              print(paste0("Gain for comm ", Cj-1, " =>", G[Cj]))
            }
            if (G[Cj] > best_increase){
              best_increase = G[Cj];
              Cnew_t = Cj; #Move to i to community of j
            }
          }
        }
        
        if (best_increase > 0){
          Cnew = Cnew_t;
          if (ddebug){
            print(paste0("Move node ", i-1, " => cluster ", Cnew-1))
          }
          Cost[i] = best_increase;
        }
        
        COM[i] = Cnew
        Ck = which(COM==Cnew); #Members of the community
        sumIn[Cnew] = sumIn[Cnew] + 2 * sum(M[i, Ck])
        sumTot[Cnew] = sumTot[Cnew] + K[i];
        if (Cnew != Ci) {
          gain = 1
        }
        
        
      }
      
      sCost = sum(Cost);
      to.reindex = reindex_com(COM) # reindex the communities
      Nco = length(unique(COM)) #Number of communities
      Nco2 = length(to.reindex$S[to.reindex$S>1]);
      mod = compute_modularity(COM, M) #computes modularity. NOTE: Need to make max_modularity compatible with this
      if (debug){
        print(paste0("Iter:", Niter, " - Modularity = ", round(mod,2), " ", Nco, " communities, (", Nco2, " non singlets)"))
      }
      Niter = Niter + 1;
      
    }
    
    
    Niter = Niter - 1;
    main.reindex = reindex_com(COM);
    COMTY = list()
    COMTY$COM[[1]] = main.reindex$C;
    COMTY$SIZE[[1]] = main.reindex$S;
    COMTY$MOD[1] = compute_modularity(main.reindex$C,M);
    COMTY$Niter[1] = Niter;
    
    
    #Perform part 2
    if (s == 1){
      Mnew = M;
      Mold = Mnew;
      COMcur = main.reindex$C; #Current community assignments
      COMfull = COMcur;
      k=2;
      
      Nco2 = length(main.reindex$S[main.reindex$S > 1])
      if (verbose){
        print(paste0("Pass number 1 - ", Nco2, " com (", Niter, " iterations)"))
      }
      
      while(1){
        Mold = Mnew;
        S2 = dim(Mold);
        Nnode = S2[1];
        
        COMu = unique(COMcur);
        Ncom = length(COMu);
        ind_com = matrix(0, nrow=Ncom, ncol=Nnode)
        ind_com_full = matrix(0, nrow=Ncom, ncol=N)
        
        for (p in 1:Ncom){
          ind = which(COMcur==p);
          ind_com[p, c(1:length(ind))] = ind; 
        }
        
        for (p in 1:Ncom){
          ind = which(COMfull==p);
          ind_com_full[p, c(1:length(ind))] = ind; 
        }
        
        Mnew = matrix(0, nrow=Ncom, ncol=Ncom)
        for (m in 1:Ncom){
          for (n in m:Ncom){
            
            ind1 = ind_com[m,]
            ind2 = ind_com[n,]
            Mnew[m,n] = sum(Mold[ind1[ind1 > 0], ind2[ind2>0]]);
            Mnew[n,m] = Mnew[m,n]
          }
        }
        
        cluster_out = cluster_jl(Mnew, s=0, self=self, debug=debug, ddebug=ddebug, verbose=verbose)
        COMTYt = cluster_out[[1]]
        endingt = cluster_out[[2]]
        
        if (endingt != 1){
          COMfull = rep(0, N)
          COMcur = COMTYt$COM[[1]];
          for (p in 1:Ncom){
            ind1 = ind_com_full[p,]
            COMfull[ind1[ind1>0]] = COMcur[p];
          }
          
          to.reindex2 = reindex_com(COMfull) # reindex the communities
          COMTY$COM[[k]] = to.reindex2$C;
          COMTY$SIZE[[k]] = to.reindex2$S;
          COMTY$MOD[k] = compute_modularity(to.reindex2$C, M);
          COMTY$Niter[k] = COMTYt$Niter;
          Nco2 = length(to.reindex2$S[to.reindex2$S > 1])
          
          if (verbose){
            print(paste0("Pass number ", k , " - ", Nco2, " communities"))
          }
          
          Ind = (COMTY$COM[[k]] == COMTY$COM[[k-1]]);
          if (sum(Ind) == length(Ind)){
            if (verbose){
              print(paste0('Identical segmentation => End'))
            }
            to.return[[1]] = COMTY
            return(to.return) 
          }
          
        } else {
          
          if (verbose){
            print(paste0("Empty matrix => End"))
          }   
          return(to.return)
          
        }
        k = k + 1;
      }
      
    } else {
      to.return[[1]] = COMTY
      return(to.return)
    }
    
  }
  
  
  #Reindex community IDs
  reindex_com = function(COMold){
    
    C = rep(0, length(COMold))
    COMu = unique(COMold)
    S = rep(0, length(COMu))
    for (l in 1:length(COMu)){
      S[l] = length(COMold[COMold == COMu[l]])
    }
    
    INDs = order(S, decreasing=TRUE)
    
    for (l in 1:length(COMu)){
      C[COMold == COMu[INDs[l]]] = l
    }
    
    to.return = list()
    to.return$C = C
    to.return$S = S[INDs]
    return(to.return)
  }
  
  #Compute modularity
  compute_modularity = function(C, Mat){
    
    m = sum(Mat);
    MOD = 0;
    COMu = unique(C);
    for (j in 1:length(COMu)){
      Cj = which(C==COMu[j])
      Ec = sum(Mat[Cj, Cj])
      Et = sum(Mat[Cj,]);
      if (Et > 0){
        MOD = MOD + Ec/m - (Et/m)^2;
      }
      
    }
    
    return(MOD)
    
  }
  
  
  #Analyzing Barcode Collisions
  
  BarcodeCollisionPlot = function(seq, title.text){
    
    freq.matrix = NucleotideFrequencyByPosition(seq)
    random.seq = MakeRandomSequences(n=length(seq),freq = freq.matrix)
    
    Lev.dist=adist(seq,seq); diag(Lev.dist) = NA
    Lev.dist.random = adist(random.seq, random.seq); diag(Lev.dist.random) = NA
    
    min.dist = apply(Lev.dist,1,function(x) min(x, na.rm=TRUE))
    min.dist.random = apply(Lev.dist.random,1,function(x) min(x, na.rm=TRUE))
    
    df = data.frame(Observed=min.dist, Random=min.dist.random)
    
    df.melt = melt(df)
    p<- ggplot(df.melt, aes(x=value)) + geom_density(aes(color=factor(variable)), size=1.5) + xlab("Min. Edit Dist.") + scale_x_discrete() + ggtitle(title.text) + xlim(c(0,5)) + 
      theme(axis.title.x = element_text(face="bold", colour="#990000", size=20), axis.text.x  = element_text(angle=0, vjust=0.5, size=16)) + 
      theme(axis.title.y = element_text(face="bold", colour="#990000", size=20), axis.text.y  = element_text(angle=0, vjust=0.5, size=16)) 
    return(p)
    
  }
  
  NucleotideFrequencyByPosition  = function(sequences){
    
    if (class(sequences) == "character")  seq.align=paste0(sequences,collapse="\n")
    freq.matrix = consensusMatrix(DNAStringSet(strsplit(seq.align, "\n")[[1]]))
    freq.matrix = scale(freq.matrix[1:4,], center=FALSE, scale=colSums(freq.matrix[1:4,]))
    return(freq.matrix)
  }
  
  
  MakeRandomSequences = function(n=3000,freq=NULL){
    freq = set.ifnull(freq, matrix(runif(16),nrow=4, dimnames = list(c("A","C","G","T"))))
    freq = scale(freq, scale=colSums(freq), center=FALSE)
    random.seq = c()
    for (i in 1:ncol(freq)){
      if (i==1){
        random.seq = sample(c("A","C","G","T"), n,replace=TRUE,prob=freq.matrix[,i])
      } else {
        random.seq = paste0(random.seq,sample(c("A","C","G","T"), n ,replace=TRUE,prob=freq[,i]))
      }
    }
    return(random.seq)
  }
  
  binompval <- function(p,N,n){
    pval   <- pbinom(n,round(N,0),p)
    pval[!is.na(pval) & pval > 0.5] <- 1-pval[!is.na(pval) & pval > 0.5]
    return(pval)
  }
  
  human_to_mouse_orth = function(query.gene.list=NULL, target.gene.list=NULL, convert=c("Human2Mouse","Mouse2Human"), orthology_file=NULL){
    
    orthology_file=set.ifnull(orthology_file,"~/Dropbox/CompResources/HMD_HumanPhenotype.rpt.txt")
    orth_table = read.table(file=orthology_file, sep="\t", header=FALSE)
    
    if (is.null(target.gene.list)){
      print("Not provided target gene list. Some entries might be missed")
    }
    convert = match.arg(convert)
    
    print(paste0("Performing conversion of genes from ", convert))
    
    df = orth_table[,c("V1","V4")]
    if (convert == "Mouse2Human"){
      colnames(df) = c("V2", "V1")
    } else {
      colnames(df) = c("V1", "V2")
    }
    
    to.return = data.frame(query=c(), target=c(), stringsAsFactors = FALSE)
    for (gene in query.gene.list){
      target = as.character(subset(df, V1==gene)$V2)
      
      if (length(target) == 0){
        
        if (is.null(target.gene.list)){
          print(paste0("No orthologs found for gene ", gene));
          next;
        } else {
          gene1 = capitalize(tolower(gene))
          gene1=gsub("rik","Rik",gene1)
          target = grep(paste0("^",gene1,"$"), target.gene.list, value=TRUE)
          if (length(target)==0){
            print(paste0("No orthologs found for gene ", gene));
            next;
          }
        }
        
      }
      
      if (length(target) >= 1){
        if (!is.null(target.gene.list)){
          target = target[target %in% target.gene.list][1]
        } else{
          target = target[1]
        }
        
      }
      to.return = rbind(to.return, data.frame(query=gene, target=target, stringsAsFactors = FALSE))
    }
    colnames(to.return) = c("query", "target")
    return(to.return)
    
  }
  
  # Plot confusion matrix
  plotConfusionMatrix = function(X,row.scale=TRUE, col.scale=FALSE, col.low="blue", col.high="red", max.size=5, ylab.use="Known", xlab.use="Predicted", order=NULL, x.lab.rot=FALSE, plot.return=TRUE, max.perc=100){
    
    if (!col.scale & row.scale){ X = t(scale(t(X), center=FALSE, scale=rowSums(X)));  X=X*100 }
    if (col.scale & !row.scale){ X = scale(X, center=FALSE, scale=colSums(X)); X = X*100 }
    if(col.scale & row.scale){
      print("Only one of row.scale or col.scale should be true. performing row scaling by default")
      X = t(scale(t(X), center=FALSE, scale=rowSums(X)))
      X=X*100
    }
    X[is.na(X)] = 0
    if (max(X) > 100){
      X=X/100
    }
    
    orig.rownames = rownames(X)
    orig.colnames = colnames(X)
    if (!is.null(order)){
      if (order == "Row"){  
        factor.levels = c()
        for (i1 in colnames(X)){
          if (max(X[,i1]) < 50) next
          ind.sort = rownames(X)[order(X[,i1], decreasing=TRUE)]
          ind.sort = ind.sort[!(ind.sort %in% factor.levels)]
          factor.levels = c(factor.levels, ind.sort[1])
        }
        factor.levels = c(factor.levels, setdiff(rownames(X), factor.levels))
        factor.levels = factor.levels[!is.na(factor.levels)]
      } 
      
      if (order == "Col") {
        factor.levels = c()
        for (i1 in rownames(X)){
          if (max(X[i1,]) < 50) next
          ind.sort = rownames(X)[order(X[i1,], decreasing=TRUE)]
          ind.sort = ind.sort[!(ind.sort %in% factor.levels)]
          factor.levels = c(factor.levels, ind.sort[1])
        }
        factor.levels = c(factor.levels, setdiff(rownames(t(X)), factor.levels))
        factor.levels = factor.levels[!is.na(factor.levels)]
      } 
    } else {
      factor.levels = rownames(t(X))
    }
    
    factor.levels = c(factor.levels, setdiff(rownames(X), factor.levels))
    X = melt(X)
    colnames(X) = c("Known", "Predicted", "Percentage")
    #X$Known = factor(X$Known, levels=rev(unique(X$Known)));
    #X$Predicted = factor(X$Predicted, levels = rev(factor.levels))
    
    if (!is.null(order)){
      if (order == "Row"){ 
        X$Known = factor(X$Known, levels=rev(factor.levels));
        X$Predicted = factor(X$Predicted, levels = orig.colnames)
        
      }
      if (order == "Col"){
        X$Predicted = factor(X$Predicted, levels = factor.levels);
        X$Known = factor(X$Known, levels=rev(orig.rownames));
      }
    } else {
      X$Known = factor(X$Known, levels=rev(unique(X$Known)));
      X$Predicted = factor(X$Predicted, levels=unique(X$Predicted));
    }
    
    #print(sum(is.na(X$Known)))
    
    
    p = ggplot(X, aes(y = Known,  x = Predicted)) + geom_point(aes(colour = Percentage,  size =Percentage)) + 
      scale_color_gradient(low =col.low,   high = col.high, limits=c(0, 100 ))+scale_size(range = c(1, max.size), limits = c(0,max.perc))+   theme_bw() #+nogrid
    p = p + xlab(xlab.use) + ylab(ylab.use) + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
      theme(axis.text.y=element_text(size=12, face="italic"))  
    
    if (x.lab.rot) {
      p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    }
    print(p)
    
    if (plot.return) return(p)
  }
  
  grab10Xmatrix = function(object_10X, collapse.genes=TRUE){
    require(cellrangerRkit)
    Count.mat = object_10X@mat
    
    #Consolidate genes
    if (collapse.genes){
      which.dup = names(table(object_10X@gene_symbols))[table(object_10X@gene_symbols)>1]
      temp.mat = data.frame()
      geneNames = object_10X@gene_symbols
      for (gene in which.dup){
        rows1 = which(geneNames %in% gene)
        geneNames = geneNames[-rows1]
        temp.mat = Matrix::colSums(Count.mat[rows1,])
        Count.mat = Count.mat[-rows1,]
        Count.mat = rbind(Count.mat, temp.mat)
        rownames(Count.mat)[nrow(Count.mat)] = gene
        geneNames = c(geneNames, gene)
      }
      rownames(Count.mat) = geneNames
    }
    
    return(Count.mat)
  }
  
  merge_dfmat = function(m1, m2){
    # merges two matrices or data frames m1 and m2 columnwise while ensuring that mutually exclusive rows are accounted for
    # dim(m)[1] = length(union(rownames(m1), rownames(m2)))
    genes_minus_m2 = rownames(m1)[!(rownames(m1) %in% rownames(m2))]
    genes_minus_m1 = rownames(m2)[!(rownames(m2) %in% rownames(m1))]
    genes_all = union(rownames(m1), rownames(m2))
    m1[genes_minus_m1,] = 0
    m2[genes_minus_m2,] = 0
    m1 = m1[genes_all,]
    m2 = m2[genes_all,]
    m = cbind(m1,m2)
    return(m)
  }
  
  
  #require(networkD3)
  plotSankey = function(ref1, map1, ref2=NULL, map2=NULL, tag1 = "ref", tag2="map1", tag3="map2", no.unmapped = TRUE){
    
    A = list()
    if (no.unmapped){
      print("Removing 0=unmapped cells")
      ref1 = ref1[ref1 != "0"]
      map1 = map1[map1 != "0"]
    }
    
    cells.common1 = intersect(names(ref1),names(map1))
    ref1 = ref1[cells.common1]; map1 = map1[cells.common1]
    ref1 = paste0(tag1,"_",ref1)
    map1 = paste0(tag2,"_", map1)
    
    res1 = reshape2::melt(table(ref1,map1))
    colnames(res1) = c("ref","map","value")
    res1$ref = as.character(res1$ref)
    res1$map = as.character(res1$map)
    res1 = res1[res1$value != 0, ]
    
    
    if (!is.null(ref2) & !is.null(map2)){
      
      if (no.unmapped){
        print("Removing 0=unmapped cells")
        ref2 = ref2[ref2 != "0"]
        map2 = map2[map2 != "0"]
      }
      cells.common2 = intersect(names(ref2),names(map2))
      ref2 = ref2[cells.common2]; map2 = map2[cells.common2]
      ref2 = paste0(tag2,"_",ref2)
      map2 = paste0(tag3,"_",map2)
      
      res2 = reshape2::melt(table(ref2,map2))
      colnames(res2) = c("ref","map","value")
      res2$ref = as.character(res2$ref)
      res2$map = as.character(res2$map)
      res2 = res2[res2$value != 0, ]
      res = rbind(res1,res2)
      
      A$nodes = data.frame(c(rownames(table(ref1,map1)), union(colnames(table(ref1,map1)), rownames(table(ref2,map2))), colnames(table(ref2,map2)))); colnames(A$nodes) = "name"
      A$links = res
      A$links$ref = match(A$links$ref, A$nodes$name)-1
      A$links$map = match(A$links$map, A$nodes$name)-1
      if (!no.unmapped){
        df = as.data.frame(matrix(c(match(paste0(tag2,"_",0), A$nodes$name), match(paste0(tag3,"_",0), A$nodes$name), 1 ), nrow=1))
        colnames(df) = colnames(A$links)
        A$links = rbind(A$links,df)
      }
      
    } else {
      res = res1
      A$nodes = data.frame(c(rownames(table(ref1,map1)), colnames(table(ref1,map1)))); colnames(A$nodes) = "name"
      A$links = res
      A$links$ref = match(A$links$ref, A$nodes$name)-1
      A$links$map = match(A$links$map, A$nodes$name)-1
    }
    
    
    sankeyNetwork(Links = A$links, Nodes = A$nodes, Source = "ref",
                  Target = "map", Value = "value", NodeID = "name",
                  units = "cells", fontSize = 12, nodeWidth = 30)
    
  }
  
  #' Samples a certain number of cells per ident
  SampleCellsPerIdent = function(object, 
                                 max.cells.per.ident = max.cells.per.ident,
                                 id=ident.use,
                                 seed.use=42){
    orig_ident = object@ident
    object = set.all.ident(object, id = ident.use)
    
    set.seed(seed.use)
    cells.use = c()
    for (ident in levels(object@ident)){
      cells.ident = which.cells(object, ident)
      cells.ident = sample(cells.ident, min(max.cells.per.ident, length(cells.ident)))
      cells.use = c(cells.use, cells.ident)
    }
    return(cells.use)
    
  }
  
  #' returns a character vector specifiying the axes for the dimensionality reduction embedding of interest
  GetDimCode <- function(reduction.use){
    if (reduction.use=="pca") dim.code = "PC"
    if (reduction.use=="tsne") dim.code = "tsne_"
    if (reduction.use=="ica") dim.code = "IC"
    if (reduction.use=="mds") dim.code="MDS"
    if (reduction.use=="cca") dim.code="CC"
    if (reduction.use=="UMAP") dim.code="UMAP"
    if (reduction.use=="ForceAtlas") dim.code="ForceAtlas"
    return(dim.code)
  }
  
  #' returns cell embeddings of interest
  GetCellEmbeddings <- function(object, reduction.use="pca", dims.use = NULL, cells.use=NULL){
    
    object.embed <- GetDimReduction(
      object = object,
      reduction.use = reduction.use,
      slot = "cell.embeddings"
    )
    
    if (length(x = object.embed) == 0) {
      stop(paste0("Cell embeddings slot for ", reduction.use, " is empty."))
    }
    
    cells.use = set.ifnull(cells.use, rownames(object.embed))
    
    if (any(! cells.use %in% rownames(x = object.embed))) {
      missing.cells <- paste0(
        cells.use[which(x = ! cells.use %in% rownames(x = object.embed))],
        collapse = ", "
      )
      warning(paste0("Could not find the following cell names: ", missing.cells))
      cells.use <- intersect(x = cells.use, y = rownames(x = object.embed))
    }
    
    dims.use <- set.ifnull(dims.use, 1:ncol(x = object.embed))
    if (any(!dims.use %in% 1:ncol(x = object.embed))) {
      missing.dims <- paste0(
        dims.use[which(x = ! dims.use %in% 1:ncol(x = object.embed))],
        collapse = ", "
      )
      stop(paste0("Could not find the following dimensions: ", missing.dims))
    }
    
    object.embed <- object.embed[cells.use, dims.use, drop = FALSE]
    object.key <- GetDimReduction(
      object = object,
      reduction.use = reduction.use,
      slot = "key"
    )
    
    if (length(x = object.key) == 0) {
      colnames(x = object.embed) <- NULL
    } else {
      colnames(x = object.embed) <- paste0(object.key, dims.use)
    }
    return(object.embed)
    
  }
  
  #' Dimensional Reduction Accessor Function (from Seurat)
  #'
  #' General accessor function for dimensional reduction objects. Pulls slot
  #' contents for specified stored dimensional reduction analysis.
  #'
  #' @param object 
  #' @param reduction.use Type of dimensional reduction to fetch (default is PCA)
  #' @param slot Specific information to pull (must be one of the following:
  #'  "cell.embeddings", "gene.loadings", "gene.loadings.full", "sdev", "key", "misc")
  #'
  #' @return Returns specified slot results from given reduction technique
  #'
  #' @export
  #'
  #' @examples
  #' pbmc_small
  #' # Get the PCA cell embeddings and print the top left corner
  #' GetDimReduction(object = pbmc_small, reduction.type = "pca",
  #'                 slot = "cell.embeddings")[1:5, 1:5]
  #' # Get the standard deviation of each PC
  #' GetDimReduction(object = pbmc_small, reduction.use = "pca", slot = "sdev")
  #'
  GetDimReduction <- function(
    object,
    reduction.use = "pca",
    slot = "gene.loadings"
  ) {
    if (! (reduction.use %in% names(object@dr))) {
      stop(paste(reduction.use, " dimensional reduction has not been computed"))
    }
    reduction <- paste0("object@dr$", reduction.use)
    reduction.slots <- slotNames(x = eval(expr = parse(text = reduction)))
    if (! (slot %in% reduction.slots)) {
      stop(paste0(slot, " slot doesn't exist"))
    }
    return(eval(expr = parse(text = paste0(reduction, "@", slot))))
  }
  
  #' Duplicate a dimensionality reduction slot with a different name
  DuplicateDimReduction <- function(
    object,
    reduction.use = "tsne",
    duplicate.name = "tsne2",
    replace.dr=FALSE
  ) {
    if (! (reduction.use %in% names(object@dr))) {
      stop(paste(reduction.use, " dimensional reduction has not been computed"))
    }
    
    if ((duplicate.name %in% names(object@dr)) & !replace.dr) {
      stop(paste(duplicate.name, " is already present. Specify replace=TRUE if you want to replace"))
    }
    
    eval(expr = parse(text = paste0("object@dr$", duplicate.name, " <- object@dr$", reduction.use)))
    return(object)
  }
  
  
  # Find the quantile of a data
  #
  # Converts a quantile in character form to a number regarding some data
  # String form for a quantile is represented as a number prefixed with 'q'
  # For example, 10th quantile is 'q10' while 2nd quantile is 'q2'
  #
  # Will only take a quantile of non-zero data values
  #
  # @param cutoff The cutoff to turn into a quantile
  # @param data The data to turn find the quantile of
  #
  # @return The numerical representation of the quantile
  #
  SetQuantile <- function(cutoff, data) {
    if (grepl(pattern = '^q[0-9]{1,2}$', x = as.character(x = cutoff), perl = TRUE)) {
      this.quantile <- as.numeric(x = sub(
        pattern = 'q',
        replacement = '',
        x = as.character(x = cutoff)
      )) / 100
      data <- unlist(x = data)
      data <- data[data > 0]
      cutoff <- quantile(x = data, probs = this.quantile)
    }
    return(as.numeric(x = cutoff))
  }
  
  #' vectorized pdist (faster than naive pdist)
  #' Pairwise distances between rows of 
  vectorized_pdist <- function(A,B=NULL){
    if (is.null(B)) B=A
    if (ncol(B) != ncol(A)) stop("Matrices A and B need to have the same number of columns")
    
    an = rowSums(A^2)
    bn = rowSums(B^2)
    
    m = nrow(A)
    n = nrow(B)
    
    tmp = matrix(rep(an, n), nrow=m) 
    tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
    return(sqrt( tmp - 2 * tcrossprod(A,B) ))
  }
  
  
  ## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
  ##   data: a data frame.
  ##   measurevar: the name of a column that contains the variable to be summariezed
  ##   groupvars: a vector containing names of columns that contain grouping variables
  ##   na.rm: a boolean that indicates whether to ignore NA's
  ##   conf.interval: the percent range of the confidence interval (default is 95%)
  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                        conf.interval=.95, .drop=TRUE) {
    library(plyr)
    
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }
    
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                     c(N    = length2(xx[[col]], na.rm=na.rm),
                       mean = mean   (xx[[col]], na.rm=na.rm),
                       sd   = sd     (xx[[col]], na.rm=na.rm)
                     )
                   },
                   measurevar
    )
    
    # Rename the "mean" column    
    datac <- plyr::rename(datac, c("mean" = measurevar))
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
  }
  
  KEst_TracyWidom = function(new_obj, p = 0.001){
    require(RMTstat)
    #p = 0.001
    # get eigenvalues from pca result
    eigenvalue = new_obj@dr$pca@eigvals
    #print(eigenvalue)
    
    # Patterson, N., Price, A.L. & Reich, D. PLoS Genet. 2, e190 (2006).
    n = ncol(new_obj@data) # number of cells
    m = nrow(new_obj@dr$pca@gene.loadings) # number of variable genes
    # get the mean value
    muTW <- ((sqrt(n - 1) + sqrt(m))^2)
    # get sigma
    sigmaTW <- ((sqrt(n - 1) + sqrt(m))) * (1/sqrt(n - 1) + 1/sqrt(m))^(1/3)
    # calculate the boundary at p level
    library(RMTstat)
    bd <- qtw(p, lower.tail = FALSE) * sigmaTW + muTW  # 3.2730 is the p=0.001 percentile point for the Tracy-Widom distribution
    print(bd*1.5)
    n_K = sum(eigenvalue>bd*1.2)
    print(paste0("Found ", n_K, " significant PCs"))
    return(n_K)
  }
  
  
  detachAllPackages <- function() {
    
    basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
    
    package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
    
    package.list <- setdiff(package.list,basic.packages)
    
    if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
    
  }
  
  topGOterms = function( fg.genes = NULL,
                         bg.genes = NULL,
                         organism = "Mouse", 
                         ontology.use = "BP",
                         stats.use = "fisher",
                         algorithm.use = "weight01",
                         topnodes.print=20,
                         num.char=100,
                         return.object = FALSE){
    
    if (is.null(fg.genes) | is.null(bg.genes)){
      stop("Error : Both gene lists are empty")
    }
    
    require(topGO)
    if (organism == "Mouse"){
      mapping.use = "org.Mm.eg.db"
      library(org.Mm.eg.db)
    } else if (organism == "Human"){
      mapping.use = "org.Hs.eg.db"
      library(org.Hs.eg.db)
    } else {
      stop("Error : Organisms other than mouse not supported currently")
    }
    
    n = length(bg.genes)
    geneList = integer(n)
    names(geneList) = bg.genes
    geneList[intersect(names(geneList), fg.genes)]=1
    print(paste0("Total ", length(geneList), " genes. ", sum(geneList), " genes in the foreground"))
    geneList = factor(geneList)
    
    if (ontology.use %in% c("BP", "CC", "MF")){
      print(paste0("Using Ontology : ", ontology.use))
    } else {
      stop("Error: Ontology not available. Should be one of BP, CC or MF")
    }
    # Make GO object
    GOdata <- new("topGOdata",
                  description = "GOanalysis",
                  ontology = ontology.use,
                  allGenes = geneList,
                  annot = annFUN.org,
                  mapping = mapping.use,
                  ID = "SYMBOL",
                  nodeSize = 10)
    print(paste0("Using the ", stats.use, " statistic with the ", algorithm.use, " algorithm"))
    res.result <- runTest(GOdata, statistic = stats.use, algorithm = algorithm.use)
    to.return = list()
    to.return$GOdata = GOdata
    to.return$res.table <- GenTable(GOdata, pval = res.result, topNodes = topnodes.print, numChar = num.char)
    if (return.object){
      return(GOdata)
    } else {
      return(to.return)
    }
    
  }
  
  # Random Forest assisted pruning
  # regress RGC
  # Random Forest classification
  ValidateClustersRF = function(object, 
                                var.genes = NULL, 
                                regress.use=TRUE, 
                                max.cells.per.ident=300, 
                                train.frac=0.7, 
                                predict.all=TRUE,
                                do.prob=FALSE,
                                margin.use=NULL,
                                margin.with.next=0.05,
                                return.new.predict=TRUE,
                                ntree.use=201){
    library(randomForest)
    
    if (is.null(var.genes)){
      var.genes = object@var.genes
    }
    
    if (regress.use){
      if (object.size(object@regress.data) < 100){ 
        predictor_Data = as.matrix(object@data[var.genes,])
      } else {
        predictor_Data = as.matrix(object@regress.data[var.genes,])
      }
      
    } else {
      predictor_Data = as.matrix(object@data[var.genes,])
    }
    
    training.set = c(); test.set=c()
    training.label = c(); test.label=c();
    print(paste0("Using mininum of ", 0.5*100, " percent cells or ", max.cells.per.ident, " cells per cluster for training"))
    for (i in as.character(levels(object@ident))){
      cells.in.clust = which.cells(object,i);
      n = min(max.cells.per.ident, round(length(cells.in.clust)*train.frac))
      train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
      test.temp = setdiff(cells.in.clust, train.temp)
      training.set = c(training.set,train.temp); test.set=c(test.set,test.temp)
      training.label = c(training.label, rep(i,length(train.temp))); test.label = c(test.label, rep(i, length(test.temp)));
    }
    
    tmp = as.vector(table(training.label))
    sampsizes = rep(min(tmp),length(tmp))
    rf =randomForest(x=t(predictor_Data[,training.set]), y=factor(training.label), importance = TRUE, ntree = ntree.use, keep.inbag=FALSE, replace=FALSE) 
    Conf_OOB0 = rf$confusion
    
    
    to.return = list()
    to.return$Conf_OOB = Conf_OOB0 # confusion matrix
    to.return$var.genes = var.genes 
    to.return$rf = rf
    
    if (predict.all){
      if (do.prob){
        to.return$predict = predict(rf,t(predictor_Data), type="prob")
        if (return.new.predict){
          if (is.null(margin.use)){
            margin.use = 5 * 1 / ncol(to.return$predict)
          }
          new.predict = apply(to.return$predict, 1, function(x){ y = sort(x,decreasing = TRUE)[c(1,2)]; if (( max(x) > margin.use) & y[1] - y[2] > margin.with.next){ 
            colnames(to.return$predict)[which.max(x)] } else {0}})
        }
        to.return$new.predict = new.predict
        
      } else {
        
        to.return$predict = predict(rf,t(predictor_Data))
      }
      
      
    }
    
    return(to.return)
    
  }
  
  # Random Forest based cleaning
  RFcleaning = function(object, 
                        rf_model = NULL,
                        var.genes = NULL, 
                        regress.use=TRUE, 
                        pmax.thresh = 0.15,
                        margin.use=NULL,
                        margin.with.next=0.03,
                        new.ident = "m_RF_clean",
                        do.plot = TRUE
  ){
    if (is.null(var.genes)){
      var.genes = object@var.genes
    }
    
    if (is.null(rf_model) | class(rf_model) != "randomForest"){
      stop("Error: provide valid random forest model")
      
    }
    
    if (regress.use){
      if (object.size(object@regress.data) < 100){ 
        predictor_Data = as.matrix(object@data[var.genes,])
      } else {
        predictor_Data = as.matrix(object@regress.data[var.genes,])
      }
      
    } else {
      predictor_Data = as.matrix(object@data[var.genes,])
    }
    
    require(randomForest)
    rf_probs = predict(rf_model, t(predictor_Data), type="prob")
    new.predict = apply(rf_probs, 1, function(x){ 
      y = sort(x,decreasing = TRUE)[c(1,2)]; 
      if (( max(x) > pmax.thresh) & (y[1] - y[2] > margin.with.next)){ 
        colnames(rf_probs)[which.max(x)] } else {"0"}
    })
    new.predict = as.numeric(new.predict); names(new.predict) = rownames(rf_probs)
    new.predict = factor(new.predict)
    Conf = table(object@ident, new.predict)
    
    # Plot the "0" cells on tSNE
    cells.0 = names(new.predict)[new.predict == 0]
    if (do.plot) plot.cluster(object, cells.plot = cells.0, nCol = 1, pt.size = 0.7, main.use = paste0("Unassigned cells ( n =", length(cells.0), ")") )
    
    # Also plot the cells that have been reassigned
    cells.reassign = names(new.predict)[((as.numeric(new.predict)-1) != as.numeric(object@ident)) & as.numeric(new.predict) != 1]
    if (do.plot) plot.cluster(object, cells.plot = cells.reassign, nCol = 1, pt.size = 0.7, main.use = paste0("Reassigned cells ( n =", length(cells.reassign), ")") )
    
    object@data.info[,new.ident] = new.predict[rownames(object@data.info)]
    cells.use = setdiff(object@cell.names, cells.0)
    object = subsetData(object, cells.use = cells.use)
    object@data.info[,new.ident] = droplevels(object@data.info[,new.ident])
    object = set.all.ident(object, id=new.ident)
    
    return(object)
    
  }
  
  #' Compare two clusterings
  #' 
  #' 
  CompareClusterings = function(old.clust.labels, new.clust.labels){
    cells.use = intersect(names(old.clust.labels), names(new.clust.labels))
    A = table(old.clust.labels[cells.use], new.clust.labels[cells.use])
    Arow = t(scale(t(A), center=FALSE, scale=rowSums(A)))
    for (x in 1:nrow(A)){
      if (max(Arow[x,]) < 0.95){
        new.clust = unname(which(Arow[x,] > 0.03))
        new.clust = new.clust[order(unname(Arow[x,new.clust]), decreasing=TRUE)]
        print(paste0("Old Cluster ", x, " split into new clust ", paste0(new.clust, collapse = ","), " ( ", paste0(round(100*Arow[x,new.clust]), collapse=","), ")"))
      }
    }
    
    Acol = scale(A, center=FALSE, scale=colSums(A))
    for (x in 1:ncol(A)){
      if (max(Acol[,x]) < 0.95){
        old.clust = unname(which(Acol[,x] > 0.03))
        old.clust = old.clust[order(unname(Acol[old.clust,x]), decreasing=TRUE)]
        print(paste0("New Cluster ", x, " came from old clust ", paste0(old.clust, collapse = ","), " ( ", paste0(round(100*Acol[old.clust,x]), collapse=","), ")"))
      } 
    }
  }
  
  compRandIndex = function(ident1, ident2){
    require(mclust)
    cells.use = intersect(names(ident1), names(ident2))
    return(adjustedRandIndex(ident1[cells.use], ident2[cells.use]))
  }
  
  compNMI= function(ident1, ident2){
    require(NMI)
    cells.use = intersect(names(ident1), names(ident2))
    df_1 = data.frame(c1 = names(ident1), c2 = ident1)
    df_2 = data.frame(c1 = names(ident2), c2 = ident2)
    NMIval = NMI(df_1[cells.use,], df_2[cells.use,])$value
    return(NMIval)
  }
  
  subsample_by_ident = function(cells.ident, downsample.frac=0.5, min.cells.per.ident=50, max.cells.per.ident = 1000){
    unique.idents = names(table(cells.ident))
    cells.subsample = c()
    for (i in unique.idents){
      cells.temp = names(cells.ident)[cells.ident == i]
      n=length(cells.temp)
      cells.samp = sample(cells.temp, 
                          max(downsample.frac*length(cells.temp),min(min.cells.per.ident, n)))
      if (length(cells.samp) > max.cells.per.ident){
        cells.samp = sample(cells.samp, max.cells.per.ident)
      }
      cells.subsample = c(cells.subsample, cells.samp)
    }
    return(cells.subsample)
    
  }
  
  
