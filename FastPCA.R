
FastPCA <- function(Df,top.k,center.use = T, scale.use = F,iterative = F, check = F) {
  #This function does PCA using a fast randomized projection method for doing the SVD
  
  #Check for issues in the matrix
  #Significantly increases the run time, only use it if you are getting errors
  if (check == T) {
    #make sure it is all numeric
    if (sum(unlist(lapply(Df,typeof))%in%c("character","logical")) > 0) {
      stop("You have logical or character columns in your data frame")
    }
    #Check for NAs
    if (sum(apply(Df,2,function(x) sum(is.na(x)))) > 0) {
      stop("Columns have NA in them, either use na.omit() to remove them, or some sort of Imputation")
    }
    #Check to make sure this is a matrix
    if (!is.matrix(Df)) {
      stop("Your data is not a matrix, use as.matrix() first")
    }
  }
  
  #First center the data
  Df.norm <- scale(Df, center = center.use,scale = scale.use)
  cen <- attr(Df.norm, "scaled:center")
  sc <- attr(Df.norm, "scaled:scale")
  
  #Do the Randomized SVD
  svd.results <- FastSVD(Df.norm,k = top.k,iterative)
  
  #Return the results in the same format as prcomp()
  dimnames(svd.results$v) <- list(colnames(Df.norm), paste0("PC", seq_len(ncol(svd.results$v))))
  r <- list(sdev = svd.results$d/sqrt(max(1, nrow(Df.norm) - 1)), rotation = svd.results$v , center = if (is.null(cen)) FALSE else cen, 
            scale = if (is.null(sc)) FALSE else sc)
  r$x <- Df.norm %*% svd.results$v
  class(r) <- "fastprcomp"
  return(r)
}

summary.fastprcomp <- function(object,Df) {
  total.var <- sum(apply(Df,2,var))
  
  vars <- object$sdev^2
  importance <- rbind(`Standard deviation` = object$sdev, `Proportion of Variance` = round(vars / sum(vars), 5), `Cumulative Proportion` = round(cumsum(vars) / total.var, 5))
  colnames(importance) <- colnames(object$rotation)
  
  print(sprintf("Total amount of variance explained is %f",sum(vars) / total.var))
  print(importance)
}

FastSVD <- function(mat.use, k, iterative = F) {
  #Based on http://arxiv.org/pdf/0909.4061v2.pdf
  
  n <- dim(mat.use)[1]
  m <- dim(mat.use)[2]
  
  Omega <- matrix(rnorm((2 * k) * m), ncol=2 * k)
  
  
  #random projection and orthonormal decomposition
  if (iterative == F) {
    Y <- mat.use %*% Omega
    Q <- qr.Q(qr(Y))
  }
  else if (iterative == T) {
    Q <- RandomizedSubspaceInteration(mat.use,Omega)
  }
  
  #Projection onto this subspace
  B <- t(Q) %*% mat.use
  # decomposing B gives us singular values and right vectors for A  
  s <- svd(B)
  U <- Q %*% s$u
  
  #Now threshold with k  
  s$v <- s$v[,1:k]
  s$d <- s$d[1:k]
  
  #Return in the same format as svd()
  return (list(u=U, v=s$v, d=s$d))
}

RandomizedSubspaceInteration <- function(A, Omega, q = 2) {
  #Algorithm 4.4 from http://arxiv.org/pdf/0909.4061v2.pdf
  #Finds Q based on an iterative approach
  #Increases accuracy
  
  Y <- A %*% Omega
  Q <- qr.Q(qr(Y))
  for (j in 1:q) {
    Yhat <- t(A) %*% Q
    Qhat <- qr.Q(qr(Yhat))
    Y <- A %*% Qhat
    Q <- qr.Q(qr(Y))
  }
  return(Q)
}
