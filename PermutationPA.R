#' Permutation Parallel Analysis
#'
#' Estimate a number of significant principal components from a permutation test.
#'
#' Adopted from sva::num.sv, and based on Buja and Eyuboglu (1992)
#'
#' @param dat a data matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param B a number (a positive integer) of resampling iterations.
#' @param threshold a numeric value between 0 and 1 to threshold p-values.
#' @param verbose a logical indicator as to whether to print the progress.
#' @param seed a seed for the random number generator.
#'
#' @return \code{permutationPA} returns
#' \item{p}{a list of p-values for significance of principal components}
#' \item{r}{an estimated number of significant principal components based on thresholding p-values at \code{threshold}}
#'
#' @references Buja A and Eyuboglu N. (1992) Remarks on parrallel analysis. Multivariate Behavioral Research, 27(4), 509-540
#' @export permutationPA
permutationPA = function (dat, B = 100, threshold = 0.05, verbose=TRUE, seed = NULL, max.pc=1000, n.cores=1) {
  ptm = proc.time()
  if (!is.null(seed)) set.seed(seed)
  n <- min(max.pc, ncol(dat))
  m <- nrow(dat)
  
  print(paste0("Considering only the top ", n, " PCs. Supply max.pc if you wish to change"))
  #uu <- fast.svd(dat, tol = 0)
  uu <- FastSVD(dat,n)
  ndf <- n - 1
  dstat <- uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
  dstat0 <- matrix(0, nrow = B, ncol = ndf)
  if(verbose==TRUE) message("Estimating a number of significant principal component: ")
  
  #permutations
  if (n.cores == 1){
        for (i in 1:B) {
          if(verbose==TRUE) cat(paste(i," "))
          dat0 <- t(apply(dat, 1, sample, replace = FALSE))
          #uu0 <- fast.svd(dat0, tol = 0)
          uu0 <- FastSVD(dat0,n)
          dstat0[i, ] <- uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
        }
  } else {
    require(parallel)
    require(foreach)
    require(doParallel)
    cl = makePSOCKcluster(n.cores, outfile="")
    registerDoParallel(cl, cores=n.cores)
    chunksize = B / n.cores
    vals = split(1:B, ceiling(seq_along(1:(B/chunksize))))
    dstat0 = foreach(run.id = 1:n.cores, .export = "FastSVD", .combine=cbind) %dopar% {
       v = vals[[run.id]]
       do.call(rbind, lapply(v, function(i) {
         if (verbose==TRUE) cat(paste(i," "))
         dat0 <- t(apply(dat, 1, sample, replace = FALSE));
         uu0 <- FastSVD(dat0,n);
         uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
       }))
    }
    cat("\nUnregistering parallel backend..")
    stopCluster(cl)
    registerDoSEQ()
    cat(" done\n");
  }
  p <- rep(1, n)
  for (i in 1:ndf) {
    p[i] <- mean(dstat0[, i] >= dstat[i])
  }
  for (i in 2:ndf) {
    p[i] <- max(p[(i - 1)], p[i])
  }
  r <- sum(p <= threshold)
  y = proc.time() - ptm
  cat(sprintf("PC permutation test completed. Runtime: %s s\n ", signif(y[["elapsed"]], 3)))
  return(list(r = r, p = p))
}
