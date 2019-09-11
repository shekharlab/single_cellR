glmDEfit = function (object, covariates = NULL,cells.use=NULL, genes.test=NULL, expression_family = "quasipoisson", reduction_method = "UMAP", cores = 1, clean_model = TRUE, scale_covariates=TRUE,
          verbose = FALSE, ...) 
{
  
  genes.test = set.ifnull(genes.test, object@var.genes)
  genes.test = set.ifnull(genes.test, rownames(object@data))
  cells.use = set.ifnull(cells.use, object@cell.names)
  covariates = set.ifnull(covariates, "nTranscripts")
  model_formula_str = paste0("~", paste0(covariates, collapse = "+"))
  
  if (!("size_factor" %in% colnames(object@data.info) )){
    object@data.info$size_factor = object@data.info$nTranscripts / median(object@data.info$nTranscripts)
  }
  
  
  
  model_form <- stats::as.formula(model_formula_str)
  coldata_df = object@data.info[cells.use,]
  tryCatch({
    coldata_df$cluster = object@ident[cells.use]
    #coldata_df$partition = partitions(cds, reduction_method)[colnames(cds)]
    #coldata_df$pseudotime = pseudotime(cds, reduction_method)
  }, error = function(e) {
  })
  
  tryCatch({
    model.frame(model_form, data = coldata_df[1, ])
  }, error = function(e) {
    stop("Error: model formula refers to something not found in object@data.info")
  })
  
  if (scale_covariates){
    for (co in covariates){
      if (is.numeric(object@data.info[,co])){
        object@data.info[,co] = as.numeric(scale(rgc_sub@data.info[,co], center = mean(object@data.info[cells.use,co]), scale = sd(object@data.info[cells.use,co])))
      } 
    }
  }
  
  disp_func <- NULL
  if (cores > 1) {
    fits <- mc_es_apply_glm(object, 1, fit_glm_helper, cells.use, genes.test, required_packages = c("BiocGenerics", 
                                                                        "Biobase", "MASS", "purrr", "pscl", "speedglm", "dplyr", 
                                                                        "Matrix"), cores = cores, reduction_method = reduction_method, 
                        model_formula_str = model_formula_str, expression_family = expression_family, 
                        disp_func = disp_func, clean_model = clean_model, 
                        verbose = verbose, ...)
    fits
  }
  else {
    fits <- smart_es_apply_glm(object, 1, fit_glm_helper, cells.use, genes.test, convert_to_dense = TRUE, 
                           model_formula_str = model_formula_str, expression_family = expression_family, 
                           reduction_method = reduction_method, disp_func = disp_func, 
                           clean_model = clean_model, verbose = verbose, ...)
    fits
  }
  return(fits)
}


########################

mc_es_apply_glm = function (object, MARGIN, FUN, cells.use,genes.test, required_packages, cores = 1, convert_to_dense = TRUE, reduction_method = "UMAP", ...) 
{
  parent <- environment(FUN)
  if (is.null(parent)) 
    parent <- emptyenv()
  e1 <- new.env(parent = parent)
  coldata_df = object@data.info[cells.use,]
  tryCatch({
    coldata_df$cluster = object@ident[cells.use]
    #coldata_df$partition = partitions(cds, reduction_method)[colnames(cds)]
    #coldata_df$pseudotime = pseudotime(cds, reduction_method)
  }, error = function(e) {
  })
  #tryCatch({
  #  coldata_df$pseudotime = pseudotime(cds)
  #}, error = function(e) {
  #})
  Biobase::multiassign(names(as.data.frame(coldata_df)), as.data.frame(coldata_df), 
                       envir = e1)
  environment(FUN) <- e1
  platform <- Sys.info()[["sysname"]]
  old_omp_num_threads = as.numeric(Sys.getenv("OMP_NUM_THREADS"))
  if (is.na(old_omp_num_threads)) {
    old_omp_num_threads = 1
  }
  RhpcBLASctl::omp_set_num_threads(1)
  old_blas_num_threads = as.numeric(Sys.getenv("OPENBLAS_NUM_THREADS"))
  if (is.na(old_omp_num_threads)) {
    old_blas_num_threads = 1
  }
  RhpcBLASctl::blas_set_num_threads(1)
  if (platform == "Windows") 
    cl <- parallel::makeCluster(cores)
  if (platform %in% c("Linux", "Darwin")) 
    cl <- parallel::makeCluster(cores, type = "FORK")
  cleanup <- function() {
    parallel::stopCluster(cl)
    RhpcBLASctl::omp_set_num_threads(old_omp_num_threads)
    RhpcBLASctl::blas_set_num_threads(old_blas_num_threads)
  }
  on.exit(cleanup)
  if (is.null(required_packages) == FALSE) {
    BiocGenerics::clusterCall(cl, function(pkgs) {
      for (req in pkgs) {
        library(req, character.only = TRUE, warn.conflicts = FALSE, 
                quietly = TRUE, verbose = FALSE)
      }
    }, required_packages)
  }
  if (MARGIN == 1) {
    suppressWarnings(res <- sparse_par_r_apply_glm(cl, object@count.data[genes.test,cells.use], 
                                               FUN, convert_to_dense, ...))
  }
  else {
    suppressWarnings(res <- sparse_par_c_apply_glm(cl, object@count.data[genes.test,cells.use], 
                                               FUN, convert_to_dense, ...))
  }
  res
}


########## TO EDIT ############

smart_es_apply_glm = function (object, MARGIN, FUN, cells.use,genes.test, convert_to_dense, reduction_method = "UMAP", 
          ...) 
{
  parent <- environment(FUN)
  if (is.null(parent)) 
    parent <- emptyenv()
  e1 <- new.env(parent = parent)
  coldata_df = object@data.info[cells.use,]
  tryCatch({
    coldata_df$cluster = object@ident[cells.use]
    #coldata_df$partition = partitions(cds, reduction_method)[colnames(cds)]
    #coldata_df$pseudotime = pseudotime(cds, reduction_method)
  }, error = function(e) {
  })
  #tryCatch({
  #  coldata_df$pseudotime = pseudotime(cds)
  #}, error = function(e) {
  #})
  Biobase::multiassign(names(as.data.frame(coldata_df)), as.data.frame(coldata_df), 
                       envir = e1)
  environment(FUN) <- e1
  if (is_sparse_matrix(object@count.data[genes.test,cells.use])) {
    res <- sparse_apply_glm(object@count.data[genes.test,cells.use], 
                        MARGIN, FUN, convert_to_dense, ...)
  }
  else {
    res <- pbapply::pbapply(object@count.data[genes.test,cells.use], 
                            MARGIN, FUN, ...)
  }
  if (MARGIN == 1) {
    names(res) <- genes.test
  }
  else {
    names(res) <- cells.use
  }
  res
}

#######################################

sparse_par_r_apply_glm = function (cl, x, FUN, convert_to_dense, ...) 
{
  par_res <- do.call(c, BiocGenerics::clusterApply(cl = cl, 
                                                   x = split_rows(x, length(cl)), fun = sparse_apply_glm, MARGIN = 1L, 
                                                   FUN = FUN, convert_to_dense = convert_to_dense, ...), 
                     quote = TRUE)
  names(par_res) <- row.names(x)
  par_res
}

# # #

sparse_par_c_apply_glm = function (cl = NULL, x, FUN, convert_to_dense, ...) 
{
  par_res <- do.call(c, BiocGenerics::clusterApply(cl = cl, 
                                                   x = split_cols(x, length(cl)), fun = sparse_apply_glm, MARGIN = 2L, 
                                                   FUN = FUN, convert_to_dense = convert_to_dense, ...), 
                     quote = TRUE)
  names(par_res) <- colnames(x)
  par_res
}

# # #

sparse_apply_glm = function (Sp_X, MARGIN, FUN, convert_to_dense, ...) 
{
  if (convert_to_dense) {
    if (MARGIN == 1) {
      Sp_X <- Matrix::t(Sp_X)
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(as.matrix(Sp_X[, i]), ...)
      }, FUN, ...)
    }
    else {
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(as.matrix(Sp_X[, i]), ...)
      }, FUN, ...)
    }
  }
  else {
    if (MARGIN == 1) {
      Sp_X <- Matrix::t(Sp_X)
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(Sp_X[, i], ...)
      }, FUN, ...)
    }
    else {
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(Sp_X[, i], ...)
      }, FUN, ...)
    }
  }
  return(res)
}

###################

fit_glm_helper = function (x, model_formula_str, expression_family, disp_func = NULL, 
                           clean_model = TRUE, verbose = FALSE, ...) 
{
  model_formula_str <- paste("f_expression", model_formula_str, 
                             sep = "")
  orig_x <- x
  if (expression_family %in% c("negbinomial", "poisson", "zinegbinomial", 
                               "zipoisson", "quasipoisson")) {
    x <- x/size_factor
    f_expression <- round(x)
  }
  else if (expression_family %in% c("binomial", "gaussian")) {
    f_expression <- x
  }
  else {
    f_expression <- log10(x)
  }
  f_expression = as.numeric(f_expression)
  model_formula = stats::as.formula(model_formula_str)
  tryCatch({
    if (verbose) 
      messageWrapper = function(expr) {
        expr
      }
    else messageWrapper = suppressWarnings
    FM_fit = messageWrapper(switch(expression_family, negbinomial = MASS::glm.nb(model_formula, 
                                                                                 epsilon = 0.001, model = FALSE, y = FALSE, ...), 
                                   poisson = speedglm::speedglm(model_formula, family = stats::poisson(), 
                                                                acc = 0.001, model = FALSE, y = FALSE, ...), 
                                   quasipoisson = speedglm::speedglm(model_formula, 
                                                                     family = stats::quasipoisson(), acc = 0.001, 
                                                                     model = FALSE, y = FALSE, ...), binomial = speedglm::speedglm(model_formula, 
                                                                                                                                   family = stats::binomial(), acc = 0.001, model = FALSE, 
                                                                                                                                   y = FALSE, ...), zipoisson = pscl::zeroinfl(model_formula, 
                                                                                                                                                                               dist = "poisson", ...), zinegbinomial = pscl::zeroinfl(model_formula, 
                                                                                                                                                                                                                                      dist = "negbin", ...)))
    FM_summary = summary(FM_fit)
    if (clean_model) 
      FM_fit = clean_glm_model_object(FM_fit)
    df = list(model = FM_fit, model_summary = FM_summary)
    df
  }, error = function(e) {
    if (verbose) {
      print(e)
    }
    list(model = NA, model_summary = NA)
  })
}


######

split_rows = function (x, ncl) 
{
  lapply(parallel::splitIndices(nrow(x), ncl), function(i) x[i, 
                                                             , drop = FALSE])
}

split_cols = function (x, ncl) 
{
  lapply(parallel::splitIndices(ncol(x), ncl), function(i) x[, 
                                                             i, drop = FALSE])
}


######

clean_glm_model_object = function (model) 
{
  if (class(model)[1] == "negbin") {
    model = clean_glm_mass_model_object(model)
  }
  else if (length(intersect(class(model), c("speedglm"))) >= 
           1) {
    model = clean_speedglm_model_object(model)
  }
  else if (class(model) == "zeroinfl") {
    model = clean_zeroinfl_model_object(model)
  }
  else {
    stop("Unrecognized model class")
  }
}


clean_glm_mass_model_object = function (cm) 
{
  cm$y <- c()
  cm$residuals <- c()
  cm$fitted.values <- c()
  cm$effects <- c()
  cm$qr$qr <- c()
  cm$linear.predictors <- c()
  cm$weights <- c()
  cm$prior.weights <- c()
  cm$data <- c()
  cm$family$variance <- c()
  cm$family$dev.resids <- c()
  cm$family$aic <- c()
  cm$family$validmu <- c()
  cm$family$simulate <- c()
  cm$family$control <- c()
  attr(cm$terms, ".Environment") <- c()
  attr(cm$formula, ".Environment") <- c()
  return(cm)
}


clean_speedglm_model_object = function (cm) 
{
  cm$y = c()
  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()
  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$aic = c()
  cm$family$validmu = c()
  cm$family$simulate = c()
  cm$family$control = c()
  attr(cm$terms, ".Environment") = c()
  attr(cm$formula, ".Environment") = c()
  return(cm)
}


clean_zeroinfl_model_object = function (cm) 
{
  cm$y = c()
  cm$model = c()
  cm$residuals = c()
  cm$fitted.values = c()
  cm$weights = c()
  cm$data = c()
  attr(cm$terms$count, ".Environment") = c()
  attr(cm$terms$zero, ".Environment") = c()
  attr(cm$terms$full, ".Environment") = c()
  attr(cm$formula, ".Environment") = c()
  return(cm)
}


extract_model_status_helper = function (model) 
{
  if (class(model)[1] == "speedglm") {
    status_str <- ifelse(model$convergence, "OK", "FAIL")
    return(status_str)
  }
  else if (class(model)[1] == "negbin") {
    status_str <- ifelse(model$converged, "OK", "FAIL")
    return(status_str)
  }
  else if (class(model) == "zeroinfl") {
    status_str <- ifelse(model$converged, "OK", "FAIL")
    return(status_str)
  }
  else {
    return("FAIL")
  }
}


Coef_table = function (model_tbl) 
{
  M_f = model_tbl %>% dplyr::mutate(terms = purrr::map2(.f = purrr::possibly(extract_coefficient_helper, 
                                                                             NA_real_), .x = model, .y = model_summary)) %>% tidyr::unnest(terms)
  M_f = M_f %>% dplyr::group_by(model_component, term) %>% 
    dplyr::mutate(q_value = stats::p.adjust(p_value)) %>% 
    dplyr::ungroup()
  return(M_f)
}


extract_coefficient_helper = function (model, model_summary, pseudo_count = 0.01) 
{
  if (class(model)[1] == "speedglm") {
    coef_mat <- model_summary$coefficients
    coef_mat <- apply(coef_mat, 2, function(x) {
      as.numeric(as.character(x))
    })
    row.names(coef_mat) = row.names(model_summary$coefficients)
    colnames(coef_mat) = c("estimate", "std_err", "test_val", 
                           "p_value")
    log_eff_over_int = log2((model$family$linkinv(coef_mat[, 
                                                           1] + coef_mat[1, 1]) + pseudo_count)/rep(model$family$linkinv(coef_mat[1, 
                                                                                                                                  1]) + pseudo_count, times = nrow(coef_mat)))
    log_eff_over_int[1] = 0
    coef_mat = tibble::as_tibble(coef_mat, rownames = "term")
    coef_mat$normalized_effect = log_eff_over_int
    coef_mat$model_component = "count"
    return(coef_mat)
  }
  else if (class(model)[1] == "negbin") {
    coef_mat = model_summary$coefficients
    coef_mat = apply(coef_mat, 2, function(x) {
      as.numeric(as.character(x))
    })
    row.names(coef_mat) = row.names(model_summary$coefficients)
    colnames(coef_mat) = c("estimate", "std_err", "test_val", 
                           "p_value")
    log_eff_over_int = log2((model$family$linkinv(coef_mat[, 
                                                           1] + coef_mat[1, 1]) + pseudo_count)/rep(model$family$linkinv(coef_mat[1, 
                                                                                                                                  1]) + pseudo_count, times = nrow(coef_mat)))
    log_eff_over_int[1] = 0
    coef_mat = tibble::as_tibble(coef_mat, rownames = "term")
    coef_mat$normalized_effect = log_eff_over_int
    coef_mat$model_component = "count"
    return(coef_mat)
  }
  else if (class(model) == "zeroinfl") {
    count_coef_mat = model_summary$coefficients$count
    colnames(count_coef_mat) = c("estimate", "std_err", "test_val", 
                                 "p_value")
    log_eff_over_int = log2((model$linkinv(count_coef_mat[, 
                                                          1] + count_coef_mat[1, 1]) + pseudo_count)/rep(model$linkinv(count_coef_mat[1, 
                                                                                                                                      1]) + pseudo_count, times = nrow(count_coef_mat)))
    log_eff_over_int[1] = 0
    count_coef_mat = tibble::as_tibble(count_coef_mat, rownames = "term")
    count_coef_mat$normalized_effect = log_eff_over_int
    count_coef_mat$model_component = "count"
    zero_coef_mat = model_summary$coefficients$zero
    colnames(zero_coef_mat) = c("estimate", "std_err", "test_val", 
                                "p_value")
    zero_coef_mat = tibble::as_tibble(zero_coef_mat, rownames = "term")
    zero_coef_mat$normalized_effect = NA
    zero_coef_mat$model_component = "zero"
    coef_mat = dplyr::bind_rows(count_coef_mat, zero_coef_mat)
    return(coef_mat)
  }
  else {
    coef_mat = matrix(NA, nrow = 1, ncol = 5)
    colnames(coef_mat) = c("estimate", "std_err", "test_val", 
                           "p_value", "normalized_effect")
    coef_mat = tibble::as_tibble(coef_mat)
    coef_mat$term = NA
    coef_mat$model_component = NA
    return(coef_mat)
  }
}
