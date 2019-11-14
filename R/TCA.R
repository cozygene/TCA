#' @title Fitting the TCA model
#'
#' @description Fits the TCA model for an input matrix of observations coming from a mixture of \code{k} sources, under the assumption that each observation is a mixture of unique source-specific values (in each feature in the data). For example, in the context of tissue-level bulk DNA methylation data coming from a mixture of cell types (i.e. the input is methylation sites by individuals), \code{tca} allows to model the methylation of each individual as a mixture of cell-type-specific methylation levels that are unique to the individual.
#'
#' @param X An \code{m} by \code{n} matrix of measurements of \code{m} features for \code{n} observations. Each column in \code{X} is assumed to be a mixture of \code{k} different sources. Note that \code{X} must include row names and column names and that NA values are currently not supported. \code{X} should not include features that are constant across all observations.
#' @param W An \code{n} by \code{k} matrix of weights - the weights of \code{k} sources for each of the \code{n} mixtures (observations). All the weights must be positive and each row, corresponding to the weights of a single observation, must sum up to 1. Note that \code{W} must include row names and column names and that NA values are currently not supported. In case where only initial estimates of \code{W} are available, \code{tca} can be set to re-estimate \code{W} (see \code{refit_W}).
#' @param C1 An \code{n} by \code{p1} design matrix of covariates that may affect the hidden source-specific values (possibly a different effect on each source). Note that \code{C1} must include row names and column names and should not include an intercept term. NA values are currently not supported.
#' @param C1.map An \code{p1} by \code{k} matrix of 0/1 values, indicating for each of the \code{p1} covariates in \code{C1} whether to consider its potential effects on the values of each of the \code{k} sources (e.g., if position \code{i,j} in \code{C1.map} is 1 then the potential effect of the \code{i}-th covariate in \code{C1} on the {j}-th source will be considered). If \code{C1.map == NULL} then effects for all covariates in \code{C1} in each of the sources will be considered.
#' @param C2 An \code{n} by \code{p2} design matrix of covariates that may affect the mixture (i.e. rather than directly the sources of the mixture; for example, variables that capture biases in the collection of the measurements). Note that \code{C2} must include row names and column names and should not include an intercept term. NA values are currently not supported.
#' @param refit_W A logical value indicating whether to re-estimate the input \code{W} under the TCA model.
#' @param refit_W.features A vector with the names of the features in \code{X} to consider when re-estimating \code{W} (i.e. when \code{refit_W == TRUE}). This is useful since oftentimes just a subset of the features in \code{X} will be informative for estimating \code{W}. If \code{refit_W.features == NULL} then the ReFACTor algorithm will be used for performing feature selection (see also \code{refit_W.sparsity, refit_W.sd_threshold}).
#' @param refit_W.sparsity A numeric value indicating the number of features to select using the ReFACTor algorithm when re-estimating \code{W} (activated only if \code{refit_W == TRUE} and \code{refit_W.features == NULL}). Note that \code{refit_W.sparsity} must be lower or equal to the number of features in \code{X}. For more information, see the argument \code{sparsity} in \link{refactor}.
#' @param refit_W.sd_threshold A numeric value indicating a standard deviation threshold to be used for excluding low-variance features in \code{X} (activated only if \code{refit_W == TRUE} and \code{refit_W.features == NULL}). For more information, see the argument \code{sd_threshold} in \link{refactor}.
#' @param tau A non-negative numeric value of the standard deviation of the measurement noise (i.e. the i.i.d. component of variation in the model). If \code{tau == NULL} then \code{tca} will estimate \code{tau}.
#' @param parallel A logical value indicating whether to use parallel computing (possible when using a multi-core machine).
#' @param num_cores A numeric value indicating the number of cores to use (activated only if \code{parallel == TRUE}). If \code{num_cores == NULL} then all available cores except for one will be used.
#' @param max_iters A numeric value indicating the maximal number of iterations to use in the optimization of the TCA model (\code{max_iters} iterations will be used as long as the optimization does not converge in earlier iterations).
#' @param log_file A path to an output log file. Note that if the file \code{log_file} already exists then logs will be appended to the end of the file. Set \code{log_file} to \code{NULL} to prevent output from being saved into a file; note that if \code{verbose == TRUE} then no output file will be generated regardless of the value of \code{log_file}.
#' @param debug A logical value indicating whether to set the logger to a more detailed debug level; please set \code{debug} to \code{TRUE} before reporting issues.
#' @param verbose A logical value indicating whether to print logs.
#'
#' @importFrom futile.logger flog.info
#' @importFrom pbapply pboptions
#'
#' @details The TCA model assumes that the hidden source-specific values are random variables. Formally, denote by \eqn{Z_{hj}^i} the source-specific value of observation \eqn{i} in feature \eqn{j} source \eqn{h}, the TCA model assumes: \deqn{Z_{hj}^i \sim N(\mu_{hj},\sigma_{hj}^2)} where \eqn{\mu_{hj},\sigma_{hj}} represent the mean and standard deviation that are specific to feature \eqn{j} source \eqn{h}. The model further assumes that the observed value of observation \eqn{i} in feature \eqn{j} is a mixture of \eqn{k} different sources: \deqn{X_{ji} =  \sum_{h=1}^k W_{ih}Z_{hj}^i + \epsilon_{ji}} where \eqn{W_{ih}} is the non-negative proportion of source \eqn{h} in the mixture of observation \eqn{i} such that \eqn{\sum_{h=1}^kW_{ih} = 1}, and \eqn{\epsilon_{ji} \sim N(0,\tau^2)} is an i.i.d. component of variation that models measurement noise. Note that the mixture proportions in \eqn{W} are, in general, unique for each individual, therefore each entry in the data matrix \eqn{X} is coming from a unique distribution (i.e. a different mean and a different variance).
#'
#' In cases where the true \code{W} is unknown, \code{tca} can be provided with initial estimates of \code{W} and then re-estimate \code{W} as part of the optimization procedure (see argument \code{refit_W}). These initial estimates should not be random but rather capture the information in \code{W} to some extent. When the argument \code{refit_W} is used, it is typically the case that only a subset of the features should be used for re-estimating \code{W}. Therefore, when re-estimating \code{W}, \code{tca} performs feature selection using the ReFACTor algorithm; alternatively, it can also be provided with a user-specified list of features to be used in the re-estimation (see argument \code{refit_W.features}).
#'
#' Factors that systematically affect the source-specific values \eqn{Z_{hj}^i} can be further considered (see argument \code{C1}). In that case, we assume: \deqn{Z_{hj}^i \sim N(\mu_{hj}+c^{(1)}_i \gamma_j^h,\sigma_{hj}^2)} where \eqn{c^{(1)}_i} is a row vector from \code{C1}, corresponding to the values of the \eqn{p_1} factors for observation \eqn{i}, and \eqn{\gamma_j^h} is a vector of \eqn{p_1} corresponding effect sizes.
#'
#' Factors that systematically affect the mixture values \eqn{X_{ji}}, such as variables that capture biases in the collection of the measurements, can also be considered (see argument \code{C2}). In that case, we assume: \deqn{X_{ji} \sim \sum_{h=1}^k W_{ih}Z_{hj}^i + c^{(2)}_i \delta_j + \epsilon_{ij}} where \eqn{c^{(2)}_i} is a row vector from \code{C2}, corresponding to the values of the \eqn{p_2} factors for observation \eqn{i}, and \eqn{\delta_j} is a vector of \eqn{p_2} corresponding effect sizes.
#'
#'
#' @return A list with the estimated parameters of the model. This list can be then used as the input to other functions such as \code{tcareg}.
#' \item{W}{An \code{n} by \code{k} matrix of weights. If \code{refit_W == TRUE} then this is the re-estimated \code{W}; otherwise this is the input \code{W}}
#' \item{mus_hat}{An \code{m} by \code{k} matrix of estimates for the mean of each source in each feature.}
#' \item{sigmas_hat}{An \code{m} by \code{k} matrix of estimates for the standard deviation of each source in each feature.}
#' \item{tau_hat}{An estimate of the standard deviation of the i.i.d. component of variation in \code{X}. If an input value was provided for \code{tau} (i.e. \code{tau != NULL}) then \code{tau_hat == tau}.}
#' \item{gammas_hat}{An \code{m} by \code{k*p1} matrix of the estimated effects of the \code{p1} factors in \code{C1} on each of the \code{m} features in \code{X}, where the first \code{p1} columns are the source-specific effects of the \code{p1} factors on the first source, the following \code{p1} columns are the source-specific effects on the second source and so on.}
#' \item{deltas_hat}{An \code{m} by \code{p2} matrix of the estimated effects of the \code{p2} factors in \code{C2} on the mixture values of each of the \code{m} features in \code{X}.}
#'
#' @examples
#' data <- test_data(100, 20, 3, 1, 1, 0.01)
#' tca.mdl <- tca(X = data$X, W = data$W, C1 = data$C1, C2 = data$C2)
#'
#' @section Note: The function \code{tca} may require a long running time when the input matrix \code{X} is very large; to alleviate this, it is strongly advised to use the \code{parallel} argument, given that a multi-core machine is available.
#'
#' @references Rahmani E, Schweiger R, Rhead B, Criswell LA, Barcellos LF, Eskin E, Rosset S, Sankararaman S, Halperin E. Cell-type-specific resolution epigenetics without the need for cell sorting or single-cell biology. Nature Communications 2018.
#' @references Rahmani E, Zaitlen N, Baran Y, Eng C, Hu D, Galanter J, Oh S, Burchard EG, Eskin E, Zou J, Halperin E. Sparse PCA corrects for cell type heterogeneity in epigenome-wide association studies. Nature Methods 2016.
#'
#' @export tca
tca <- function(X, W, C1 = NULL, C1.map = NULL, C2 = NULL, refit_W = FALSE, refit_W.features = NULL, refit_W.sparsity = 500, refit_W.sd_threshold = 0.02, tau = NULL, parallel = FALSE, num_cores = NULL, max_iters = 10, log_file = "TCA.log", debug = FALSE, verbose = TRUE){

  start_logger(log_file, debug, verbose)

  flog.info("Starting tca...")

  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)

  # progress bar options
  op <- if (verbose) pboptions(nout = config$nout, type = config[["type"]]) else pboptions(nout = config$nout, type = "none")

  X <- if (is.matrix(X)) X else as.matrix(X)
  W <- if (is.matrix(W)) W else as.matrix(W)
  C1 <- if (is.matrix(C1) | is.null(C1)) C1 else as.matrix(C1)
  C2 <- if (is.matrix(C2) | is.null(C2)) C2 else as.matrix(C2)

  flog.info("Validating input...")
  tca.validate_input(X, W, C1, C1.map, C2, refit_W, refit_W.features, refit_W.sparsity, refit_W.sd_threshold, tau, parallel, num_cores, max_iters, log_file, debug)
  if (is.null(C1)) C1 <- matrix(0, nrow=ncol(X), ncol=0)
  if (is.null(C2)) C2 <- matrix(0, nrow=ncol(X), ncol=0)

  C1.map <- if (is.null(C1.map)) matrix(1,ncol(C1),ncol(W)) else C1.map

  msg <- "Fitting the TCA model..."
  if (refit_W){
    flog.info("Starting re-estimation of W...")
    if (is.null(refit_W.features)){
      flog.info("Performing feature selection using refactor...")
      #ref <- refactor(X, ncol(W), sparsity = refit_W.sparsity, C = if (ncol(cbind(C1,C2))) cbind(C1,C2) else NULL, sd_threshold = refit_W.sd_threshold, rand_svd = config$rand_svd, log_file = NULL)
      ref <- refactor.run(X = X, k = ncol(W), sparsity = refit_W.sparsity, C = if (ncol(cbind(C1,C2))) cbind(C1,C2) else NULL, C.remove = FALSE, sd_threshold = refit_W.sd_threshold, num_comp = NULL, rand_svd = config$rand_svd, log_file = NULL, logger_on = FALSE, verbose = verbose)
      refit_W.features <- ref$ranked_list[1:refit_W.sparsity]
    }
    X_sub <- subset(X, subset = rownames(X) %in% refit_W.features)
    flog.info("Fitting the TCA model using the selected features for re-estimating W...")
    mdl0 <- tca.fit(X = t(X_sub), W = W, C1 = C1, C1.map = C1.map, C2 = C2, refit_W = TRUE, tau = tau, parallel = parallel, num_cores = num_cores, max_iters = max_iters)
    W <- mdl0[["W"]]
    msg <- "Fitting the TCA model given the updated W..."
  }
  X <- t(X)
  flog.info(msg)
  mdl <- tca.fit(X = X, W = W, C1 = C1, C1.map = C1.map, C2 = C2, refit_W = FALSE, tau = tau, parallel = parallel, num_cores = num_cores, max_iters = max_iters)
  W <- mdl[["W"]]
  mus_hat <- mdl[["mus_hat"]]
  sigmas_hat <- mdl[["sigmas_hat"]]
  tau_hat <- mdl[["tau_hat"]]
  deltas_hat <- mdl[["deltas_hat"]]
  gammas_hat <- mdl[["gammas_hat"]]
  flog.info("Finished tca.")
  return( list("W" = W, "mus_hat" = mus_hat, "sigmas_hat" = sigmas_hat, "tau_hat" = tau_hat, "deltas_hat" = deltas_hat, "gammas_hat" = gammas_hat, "C1" = C1, "C2" = C2) )
}


#' @title Subsetting features from a TCA model
#'
#' @description Extracts from a fitted TCA model (i.e. a value returned by the function \code{tca}) a subset of the features.
#'
#' @param tca.mdl The value returned by applying the function \code{tca} to some data matrix\code{X}.
#' @param features A vector with the identifiers of the features to extract (as they appear in the rows of \code{X}).
#' @param log_file A path to an output log file. Note that if the file \code{log_file} already exists then logs will be appended to the end of the file. Set \code{log_file} to \code{NULL} to prevent output from being saved into a file; note that if \code{verbose == TRUE} then no output file will be generated regardless of the value of \code{log_file}.
#' @param debug A logical value indicating whether to set the logger to a more detailed debug level; please set \code{debug} to \code{TRUE} before reporting issues.
#' @param verbose A logical value indicating whether to print logs.
#'
#' @importFrom futile.logger flog.info
#'
#' @details This function allows to extract a subset of the features from a fitted TCA model (i.e. from a value returned by the function \code{tca}). This allows, for example, to extract and then perform post-hoc tests on only a small set of candidate features (e.g., using the function \code{tcareg}), without the need to run \code{tca} again for fitting the model to the candidate features.
#'
#' @return A list with the estimated parameters of the model for the given set of features.
#' \item{W}{Equals to \code{tca.mdl$W}}
#' \item{mus_hat}{A \code{q} by \code{k} matrix which is a subset of the matrix \code{tca.mdl$mus_hat}, where \code{q} is the number of features in the argument \code{features}.}
#' \item{sigmas_hat}{A \code{q} by \code{k} matrix which is a subset of the matrix \code{tca.mdl$sigmas_hat}, where \code{q} is the number of features in the argument \code{features}.}
#' \item{tau_hat}{Equals to \code{tca.mdl$tau_hat}}
#' \item{gammas_hat}{A \code{q} by \code{k*p1}  matrix which is a subset of the matrix \code{tca.mdl$gammas_hat}, where \code{q} is the number of features in the argument \code{features}.}
#' \item{deltas_hat}{A \code{q} by \code{p2} matrix which is a subset of the matrix \code{tca.mdl$deltas_hat}, where \code{q} is the number of features in the argument \code{features}.}
#'
#' @examples
#' data <- test_data(50, 20, 3, 0, 0, 0.01)
#' tca.mdl <- tca(X = data$X, W = data$W)
#' tca.mdl.subset <- tcasub(tca.mdl, rownames(data$X)[1:10])
#' y <- matrix(rexp(50, rate=.1), ncol=1)
#' # run tcareg test with an outcome y:
#' res <- tcareg(data$X[1:10,], tca.mdl.subset, y, test = "joint", save_results = FALSE)
#'
#' @export tcasub
tcasub <- function(tca.mdl, features, log_file = "TCA.log", debug = FALSE, verbose = TRUE){

  start_logger(log_file, debug, verbose)
  flog.info("Starting tcasub...")
  tcasub.validate_input(tca.mdl, features, log_file, debug)
  flog.info("Finished tcasub.")
  if (length(features) == 1){
    mus_hat <- t(as.matrix(tca.mdl[["mus_hat"]][features,]))
    sigmas_hat <- t(as.matrix(tca.mdl[["sigmas_hat"]][features,]))
    deltas_hat <- t(as.matrix(tca.mdl[["deltas_hat"]][features,]))
    gammas_hat <- t(as.matrix(tca.mdl[["gammas_hat"]][features,]))
  }else{
    mus_hat <- as.matrix(tca.mdl[["mus_hat"]][features,])
    sigmas_hat <- as.matrix(tca.mdl[["sigmas_hat"]][features,])
    deltas_hat <- as.matrix(tca.mdl[["deltas_hat"]][features,])
    gammas_hat <- as.matrix(tca.mdl[["gammas_hat"]][features,])
  }
  return( list("W" = tca.mdl[["W"]], "mus_hat" = mus_hat, "sigmas_hat" = sigmas_hat, "tau_hat" = tca.mdl[["tau_hat"]], "deltas_hat" = deltas_hat, "gammas_hat" = gammas_hat, "C1" = as.matrix(tca.mdl[["C1"]]), "C2" = as.matrix(tca.mdl[["C2"]])) )

}


#' @title Extracting hidden 3D signals from 2D input
#'
#' @description Estimates 3-dimensional signals (features by observations by sources) from input of mixtures (features by observations), under the assumption of the TCA model that each observation is a mixture of unique source-specific values (in each feature in the data). For example, in the context of  tissue-level bulk DNA methylation data coming from a mixture of cell types (i.e. the input is methylation sites by individuals), \code{tensor} allows to estimate a tensor of cell-type-specific levels for each individual in each methylation site (i.e. a tensor of methylation sites by individuals by cell types).
#'
#' @param X An \code{m} by \code{n} matrix of measurements of \code{m} features for \code{n} observations. Each column in \code{X} is assumed to be a mixture of \code{k} different sources. Note that \code{X} must include row names and column names and that NA values are currently not supported. \code{X} should not include features that are constant across all observations.
#' @param tca.mdl The value returned by applying the function \code{tca} to \code{X}.
#' @param parallel A logical value indicating whether to use parallel computing (possible when using a multi-core machine).
#' @param num_cores A numeric value indicating the number of cores to use (activated only if \code{parallel == TRUE}). If \code{num_cores == NULL} then all available cores except for one will be used.
#' @param log_file A path to an output log file. Note that if the file \code{log_file} already exists then logs will be appended to the end of the file. Set \code{log_file} to \code{NULL} to prevent output from being saved into a file; note that if \code{verbose == TRUE} then no output file will be generated regardless of the value of \code{log_file}.
#' @param debug A logical value indicating whether to set the logger to a more detailed debug level; please set \code{debug} to \code{TRUE} before reporting issues.
#' @param verbose A logical value indicating whether to print logs.
#'
#' @details See \link{tca} for notations and details about the TCA model. Given estimates of the parameters of the model (given by \link{tca}), \code{tensor} uses the conditional distribution \eqn{Z_{hj}^i|X_{ji}} for estimating the \eqn{k} source-specific levels of each observation \eqn{i} in each feature \eqn{j}.
#
#' @return A list with the estimated source-specific values. The first element in the list is an \code{m} by \code{n} matrix (features by observations) corresponding to the estimated values coming from the first source, the second element in the list is another \code{m} by \code{n} matrix (features by observations) corresponding to the estimated values coming from the second source and so on.
#'
#' @examples
#' data <- test_data(50, 20, 3, 2, 2, 0.01)
#' tca.mdl <- tca(X = data$X, W = data$W, C1 = data$C1, C2 = data$C2)
#' Z_hat <- tensor(data$X, tca.mdl)
#'
#' @references Rahmani E, Schweiger R, Rhead B, Criswell LA, Barcellos LF, Eskin E, Rosset S, Sankararaman S, Halperin E. Cell-type-specific resolution epigenetics without the need for cell sorting or single-cell biology. Nature Communications 2018.
#'
#' @export tensor
tensor <- function(X, tca.mdl, parallel = FALSE, num_cores = NULL, log_file = "TCA.log", debug = FALSE, verbose = TRUE){

  start_logger(log_file, debug, verbose)

  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)

  # progress bar options
  op <- if (verbose) pboptions(nout = config$nout, type = config[["type"]]) else pboptions(nout = config$nout, type = "none")

  flog.info("Starting tensor for estimating Z...")

  X <- if (is.matrix(X)) X else as.matrix(X)

  W <- tca.mdl[["W"]]
  mus_hat <- tca.mdl[["mus_hat"]]
  sigmas_hat <- tca.mdl[["sigmas_hat"]]
  tau_hat <- tca.mdl[["tau_hat"]]
  deltas_hat <- tca.mdl[["deltas_hat"]]
  gammas_hat <- tca.mdl[["gammas_hat"]]
  C1 <- tca.mdl[["C1"]]
  C2 <- tca.mdl[["C2"]]

  return( estimate_Z(t(X), W, mus_hat, sigmas_hat, tau_hat, C2, deltas_hat, C1, gammas_hat, parallel, num_cores) )
}



#' @title Fitting a TCA regression model
#'
#' @description TCA regression allows to test for several types of statistical relations between source-specific values and an outcome of interest. For example, in the context of tissue-level bulk DNA methylation data coming from a mixture of cell types (i.e. the input is methylation sites by individuals), \code{tcareg} allows to test for cell-type-specific effects of methylation on an outcome of interest.
#'
#' @param X An \code{m} by \code{n} matrix of measurements of \code{m} features for \code{n} observations. Each column in \code{X} is assumed to be a mixture of \code{k} different sources. Note that \code{X} must include row names and column names and that NA values are currently not supported. \code{X} should not include features that are constant across all observations.
#' @param tca.mdl The value returned by applying the function \code{tca} to \code{X}.
#' @param y An \code{n} by 1 matrix of an outcome of interest for each of the \code{n} observations in \code{X}. Note that \code{y} must include row names and column names and that NA values are currently not supported.
#' @param C3 An \code{n} by \code{p3} design matrix of covariates that may affect \code{y}. Note that \code{C3} must include row names and column names and should not include an intercept term. NA values are currently not supported.
#' @param test A character vector with the type of test to perform on each of the features in \code{X}; one of the following options: \code{'marginal'}, \code{'marginal_conditional'}, \code{'joint'}, \code{'single_effect'}, or \code{'custom'}. Setting \code{'marginal'} or \code{'marginal_conditional'} corresponds to testing each feature in \code{X} for a statistical relation between \code{y} and each of the \code{k} sources separately; for any particular source under test, the \code{marginal_conditional} option further accounts for possible effects of the rest of the \code{k-1} sources (\code{'marginal'} will therefore tend to be more powerful in discovering truly related features, but at the same time more prone to falsely tagging the correct related sources if sources are highly correlated). Setting \code{'joint'} or \code{'single_effect'} corresponds to testing each feature for an overall statistical relation with \code{y}, while modeling source-specific effects; the latter option further assumes that the source-specific effects are the same within each feature (\code{'single_effect'} means only one degree of freedom and will therefore be more powerful when the assumption of a single effect within a feature holds). Finally, \code{'custom'} corresponds to testing each feature in \code{X} for a statistical relation with \code{y} under a user-specified model (alternative model) with respect to a null model (null model); for example, for testing for relation of the combined (potentially different) effects of features 1 and 2 while accounting for the (potentially different) effects of 3 and 4, set the null model to be sources 3, 4 and the alternative model to be sources 1, 2, 3, 4. Indicating that \code{null_model} assumes no effect for any of the sources can be done by setting it to \code{NULL}.
#' @param null_model A vector with a subset of the names of the sources in \code{tca.mdl$W} to be used as a null model (activated only if \code{test == 'custom'}). Note that the null model must be nested within the alternative model; set \code{null_model} to be \code{NULL} for indicating no effect for any of the sources under the null model.
#' @param alternative_model A vector with a subset (or all) of the names of the sources in \code{tca.mdl$W} to be used as an alternative model (activated only if \code{test == 'custom'}).
#' @param save_results A logical value indicating whether to save the returned results in a file. If \code{TRUE} and \code{test == 'marginal'} or \code{test == 'marginal_conditional'} then \code{k} files will be saved (one for the results of each source).
#' @param output Prefix for output files (activated only if \code{save_results == TRUE}).
#' @param sort_results A logical value indicating whether to sort the results by their p-value (i.e. features with lower p-value will appear first in the results).
#' @param parallel A logical value indicating whether to use parallel computing (possible when using a multi-core machine).
#' @param num_cores A numeric value indicating the number of cores to use (activated only if \code{parallel == TRUE}). If \code{num_cores == NULL} then all available cores except for one will be used.
#' @param log_file A path to an output log file. Note that if the file \code{log_file} already exists then logs will be appended to the end of the file. Set \code{log_file} to \code{NULL} to prevent output from being saved into a file; note that if \code{verbose == TRUE} then no output file will be generated regardless of the value of \code{log_file}.
#' @param features_metadata A path to a csv file containing metadata about the features in \code{X} that will be added to the output files (activated only if \code{save_results == TRUE}). Each row in the metadata file should correspond to one feature (with the row name being the feature identifier, as it appears in the rows of \code{X}) and each column should correspond to one metadata descriptor (with an appropriate column name). Features that do not exist in \code{X} will be ignored and features in \code{X} with missing metadata information will show missing values.
#' @param debug A logical value indicating whether to set the logger to a more detailed debug level; please set \code{debug} to \code{TRUE} before reporting issues.
#' @param verbose A logical value indicating whether to print logs.
#'
#' @details TCA models \eqn{Z_{hj}^i} as the source-specific value of observation \eqn{i} in feature \eqn{j} coming from source \eqn{h} (see \link{tca} for more details). A TCA regression model tests an outcome \eqn{Y} for a linear statistical relationwith the source-specific values of a feature \eqn{j} by assuming: \deqn{Y_i = \sum_{h=1}^k \beta_{hj} Z_{hj}^i + e_i} where \eqn{e_i \sim N(0,\phi^2)}. In practice, \code{tcareg} fits this model using the conditional distribution \eqn{Y|X}, which, effectively, integrates over the latent \eqn{Z_{hj}^i} parameters. Statistical significance is then calculated using a likelihood ratio test (LRT). Note that the null and alternative models will be set automatically, except when \code{test == 'custom'}, in which case they will be set according to the user-specified null and alternative hypotheses.
#'
#' Under the TCA regression model, several statistical tests can be performed by setting the argument \code{test} according to one of the following options.
#'
#' 1. If \code{test == 'marginal'}, \code{tcareg} will perform the following for each source \eqn{l}. For each feature \eqn{j}, \eqn{\beta_{lj}} will be estimated and tested for a non-zero effect, while assuming \eqn{\beta_{hj}=0} for all other sources \eqn{h\neq l}.
#'
#' 2. If \code{test == 'marginal_conditional'}, \code{tcareg} will perform the following for each source \eqn{l}. For each feature \eqn{j}, \eqn{\beta_{lj}} will be estimated and tested for a non-zero effect, while also estimating the effect sizes \eqn{\beta_{hj}} for all other sources \eqn{h\neq l}.
#'
#' 3. If \code{test == 'joint'}, \code{tcareg} will estimate for each feature \eqn{j} the effect sizes of all \eqn{k} sources \eqn{\beta_{1j},…,\beta_{kj}} and then test the set of \eqn{k} estimates of each feature \code{j} for a joint effect.
#'
#' 4. If \code{test == 'single_effect'}, \code{tcareg} will estimate for each feature \eqn{j} the effect sizes of all \eqn{k} sources \eqn{\beta_{1j},…,\beta_{kj}}, under the assumption that \eqn{\beta_{1j} = … = \beta_{kj}}, and then test the set of \eqn{k} estimates of each feature \code{j} for a joint effect.
#'
#' 5. If \code{test == 'custom'}, \code{tcareg} will estimate for each feature \eqn{j} the effect sizes of a predefined set of sources (defined by a user-specified alternative model) and then test their estimates for a joint effect, while accounting for a nested predefined set of sources (defined by a user-specified null model).
#'
#' @return A list with the results of applying the TCA regression model to each of the features in the data. If \code{test == 'marginal'} or \code{test == 'marginal_conditional'} then a list of \code{k} such lists of results are returned, one for the results of each source.
#' \item{phi}{An estimate of the standard deviation of the i.i.d. component of variation in the TCA regression model.}
#' \item{beta}{A matrix of effect size estimates for the source-specific effects, such that each row corresponds to the estimated effect sizes of one feature. The number of columns corresponds to the number of estimated effects (e.g., if \code{test} is set to \code{marginal} then \code{beta} will include a single column, if \code{test} is set to \code{joint} then \code{beta} will include \code{k} columns and so on).}
#' \item{intercept}{An \code{m} by \code{1} matrix of estimates for the intercept of each feature.}
#' \item{alpha}{An \code{m} by \code{p3} matrix of effect size estimates for the \code{p3} covariates in \code{C3}, such that each row corresponds to the estimated effect sizes of one feature.}
#' \item{null_ll}{An \code{m} by \code{1} matrix of the log-likelihood of the model under the null hypothesis.}
#' \item{alternative_ll}{An \code{m} by \code{1} matrix of the log-likelihood of the model under the alternative hypothesis.}
#' \item{stats}{An \code{m} by \code{1} matrix of the LRT statistic for each feature in the data.}
#' \item{df}{The degrees of freedom for deriving p-values using LRT.}
#' \item{pvals}{An \code{m} by \code{1} matrix of the p-value for each feature in the data.}
#' \item{qvals}{An \code{m} by \code{1} matrix of the q-value (FDR-adjusted p-values) for each feature in the data.}
#'
#' @section Note: The function \code{tcareg} may require a long running time when the input matrix \code{X} is very large; to alleviate this, it is strongly advised to use the \code{parallel} argument, given that a multi-core machine is available.
#'
#' @examples
#' n <- 50
#' m <- 10
#' k <- 3
#' p1 <- 1
#' p2 <- 1
#' data <- test_data(n, m, k, p1, p2, 0.01)
#' tca.mdl <- tca(X = data$X, W = data$W, C1 = data$C1, C2 = data$C2)
#' y <- matrix(rexp(n, rate=.1), ncol=1)
#' # joint test:
#' res1 <- tcareg(data$X, tca.mdl, y, test = "joint", save_results = FALSE)
#' # custom test, testing for a joint effect of sources 1,2 while accounting for source 3
#' res2 <- tcareg(data$X, tca.mdl, y, test = "custom", null_model = c("3"),
#' alternative_model = c("1","2","3"), save_results = FALSE)
#' # custom test, testing for a joint effect of sources 1,2 assuming no effects under the null
#' res3 <- tcareg(data$X, tca.mdl, y, test = "custom", null_model = NULL,
#' alternative_model = c("1","2"), save_results = FALSE)
#'
#' @references Rahmani E, Schweiger R, Rhead B, Criswell LA, Barcellos LF, Eskin E, Rosset S, Sankararaman S, Halperin E. Cell-type-specific resolution epigenetics without the need for cell sorting or single-cell biology. Nature Communications 2018.
#'
#' @export tcareg
#'
tcareg <- function(X, tca.mdl, y, C3 = NULL, test = "marginal", null_model = NULL, alternative_model = NULL, save_results = TRUE, output = "TCA", sort_results = TRUE, parallel = FALSE, num_cores = NULL, log_file = "TCA.log", features_metadata = NULL, debug = FALSE, verbose = TRUE){

  start_logger(log_file, debug, verbose)

  flog.info("Starting tcareg...")

  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)

  # set options for the progress bars
  op <- if (verbose) pboptions(nout = config$nout, type = config[["type"]]) else pboptions(nout = config$nout, type = "none")

  X <- if (is.matrix(X)) X else as.matrix(X)
  y <- if (is.matrix(y)) y else as.matrix(y)
  C3 <- if (is.matrix(C3) | is.null(C3)) C3 else as.matrix(C3)

  W <- tca.mdl[['W']]
  flog.info("Validating input...")
  tcareg.validate_input(X, W, y, C3, test, null_model, alternative_model, save_results, output, sort_results, parallel, num_cores, log_file, features_metadata, debug)

  flog.info("Preparing data...")
  X <- t(X)
  mus_hat <- tca.mdl[["mus_hat"]]
  sigmas_hat <- tca.mdl[["sigmas_hat"]]
  tau_hat <- tca.mdl[["tau_hat"]]
  C2 <- tca.mdl[["C2"]]
  deltas_hat <- tca.mdl[["deltas_hat"]]
  C1 <- tca.mdl[["C1"]]
  gammas_hat <- tca.mdl[["gammas_hat"]]
  C3 <- if(!is.null(C3)) C3 else matrix(0, nrow=nrow(X), ncol=0)
  #null_model <- match(null_model,colnames(W))
  null_model <- if (is.null(null_model)) NULL else match(null_model,colnames(W))
  alternative_model <- match(alternative_model,colnames(W))
  k <- ncol(W)

  flog.info("Running TCA regression...")
  res <- tcareg.fit(X, y, W, mus_hat, sigmas_hat, C2, deltas_hat, C1, gammas_hat, tau_hat, C3, test, null_model, alternative_model, parallel, num_cores)

  if (save_results){
    flog.info("Saving results...")
    save_association_results(res, output, test, alternative_model, colnames(X), colnames(W), colnames(C3), sort_results, features_metadata)
  }

  flog.info("Finished tcareg.")
  return(res)
}



#' @title Sparse principal component analysis using ReFACTor
#'
#' @description Performs unsupervised feature selection followed by principal component analysis (PCA) under a row-sparse model using the ReFACTor algorithm. For example, in the context of tissue-level bulk DNA methylation data coming from a mixture of cell types (i.e. the input is methylation sites by individuals), \code{refactor} allows to capture the variation in cell-type composition, which was shown to be a dominant sparse signal in methylation data.
#'
#' @param X An \code{m} by \code{n} matrix of measurements of \code{m} features for \code{n} observations. Each column in \code{X} is assumed to be a mixture of \code{k} different sources. Note that \code{X} must include row names and column names and that NA values are currently not supported. \code{X} should not include features that are constant across all observations.
#' @param k A numeric value indicating the dimension of the signal in the data (i.e. the number of sources).
#' @param sparsity A numeric value indicating the sparsity of the signal in the data (the number of signal rows).
#' @param C An \code{n} by \code{p} design matrix of covariates that will be accounted for in the feature selection step. Note that \code{C} must include row names and column names and that NA values are currently not supported; ; set \code{C} to be \code{NULL} if there are no such covariates.
#' @param C.remove A logical value indicating whether the covariates in X should be accounted for not only in the feature selection step, but also in the final calculation of the principal components (i.e. if \code{C.remove == TRUE} then the selected features will be adjusted for the covariates in \code{C} prior to calculating principal components). Note that setting \code{C.remove} to be \code{TRUE} is desired when ReFACTor is intended to be used for correction in downstream analysis, whereas setting \code{C.remove} to be \code{FALSE} is desired when ReFACTor is merely used for capturing the sparse signals in the data (i.e. regardless of correction).
#' @param sd_threshold A numeric value indicating a standard deviation threshold to be used for excluding low-variance features in \code{X} (i.e. features with standard deviation lower than \code{sd_threshold} will be excluded). Set \code{sd_threshold} to be \code{NULL} for turning off this filter. Note that removing features with very low variability tends to improve speed and performance.
#' @param num_comp A numeric value indicating the number of ReFACTor components to return.
#' @param rand_svd A logical value indicating whether to use random svd for estimating the low-rank structure of the data in the first step of the algorithm; random svd can result in a substantial speedup for large data.
#' @param log_file A path to an output log file. Note that if the file \code{log_file} already exists then logs will be appended to the end of the file. Set \code{log_file} to \code{NULL} to prevent output from being saved into a file; note that if \code{verbose == TRUE} then no output file will be generated regardless of the value of \code{log_file}.
#' @param debug A logical value indicating whether to set the logger to a more detailed debug level; please set \code{debug} to \code{TRUE} before reporting issues.
#' @param verbose A logical value indicating whether to print logs.
#'
#' @importFrom rsvd rpca
#' @importFrom gmodels fast.prcomp
#'
#' @details ReFACTor is a two-step algorithm for sparse principal component analysis (PCA) under a row-sparse model. The algorithm performs an unsupervised feature selection by ranking the features based on their correlation with their values under a low-rank representation of the data, followed by a calculation of principal components using the top ranking features (ReFACTor components).
#'
#' Note that ReFACTor is tuned towards capturing sparse signals of the dominant sources of variation in the data. Therefore, in the presence of other potentially dominant factors in the data (i.e. beyond the variation of interest), these factors should be accounted for by including them as covariates (see argument \code{C}). In cases where the ReFACTor components are designated to be used as covariates in a downstream analysis alongside the covariates in \code{C} (e.g., in a standard regression analysis or in a TCA regression), it is advised to set the argument \code{C.remove} to be \code{TRUE}. This will adjust the selected features for the information in \code{C} prior to the calculation of the ReFACTor components, which will therefore capture only signals that is not present in \code{C} (and as a result may benefit the downstream analysis by potentially capturing more signals beyond the information in \code{C}).
#'
#' @return A list with the estimated components of the ReFACTor model.
#' \item{scores}{An \code{n} by \code{num_comp} matrix of the ReFACTor components (the projection scores).}
#' \item{coeffs}{A \code{sparsity} by \code{num_comp} matrix of the coefficients of the ReFACTor components (the projection loadings).}
#' \item{ranked_list}{A vector with the features in the data, ranked by their scores in the feature selection step of the algorithm; the top scoring features (set according to the argument \code{sparsity}) are used for calculating the ReFACTor components. Note that features that were excluded according to \code{sd_threshold} will not appear in this \code{ranked_list}.}
#'
#' @section Note: For very large input matrices it is advised to use random svd for speeding up the feature selection step (see argument \code{rand_svd}).
#'
#' @examples
#' data <- test_data(100, 200, 3, 0, 0, 0.01)
#' ref <- refactor(data$X, k = 3, sparsity = 50)
#'
#' @references Rahmani E, Zaitlen N, Baran Y, Eng C, Hu D, Galanter J, Oh S, Burchard EG, Eskin E, Zou J, Halperin E. Sparse PCA corrects for cell type heterogeneity in epigenome-wide association studies. Nature Methods 2016.
#' @references Rahmani E, Zaitlen N, Baran Y, Eng C, Hu D, Galanter J, Oh S, Burchard EG, Eskin E, Zou J, Halperin E. Correcting for cell-type heterogeneity in DNA methylation: a comprehensive evaluation. Nature Methods 2017.
#'
#' @export refactor
refactor <- function(X, k, sparsity = 500, C = NULL, C.remove = FALSE, sd_threshold = 0.02, num_comp = NULL, rand_svd = FALSE, log_file = "TCA.log", debug = FALSE, verbose = TRUE){

  return (refactor.run(X = X, k = k, sparsity = sparsity, C = C, C.remove = C.remove, sd_threshold = sd_threshold, num_comp = num_comp, rand_svd = rand_svd, log_file = log_file, debug = debug, logger_on = TRUE, verbose))

}



#' @title Generate test data
#'
#' @description Generates simple test data following the TCA model.
#'
#' @param n The number of observations to simulate.
#' @param m The number of features to simulate.
#' @param k The number of sources to simulate.
#' @param p1 The number of covariates that affect the source-specific values to simulate.
#' @param p2 The number of covariates that affect the mixture values to simulate.
#' @param tau The variance of the i.i.d. component of variation to add on top of the simulated mixture values.
#' @param log_file A path to an output log file. Note that if the file \code{log_file} already exists then logs will be appended to the end of the file. Set \code{log_file} to \code{NULL} to prevent output from being saved into a file; note that if \code{verbose == TRUE} then no output file will be generated regardless of the value of \code{log_file}.
#' @param verbose A logical value indicating whether to print logs.
#'
#' @importFrom futile.logger flog.info
#' @importFrom pracma Reshape
#'
#' @details See \link{tca} for details about the TCA model.
#'
#' @return A list with the simulated data and parameters.
#' \item{X}{An \code{m} by \code{n} matrix of simulated data with \code{m} features for \code{n} observations.}
#' \item{Z}{A list with the simulated source-specific values, where the first element in the list is an \code{m} by \code{n} matrix (features by observations) corresponding to the values coming from the first source, the second element in the list is another \code{m} by \code{n} matrix (features by observations) corresponding to the values coming from the second source and so on.}
#' \item{W}{An \code{n} by \code{k} matrix of simulated weights - the weights of the \code{k} sources for each of the \code{n} mixtures (observations).}
#' \item{mus}{An \code{m} by \code{k} matrix of the mean of each of the \code{m} features for each of the \code{k} sources.}
#' \item{sigmas}{An \code{m} by \code{k} matrix of the standard variation of each of the \code{m} features for each of the \code{k} sources.}
#' \item{C1}{ An \code{n} by \code{p1} design matrix of simulated covariates that affect the hidden source-specific values.}
#' \item{C2}{ An \code{n} by \code{p2} design matrix of simulated covariates that affect the mixture.}
#' \item{gammas}{An \code{m} by \code{k*p1} matrix of the effects of the \code{p1} factors in \code{C1} on each of the \code{m} features in \code{X}, where the first \code{p1} columns are the source-specific effects of the \code{p1} factors on the first source, the following \code{p1} columns are the source-specific effects on the second source and so on.}
#' \item{deltas}{An \code{m} by \code{p2} matrix of the effects of the \code{p2} factors in \code{C2} on the mixture values of each of the \code{m} features in \code{X}.}
#'
#' @examples
#' data <- test_data(100, 50, 3, 2, 2, 0.01)
#'
#' @export test_data
#'
test_data <- function(n, m, k, p1, p2, tau, log_file = "TCA.log", verbose = TRUE){

  start_logger(log_file, FALSE, verbose)

  flog.info("Starting test_data...")

  mus <- Reshape(runif(m*k, min = 0, max = 1),m,k)
  sigmas <- Reshape(runif(m*k, min = 0, max = 0.2),m,k)
  W <- Reshape(runif(n*k, min = 0, max = 1),n,k)
  W <- W/t(repmat(rowSums(W),k,1)) # normalize so the proportions are in the range [0,1] and sum up to 1 for each observation.
  C1 <- Reshape(rnorm(n*p1, mean = 0, sd = 1), n, p1)
  C2 <- Reshape(rnorm(n*p2, mean = 0, sd = 1), n, p2)
  gammas <- Reshape(rnorm(m*p1*k, mean = 0, sd = 0.1), m, p1*k)
  deltas <- Reshape(rnorm(m*p2, mean = 0, sd = 0.1), m, p2)
  Z <- list()
  for (h in 1:k){
    Z[[h]] <- matrix(0,n,m)
  }
  for (i in 1:n){
    for (h in 1:k){
      Z[[h]][i,] <- rnorm(m, mean = mus[,h], sd = sigmas[,h])
      if (p1){
        Z[[h]][i,] <- Z[[h]][i,] + C1[i,] %*% t(gammas[,((h-1)*p1+1):(p1*h)])
      }
    }
  }
  X = matrix(0,n,m)
  for (h in 1:k){
    X <- X + Z[[h]]*t(repmat(W[,h],m,1))
  }
  if (p2){
    X <- X + C2 %*% t(deltas)
  }
  # add i.i.d. component of variation
  X <- X + Reshape(rnorm(n*m, mean = 0, sd = tau), n, m)
  X<-t(X)
  for (h in 1:k){
    Z[[h]] <- t(Z[[h]])
    rownames(Z[[h]]) <- 1:m
    colnames(Z[[h]]) <- 1:n
  }
  colnames(X) <- 1:n; rownames(X) <- 1:m
  colnames(W) <- 1:k; rownames(W) <- 1:n
  colnames(mus) <- colnames(W); colnames(sigmas) <- colnames(W)
  rownames(mus) <- 1:m; rownames(sigmas) <- 1:m
  if (p1){
    rownames(C1) <- 1:n
    colnames(C1) <- 1:dim(C1)[2]
    rownames(gammas) <- 1:m
    colnames(gammas) <- 1:dim(gammas)[2]
  }
  if (p2){
    rownames(C2) <- 1:n
    colnames(C2) <- 1:dim(C2)[2]
    rownames(deltas) <- 1:m
    colnames(deltas) <- 1:dim(deltas)[2]
  }

  flog.info("Finished test_data.")
  return( list("X" = X, "Z" = Z, "W" = W, "mus" = mus, "sigmas" = sigmas, "C1" = C1, "C2" = C2, "gammas" = gammas, "deltas" = deltas) )
}
