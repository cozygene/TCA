#' @importFrom rsvd rpca
#' @importFrom gmodels fast.prcomp

refactor.run <- function(X, k, sparsity, C, C.remove, sd_threshold, num_comp, rand_svd, log_file, debug, logger_on = TRUE, verbose){

  if (logger_on) start_logger(log_file, debug, verbose)

  flog.info("Starting refactor...")

  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)

  X <- if (is.matrix(X)) X else as.matrix(X)
  C <- if (is.matrix(C) | is.null(C)) C else as.matrix(C)

  if (!is.null(sd_threshold)){
    flog.debug("Excluding features with low variance (sd_threshold < %s)...", sd_threshold)
    sds <- apply(t(X), 2, sd)
    m <- length(sds)
    include <- which(sds >= sd_threshold)
    X <- X[include,]
    flog.debug("%s features were excluded due to low variance", m-length(include))
  }

  num_comp <- if(is.null(num_comp)) k else num_comp
  sparsity <- if(sparsity > nrow(X)) nrow(X) else sparsity

  X_adj <- X
  if (!is.null(C)){
    flog.info("Adjusting the data for covariates...")
    X_adj <- X - tcrossprod(tcrossprod(X,tcrossprod(matrix.inverse(crossprod(C,C)),C)),C)
  }

  flog.info("Running PCA on X using rand_svd == %s...", rand_svd)
  pcs <- if (rand_svd) rpca(t(X_adj), scale = TRUE, k=num_comp) else fast.prcomp(scale(t(X_adj)))
  coeff <- pcs$rotation
  score <- pcs$x

  flog.info("Computing a low rank approximation of X...")
  X_hat <- score[,1:k]%*%t(coeff[,1:k])
  An <- scale(t(X_adj),center=T,scale=F)
  Bn <- scale(X_hat,center=T,scale=F)
  An <- t(t(An)*(1/sqrt(apply(An^2,2,sum))))
  Bn <- t(t(Bn)*(1/sqrt(apply(Bn^2,2,sum))))

  flog.info("Calculating the distance of each feature from its low rank approximation...")
  distances <- apply((An-Bn)^2,2,sum)^0.5
  dsort <- sort(distances,index.return=T)
  ranked_list <- rownames(X)[dsort$ix]

  flog.info("Computing the ReFACTor components based on the top %s features with lowest distances...", sparsity)
  features <- ranked_list[1:sparsity]

  if (C.remove){
    refactor_pcs <- fast.prcomp(scale(t(X_adj[features,])))
  }else{
    refactor_pcs <- fast.prcomp(scale(t(X[features,])))
  }

  flog.info("Finished refactor.")

  return( list("scores" = refactor_pcs$x[,1:num_comp], "coeffs" = refactor_pcs$rotation[,1:num_comp], "ranked_list" = ranked_list) )
}
