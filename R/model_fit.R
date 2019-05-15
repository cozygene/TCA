
# Fits the TCA model
tca.fit <- function(X, W, C1, C2, refit_W, parallel, num_cores, max_iters){

  flog.debug("Starting function 'tca.fit'...")

  #config <- config::get()
  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)

  n <- nrow(X)
  m <- ncol(X)
  k <- ncol(W)
  p1 <- ncol(C1)
  p2 <- ncol(C2)

  # init estimates
  init <- init_means_vars(colnames(C1), colnames(C2), colnames(X), colnames(W))
  mus_hat <- init[["mus_hat"]]
  sigmas_hat <- init[["sigmas_hat"]]
  deltas_hat <- init[["deltas_hat"]]
  gammas_hat <- init[["gammas_hat"]]
  tau_hat <- init[["tau_hat"]]

  const <- -n*log(2*pi)
  ll_prev <- -Inf
  for (iter in 1:max_iters){
    if (refit_W) flog.info("Iteration %s out of %s external iterations (fitting all parameters including W)", iter, max_iters)

    flog.info("Fitting means and variances...")
    res <- tca.fit_means_vars(X, W, mus_hat, sigmas_hat, tau_hat, C2, deltas_hat, C1, gammas_hat, max_iters, parallel, num_cores)
    mus_hat <- res[["mus_hat"]]
    sigmas_hat <- res[["sigmas_hat"]]
    deltas_hat <- res[["deltas_hat"]]
    gammas_hat <- res[["gammas_hat"]]
    tau_hat <- res[["tau_hat"]]

    if (!refit_W) break
    flog.info("Fitting W...")
    W <- tca.fit_W(X, W, mus_hat, sigmas_hat, tau_hat, C2, deltas_hat, C1, gammas_hat, max_iters, parallel, num_cores)

    C1_ <- calc_C1_W_interactions(W,C1)
    flog.debug("Test for convergence.")
    U <- (tcrossprod(W,mus_hat) + tcrossprod(C2,deltas_hat) + tcrossprod(C1_,gammas_hat) - X)**2
    W_squared <- W**2
    ll_new <- -minus_log_likelihood_tau(U,W_squared,sigmas_hat,const,tau_hat)[[1]]
    flog.debug("~~Main loop ll=%s",ll_new)
    # Test for convergence
    ll_diff = ll_new-ll_prev
    flog.debug("ll_new=%s, ll_prev=%s, ll_diff=%s, threshold=%s",ll_new,ll_prev,ll_diff,config[["epsilon"]]*abs(ll_new))
    if (ll_diff < config[["epsilon"]]*abs(ll_new)){
      flog.debug("break")
      break
    }
    ll_prev <- ll_new

  }

  return (list("W" = W, "mus_hat" = mus_hat, "sigmas_hat" = sigmas_hat, "tau_hat" = tau_hat, "deltas_hat" = deltas_hat, "gammas_hat" = gammas_hat))

}


init_means_vars <- function(C1_names, C2_names, feature_ids, source_ids){
  #config <- config::get()
  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)
  k <- length(source_ids)
  m <- length(feature_ids)
  p1 <- length(C1_names)
  p2 <- length(C2_names)
  C1_names_ <- matrix("0",k*p1,1)
  if (p1){
    for (j in 1:k){
      C1_names_[((j-1)*p1+1):(p1*j)] <- unlist(lapply(C1_names,function(i) paste(source_ids[j],".",i,sep="")))
    }
  }
  # init
  mus_hat <- matrix(0, nrow=m, ncol=k)
  sigmas_hat <- matrix(0, nrow=m, ncol=k)
  deltas_hat <- matrix(0, nrow=m, ncol=p2)
  gammas_hat <- matrix(0, nrow=m, ncol=p1*k)
  tau_hat <- config[["tau_hat_init"]]
  # set row and column names
  rownames(mus_hat) <- feature_ids
  colnames(mus_hat) <- source_ids
  rownames(sigmas_hat) <- feature_ids
  colnames(sigmas_hat) <- source_ids
  rownames(deltas_hat) <- feature_ids
  colnames(deltas_hat) <- C2_names
  rownames(gammas_hat) <- feature_ids
  colnames(gammas_hat) <- C1_names_

  return( list("mus_hat" = mus_hat, "sigmas_hat" = sigmas_hat, "gammas_hat" = gammas_hat, "deltas_hat"= deltas_hat, "tau_hat" = tau_hat) )
}


#' @importFrom parallel clusterExport
#' @importFrom pbapply pblapply
#' @importFrom pracma repmat
#' @importFrom pracma lsqlincon
#' @importFrom nloptr nloptr
#' @importFrom matrixStats colVars
tca.fit_means_vars <- function(X, W, mus_hat, sigmas_hat, tau_hat, C2, deltas_hat, C1, gammas_hat, max_iters, parallel, num_cores){

  flog.debug("Starting function 'tca.fit_means_vars'...")

  #config <- config::get()
  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)
  nloptr_opts = list("algorithm"=config[["nloptr_opts_algorithm"]], "xtol_rel"=config[["nloptr_opts_xtol_rel"]], "print_level" = config[["nloptr_opts_print_level"]], "check_derivatives" = config[["nloptr_opts_check_derivatives"]])

  n <- nrow(X)
  m <- ncol(X)
  k <- ncol(W)
  p1 <- ncol(C1)
  p2 <- ncol(C2)

  C1_ <- calc_C1_W_interactions(W,C1)

  cl <- if (parallel) init_cluster(num_cores) else NULL

  const <- -n*log(2*pi)
  W_squared <- W**2

  ll_prev <- -Inf
  # Perform an alternative optimization of the means (mus, deltas, gammas) and variances (sigmas and tau)
  for (iter in 1:max_iters){

    flog.info("Iteration %s out of %s internal iterations", iter, max_iters)

    # (1) Estimate the means (mus, deltas, gammas)

    ub <- numeric(k+p2+p1*k)+config[["lsqlincon_inf"]]
    ub[1:k] <- if(is.null(config[["mu_max"]])) max(X) else config[["mu_max"]]
    ub[1:k] <- ub[1:k] - config[["mu_epsilon"]]
    lb <- numeric(k+p2+p1*k)-config[["lsqlincon_inf"]]
    lb[1:k] <- if(is.null(config[["mu_min"]])) min(X) else config[["mu_min"]]
    lb[1:k] <- lb[1:k] + config[["mu_epsilon"]]

    flog.debug("Calculate W_norms and related quantities")
    if (sum(colSums(mus_hat)) == 0){
      # Use the following for getting initial estimates of mus, deltas, gammas; under the assumptions that tau=0 and sigmas_{1j},...,sigmas_{kj} for each j.
      W_norms <- rowSums(W**2)**0.5
      # Since W_norms is the same for all features in this case can already calculate thesethe following quantities
      W_tilde <- W/t(repmat(W_norms,k,1))
      C1_tilde <- if (p1>0) C1_/t(repmat(W_norms,k*p1,1)) else C1_
      C2_tilde <- if (p2>0) C2/t(repmat(W_norms,p2,1)) else C2
    }else{
      flog.debug("Calculate W_norms")
      W_norms <- (tcrossprod(W**2,sigmas_hat**2)+tau_hat**2 )**0.5
    }
    X_tilde <- X/W_norms

    flog.debug("Estimate mus, deltas, gammas.")

    if (sum(colSums(mus_hat)) == 0){
      if (parallel) clusterExport(cl, varlist = c("W_tilde","C2_tilde","C1_tilde","X_tilde","lb","ub","k","p1","p2","lsqlincon"), envir=environment())
      res <- pblapply(1:m,function(j) lsqlincon(cbind(W_tilde,C2_tilde,C1_tilde), X_tilde[,j], lb = lb,  ub = ub), cl = cl )
    }else{
      if (parallel) clusterExport(cl, c("W_norms","W","C2","C1_","X_tilde","lb","ub","k","p1","p2","lsqlincon"), envir=environment())
      res <- pblapply(1:m,function(j) lsqlincon(cbind(W/t(repmat(W_norms[,j],k,1)),
                                                        if (p2>0) C2/t(repmat(W_norms[,j],p2,1)) else C2,
                                                        if (p1>0) C1_/t(repmat(W_norms[,j],k*p1,1)) else C1_),
                                                  X_tilde[,j], lb = lb,  ub = ub), cl = cl )
    }

    # Update estimtes
    for (j in 1:m){
      mus_hat[j,] = res[[j]][1:k]
      deltas_hat[j,seq(1,p2,length=p2)] = res[[j]][seq(k+1,k+p2,length=p2)]
      gammas_hat[j,seq(1,k*p1,length=k*p1)] = res[[j]][seq(k+p2+1,k+p2+p1*k,length=p1*k)]
    }


    # (2) Estimate the variances (sigmas, tau)

    # Calculate some quantities that will be repeatedly used throughout the optimization in this step
    U <- (tcrossprod(W,mus_hat) + tcrossprod(C2,deltas_hat) + tcrossprod(C1_,gammas_hat) - X)**2

    if (sum(colSums(sigmas_hat)) == 0){
      # Set a starting point for the optimization
      flog.debug("Get initial estimates of sigmas")
      row_names <- rownames(sigmas_hat)
      col_names <- colnames(sigmas_hat)
      sigmas_hat <- t(repmat((colVars(X)/k)**0.5,k,1))
      rownames(sigmas_hat) <- row_names
      colnames(sigmas_hat) <- col_names
    }

    # # (2.1) Get initial estimate of sigmas, tau under the assumption that tau is feature-specific
    #
    # flog.debug("Estimate sigmas, tau under the assumption that tau is feature-specific.")
    # lb <- numeric(k+1)+config$min_var
    # ub <- numeric(k+1)+Inf
    # if (parallel){
    #   clusterExport(cl, c("W_tilde","C2_tilde","C1_tilde","X_tilde","lb","ub","n","k","p1","p2","U","const","W_squared","sigmas_hat","tau_hat","nloptr_opts"))
    #   res <- parLapply(cl, 1:m,function(j) nloptr( x0=c(sigmas_hat[j,],tau_hat),
    #                                          eval_f=function(x,U_j,W_squared,const)
    #                                          {k<-length(x)-1;
    #                                          V <- tcrossprod(W_squared,t(x[1:k]**2))+x[k+1]**2;
    #                                          V_squared <- V**2;
    #                                          g<-numeric(k+1)
    #                                          g[1:k] <- -(colSums(W_squared*repmat(x[1:k],n,1)*t(repmat(U_j,k,1))/repmat(V_squared,1,k)) - colSums(W_squared*repmat(x[1:k],n,1)/repmat(V,1,k)))
    #                                          g[k+1] <- -(x[k+1]*(sum(U_j/V_squared)-sum(1./V)))
    #                                          return(list("objective"= -0.5*(const-sum(log(V))-sum(U_j/V)),
    #                                                      "gradient" = g))},
    #                                          lb = lb,
    #                                          ub = ub,
    #                                          opts = nloptr_opts,
    #                                          U_j = U[,j],
    #                                          W_squared = W_squared,
    #                                          const = const)$solution[1:k] )
    # }else{
    #   res <- lapply(1:m,function(j) nloptr( x0=c(sigmas_hat[j,],tau_hat),
    #                                   eval_f=function(x,U_j,W_squared,const)
    #                                   {k<-length(x)-1;
    #                                   V <- tcrossprod(W_squared,t(x[1:k]**2))+x[k+1]**2;
    #                                   V_squared <- V**2;
    #                                   g<-numeric(k+1)
    #                                   g[1:k] <- -(colSums(W_squared*repmat(x[1:k],n,1)*t(repmat(U_j,k,1))/repmat(V_squared,1,k)) - colSums(W_squared*repmat(x[1:k],n,1)/repmat(V,1,k)))
    #                                   g[k+1] <- -(x[k+1]*(sum(U_j/V_squared)-sum(1./V)))
    #                                   return(list("objective"= -0.5*(const-sum(log(V))-sum(U_j/V)),
    #                                               "gradient" = g))},
    #                                   lb = lb,
    #                                   ub = ub,
    #                                   opts = nloptr_opts,
    #                                   U_j = U[,j],
    #                                   W_squared = W_squared,
    #                                   const = const)$solution[1:k] )
    # }
    # for (j in 1:m){
    #   sigmas_hat[j,] = res[[j]]
    # }

    # (2.2) Estimate sigmas
    flog.debug("Estimate sigmas.")
    lb <- numeric(k)+config[["min_sd"]]
    ub <- numeric(k)+Inf
    if (parallel) clusterExport(cl, c("lb","ub","n","k","U","const","W_squared","sigmas_hat","tau_hat","nloptr_opts","minus_log_likelihood_sigmas"), envir=environment())
    res <- pblapply(1:m,function(j) nloptr( x0=sigmas_hat[j,],
                                          eval_f = function(x,U_j,W_squared,const,tau_hat) minus_log_likelihood_sigmas(x,U_j,W_squared,const,tau_hat),
                                          lb = lb,
                                          ub = ub,
                                          opts = nloptr_opts,
                                          U_j = U[,j],
                                          W_squared = W_squared,
                                          const = const,
                                          tau_hat = tau_hat)$solution, cl = cl )

    for (j in 1:m){
      sigmas_hat[j,] = res[[j]]
    }

    # (2.3) Estimate tau
    flog.debug("Estimate tau.")
    lb <- config[["min_sd"]]
    ub <- Inf
    tau_hat = nloptr(x0=tau_hat,
                     eval_f = function(x,U,W_squared,sigmas_hat,const) minus_log_likelihood_tau(U,W_squared,sigmas_hat,const,x),
                     lb = lb,
                     ub = ub,
                     opts = nloptr_opts,
                     U = U,
                     W_squared = W_squared,
                     sigmas_hat = sigmas_hat,
                     const = const)$solution

    flog.debug("Test for convergence.")
    ll_new <- -minus_log_likelihood_tau(U,W_squared,sigmas_hat,const,tau_hat)[[1]]
    # Test for convergence
    ll_diff = ll_new-ll_prev
    flog.debug("ll_new=%s, ll_prev=%s, ll_diff=%s, diff_threshold=%s",ll_new,ll_prev,ll_diff,config[["epsilon"]]*abs(ll_new))
    if (ll_diff < config[["epsilon"]]*abs(ll_new)){
      flog.debug("break")
      break
    }
    ll_prev <- ll_new

  }

  if (parallel) stop_cluster(cl)

  return(list("mus_hat"=mus_hat, "sigmas_hat"=sigmas_hat, "tau_hat"=tau_hat, "deltas_hat"=deltas_hat, "gammas_hat"=gammas_hat, "tau_hat"=tau_hat))

}


tca.fit_W <- function(X, W, mus_hat, sigmas_hat, tau_hat, C2, deltas_hat, C1, gammas_hat, max_iters, parallel, num_cores){

  flog.debug("Starting function 'tca.fit_W'...")

  #config <- config::get()
  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)
  nloptr_opts = list("algorithm"=config[["nloptr_opts_fit_W_algorithm"]], "xtol_rel"=config[["nloptr_opts_xtol_rel"]], "print_level" = config[["nloptr_opts_print_level"]], "check_derivatives" = config[["nloptr_opts_check_derivatives"]])

  n <- nrow(X)
  m <- ncol(X)
  k <- ncol(W)
  p1 <- ncol(C1)
  lb <- numeric(k)
  ub <- numeric(k)+1
  ones <- numeric(k)+1
  sigmas_squared <- sigmas_hat**2
  W_hat <- matrix(0, n, k)
  colnames(W_hat) <- colnames(W)
  rownames(W_hat) <- rownames(W)
  const <- -m*log(2*pi)

  cl <- if (parallel) init_cluster(num_cores) else NULL
  if (parallel) clusterExport(cl, c("lb","ub","ones","nloptr_opts","X","W","C1","C2","n","k","p1","m","const","tau_hat","mus_hat","gammas_hat","deltas_hat","sigmas_squared","minus_log_likelihood_w","calc_C1_W_interactions"), envir=environment())
  res <- pblapply(1:n,function(i) nloptr( x0=W[i,],
                                              eval_f = function(x, x_i, c1_i, mus, tau, gammas, const, C_tilde, sigmas_squared, crossprod_deltas_c2_i) minus_log_likelihood_w(t(x), x_i, c1_i, mus, tau, gammas, const, C_tilde, sigmas_squared, crossprod_deltas_c2_i),
                                              eval_g_eq = function(x, x_i, c1_i, mus, tau, gammas, const, C_tilde, sigmas_squared, crossprod_deltas_c2_i) list("constraints" = crossprod(x,ones)-1, "jacobian" = ones),
                                              lb = lb,
                                              ub = ub,
                                              opts = nloptr_opts,
                                              x_i = t(X[i,]),
                                              c1_i = t(C1[i,]),
                                              mus = mus_hat,
                                              tau = tau_hat,
                                              gammas = gammas_hat,
                                              const = const,
                                              C_tilde = if (p1>0) apply(as.matrix(1:k),1,function(h) tcrossprod(gammas_hat[,(1+(h-1)*p1):(h*p1)],t(C1[i,]))) else matrix(0,m,k),
                                              sigmas_squared = sigmas_squared,
                                              crossprod_deltas_c2_i = tcrossprod(deltas_hat,t(C2[i,])) )$solution, cl = cl )
  if (parallel) stop_cluster(cl)
  for (i in 1:n){
    W_hat[i,] = res[[i]]
  }
  return(W_hat)
}


# Creates a p1*k matrix of the interaction terms betweem each covariate in C1 and each column in W
# Both W, C1 have to be of type matrix
#' @importFrom matrixcalc hadamard.prod
calc_C1_W_interactions <- function(W,C1){
  n <- nrow(W)
  k <- ncol(W)
  p1 <- ncol(C1)
  if (p1){
    return( hadamard.prod(Reshape(Reshape(apply(W, 2, function(v) repmat(v,1,p1)), n*p1*k,1), n,p1*k), repmat(C1, 1, k)) )
  }else{
    return(matrix(0,n,0))
  }
}

# Returns the (minus) log likelihood of the entire model in a list together with the derivative of tau
# Input:
#  U <- (tcrossprod(W,mus_hat) + tcrossprod(C2,deltas_hat) + tcrossprod(C1_,gammas_hat) - X)**2
#  const <- -n*log(2*pi)
#  W_squared <- W**2
#  tau
minus_log_likelihood_tau <- function(U,W_squared,sigmas,const,tau){
  m <- ncol(U)
  res <- matrix(0,m,2)
  tmp <- lapply(1:m, function(j)
    {V <- tcrossprod(W_squared,t(sigmas[j,]**2))+tau**2;
    V_squared <- V**2;
    return (c(-0.5*(const-sum(log(V))-sum(U[,j]/V)), -(tau*(sum(U[,j]/V_squared)-sum(1./V))))) } )
  for (j in 1:m){
    res[j,] = tmp[[j]]
  }
  res <- colSums(res)
  return(list( "objective" = res[1], "gradient" = res[2] ))
}

# Returns the (minus) log likelihood of the model for a particular feature j in a list together with the gradient with respect to sigmas of feature j
# Input:
#  sigmas_hat (of a particular feature j)
#  U_j is the j-th columns of U, where U = (tcrossprod(W,mus_hat) + tcrossprod(C2,deltas_hat) + tcrossprod(C1_,gammas_hat) - X)**2
#  const <- -n*log(2*pi)
#  W_squared <- W**2
#  tau
minus_log_likelihood_sigmas <- function(sigmas,U_j,W_squared,const,tau){
  k <- length(sigmas)
  n <- nrow(W_squared)
  V <- tcrossprod(W_squared,t(sigmas**2))+tau**2
  V_squared <- V**2
  return(list("objective"= -0.5*(const-sum(log(V))-sum(U_j/V)),
              "gradient" = -(colSums(W_squared*repmat(sigmas,n,1)*t(repmat(U_j,k,1))/repmat(V_squared,1,k)) - colSums(W_squared*repmat(sigmas,n,1)/repmat(V,1,k)))))
}


# Returns the (minus) log likelihood of the model for a particular observation i in a list together with the gradient with respect to w_i
# Input:
# C_tilde - if p1=0 then C_tilde = matrix(0,m,k), otherwise it is calculated as follows:
# for (h in 1:k){
#    C_tilde[,h] = tcrossprod(gammas[,(1+(h-1)*p1):(h*p1)],c1_i)
# }
# sigmas_squared = sigmas**2
# crossprod_deltas_c2_i = tcrossprod(deltas,t(C2[i,]))
minus_log_likelihood_w <- function(w_i, x_i, c1_i, mus, tau, gammas, const, C_tilde, sigmas_squared, crossprod_deltas_c2_i){
  k <- length(w_i)
  m <- dim(mus)[1]
  c1_i_ <- calc_C1_W_interactions(w_i,c1_i)
  V <- tcrossprod(sigmas_squared,w_i**2)+tau**2
  V_rep <- repmat(V,1,k)
  U_i <- tcrossprod(mus,w_i) + crossprod_deltas_c2_i + tcrossprod(gammas,c1_i_) - t(x_i)
  U_i_squared <- U_i**2
  w_i_rep <- repmat(w_i,m,1)
  fval <- -0.5*(const-sum(log(V))-sum(U_i_squared/V))
  gval <- colSums(w_i_rep*sigmas_squared/V_rep) + colSums(( (mus+C_tilde)*repmat(U_i,1,k)*V_rep - w_i_rep*sigmas_squared*repmat(U_i_squared,1,k) ) / repmat(V**2,1,k))
  return(list("objective"= fval, "gradient" = gval))
}


# X is n by m
estimate_Z <- function(X, W, mus_hat, sigmas_hat, tau_hat, C2, deltas_hat, C1, gammas_hat, parallel, num_cores){
  flog.info("Estimate tensor...")
  # init
  n <- dim(X)[1]
  m <- dim(X)[2]
  k <- ncol(W)
  Z_hat <- list()
  for (h in 1:k){
    Z_hat[[h]] <- matrix(0,n,m)
    colnames(Z_hat[[h]]) <- colnames(X)
    rownames(Z_hat[[h]]) <- rownames(X)
  }
  # Calculate quantities that can be calculated only once
  W_prime <- replicate(n, matrix(0,k,k), simplify=F)
  for (i in 1:n){
    W_prime[[i]] = tcrossprod(W[i,],W[i,])/(tau_hat**2)
  }
  C2_prime <- (tcrossprod(C2,deltas_hat) + X)/(tau_hat**2)
  cl <- if (parallel) init_cluster(num_cores) else NULL
  if (parallel) clusterExport(cl, c("W","mus_hat","sigmas_hat","tau_hat","C1","gammas_hat","W_prime","C2_prime","estimate_Z_j"), envir=environment())
  # Estimate Z
  res <- pblapply(1:m, function(j) estimate_Z_j(W, mus_hat[j,], sigmas_hat[j,], tau_hat, C1, gammas_hat[j,], W_prime, C2_prime[,j]), cl = cl )
  if (parallel) stop_cluster(cl)
  for (j in 1:m){
    for (h in 1:k){
      Z_hat[[h]][,j] = res[[j]][,h]
    }
  }
  # add rownames and colnames and transpose matrices
  for (h in 1:k){
    rownames(Z_hat[[h]]) <- rownames(X)
    colnames(Z_hat[[h]]) <- colnames(X)
    Z_hat[[h]] <- t(Z_hat[[h]])
  }
  return(Z_hat)
}


# Estimate Z for one feature j
#' @importFrom matrixcalc matrix.inverse
estimate_Z_j <- function(W, mus_hat_j, sigmas_hat_j, tau_hat, C1, gammas_hat_j, W_prime, C2_prime_j){
  n <- nrow(W)
  k <- ncol(W)
  p1 <- ncol(C1)
  Z_j_hat <- matrix(0,n,k)
  Sig_j <- matrix.inverse(diag(sigmas_hat_j**2))
  C1_prime <- tcrossprod(C1, t(Reshape(gammas_hat_j,p1,k)))
  for (i in 1:n){
    Z_j_hat[i,] = crossprod(matrix.inverse(W_prime[[i]]+Sig_j), ( tcrossprod(Sig_j,t(mus_hat_j+C1_prime[i,])) + W[i,]*C2_prime_j[i] ) )
  }
  return(Z_j_hat)
}


tcareg.fit <- function(X, y, W, mus_hat, sigmas_hat, C2, deltas_hat, C1, gammas_hat, tau_hat, C3, test, null_model, alternative_model, parallel, num_cores){

  flog.debug("Starting function 'tcareg.fit'...")

  single_effect <- if (test == "single_effect") TRUE else FALSE

  if (test == "joint" | test == "single_effect") results <- tcareg.fit_joint(X, y, W, mus_hat, sigmas_hat, C2, deltas_hat, C1, gammas_hat, tau_hat, C3, parallel, num_cores, single_effect)
  if (test == "marginal") results <- tcareg.fit_marginal(X, y, W, mus_hat, sigmas_hat, C2, deltas_hat, C1, gammas_hat, tau_hat, C3, parallel, num_cores)
  if (test == "marginal_conditional") results <- tcareg.fit_marginal_conditional(X, y, W, mus_hat, sigmas_hat, C2, deltas_hat, C1, gammas_hat, tau_hat, C3, parallel, num_cores)
  if (test == "custom") results <- tcareg.fit_custom(X, y, W, mus_hat, sigmas_hat, C2, deltas_hat, C1, gammas_hat, tau_hat, C3, parallel, num_cores, null_model, alternative_model)

  return(results)

}


tcareg.fit_joint <- function(X, y, W, mus_hat, sigmas_hat, C2, deltas_hat, C1, gammas_hat, tau_hat, C3, parallel, num_cores, single_effect){
  k <- ncol(W)
  mdl1 <- tcareg.optimize(X, y, W, mus_hat, sigmas_hat, C2, deltas_hat, C1, gammas_hat, tau_hat, C3, 1:k, parallel, num_cores, single_effect)
  ll1 <- mdl1[["ll"]]
  d <- data.frame(y,C3)
  mdl0 <- lm(y~.,data = d) # the null log-likelihood under the model y ~ C3 (i.e. the null under the case of no source-specific effects)
  ll0 <-numeric(length(mdl1[["ll"]])) + as.numeric(logLik(mdl0))
  df <- if(single_effect) 1 else ncol(W)
  lrt <- lrt.test(ll0, mdl1[["ll"]], df = df)
  qvals <- p.adjust(lrt[["pvals"]], method = "BH")
  return( list("phi" = mdl1[["phi"]], "beta" = mdl1[["beta"]], "intercept" = mdl1[["intercept"]], "alpha" = mdl1[["alpha"]], "null_ll" = ll0, "alternative_ll" = ll1, "pvals" = lrt[["pvals"]], "qvals" = qvals, "stats" = lrt[["stats"]], "df" = df) )
}


tcareg.fit_marginal <- function(X, y, W, mus_hat, sigmas_hat, C2, deltas_hat, C1, gammas_hat, tau_hat, C3, parallel, num_cores){
  m <- ncol(X)
  k <- ncol(W)
  results <- list()
  d <- data.frame(y,C3)
  mdl0 <- lm(y~.,data = d) # the null log-likelihood under the model y ~ C3 (i.e. the null under the case of no source-specific effects)
  for (h in 1:k){
    mdl1 <- tcareg.optimize(X, y, W, mus_hat, sigmas_hat, C2, deltas_hat, C1, gammas_hat, tau_hat, C3, h, parallel, num_cores, FALSE)
    ll1 <- mdl1[["ll"]]
    ll0 <-numeric(length(mdl1[["ll"]])) + as.numeric(logLik(mdl0))
    df <- 1
    lrt <- lrt.test(ll0, ll1, df = df)
    qvals <- p.adjust(lrt[["pvals"]], method = "BH")
    results[[length(results)+1]] = list("phi" = mdl1[["phi"]], "beta" = mdl1[["beta"]], "intercept" = mdl1[["intercept"]], "alpha" = mdl1[["alpha"]], "null_ll" = ll0, "alternative_ll" = ll1, "pvals" = lrt[["pvals"]], "qvals" = qvals, "stats" = lrt[["stats"]], "df" = df)
  }
  return(results)
}


tcareg.fit_marginal_conditional <- function(X, y, W, mus_hat, sigmas_hat, C2, deltas_hat, C1, gammas_hat, tau_hat, C3, parallel, num_cores){
  m <- ncol(X)
  k <- ncol(W)
  results <- list()
  mdl1 <- tcareg.optimize(X, y, W, mus_hat, sigmas_hat, C2, deltas_hat, C1, gammas_hat, tau_hat, C3, 1:k, parallel, num_cores, FALSE)
  ll1 <- mdl1[["ll"]]
  for (h in 1:k){
    mdl0 <- tcareg.optimize(X, y, W, mus_hat, sigmas_hat, C2, deltas_hat, C1, gammas_hat, tau_hat, C3, setdiff(1:k,h), parallel, num_cores, FALSE)
    ll0 <- mdl0[["ll"]]
    df <- 1
    lrt <- lrt.test(ll0, ll1, df = df)
    qvals <- p.adjust(lrt[["pvals"]], method = "BH")
    results[[length(results)+1]] = list("phi" = mdl1[["phi"]], "beta" = mdl1[["beta"]], "intercept" = mdl1[["intercept"]], "alpha" = mdl1[["alpha"]], "null_ll" = ll0, "alternative_ll" = ll1, "pvals" = lrt[["pvals"]], "qvals" = qvals, "stats" = lrt[["stats"]], "df" = df)
  }
  return(results)
}


tcareg.fit_custom <- function(X, y, W, mus_hat, sigmas_hat, C2, deltas_hat, C1, gammas_hat, tau_hat, C3, parallel, num_cores, null_model, alternative_model){
  mdl0 <- tcareg.optimize(X, y, W, mus_hat, sigmas_hat, C2, deltas_hat, C1, gammas_hat, tau_hat, C3, null_model, parallel, num_cores, FALSE)
  ll0 <- mdl0[["ll"]]
  mdl1 <- tcareg.optimize(X, y, W, mus_hat, sigmas_hat, C2, deltas_hat, C1, gammas_hat, tau_hat, C3, alternative_model, parallel, num_cores, FALSE)
  ll1 <- mdl1[["ll"]]
  df <- length(alternative_model) - length(null_model)
  lrt <- lrt.test(ll0, ll1, df = df)
  qvals <- p.adjust(lrt[["pvals"]], method = "BH")
  return( list("phi" = mdl1[["phi"]], "beta" = mdl1[["beta"]], "intercept" = mdl1[["intercept"]], "alpha" = mdl1[["alpha"]], "null_ll" = ll0, "alternative_ll" = ll1, "pvals" = lrt[["pvals"]], "qvals" = qvals, "stats" = lrt[["stats"]], "df" = df) )
}


tcareg.optimize <- function(X, y, W, mus_hat, sigmas_hat, C2, deltas_hat, C1, gammas_hat, tau_hat, C3, mdl, parallel, num_cores, single_effect){

  flog.debug("...")

  #config <- config::get()
  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)

  nloptr_opts = list("algorithm"=config[["nloptr_opts_fit_conditional_algorithm"]], "xtol_rel"=config[["nloptr_opts_xtol_rel"]], "print_level" = config[["nloptr_opts_print_level"]], "check_derivatives" = config[["nloptr_opts_check_derivatives"]])
  lambda <- config[["lambda"]]

  l <- length(mdl)
  k <- ncol(W)
  n <- nrow(W)
  m <- ncol(X)
  p1 <- ncol(C1)
  p3 <- ncol(C3)

  # define lower of upper bounds for the optimization
  num_betas <- if(single_effect) 1 else l
  ub <- numeric(1+num_betas+p3+1)+Inf # adding one for intercept
  lb <- numeric(1+num_betas+p3+1)-Inf # adding one for intercept
  lb[1] = config[["min_sd"]]

  # # find initial estimates for phi, beta, alpha for each feature j using the model y~Z_hat[,j]+C3; note that an intercept term is implicitly added to C3; In practice, this didn't work well.
  # X0 <- get_initial_estimates(Z_hat, y, C3, mdl, lambda, single_effect, parallel, num_cores)

  # Set initial estimates for phi and alpha for each feature j using the model y~C3; initialize all betas with 0
  d <- data.frame(y,C3)
  mdl0 <- lm(y~.,data = d)
  X0 <- list("phi0" = matrix(sqrt(sum(mdl0$residuals**2)/(n-1)),m,1),
             "beta0" = matrix(0,m,num_betas),
             "alpha0" = repmat(t(as.matrix(mdl0$coefficients)),m,1) )

  # include an intercept term in C3; required for TCA
  C3 <- cbind(as.matrix(numeric(n)+1), C3)

  # calculate quantaties that can be calculated only once
  const <- -n*log(2*pi)
  sigmas_squared <- sigmas_hat**2
  C1_ <- calc_C1_W_interactions(W,C1)
  X_tilde = X - (tcrossprod(C2,deltas_hat) + tcrossprod(C1_,gammas_hat) + tcrossprod(W,mus_hat))
  V <- tau_hat**2 + tcrossprod(W**2,sigmas_squared)

  # Fit the model for each feature
  cl <- if (parallel) init_cluster(num_cores) else NULL
  if (parallel) clusterExport(cl, c("lb","ub","p1","k","W","y","C1","C3","n","l","X0","const","sigmas_squared","X_tilde","V","mdl","gammas_hat","sigmas_squared","mus_hat","lambda","single_effect","nloptr_opts","conditional_model_minus_log_likelihood","tcareg.optimize_j"), envir=environment())
  res <- pblapply(1:m, function(j) tcareg.optimize_j(j, y, W, C1, C3, mus_hat, gammas_hat, const, sigmas_squared, V, X_tilde, X0, mdl, lb, ub, nloptr_opts, single_effect, lambda), cl = cl)
  if (parallel) stop_cluster(cl)
  intercept <- matrix(0, m, 1)
  alpha <- matrix(0, m, p3)
  beta <- matrix(0, m, num_betas)
  phi <- matrix(0, m, 1)
  ll <- matrix(0, m, 1)
  for (j in 1:m){
    phi[j,] <- res[[j]][["solution"]][1]
    beta[j,] <- res[[j]][["solution"]][2:(1+num_betas)]
    intercept[j,] <- res[[j]][["solution"]][2+num_betas]
    alpha[j,] <- res[[j]][["solution"]][-1:-(2+num_betas)]
    ll[j,] <- -res[[j]][["objective"]] # take the negative to get the log likelihood
  }
  return( list("phi" = phi, "beta" = beta, "intercept" = intercept, "alpha" = alpha, "ll" = ll) )
}


# fit the model for a particular feature
tcareg.optimize_j <- function(j, y, W, C1, C3, mus_hat, gammas_hat, const, sigmas_squared, V, X_tilde, X0, mdl, lb, ub, nloptr_opts, single_effect, lambda){

  n <- nrow(W)
  k <- ncol(W)
  l <- length(mdl)
  p1 <- ncol(C1)
  p3 <- ncol(C3)

  # calculate external quantities (e1,e2,...) that need to be calculated only once per feature j
  gammas_hat_j_tilde <- t(Reshape(gammas_hat[j,],p1,k)[,mdl])
  gammas_hat_j_tilde <- if (p1 == 1) t(gammas_hat_j_tilde) else gammas_hat_j_tilde
  sigmas_squared_j <- sigmas_squared[j,]
  e1 <- V[,j]
  e2 <- W[,mdl]*repmat(sigmas_squared[j,mdl],n,1)
  e3 <- t(repmat(e1,l,1))
  e4 <- repmat(mus_hat[j,mdl],n,1) + tcrossprod(C1, gammas_hat_j_tilde)  + (e2 * t(repmat(X_tilde[,j],l,1)) )/e1
  e5 <- repmat(sigmas_squared[j,mdl],n,1)
  e6 <- e2/e3
  if (single_effect){
    nloptr_res <- nloptr( x0=c(X0[["phi0"]][j,1], X0[["beta0"]][j,], X0[["alpha0"]][j,]), eval_f = function(x, y, C3, const, sigmas_squared_j, e1, e2, e3, e4, e5, e6, mdl, lambda){
          res <- conditional_model_minus_log_likelihood(x = c(x[1],repmat(x[2],1,length(mdl)),x[-1:-2]), y, C3, const, sigmas_squared_j, e1, e2, e3, e4, e5, e6, mdl, lambda)
          return( list( "objective" = res[["objective"]], "gradient" = c(res[["gradient"]][1], sum(res[["gradient"]][2:(1+length(mdl))]), res[["gradient"]][-1:-(1+length(mdl))]) ) )
          },
        lb = lb, ub = ub, opts = nloptr_opts,
        y = y, C3 = C3, const = const, sigmas_squared_j = sigmas_squared_j, e1 = e1, e2 = e2, e3 = e3, e4 = e4, e5 = e5, e6 = e6, mdl = mdl, lambda = lambda)
  }else{
    nloptr_res <- nloptr( x0=c(X0[["phi0"]][j,1], X0[["beta0"]][j,], X0[["alpha0"]][j,]), eval_f = function(x, y, C3, const, sigmas_squared_j, e1, e2, e3, e4, e5, e6, mdl, lambda)
      conditional_model_minus_log_likelihood(x, y, C3, const, sigmas_squared_j, e1, e2, e3, e4, e5, e6, mdl, lambda),
      lb = lb, ub = ub, opts = nloptr_opts,
      y = y, C3 = C3, const = const, sigmas_squared_j = sigmas_squared_j, e1 = e1, e2 = e2, e3 = e3, e4 = e4, e5 = e5, e6 = e6, mdl = mdl, lambda = lambda)
  }
  return( list("solution" = nloptr_res[["solution"]], "objective" = nloptr_res[["objective"]]) )
}


conditional_model_minus_log_likelihood <- function(x, y, C3, const, sigmas_squared_j, e1, e2, e3, e4, e5, e6, mdl, lambda){
  l <- length(mdl)
  p3 <- ncol(C3)
  phi = x[1]
  beta_j = x[2:(l+1)]
  alpha_j = x[-1:-(1+l)]
  s1 <- beta_j**2
  s2 <- phi**2 + sum((s1)*(sigmas_squared_j[mdl])) - (tcrossprod(e2,t(beta_j))**2)/e1 # sigma_ij_tilde_squared
  s3 <- tcrossprod(C3,t(alpha_j)) + tcrossprod(e4,t(beta_j)) - y  # mu_ij_tilde - y
  s4 <- s3/s2
  s5 <- s3*s4
  s6 <- 1/s2
  s7 <- (e5 - e6*repmat(tcrossprod(e2,t(beta_j)),1,l) )/repmat(s2,1,l)
  fval <- -0.5*(const - sum(log(s2)) - sum(s5)) + 0.5*lambda*sum(s1)
  gval = numeric(length(x))
  gval[1] <- phi*(sum(s6) - sum(s5*s6))
  gval[2:(l+1)] <- beta_j * colSums(s7) + colSums(e4*repmat(s4,1,l)) - beta_j * colSums(s7*repmat(s5,1,l)) + lambda*beta_j
  gval[-1:-(1+l)] <- colSums(C3 * repmat(s4,1,p3))
  return(list("objective"= fval, "gradient" = gval))
}


## The following can be used for getting initial estimates for the optimization using Z_hat; in practice, this worked worse than a naive initialization.
#'
#' # the returned alpha0 contains intercept term at the beginning
#' get_initial_estimates <- function(Z_hat, y, C3, mdl, lambda, single_effect, parallel, num_cores){
#'   l <- length(mdl)
#'   num_betas <- if(single_effect) 1 else l
#'   n <- length(y)
#'   m <- ncol(Z_hat[[1]])
#'   p3 <- ncol(C3)
#'   phi0 <- matrix(0,m,1)
#'   beta0 <- matrix(0,m,num_betas)
#'   alpha0 <- matrix(0,m,p3+1) # include intercept terms
#'   cl <- if (parallel) init_cluster(num_cores) else NULL
#'   if (parallel) clusterExport(cl, c("n","l","p3","Z_hat","mdl","y","lambda","C3","get_initial_estimates_j","single_effect"), envir=environment())
#'   res <- pblapply(1:m, function(j) get_initial_estimates_j(j, n, l, p3, Z_hat, mdl, y, lambda, C3, single_effect), cl = cl )
#'   if (parallel) stop_cluster(cl)
#'   for (j in 1:m){
#'     beta0[j,] = res[[j]][["beta0"]]
#'     alpha0[j,] = res[[j]][["alpha0"]]
#'     phi0[j,1] = res[[j]][["phi0"]]
#'   }
#'   return( list("phi0" = phi0, "beta0" = beta0, "alpha0" = alpha0) )
#' }
#'
#' #' @importFrom glmnet glmnet
#' get_initial_estimates_j <- function(j, n, l, p3, Z_hat, mdl, y, lambda, C3, single_effect){
#'   if (single_effect){
#'     U <- matrix(0,n,1+p3)
#'     for (h in mdl){
#'       U[,1] = U[,1] + Z_hat[[h]][,j]
#'     }
#'     U[,-1] <- C3
#'   }else{
#'     U <- matrix(0,n,l+p3)
#'     for (h in mdl){
#'       U[,h] = Z_hat[[h]][,j]
#'     }
#'     U[,-1:-l] <- C3
#'   }
#'   num_betas <- if (single_effect) 1 else l
#'   p <- numeric(num_betas+p3)
#'   p[1:num_betas] <- 1
#'   coeffs <- as.numeric(coefficients(glmnet(U, y, family="gaussian", alpha=0, lambda = lambda, penalty.factor = p, intercept = TRUE)))
#'   phi0 <- sqrt(sum((y-tcrossprod(cbind(numeric(n), U),t(coeffs)))**2)/n)
#'   beta0 <- coeffs[2:(num_betas+1)]
#'   alpha0 <- c(coeffs[1], coeffs[-1:-(1+num_betas)])
#'   return( list("beta0" = beta0, "alpha0" = coeffs[-1:-num_betas], "phi0" =  phi0) )
#' }

