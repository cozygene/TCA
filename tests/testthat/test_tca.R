library("TCA")

context("Test TCA")

# # Simulated data from matlab was simulated with using n = 200; m = 100; k = 3; tau = 0.01; p1 = 2; p2 = 2;
# test_that("Results of tca", {
#   data <- test_data(50, 100, 3, 2, 2, 0.01, log_file = NULL)
#   expect_is(data$W, "matrix")
# })


test_that("Comparison with the results given by the matlab version", {

  skip_on_cran()

  basedir <- "../assets/"
  # load data
  X <- as.matrix(read.table(paste(basedir,"X.txt",sep=""), header = FALSE, sep=","))
  W <- as.matrix(read.table(paste(basedir,"W.txt",sep=""), header = FALSE, sep=","))
  C1 <- as.matrix(read.table(paste(basedir,"C1.txt",sep=""), header = FALSE, sep=","))
  C2 <- as.matrix(read.table(paste(basedir,"C2.txt",sep=""), header = FALSE, sep=","))
  rownames(X) <- 1:nrow(X)
  colnames(X) <- 1:ncol(X)
  rownames(W) <- rownames(X)
  colnames(W) <- 1:ncol(W)
  rownames(C1) <- rownames(X)
  colnames(C1) <- 1:ncol(C1)
  rownames(C2) <- rownames(X)
  colnames(C2) <- 1:ncol(C2)

  # load real parameters
  gammas <- as.matrix(read.table(paste(basedir,"gammas.txt",sep=""), header = FALSE, sep=","))
  deltas <- as.matrix(read.table(paste(basedir,"deltas.txt",sep=""), header = FALSE, sep=","))
  mus <- as.matrix(read.table(paste(basedir,"mus.txt",sep=""), header = FALSE, sep=","))
  sigmas <- as.matrix(read.table(paste(basedir,"sigmas.txt",sep=""), header = FALSE, sep=","))

  # load matlab estimates
  matlab.mus_hat <- as.matrix(read.table(paste(basedir,"mus_hat.txt",sep=""), header = FALSE, sep=","))
  matlab.sigmas_hat <- as.matrix(read.table(paste(basedir,"sigmas_hat.txt",sep=""), header = FALSE, sep=","))
  matlab.gammas_hat <- as.matrix(read.table(paste(basedir,"gammas_hat.txt",sep=""), header = FALSE, sep=","))
  matlab.deltas_hat <- as.matrix(read.table(paste(basedir,"deltas_hat.txt",sep=""), header = FALSE, sep=","))
  matlab.W_hat <- as.matrix(read.table(paste(basedir,"W_hat.txt",sep=""), header = FALSE, sep=","))

  # run tca
  res <- tca(X = t(X), W = W, C1 = C1, C2 = C2, refit_W = TRUE, refit_W.sparsity = ncol(X), parallel = FALSE, log_file = NULL, constrain_mu = TRUE, vars.mle = TRUE)

  # compare parameter estimates
  for (h in 1:ncol(res$mus_hat)){
    expect_equal(cor(res$mus_hat[,h],matlab.mus_hat[,h]) > 0.99, TRUE)
    expect_equal(cor(res$sigmas_hat[,h],matlab.sigmas_hat[,h]) > 0.95, TRUE)
  }
  for (h in 1:ncol(res$gammas_hat)){
    expect_equal(cor(res$gammas_hat[,h],matlab.gammas_hat[,h]) > 0.95, TRUE)
  }
  for (h in 1:ncol(res$deltas_hat)){
    expect_equal(cor(res$deltas_hat[,h],matlab.deltas_hat[,h]) > 0.99, TRUE)
  }
  for (h in 1:ncol(res$W)){
    expect_equal(cor(res$W[,h],matlab.W_hat[,h]) > 0.99, TRUE)
  }

  # use noisy initial estimates of W
  W_noisy <- as.matrix(read.table(paste(basedir,"W_noisy.txt",sep=""), header = FALSE, sep=","))
  colnames(W_noisy) <- 1:ncol(W_noisy)
  rownames(W_noisy) <- 1:nrow(W_noisy)
  res2 <- tca(X = t(X), W = W_noisy, C1 = C1, C2 = C2, refit_W = TRUE, refit_W.sparsity = ncol(X), parallel = FALSE, log_file = NULL, constrain_mu = TRUE, vars.mle = TRUE)
  matlab.W_noisy_hat <- as.matrix(read.table(paste(basedir,"W_noisy_hat.txt",sep=""), header = FALSE, sep=","))
  for (h in 1:ncol(res2$W)){
    expect_equal(cor(res2$W[,h],matlab.W_noisy_hat[,h]) > 0.99, TRUE)
  }

})

#
test_that("Comparison with the results given by the matlab version for tensor", {

  skip_on_cran()

  basedir <- "../assets/"
  # load data
  X <- as.matrix(read.table(paste(basedir,"X.txt",sep=""), header = FALSE, sep=","))
  colnames(X) <- 1:ncol(X)
  rownames(X) <- 1:nrow(X)
  tca.mdl <- list()
  tca.mdl[["W"]] <- as.matrix(read.table(paste(basedir,"W.txt",sep=""), header = FALSE, sep=","))
  tca.mdl[["mus_hat"]] <- t(as.matrix(read.table(paste(basedir,"mus.txt",sep=""), header = FALSE, sep=",")))
  tca.mdl[["sigmas_hat"]] <- t(as.matrix(read.table(paste(basedir,"sigmas.txt",sep=""), header = FALSE, sep=",")))
  tca.mdl[["gammas_hat"]] <- as.matrix(read.table(paste(basedir,"gammas.txt",sep=""), header = FALSE, sep=","))
  tca.mdl[["deltas_hat"]] <- as.matrix(read.table(paste(basedir,"deltas.txt",sep=""), header = FALSE, sep=","))
  tca.mdl[["tau_hat"]] <- 0.01
  tca.mdl[["C1"]] <- as.matrix(read.table(paste(basedir,"C1.txt",sep=""), header = FALSE, sep=","))
  tca.mdl[["C2"]] <- as.matrix(read.table(paste(basedir,"C2.txt",sep=""), header = FALSE, sep=","))

  # load matlab's estimates of Z
  matlab.Z1_hat <- t(as.matrix(read.table(paste(basedir,"Z1_hat.txt",sep=""), header = FALSE, sep=",")))
  matlab.Z2_hat <- t(as.matrix(read.table(paste(basedir,"Z2_hat.txt",sep=""), header = FALSE, sep=",")))
  matlab.Z3_hat <- t(as.matrix(read.table(paste(basedir,"Z3_hat.txt",sep=""), header = FALSE, sep=",")))

  # run tensor
  Z_hat <- tensor(t(X), tca.mdl, log_file = NULL)

  # compare Z estimates
  th = 0.001
  expect_equal(mean((Z_hat[[1]]-matlab.Z1_hat)**2) < th, TRUE)
  expect_equal(mean((Z_hat[[2]]-matlab.Z2_hat)**2) < th, TRUE)
  expect_equal(mean((Z_hat[[3]]-matlab.Z3_hat)**2) < th, TRUE)
  for (h in 1:ncol(matlab.Z1_hat)){
    expect_equal(cor(matlab.Z1_hat[,h],Z_hat[[1]][,h]) > 0.99, TRUE)
  }
  for (h in 1:ncol(matlab.Z1_hat)){
    expect_equal(cor(matlab.Z2_hat[,h],Z_hat[[2]][,h]) > 0.99, TRUE)
  }
  for (h in 1:ncol(matlab.Z1_hat)){
    expect_equal(cor(matlab.Z3_hat[,h],Z_hat[[3]][,h]) > 0.99, TRUE)
  }

})



run_ewas <- function(X, W, tca.mdl, basedir, yfile, test, fast_mode = FALSE){

  Y <- as.matrix(read.table(paste(basedir,yfile,sep=""), header = FALSE, sep=","))
  rownames(Y) <- 1:nrow(Y)
  m <- ncol(Y)
  k <- ncol(W)
  tca.mdl2 <- tca.mdl
  pvals <- matrix(0,m,1)
  if (test == "marginal" | test == "marginal_conditional") pvals <- matrix(0,m,k)
  for (j in 1:m){
    tca.mdl2[["mus_hat"]] <- t(as.matrix(tca.mdl[["mus_hat"]][j,]))
    tca.mdl2[["sigmas_hat"]] <- t(as.matrix(tca.mdl[["sigmas_hat"]][j,]))
    tca.mdl2[["deltas_hat"]] <- t(as.matrix(tca.mdl[["deltas_hat"]][j,]))
    tca.mdl2[["gammas_hat"]] <- t(as.matrix(tca.mdl[["gammas_hat"]][j,]))
    y <- Y[,j,drop=F]
    x <- t(X[,j,drop=F])
    rownames(y) <- 1:length(y)
    if (test == "custom"){
      res <- tcareg(X = x, tca.mdl = tca.mdl2, y = y, test = test, save_results = FALSE, null_model = c("3"), alternative_model = c("1","2","3"), log_file = NULL, fast_mode = fast_mode, verbose = FALSE)
    }else{
      res <- tcareg(X = x, tca.mdl = tca.mdl2, y = y, test = test, save_results = FALSE, log_file = NULL, fast_mode = fast_mode, verbose = FALSE)
    }
    if (test == "single_effect" | test == "custom" | test == "joint") pvals[j,] <- res$pval[1]
    if (test == "marginal" | (fast_mode == FALSE & test == "marginal_conditional")) for (h in 1:k) pvals[j,h] <- res[[h]]$pval[1]
    if (fast_mode == TRUE & test == "marginal_conditional") pvals[j,] <- res$pval
  }
  return(pvals)
}


test_that("Evaluate tcareg with the results of the matlab version", {

  skip_on_cran()

  basedir <- "../assets/ewas/"
  # load data
  X <- as.matrix(read.table(paste(basedir,"TCA_methylation_exp2_sim1_k_6_n_500_data.txt",sep=""), header = FALSE, sep=","))
  colnames(X) <- 1:ncol(X)
  rownames(X) <- 1:nrow(X)
  W <- as.matrix(read.table(paste(basedir,"TCA_methylation_exp2_sim1_k_6_n_500_W.txt",sep=""), header = FALSE, sep=","))
  colnames(W) <- 1:ncol(W)
  rownames(W) <- 1:nrow(W)
  m <- ncol(X)
  k <- ncol(W)

  tca.mdl <- tca(X = t(X), W = W, refit_W = FALSE, parallel = FALSE, log_file = NULL, constrain_mu = TRUE, vars.mle = TRUE)

  # (1) joint test
  # power
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_joint_effect_size_10_W_noise_level_0.txt", test="joint", fast_mode = FALSE)
  expect_equal(sum(pvals < 0.05/m)/m > 0.5, TRUE)
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_joint_effect_size_10_W_noise_level_0.txt", test="joint", fast_mode = TRUE)
  expect_equal(sum(pvals < 0.05/m)/m > 0.5, TRUE)
  # fps rate
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_fps_effect_size_1_W_noise_level_0.txt", test="joint", fast_mode = FALSE)
  expect_equal(sum(pvals < 0.05/m)/m < 0.05, TRUE)
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_fps_effect_size_1_W_noise_level_0.txt", test="joint", fast_mode = TRUE)
  expect_equal(sum(pvals < 0.05/m)/m < 0.05, TRUE)

  # (2) marginal test
  # power
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_marginal_ct_1_effect_size_5_W_noise_level_0.txt", test="marginal", fast_mode = FALSE)
  expect_equal(sum(pvals[,1] < 0.05/m)/m > 0.9, TRUE)
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_marginal_ct_1_effect_size_5_W_noise_level_0.txt", test="marginal", fast_mode = TRUE)
  expect_equal(sum(pvals[,1] < 0.05/m)/m > 0.9, TRUE)
  # fps rate
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_fps_effect_size_1_W_noise_level_0.txt", test="marginal", fast_mode = FALSE)
  for (h in 1:k){
    expect_equal(sum(pvals[,h] < 0.05)/m < 0.1, TRUE)
  }
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_fps_effect_size_1_W_noise_level_0.txt", test="marginal", fast_mode = TRUE)
  for (h in 1:k){
    expect_equal(sum(pvals[,h] < 0.05)/m < 0.1, TRUE)
  }

  # (3) conditional marginal test
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_marginal_ct_1_effect_size_5_W_noise_level_0.txt", test="marginal_conditional", fast_mode = FALSE)
  expect_equal(sum(pvals[,1] < 0.05)/m > 0.9, TRUE)
  for (h in 2:k){
    expect_equal(sum(pvals[,h] < 0.05)/m > 0.1, FALSE)
  }
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_marginal_ct_1_effect_size_5_W_noise_level_0.txt", test="marginal_conditional", fast_mode = TRUE)
  expect_equal(sum(pvals[,1] < 0.05)/m > 0.9, TRUE)
  for (h in 2:k){
    expect_equal(sum(pvals[,h] < 0.05)/m > 0.1, FALSE)
  }

  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_marginal_ct_3_effect_size_10_W_noise_level_0.txt", test="marginal_conditional", fast_mode = FALSE)
  l <- 3
  expect_equal(sum(pvals[,l] < 0.05)/m > 0.1, TRUE)
  for (h in setdiff(1:k,l)){
    expect_equal(sum(pvals[,h] < 0.05)/m > 0.1, FALSE)
  }
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_marginal_ct_3_effect_size_10_W_noise_level_0.txt", test="marginal_conditional", fast_mode = TRUE)
  l <- 3
  expect_equal(sum(pvals[,l] < 0.05)/m > 0.1, TRUE)
  for (h in setdiff(1:k,l)){
    expect_equal(sum(pvals[,h] < 0.05)/m > 0.1, FALSE)
  }

  # (4) joint test - single beta
  # power
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_joint_single_beta_effect_size_10_W_noise_level_0.txt", test="single_effect", fast_mode = FALSE)
  expect_equal(sum(pvals < 0.05/m)/m > 0.8, TRUE)
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_joint_single_beta_effect_size_10_W_noise_level_0.txt", test="single_effect", fast_mode = TRUE)
  expect_equal(sum(pvals < 0.05/m)/m > 0.8, TRUE)
  # fps rate
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_fps_effect_size_1_W_noise_level_0.txt", test="single_effect", fast_mode = FALSE)
  expect_equal(sum(pvals < 0.05/m)/m > 0.05, FALSE)
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_fps_effect_size_1_W_noise_level_0.txt", test="single_effect", fast_mode = TRUE)
  expect_equal(sum(pvals < 0.05/m)/m > 0.05, FALSE)

  # (5) custom test
  # power
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_joint_effect_size_10_W_noise_level_0.txt", test="custom", fast_mode = FALSE)
  expect_equal(sum(pvals < 0.05/m)/m > 0.4, TRUE)
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_joint_effect_size_10_W_noise_level_0.txt", test="custom", fast_mode = TRUE)
  expect_equal(sum(pvals < 0.05/m)/m > 0.4, TRUE)
  # fps rate
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_fps_effect_size_1_W_noise_level_0.txt", test="custom", fast_mode = FALSE)
  expect_equal(sum(pvals < 0.05/m)/m > 0.1, FALSE)
  pvals <- run_ewas(X, W, tca.mdl, basedir, "TCA_methylation_exp2_sim1_k_6_n_500_Y_fps_effect_size_1_W_noise_level_0.txt", test="custom", fast_mode = TRUE)
  expect_equal(sum(pvals < 0.05/m)/m > 0.1, FALSE)

})


test_that("Verify that tcareg with save_results set to TRUE functions well ", {

  skip_on_cran()

  outfile.prefix <- "test_tca.tmp_file"
  n <- 50
  m <- 5
  k <- 3
  p1 <- 2
  p2 <- 2
  data <- test_data(n, m, k, p1, p2, 0.01)
  tca.mdl <- tca(X = data$X, W = data$W, C1 = data$C1, C2 = data$C2, log_file = NULL)
  y <- matrix(rexp(n, rate=.1), ncol=1)
  C3 <- data$C2
  colnames(C3) <- c("C1_1","C1_2")
  rownames(y) <- rownames(data$W)

  pass <- FALSE
  result = tryCatch({
    for (fast_mode in c(TRUE,FALSE)){
      res <- tcareg(data$X, tca.mdl, y, C3 = C3, test = "joint", save_results = TRUE, fast_mode = fast_mode, output = outfile.prefix, log_file = NULL)
      res <- tcareg(data$X, tca.mdl, y, C3 = C3, test = "marginal", save_results = TRUE, fast_mode = fast_mode, output = outfile.prefix, log_file = NULL)
      res <- tcareg(data$X, tca.mdl, y, C3 = C3, test = "marginal_conditional", save_results = TRUE, fast_mode = fast_mode, output = outfile.prefix, log_file = NULL)
      res <- tcareg(data$X, tca.mdl, y, C3 = C3, test = "single_effect", save_results = TRUE, fast_mode = fast_mode, output = outfile.prefix, log_file = NULL)
      res <- tcareg(data$X, tca.mdl, y, C3 = C3, test = "custom", save_results = TRUE, null_model = c("3"), alternative_model = c("1","2","3"), fast_mode = fast_mode, output = outfile.prefix, log_file = NULL)
    }
    pass <- TRUE
  }, finally = {
    for (f in list.files(path=".", pattern=glob2rx(paste(outfile.prefix,"*",sep="")), full.names=TRUE)){
      if (file.exists(f)) file.remove(f)
    }
  })
  expect_equal(pass, TRUE)

})



test_that("Verify that tcareg with parallel set to TRUE functions well ", {

  skip_on_cran()

  n <- 50
  m <- 5
  k <- 3
  p1 <- 2
  p2 <- 2
  data <- test_data(n, m, k, p1, p2, 0.01)
  tca.mdl <- tca(X = data$X, W = data$W, C1 = data$C1, C2 = data$C2, log_file = NULL)
  y <- matrix(rexp(n, rate=.1), ncol=1)
  C3 <- data$C2
  colnames(C3) <- c("C1_1","C1_2")
  rownames(y) <- rownames(data$W)

  pass <- FALSE
  result = tryCatch({
    for (fast_mode in c(TRUE,FALSE)){
      res <- tcareg(data$X, tca.mdl, y, C3 = C3, test = "joint", save_results = FALSE, fast_mode = fast_mode, log_file = NULL, parallel = TRUE)
      res <- tcareg(data$X, tca.mdl, y, C3 = C3, test = "marginal", save_results = FALSE, fast_mode = fast_mode, log_file = NULL, parallel = TRUE)
      res <- tcareg(data$X, tca.mdl, y, C3 = C3, test = "marginal_conditional", save_results = FALSE, fast_mode = fast_mode, log_file = NULL, parallel = TRUE)
      res <- tcareg(data$X, tca.mdl, y, C3 = C3, test = "single_effect", save_results = FALSE, fast_mode = fast_mode, log_file = NULL, parallel = TRUE)
      res <- tcareg(data$X, tca.mdl, y, C3 = C3, test = "custom", save_results = FALSE, null_model = c("2","3"), alternative_model = c("1","2","3"), fast_mode = fast_mode, log_file = NULL, parallel = TRUE)
    }
    pass <- TRUE
  })
  expect_equal(pass, TRUE)

})
