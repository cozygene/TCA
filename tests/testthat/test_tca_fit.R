library("TCA")

context("Test tca fit")

test_that("Verify that p-values are not returned for gamma and delta of constrain_mu == TRUE", {

  skip_on_cran()

  n = 20; m = 30; k = 3; tau = 0.01; p1 = 2; p2 = 2;
  res <- test_data(n, m, k, p1, p2, tau, log_file = NULL)

  tca.mdl <- tca(X = res$X, W = res$W, C1 = res$C1, C2 = res$C2, refit_W = FALSE, parallel = FALSE, log_file = NULL, constrain_mu = TRUE)
  print(tca.mdl$gammas_hat_pvals)
  print(tca.mdl$deltas_hat_pvals)
  expect_equal(is.null(tca.mdl$gammas_hat_pvals), TRUE)
  expect_equal(is.null(tca.mdl$deltas_hat_pvals), TRUE)

  tca.mdl <- tca(X = res$X, W = res$W, C1 = res$C1, C2 = res$C2, refit_W = FALSE, parallel = FALSE, log_file = NULL, constrain_mu = FALSE)
  expect_equal(nrow(tca.mdl$gammas_hat_pvals) == m & ncol(tca.mdl$gammas_hat_pvals) == p1*k, TRUE)
  expect_equal(nrow(tca.mdl$deltas_hat_pvals) == m & ncol(tca.mdl$deltas_hat_pvals) == p2, TRUE)

})



test_that("Verify that setting vars.mle to either TRUE or FALSE provides similar results", {

  skip_on_cran()

  n = 500; m = 30; k = 3; tau = 0.01; p1 = 2; p2 = 2;
  res <- test_data(n, m, k, p1, p2, tau, log_file = NULL)

  tca.mdl1 <- tca(X = res$X, W = res$W, C1 = res$C1, C2 = res$C2, refit_W = FALSE, parallel = FALSE, log_file = NULL, constrain_mu = TRUE, vars.mle = TRUE)
  tca.mdl2 <- tca(X = res$X, W = res$W, C1 = res$C1, C2 = res$C2, refit_W = FALSE, parallel = FALSE, log_file = NULL, constrain_mu = TRUE, vars.mle = FALSE)

  parameters <- c("mus_hat","sigmas_hat","deltas_hat","gammas_hat")
  for (params in parameters){
    expect_equal( cor(Reshape(tca.mdl1[[params]],ncol(tca.mdl1[[params]])*nrow(tca.mdl1[[params]])), Reshape(tca.mdl2[[params]],ncol(tca.mdl2[[params]])*nrow(tca.mdl2[[params]])) )[1] > 0.9, TRUE)
  }

  # do not constrain mu this time
  tca.mdl1 <- tca(X = res$X, W = res$W, C1 = res$C1, C2 = res$C2, refit_W = FALSE, parallel = FALSE, log_file = NULL, constrain_mu = FALSE, vars.mle = TRUE)
  tca.mdl2 <- tca(X = res$X, W = res$W, C1 = res$C1, C2 = res$C2, refit_W = FALSE, parallel = FALSE, log_file = NULL, constrain_mu = FALSE, vars.mle = FALSE)

  for (params in parameters){
    expect_equal( cor(Reshape(tca.mdl1[[params]],ncol(tca.mdl1[[params]])*nrow(tca.mdl1[[params]])), Reshape(tca.mdl2[[params]],ncol(tca.mdl2[[params]])*nrow(tca.mdl2[[params]])) )[1] > 0.9, TRUE)
  }

})

