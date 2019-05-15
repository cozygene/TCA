library("TCA")

context("Test test_data")

test_that("generate test data", {

  n = 200; m = 100; k = 3; tau = 0.01; p1 = 2; p2 = 2;
  res <- test_data(n, m, k, p1, p2, tau)

  X <- res$X
  W <- res$W
  mus <- res$mus
  sigmas <- res$sigmas
  C1 <- res$C1
  C2 <- res$C2
  gammas <- res$gammas
  deltas <- res$deltas

  tca.mdl <- tca(X, W, C1 = C1, C2 = C2, refit_W = TRUE, refit_t = nrow(X), parallel = FALSE, log_file = NULL)
  Z_hat <- tensor(X, tca.mdl)
  expect_is(Z_hat[[1]], "matrix")
  expect_is(Z_hat[[2]], "matrix")
  expect_is(Z_hat[[3]], "matrix")

})
