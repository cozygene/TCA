library("TCA")

context("Test tcasub")

test_that("subset the results of tca", {
  n = 20; m = 30; k = 3; tau = 0.01; p1 = 0; p2 = 0;
  res <- test_data(n, m, k, p1, p2, tau, log_file = NULL)
  X <- res$X
  W <- res$W
  mus <- res$mus
  sigmas <- res$sigmas
  tca.mdl <- tca(X, W, refit_W = FALSE, parallel = FALSE, log_file = NULL)
  remove <- c("1","3","5")
  features <- setdiff(as.character(1:m),remove)
  tca.mdl.subset <- tcasub(tca.mdl, features, log_file = NULL)
  expect_equal(sum(rownames(tca.mdl.subset$mus_hat) == features) == length(features), TRUE)
})
