library("TCA")

context("Test ReFACTor")

test_that("Compare the result of refactor with those of the matlab version of refactor", {

  skip_on_cran()

  basedir <- "../assets/"

  # load data and the matlab results
  X <- as.matrix(read.table(paste(basedir,"X.txt",sep=""), header = FALSE, sep=","))
  rownames(X) <- 1:nrow(X)
  C <- as.matrix(read.table(paste(basedir,"C2.txt",sep=""), header = FALSE, sep=","))
  t <- 50
  k <- 3
  matlab.ref_comp <- as.matrix(read.table(paste(basedir,"refactor.R_est.txt",sep=""), header = FALSE, sep=","))
  matlab.ranked_list <- as.matrix(read.table(paste(basedir,"refactor.ranked_list.txt",sep=""), header = FALSE, sep=","))
  matlab.ref_comp.adj <- as.matrix(read.table(paste(basedir,"refactor.R_est.adjusted.txt",sep=""), header = FALSE, sep=","))
  matlab.ref_comp.adj.C.remove <- as.matrix(read.table(paste(basedir,"refactor.R_est.adjusted.C.remove.txt",sep=""), header = FALSE, sep=","))
  matlab.ranked_list.adj <- as.matrix(read.table(paste(basedir,"refactor.ranked_list.adjusted.txt",sep=""), header = FALSE, sep=","))

  # run refactor
  ref1 <- refactor(t(X), k, sparsity = t, C = NULL, sd_threshold = 0, num_comp = ncol(matlab.ref_comp), debug = TRUE, log_file = NULL)
  ref2 <- refactor(t(X), k, sparsity = t, C = C, C.remove = FALSE, sd_threshold = 0, num_comp = ncol(matlab.ref_comp), debug = TRUE, log_file = NULL)
  ref3 <- refactor(t(X), k, sparsity = t, C = C, C.remove = TRUE, sd_threshold = 0, num_comp = ncol(matlab.ref_comp), debug = TRUE, log_file = NULL)

  # evaluating the feature ranking
  expect_equal(sum(as.numeric(substring(ref1[["ranked_list"]], 2)) == matlab.ranked_list), ncol(X))
  expect_equal(sum(as.numeric(substring(ref2[["ranked_list"]], 2)) == matlab.ranked_list.adj), ncol(X))
  expect_equal(sum(as.numeric(substring(ref3[["ranked_list"]], 2)) == matlab.ranked_list.adj), ncol(X))

  # evaluate correlation of the refactor components
  for (h in 1:ncol(matlab.ref_comp)){
    expect_equal(abs(cor(matlab.ref_comp[,h],ref1[["scores"]][,h])) > 0.99, TRUE)

    expect_equal(abs(cor(matlab.ref_comp.adj.C.remove[,h],ref3[["scores"]][,h])) > 0.99, TRUE)
  }
  # the version with C.remove=TRUE has more differences compared with the matlab version - evaluate the the first two components only.
  for (h in 1:2){
    expect_equal(abs(cor(matlab.ref_comp.adj[,h],ref2[["scores"]][,h])) > 0.99, TRUE)
  }

})


test_that("Verify the C.remove argument", {

  skip_on_cran()

  n = 100; m = 50; k = 3; tau = 0.01; p1 = 0; p2 = 2;
  res <- test_data(n, m, k, p1, p2, tau, log_file = NULL)
  X <- res$X
  C2 <- res$C2

  ref1 <- refactor(X, k, sparsity = 10, C = C2, C.remove = TRUE, debug = TRUE, log_file = NULL)
  expect_equal(sum(abs(sort(cor(ref1$scores,C2))) > 10e-14) == 0, TRUE)

  ref2 <- refactor(X, k, sparsity = 10, C = C2, C.remove = FALSE, debug = TRUE, log_file = NULL)
  expect_equal(sum(abs(sort(cor(ref2$scores,C2))) > 10e-14) > 0, TRUE)

})

