library("TCA")

context("Test ReFACTor")

# test_that("Compare the result of refactor with those of the matlab version of refactor", {
#
#   skip_on_cran()
#
#   basedir <- "../assets/"
#
#   # load data and the matlab results
#   X <- as.matrix(read.table(paste(basedir,"X.txt",sep=""), header = FALSE, sep=","))
#   rownames(X) <- 1:nrow(X)
#   C <- as.matrix(read.table(paste(basedir,"C2.txt",sep=""), header = FALSE, sep=","))
#   t <- 50
#   k <- 3
#   matlab.ref_comp <- as.matrix(read.table(paste(basedir,"refactor.R_est.txt",sep=""), header = FALSE, sep=","))
#   matlab.ranked_list <- as.matrix(read.table(paste(basedir,"refactor.ranked_list.txt",sep=""), header = FALSE, sep=","))
#
#   # run refactor
#   ref <- refactor(t(X), k, sparsity = t, C = C, sd_threshold = 0, num_comp = ncol(matlab.ref_comp), debug = TRUE, log_file = NULL)
#
#   # evaluating the feature ranking
#   expect_equal(sum(as.numeric(substring(ref[["ranked_list"]], 2)) == matlab.ranked_list), ncol(X))
#
#   # evaluate correlation of the refactor components
#   for (h in 1:2){
#     expect_equal(abs(cor(matlab.ref_comp[,h],ref[["scores"]][,h])) > 0.99, TRUE)
#   }
#
# })
