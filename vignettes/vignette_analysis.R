library("TCA")
library("ggplot2")
library("ggpubr")
library("pracma")
library("matrixStats")

prep_data <- function(data_path){

  file_name1 <- paste(data_path,.Platform$file.sep,"hannum.chr22.RData",sep="")
  file_name2 <- paste(data_path,.Platform$file.sep,"liu.cd4.chr22.RData",sep="")
  file_name3 <- paste(data_path,.Platform$file.sep,"paquette.chr22.RData",sep="")

  file_names <- c(file_name1, file_name2, file_name3)

  if(any(!file.exists(file_name1, file_name2, file_name3))){
    library("GEOquery")
    library("data.table")
    library("EpiDISH")
  }

  if (!file.exists(file_name1)){

    # Download the Hannum et al. data
    gse <- GEOquery::getGEO("GSE40279", destdir = data_path, GSEMatrix = TRUE)

    # Extract methylation data
    X.hannum <- Biobase::exprs(gse[[1]])

    # covariates; note that there's also 'ethnicity' covariate in the data, however, it is perfectly captured by the plateinformation
    plate.hannum <- as.numeric(as.factor((Biobase::pData(gse[[1]])[["plate:ch1"]])))
    cov.hannum <- cbind(as.numeric((Biobase::pData(gse[[1]])[["age (y):ch1"]])), as.numeric(as.factor((Biobase::pData(gse[[1]])[["gender:ch1"]]))), indicator_vars(plate.hannum))
    rownames(cov.hannum) <- colnames(X.hannum)
    colnames(cov.hannum) <- c("age", "gender", "plate1", "plate2", "plate3", "plate4", "plate5", "plate6", "plate7", "plate8")

    # Calculate principal components from control probes that reflect technical variability; since we are not working with IDAT files we don't have "real" control probes - instead, we use low variance sites that are not expected to capture any true biological variability.
    low_var_pcs.hannum <- low_var_pcs(X.hannum, rank = 10, p = 1000)
    cov.hannum <- cbind(cov.hannum, low_var_pcs.hannum)

    # estimate cell-type proportions usign a reference-based approach
    ref <- as.matrix(EpiDISH::centDHSbloodDMC.m[,c("Neutro","Eosino","CD4T","CD8T","Mono","B","NK")])
    W.hannum <- EpiDISH::epidish(X.hannum, ref)$estF
    # merge Neutro and Eosino
    W.hannum.gran <- W.hannum[,"Neutro",drop=F] + W.hannum[,"Eosino",drop=F]
    colnames(W.hannum.gran) <- "Gran"
    W.hannum <- cbind(W.hannum.gran,W.hannum[,c("CD4T","CD8T","Mono","B","NK")])

    # remove polymorphic sites, cross-reactive sites, and non-autosomal sites according to Chen et al.; in addition, remove low variance sites, as those are unlikely to demonstrate biological signal.
    X.hannum.processed <- filter_data(X.hannum)

    # keep only sites on chromosome 22 (so that tutorial can run quickly)
    map <- read.table("https://github.com/cozygene/glint/blob/master/utils/assets/HumanMethylationSites?raw=true",sep=",",row.names = 1)
    chrs <- map[rownames(X.hannum.processed),1]
    chr22sites <- which(chrs == 22)
    X.hannum.processed.22 <- X.hannum.processed[chr22sites,]

    hannum <- list(X = X.hannum.processed.22, cov = cov.hannum, W = W.hannum)
    save(hannum, file = file_name1, compress = "bzip2")

    rm(hannum, gse, X.hannum, X.hannum.processed)
    file.remove(paste(data_path,.Platform$file.sep,"GPL13534.soft",sep=""))
    file.remove(paste(data_path,.Platform$file.sep,"GSE40279_series_matrix.txt.gz",sep=""))

  }

  if (!file.exists(file_name2)){

    # Download the Liu et al. CD4  data

    # covariates
    gse <- GEOquery::getGEO("GSE56581", destdir = data_path, GSEMatrix = TRUE)
    liu.cd4.age <- as.matrix(as.numeric(Biobase::pData(gse[[1]])[["age (yrs):ch1"]]))
    colnames(liu.cd4.age) <- "age"

    # methylation data
    liu.cd4.methfile <- paste(data_path,.Platform$file.sep,"GSE56581_methylome_normalized.txt",sep="")
    download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56581/suppl/GSE56581_methylome_normalized.txt.gz", paste(liu.cd4.methfile,".gz",sep=""))
    gunzip(paste(liu.cd4.methfile,".gz",sep=""))
    X.liu.cd4 <- fread(liu.cd4.methfile)

    X.liu.cd4.rownames <- X.liu.cd4$ID_REF
    X.liu.cd4 <- as.matrix(X.liu.cd4[,2:ncol(X.liu.cd4)])
    rownames(X.liu.cd4) <- X.liu.cd4.rownames
    X.liu.cd4 <- X.liu.cd4[,seq(1,ncol(X.liu.cd4),2)]

    # Calculate principal components from control probes that reflect technical variability; since we are not working with IDAT files we don't have "real" control probes - instead, we use low variance sites that are not expected to capture any true biological variability.
    low_var_pcs.liu.cd4 <- low_var_pcs(X.liu.cd4, rank = 10)
    cov.liu.cd4 <- cbind(liu.cd4.age, low_var_pcs.liu.cd4)

    # keep only sites on chromosome 22 (so that tutorial can run quickly)
    map <- read.table("https://github.com/cozygene/glint/blob/master/utils/assets/HumanMethylationSites?raw=true",sep=",",row.names = 1)
    chrs <- map[rownames(X.liu.cd4),1]
    chr22sites <- which(chrs == 22)
    X.liu.cd4.22 <- X.liu.cd4[chr22sites,]

    liu.cd4 <- list(X = X.liu.cd4.22, cov = cov.liu.cd4)
    save(liu.cd4, file = file_name2)

    rm(liu.cd4, gse, X.liu.cd4.22, X.liu.cd4)
    file.remove(paste(data_path,.Platform$file.sep,"GSE56581_methylome_normalized.txt",sep=""))
    file.remove(paste(data_path,.Platform$file.sep,"GSE56581_series_matrix.txt.gz",sep=""))

  }

  if (!file.exists(file_name3)){
    # Paquette et al. (Epigenetics 2016)
    gse <- GEOquery::getGEO("GSE75248", destdir = data_path, GSEMatrix = TRUE)

    # methylation data
    X.paquette <- Biobase::exprs(gse[[1]])

    # Calculate principal components from control probes that reflect technical variability; since we are not working with IDAT files we don't have "real" control probes - instead, we use low variance sites that are not expected to capture any true biological variability.
    low_var_pcs.paquette <- low_var_pcs(X.paquette, rank = 10, p = 1000)

    # remove polymorphic sites, cross-reactive sites, and non-autosomal sites according to Chen et al.; in addition, remove low variance sites, as those are unlikely to demonstrate biological signal.
    X.paquette.processed <- filter_data(X.paquette)

    # keep only sites on chromosome 22 (so that tutorial can run quickly); Paquette et al. reported some signal in chromosome 16
    map <- read.table("https://github.com/cozygene/glint/blob/master/utils/assets/HumanMethylationSites?raw=true",sep=",",row.names = 1)
    chrs <- as.character(map[rownames(X.paquette.processed),1])
    chr22sites <- which(chrs == "22")
    X.paquette.processed.22 <- X.paquette.processed[chr22sites,]

    # Extract covaraites
    gestational_age <-as.numeric(unlist(lapply(Biobase::pData(gse[[1]])[["gestational age:ch1"]], function(x) strsplit(x,";")[[1]][1])))
    arousal <- as.numeric(Biobase::pData(gse[[1]])[["arousal:ch1"]])
    batch <- as.numeric(as.factor((Biobase::pData(gse[[1]])[["batch:ch1"]])))
    batch[is.na(batch)] <- 3
    batch <- indicator_vars(batch)
    # since gender information is not available on the GEO record, get gender information that was inferred from the IDAR files based on X chromosome methylation
    gender <- read.table("https://raw.githubusercontent.com/cozygene/TCA/master/vignettes/gender.paquette.txt", row.names = 1)
    cov.paquette <- data.frame(gender,gestational_age,arousal,batch)
    colnames(cov.paquette) <- c("gender","gestational_age","arousal","batch1","batch2")
    rownames(cov.paquette) <- colnames(X.paquette.processed.22)
    # remove na values
    keep <- setdiff(1:nrow(cov.paquette) ,which(rowSums(is.na(cov.paquette)) > 0))
    cov.paquette <- cbind(cov.paquette[keep,], low_var_pcs.paquette[keep,])

    X.paquette.processed.22 <- X.paquette.processed.22[,keep]

    # Load cell-type proportions usign a reference-based approach; these were estimated using the raw IDAT files of the data (available on GEO) with the package ENmix.
    W.paquette <- read.table("https://raw.githubusercontent.com/cozygene/TCA/master/vignettes/W.paquette.txt", row.names = 1)
    W.paquette <- W.paquette[,2:ncol(W.paquette)]
    W.paquette[W.paquette < 0] <- 0
    # remove lowly abundant cell types
    W.paquette <- W.paquette[keep, colMeans(W.paquette)>=0.01]
    rownames(W.paquette) <- rownames(cov.paquette)
    # normalize cell type proportions (output from the ENmix package didn't require the proportions of an individual to sum up to 1)
    W.paquette <- W.paquette/t(repmat(rowSums(W.paquette),ncol(W.paquette),1))

    # get the set of cord blood reference CpGs that were used for estimating cell-type proportions; we will use these CpGs for re-estimating W
    library("FlowSorted.CordBlood.450k")
    ref.cord <- get("FlowSorted.CordBlood.450k")
    ref.cord <- preprocessRaw(ref.cord)
    ref.cord.cpgs <- rownames(FlowSorted.CordBlood.450k.ModelPars)
    ref.cord.cpgs.intersect <- intersect(ref.cord.cpgs, rownames(X.paquette))

    paquette <- list(X = X.paquette.processed.22, cov = cov.paquette, W = W.paquette, X.ref_cpgs = X.paquette[ref.cord.cpgs.intersect,keep])
    save(paquette, file = file_name3)

    rm(paquette, gse, X.paquette, X.paquette.processed, X.paquette.processed.22)
    file.remove(paste(data_path,.Platform$file.sep,"GSE75248_series_matrix.txt.gz",sep=""))

    ## note - the files W.paquette.txt and gender.paquette.txt were generated using the original IDAT files of the Paquette data.
    # Gender information was inferred from methylation levels frmo the X chromosome as it is not available on Tthe GEO record.
    # Code attched below; in order to run it, download the IDAT files from GEO accession number GSE75248 and set 'path' in 'readidat' below to the location of the files.
    # require("ENmix")
    # rgdat <- readidat(path = ".", manifestfile=NULL, recursive = FALSE, verbose = TRUE)
    # meth <- getmeth(rgdat); rm(rgdat)
    # W <- estimateCellProp(meth, refdata="FlowSorted.CordBlood.450k", cellTypes=NULL, nonnegative = TRUE, nProbes=50, normalize=TRUE, refplot=FALSE)
    # write.table(W, file = "W.paquette.txt", quote = FALSE, sep = " ")
    # beta <- getB(meth)
    # map <- read.table("https://github.com/cozygene/glint/blob/master/utils/assets/HumanMethylationSites?raw=true",sep=",",row.names = 1)
    # chrs <- as.character(map[rownames(beta),1])
    # sites.X <- which(chrs == "X")
    # sex_cpgs.pca.paquette <- prcomp(t(beta[sites.X,]), center=TRUE, scale=TRUE, rank = 2)
    # write.table(sex_cpgs.pca.paquette$x[,1:2], file = "sex_cpgs.pca.paquette.txt", quote = FALSE, sep = " ")
    # plot(sex_cpgs.pca.paquette$x[,1],sex_cpgs.pca.paquette$x[,2]) # shows a separation between males and females; note that the males/females ratio here perfectly match the one in the Paquette et al. paper.
    # df <-data.frame(as.numeric(x[,1] > 0))
    # colnames(df) <-"gender"; rownames(df) <- rownames(x)
    # write.table(df,file = "gender.paquette.txt", quote = FALSE, sep = " ")

  }

  return(file_names)

}


indicator_vars <- function(x){
  x.indicators <- matrix(0,length(x),length(unique(x))-1)
  counter <- 1
  u <- unique(x)
  for (i in setdiff(u,u[length(u)])){
    x.indicators[,counter] <- as.numeric(x == i)
    counter <- counter + 1
  }
  return(x.indicators)
}

low_var_pcs <- function(X, rank = 10, p = 1000){
  site.variances <- matrixStats::rowVars(X)
  names(site.variances) <- rownames(X)
  low.var.sites <- names(head(sort(site.variances), p))
  low.var.pca <- prcomp(t(X[low.var.sites,]), center=TRUE, scale=TRUE, rank = rank)
  return(low.var.pca$x)
}

filter_data <- function(X, var_th = 0.0001){
  nonspecific_probes <- read.table("https://raw.githubusercontent.com/cozygene/glint/master/parsers/assets/nonspecific_probes.txt")[,1]
  XY_probes <- read.table("https://raw.githubusercontent.com/cozygene/glint/master/parsers/assets/HumanMethylationSites_X_Y.txt")[,1]
  polymorphic_probes <- read.table("https://raw.githubusercontent.com/cozygene/glint/master/parsers/assets/polymorphic_cpgs.txt")[,1]
  # remove sites with very low variance that are unlikely to exhibit biological signal
  site.variances <- matrixStats::rowVars(X)
  names(site.variances) <- rownames(X)
  low_var_sites <- names(which(site.variances < var_th))
  exclude <- union(nonspecific_probes,union(XY_probes,union(low_var_sites,polymorphic_probes)))
  return(X[setdiff(rownames(X),exclude),])
}

plot_qq <- function(pvals, labels, ggarrange.nrow = 1, ggarrange.ncol = 1, alpha = 0.05, experiment_wide_line = TRUE){
  significance_th <- list(alpha/length(pvals[[1]]))
  if(length(pvals)-1) significance_th[[2]] <- alpha/(length(pvals)*length(pvals[[1]]))
  qqplots <- lapply(1:length(pvals), function(p){
    df <- data.frame(pvals.obs = -log10(sort(pvals[[p]])), pvals.exp = -log10(sort((1:length(pvals[[1]]))/length(pvals[[1]]))));
    qqplot <- ggplot(df, aes(x = pvals.exp, y = pvals.obs)) +
      stat_binhex(geom = "point", bins=1000, size=1) +
      geom_abline() +
      ggtitle(labels[p]) +
      xlab(expression(Expected~-log[10](P))) + ylab(expression(Observed~-log[10](P))) +
      theme_bw() +
      guides(fill="none") +
      geom_hline(yintercept=-log10(significance_th[[1]]), linetype="dashed", color = "red", size=1)
    if (length(significance_th)-1 & experiment_wide_line) qqplot <- qqplot + geom_hline(yintercept=-log10(significance_th[[2]]), linetype="dashed", color = "red", size=0.5)
    return(qqplot)
  })
  ggarrange(plotlist = qqplots, ncol = ggarrange.ncol, nrow = ggarrange.nrow)
}

plot_scatter <- function(dfs, ggarrange.ncol, ggarrange.nrow, xlab, ylab, titles){
  plots <- vector("list", length = length(dfs))
  for (i in 1:length(dfs)){
    df <- data.frame(y = dfs[[i]]$y, x = dfs[[i]]$x)
    plots[[i]] <- ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method=lm) +
      stat_cor(method = "pearson", colour = "blue") +
      xlab(xlab) +
      ylab(ylab) +
      ggtitle(titles[i])
  }
  ggarrange(plotlist = plots, ncol = ggarrange.ncol, nrow = ggarrange.nrow)
}


if(FALSE){

  # Set a path for storing data
  data_path <- "./"

  # Load the data
  filenames <- prep_data(data_path)
  for (filename in filenames) load(filename)

  ## Experiment #1: detecting CD4 differential methylation with age; working under the assumption X|Y (i.e. age affects methylation levels)

  # Apply the TCA model to the Hannum whole-blood data
  tca.mdl.hannum <- tca(X = hannum$X,
                        W = hannum$W,
                        C1 = hannum$cov[,c("gender","age")],
                        C2 = hannum$cov[,3:ncol(hannum$cov)])

  # Extract p-values of a joint test
  tca.mdl.hannum.pvals.joint <- tca.mdl.hannum$gammas_hat_pvals.joint[,"age"]

  # Extract p-values for each cell type, under a marginal conditional test
  tca.mdl.hannum.pvals.marg_cond <- tca.mdl.hannum$gammas_hat_pvals[,paste(colnames(hannum$W),".age",sep="")]

  # qq-plots - for the p-values of the joint test, and for the p-values in CD4, under a marginal conditional test
  plot_qq(list(tca.mdl.hannum.pvals.joint, tca.mdl.hannum.pvals.marg_cond[,"CD4T.age"]),
          labels = c("Joint test with age", "CD4 marginal conditional test with age"),
          ggarrange.nrow = 1,
          ggarrange.ncol = 2,
          experiment_wide_line = FALSE)

  # Run ReFACTor to capture more of the variation of cell-type composition
  refactor.mdl.hannum <- refactor(X = hannum$X,
                                  k = 6,
                                  C = hannum$cov[,3:ncol(hannum$cov)])

  # Rerun TCA, this time include the ReFACTor components as additional covariates
  tca.mdl.hannum.2 <- tca(X = hannum$X,
                          W = hannum$W,
                          C1 = hannum$cov[,c("gender","age")],
                          C2 = cbind(hannum$cov[,3:ncol(hannum$cov)],refactor.mdl.hannum$scores))

  # Extract the updated p-values of a joint test
  tca.mdl.hannum.2.pvals.joint <- tca.mdl.hannum.2$gammas_hat_pvals.joint[,"age"]

  # Extract the updated marginal conditional p-values
  tca.mdl.hannum.2.pvals.marg_cond <- tca.mdl.hannum.2$gammas_hat_pvals[,paste(colnames(hannum$W),".age",sep="")]

  # qq-plots - for the p-values of the joint test, and for the p-values in CD4, under a marginal conditional test
  plot_qq(list(tca.mdl.hannum.2.pvals.joint, tca.mdl.hannum.2.pvals.marg_cond[,"CD4T.age"]),
          labels = c("Joint test with age", "CD4 marginal conditional test with age"),
          ggarrange.nrow = 1,
          ggarrange.ncol = 2,
          experiment_wide_line = FALSE)

  # Extract the hits based on the joint test
  hits.joint <- names(which(tca.mdl.hannum.2.pvals.joint < 0.05/nrow(hannum$X)))

  # Extract the hits from hits.joint where CD4 cells demonstrate the lowest p-value across all cell types
  cd4.hits <- names(which(tca.mdl.hannum.2.pvals.marg_cond[hits.joint,"CD4T.age"] == rowMins(tca.mdl.hannum.2.pvals.marg_cond[hits.joint,])))

  sprintf("Detected %d associations using a joint test, %d associations using a marginal conditional test, and %d associations in CD4 cells using a marginal conditional test.",
          sum(tca.mdl.hannum.2.pvals.joint <= 0.05/nrow(hannum$X)),
          sum(tca.mdl.hannum.2.pvals.marg_cond <= 0.05/(nrow(hannum$X)*nrow(hannum$W))),
          sum(tca.mdl.hannum.2.pvals.marg_cond[,"CD4T.age"] <= 0.05/(nrow(hannum$X)*nrow(hannum$W))))

  # Replicate the CD4 hits in the Liu purified CD4 data
  cd4.hits.liu.pvals <- unlist(lapply(1:length(cd4.hits),
                                      function(x) summary(lm(y~., data.frame(y=liu.cd4$X[cd4.hits[x],], liu.cd4$cov)))$coefficients["age","Pr(>|t|)"]))

  # Plot adjusted methylation (i.e. adjusted for covariates) as a function of age in all four replicated CD4 associations
  dfs <-  vector("list", length = 4)
  for (i in 1:4){
    r <- scale(residuals(lm(y~., data.frame(y=liu.cd4$X[cd4.hits[i],], liu.cd4$cov[,2:ncol(liu.cd4$cov)]))))
    dfs[[i]] <- data.frame(x = liu.cd4$cov[,"age"], y = r)
  }
  plot_scatter(dfs = dfs,
               ggarrange.ncol = 2,
               ggarrange.nrow = 2,
               xlab = "Age",
               ylab = "Adjusted methylation level",
               titles = paste("CD4 methylation in ",cd4.hits,sep=""))

  # Verify that p-values are calibrated in the Liu data
  liu.cd4.regression.pvals <- unlist(lapply(1:nrow(liu.cd4$X),
                                            function(x) summary(lm(y~.,data.frame(y = liu.cd4$X[x,],liu.cd4$cov)))$coefficients["age","Pr(>|t|)"]))
  plot_qq(list(liu.cd4.regression.pvals), labels = "Linear regression (sorted CD4)")





  ## Experiment #2: detecting differential methylation with infant arousal in granulocytes; working under the assumption Y|X (i.e. methylation levels affect or mediate components affecting infant arousal).

  # Apply the TCA model to the Paquette data
  tca.mdl.paquette <- tca(X = paquette$X,
                          W = paquette$W,
                          C1 = paquette$cov[,c("gender","gestational_age")],
                          C2 = paquette$cov[,4:ncol(paquette$cov)],
                          constrain_mu = TRUE)

  # Run tcareg with a joint test and extract p-values; generate a qq-plot
  C3_names <- c("gender","gestational_age","batch1","batch2")
  tcareg.mdl.paquette.joint <- tcareg(X = paquette$X,
                             tca.mdl = tca.mdl.paquette,
                             y = paquette$cov[,"arousal",drop=F],
                             C3 = paquette$cov[,C3_names],
                             test = "joint")
  plot_qq(list(tcareg.mdl.paquette.joint$pvals), labels = "Joint test with infant arousal")

  # Since there is an inflation in the qq-plot we try to re-estimate the cell-type proportions
  # First, run TCA using only the reference sites for getting a new estimate of W
  tca.mdl.paquette.refit_W <- tca(X = paquette$X.ref_cpgs,
                                  W = paquette$W,
                                  C1 = paquette$cov[,c("gender","gestational_age")],
                                  C2 = paquette$cov[,4:ncol(paquette$cov)],
                                  constrain_mu = TRUE,
                                  refit_W = TRUE,
                                  refit_W.features = rownames(paquette$X.ref_cpgs))
  # Use the updated estimate of W in a new execution of TCA on the data
  tca.mdl.paquette.2 <- tca(X = paquette$X,
                            W = tca.mdl.paquette.refit_W$W,
                            C1 = paquette$cov[,c("gender","gestational_age")],
                            C2 = paquette$cov[,4:ncol(paquette$cov)],
                            constrain_mu = TRUE)

  # Run tcareg with a joint test and extract p-values; generate a qq-plot
  tcareg.mdl.paquette.joint.2 <- tcareg(X = paquette$X,
                               tca.mdl = tca.mdl.paquette.2,
                               y = paquette$cov[,"arousal",drop=F],
                               C3 = paquette$cov[,C3_names],
                               test = "joint")
  plot_qq(list(tcareg.mdl.paquette.joint.2$pvals), labels = "Joint test with infant arousal")

  # Run tcareg with a marginal conditional test; generate qq plots
  tcareg.mdl.paquette.2.marg_cond <- tcareg(X = paquette$X,
                                         tca.mdl = tca.mdl.paquette.2,
                                         y = paquette$cov[,"arousal",drop=F],
                                         C3 = paquette$cov[,C3_names],
                                         test = "marginal_conditional")
  plot_qq(split(tcareg.mdl.paquette.2.marg_cond$pvals,rep(1:ncol(paquette$W), each = nrow(paquette$X))),
          labels = paste(colnames(paquette$W)," marginal conditional with arousal", sep=""),
          ggarrange.nrow = 2,
          ggarrange.ncol = 2)

  # Exctract the hit we got from the joint test
  hit.joint <- rownames(paquette$X)[which(tcareg.mdl.paquette.joint.2$pvals < 0.05/nrow(paquette$X))]

  # The p-values of the marginal conditional test in the hit we found with the joint test suggest that the association is in granulocytes
  tcareg.mdl.paquette.2.marg_cond$pvals[hit.joint,]

  # Extract from the tca model the part that is relevant for the hit found
  tcasub.mdl.paquette.2 <- tcasub(tca.mdl = tca.mdl.paquette.2,
                                  features = hit.joint)

  # Calculate cell-type-specific methylation for the samples in the detected CpG
  tensor.mdl.hit <- tensor(tca.mdl = tcasub.mdl.paquette.2,
                           X = paquette$X[hit.joint,,drop=F])

  # Plot the estimated granulocyte-specific methylation with arasual; adjust methylation levels to covariates
  # Compare with using the bulk data
  r.gran <- scale(residuals(lm(y~., data.frame(y=tensor.mdl.hit[[2]][1,], paquette$cov[,setdiff(1:ncol(paquette$cov),3)]))))
  df.gran <- data.frame(y = r.gran, x = paquette$cov[,"arousal"])
  # The bulk data should be further adjusted for cell-type composition (i.e. tca.mdl.paquette.2$W)
  r.bulk <- scale(residuals(lm(y~., data.frame(y=paquette$X[hit.joint,], tca.mdl.paquette.2$W, paquette$cov[,setdiff(1:ncol(paquette$cov),3)]))))
  df.bulk <- data.frame(y = r.bulk, x = paquette$cov[,"arousal"])
  plot_scatter(dfs = list(df.bulk, df.gran),
               ggarrange.ncol = 2,
               ggarrange.nrow = 1,
               xlab = "Infant arousal",
               ylab = "Adjusted methylation level",
               titles = paste(c("Whole-blood methylation in ", "Granulocyte methylation in "),hit.joint,sep=""))

}
