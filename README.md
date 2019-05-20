# Tensor Composition Analysis (TCA)

Tensor Composition Analysis (TCA) allows the deconvolution of two-dimensional data (features by observations) coming from a mixture of sources into a three-dimensional matrix of signals (features by observations by sources). TCA further allows to test the features in the data for different statistical relations with an outcome of interest while modeling source-specific effects (TCA regression); particularly, it allows to look for statistical relations between source-specific signals and an outcome.

In the context of DNA methylation data, TCA can deconvolve tissue-level bulk methylation (methylation sites by individuals) into a tensor of cell-type-specific methylation levels for each individual (methylation sites by individuals by cell types) and it allows to detect cell-type-specific relations (associations) with an outcome of interest. For more details see Rahmani et al. (2019)<sup>[1](#myfootnote1)</sup>.

TCA is available in both R and Matlab. Note that the Matlab version was used for deriving the results in the publication describing TCA<sup>[1](#myfootnote1)</sup>. 

## R version

[![Travis build status](https://travis-ci.com/cozygene/TCA.svg?branch=master)](https://travis-ci.org/cozygene/TCA)
 
**The R package of TCA will soon be available on CRAN**

The full documentation of the TCA R package can be found <a href="https://github.com/cozygene/TCA/blob/master/manual.pdf" target="_blank">here</a>.

You can also find a full working example of TCA in this vignette about cell-type-specific resolution epigenetics using TCA in R.

<!--describe the config file.-->

## Matlab version

implemented and tested using Matlab R2015b

### Usage

TCA requires cell-type proportion estimates for the samples in the data. These can be obtained by either using the reference-based model by Houseman et al. 2012<sup>[2](#myfootnote2)</sup> (see an implementation <a href="http://glint-epigenetics.readthedocs.io/" target="_blank">here</a>) or using the semi-supervised model by Rahmani et al. 2018<sup>[3](#myfootnote3)</sup> (does not require reference data; see an implementation <a href="https://github.com/cozygene/BayesCCE" target="_blank">here</a>).

There are two main functions in this distribution. A full documentation of the input arguments and output values of these functions is provided in the headers of these function.
* **TCA_EWAS.m** - for performing cell-type-specific EWAS under the TCA model.

* **TCA.m** - for estimating cell-type-specific methylation levels (in case only these estimates are desired rather than performing a cell-type-specific EWAS). The required input arguments for this function can be estimated using the function TCA_fit_model.m.

### Demo

We provide small simulated demo data, wherein the phenotype is associated with the last site in the data matrix. For performing EWAS on the demo files following the TCA model, execute in matlab the following commands from the 'demo' directory:
```matlab
% <matlab code>
% Add the .m files to the path
addpath '../'
% Read the data files
y = dlmread('demo_y.txt');  % phenotype
X = dlmread('demo_X.txt');  % methylation matrix (individuals by sites)
W = dlmread('demo_W.txt');  % proportions matrix (individuals by cell types)
% Fit the parameters of the TCA model
[W,mus_hat,sigmas_hat,tau_hat] = TCA_fit_model(X, W);
% Perform EWAS under the TCA model with cell-type-specific effects
pvals = TCA_EWAS(y, X, W, mus_hat, sigmas_hat, tau_hat);
```

---

### License

Both the R and Matlab version of TCA are available under the <a href="https://opensource.org/licenses/GPL-3.0" target="_blank">GPL-3 license</a>.

<!---
#### Citing TCA

If you use TCA in any published work, please cite the manuscript describing the method:

Elior Rahmani, Regev Schweiger, Brooke Rhead, Lindsey A. Criswell, Lisa F. Barcellos, Eleazar Eskin, Saharon Rosset, Sriram Sankararaman, and Eran Halperin. *bioRxiv*, 2018.
-->


### Author

This software was developed by Elior Rahmani (elior.rahmani@gmail.com).

### Bug reports
Please <a href="https://github.com/cozygene/TCA/issues/" target="_blank">open an issue</a>) for reporting bugs. If you are reporting bugs with the R version, please make sure to set the argument debug to TRUE and attach your log.


___

<a name="myfootnote1">1</a>: Rahmani et al. "Cell-type-specific resolution epigenetics without the need for cell sorting or single-cell biology." Nature Communications, in press (2019).

<a name="myfootnote2">2</a>: Houseman et al. "DNA methylation arrays as surrogate measures of cell mixture distribution." BMC bioinformatics (2012).

<a name="myfootnote3">3</a>: Rahmani et al. "BayesCCE: a Bayesian framework for estimating cell-type composition from DNA methylation without the need for methylation reference." Genome biology (2018).
