# Tensor Composition Analysis (TCA)

Tensor Composition Analysis (TCA)<sup>[1](#myfootnote1)</sup> is a method for estimating cell-type-specific methylation levels and performing cell-type-specific epigenetic association studies (EWAS) using bulk methylation data collected from heterogeneous source.

Currently, we only provide a Matlab implementation of the method (implemented and tested using Matlab R2015b). **We are currently working on an additional implementation in R**.

## Usage

TCA requires cell-type proportion estimates for the samples in the data. These can be obtained by either using the reference-based model by Houseman et al. 2012<sup>[2](#myfootnote2)</sup> (see an implementation <a href="http://glint-epigenetics.readthedocs.io/" target="_blank">here</a>) or using the semi-supervised model by Rahmani et al. 2018<sup>[3](#myfootnote3)</sup> (does not require reference data; see an implementation <a href="https://github.com/cozygene/BayesCCE" target="_blank">here</a>).

There are two main functions in this distribution. A full documentation of the input arguments and output values of these functions is provided in the headers of these function.
* **TCA_EWAS.m** - for performing cell-type-specific EWAS under the TCA model.

* **TCA.m** - for estimating cell-type-specific methylation levels (in case only these estimates are desired rather than performing a cell-type-specific EWAS). The required input arguments for this function can be estimated using the function TCA_fit_model.m.

## Demo

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


<!---
### Citing TCA

If you use TCA in any published work, please cite the manuscript describing the method:

Elior Rahmani, Regev Schweiger, Brooke Rhead, Lindsey A. Criswell, Lisa F. Barcellos, Eleazar Eskin, Saharon Rosset, Sriram Sankararaman, and Eran Halperin. *bioRxiv*, 2018.
-->

### License

TCA is available under the <a href="https://opensource.org/licenses/GPL-3.0" target="_blank">GPL-3 license</a>.

### Author

This software was developed by Elior Rahmani (elior.rahmani@gmail.com).

 [![Travis build status](https://travis-ci.org/cozygene/TCA.svg?branch=master)](https://travis-ci.org/cozygene/TCA)
___

<a name="myfootnote1">1</a>: Rahmani et al. "Cell-type-specific resolution epigenetics without the need for cell sorting or single-cell biology." bioRxiv (2018).

<a name="myfootnote2">2</a>: Houseman et al. "DNA methylation arrays as surrogate measures of cell mixture distribution." BMC bioinformatics (2012).

<a name="myfootnote3">3</a>: Rahmani et al. "BayesCCE: a Bayesian framework for estimating cell-type composition from DNA methylation without the need for methylation reference." Genome biology (2018).
