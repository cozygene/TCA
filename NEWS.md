TCA v1.2.1 (Release date: 02-13-2021)
==============

*NOTE* v1.2.1 is not backward compatible with earlier versions

Changes:

* A new argument 'vars.mle' in 'tca' which controls the method used for optimization; setting to FALSE allows a much faster optimization.

* A new argument 'fast_mode' in 'tcareg' allows a much faster optimization.

* A new argument 'constrain_mu' in 'tca' allows to determine whether to constrain the mean parameters in the optimization; setting to false now allows to obtain p-values for the covariates in 'C1' and 'C2'.

* A new 'scale' argument in 'tensor' for dividing estimates by estimated standard deviations.

* Default values of several arguments have been changed.

* Switched the default FDR method used in ‘tcareg’ back from “BY” to “BH” (can be changed by editing the config file).


TCA v1.1.0 (Release date: 11-14-2019)
==============

Changes:

* Introducing the argument 'C1.map' in the function 'tca' for allowing source-specific covariates 'C1' to affect only a user-defined subset of the sources.

* Now allowing to provide 'tau' as an input to the 'tca' function.

* Switched the FDR method used in 'tcareg' from "BH" to "BY".

* Added a 'verbose' argument to all functions.

Bug fixes:

* Added glmnet as a required dependency to the description file.

* Fixed a formatting error in 'tcareg' in some cases where both 'C1' and 'C3' were  used.

* Fixed a problem in 'tca' in cases where 'refit_W=TRUE', 'C1=NULL', 'C2=NULL'.

* Fixed logging issue in cases where 'debug=TRUE', 'refit_W=TRUE'.

* Fixed a bug in the calculation of the returned value of 'tensor' in cases where 'C2' were provided.


TCA v1.0.0 (Release date: 05-22-2019)
==============
