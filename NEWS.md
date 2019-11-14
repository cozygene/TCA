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
