# IndeterminateR

This R-package verifies the claims in the article 

**The Full Set of Solutions of Linear Rational Expectation Models**, 
published in Economics Letters, Volume 161, December 2017, Pages 47–51, 
https://doi.org/10.1016/j.econlet.2017.09.021

In particular, it shows that some equilibria are classified as indeterminate by Lubik and Schorfheide even though they have the same observable outcome.

## Installation

Run `devtools::install_github("bfunovits/indeterminateR", build_vignettes = TRUE)`.
The package `devtools` needs to be installed in order that the above line of code work.
By default, vignettes are not built when using `devtools::install_github()`, thus the argument `build_vignettes` must be set to `TRUE`.
