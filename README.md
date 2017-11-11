# IndeterminateR

This R-package verifies the claims in the article "The Full Set of Solutions of Linear Rational Expectation Models". In particular, it shows that some equilibria are classified as indeterminate by Lubik and Schorfheide even though they have the same observable outcome.

## Installation

Run `devtools::install_github("bfunovits/indeterminateR", build_vignettes = TRUE)`.
The package `devtools` needs to be installed in order that the above line of code work.
By default, vignettes are not built when using `install_github()`, thus the argument `build_vignettes` must be set `TRUE`.
