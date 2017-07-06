# Lubik Marzo (26_0) ##########################################

#' Obtain linear parameters of model by Lubik and Marzo (2007)
#'
#' This function returns a verison of the New Keynesian monetary model with a monetary rule which is similar to the one in \href{https://doi.org/10.1016/j.iref.2004.09.008}{Lubik Marzo (2007)}, see the vignette for more detail.
#' In particular, it takes the deep parameters as input and returns the matrices occurring in Sims' canonical form.
#'
#' @param s Elasticity of intertemporal substitution and the inverse of the coefficient of relative risk aversion.
#' @param b Subjective time preference rate
#' @param k Inverse of elasticity of aggregate supply with respect to inflation
#' @param f_y Measures the elasticity of interest rate response w.r.t. output.
#' @param f_p Measures the elasticity of interest rate response w.r.t. inflation.
#' @param f_r Interest smoothing coefficient (page 27 in LM07)
#'
#' @return List object containing the parameter matrices \deqn{\Gamma_0, \Gamma_1, \Psi, \Pi}.
#' @export
#'
linear_parameters_LM_26_0 <- function(s = 1,
                                      b = 0.99,
                                      k = 0.5,
                                      f_y = 0.8,
                                      f_p = 0.7,
                                      f_r = 1.2){
  Gamma_0 <- matrix(c(1, 0, -s,
                      0, b, 0,
                      0, 0, 1),
                    byrow = TRUE,
                    nrow = 3)
  Gamma_1 <- matrix(c(1, -s, 0,
                      -k, 1, 0,
                      0, 0, f_r),
                    byrow = TRUE,
                    nrow = 3)
  Psi <- diag(c(1, 1, 1))
  Pi <- matrix(c(1, -s,
                 -k, 1,
                 f_y, f_p), # there might be f_y and f_p missing
               byrow = TRUE,
               nrow = 3)
  return(list(
    Gamma_0 = Gamma_0,
    Gamma_1 = Gamma_1,
    Psi = Psi,
    Pi = Pi)
  )
}

