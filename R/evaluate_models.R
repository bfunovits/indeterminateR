#' Determine existence and uniqueness properties of the New Keynesian Monetary Model
#'
#' Function for evaluating the existence and uniqueness properties of a system in Sims' canonical form.
#'
#' @param deep_params List object containing the deep parameters for the respective model.
#'   Default is the model by Lubik and Marzo (2007), equation (26) with \code{i = 0}.
#'   For more detail on the deep parameters see \link{linear_parameters_LM_26_0}.
#' @param threshold Numeric value. Generalized eigenvalues larger than "threshold" are considered unstable.
#' @param my_model String indicating which model should be evaluated.
#'   Default is the model by Lubik and Marzo (2007), equation (26) with \code{i = 0}.
#' @param tol Threshold to determine what is considered practically zero.
#'   Default set to square root of machine precision.
#' @param verbose Boolean variable to indicate whether matrices in Sims' canonical form,
#'   the QZ decomposition, and the matrices determining existence and uniqueness should be printed.
#'   Default is set to FALSE.

#' @return Object of type "data_frame". Contains information about the existence and uniqueness properties.
#'   In particular, the variables are
#'   \itemize{
#'     \item existence_boolean TRUE if a stable and causal solution exists
#'     \item existence_dim_kernel Dimension of the kernel of the matrix pertaining to the endogenous forecast errors in the existence cóndition.
#'     \item uniqueness_boolean TRUE if the uniqueness condition is satisfied
#'     \item indeterminacy_dim Dimension of indeterminacy as characterized in the paper accompanying this package
#'     \item indeterminacy_smallest_sv Smallest eigenvalue of the intersection of the spaces in the uniqueness condition. The smaller, the "more unique".
#'   }
#'
#' @importFrom purrr %>%
#' @export
#'
evaluate_my_model <- function(deep_params = NULL, threshold = 1.01, my_model = "Lubik_Marzo_26_0", tol = sqrt(.Machine$double.eps), verbose = FALSE){

  # Choose a model and obtain its Sims canonical form as a function of the deep parameters ----

  if (my_model == "Lubik_Marzo_26_0"){
    if (is.null(deep_params)){
      deep_params <- list(s = 1, b = 0.95, f_y = 0.8, f_p = 0.7, k = 0.5, f_r = 1.2)
    }
    list_lin_params <- linear_parameters_LM_26_0(s = deep_params$s, b = deep_params$b, k = deep_params$k, f_y = deep_params$f_y, f_p = deep_params$f_p, f_r = deep_params$f_r)
  } else {
    cat("Hi. \n
        Your model does not exist (yet).")
  }
  my_Gamma_0 <- list_lin_params$Gamma_0
  my_Gamma_1 <- list_lin_params$Gamma_1
  my_Pi <- list_lin_params$Pi
  my_Psi <- list_lin_params$Psi

  # Dimensions in Sims canoncial form ----
  n <- nrow(my_Gamma_0)
  n_endog_fc_err <- ncol(my_Pi)
  n_exog <- ncol(my_Psi)

  # Calculate the QZ decomposition and define stable and unstable part of the orthogonal matrix Q ----
  list_decoupling_orth_mat <- qz_ordered(my_Gamma_0 = my_Gamma_0, my_Gamma_1 = my_Gamma_1)

  if (!is.null(list_decoupling_orth_mat$Q_stable)){
    Q_s <- list_decoupling_orth_mat$Q_stable
    n_stable <- nrow(Q_s)
  } else {
    Q_s <- matrix(0, nrow = 1, ncol = n)
    n_stable <- 0
  }

  if (!is.null(list_decoupling_orth_mat$Q_unstable)){
    Q_u <- list_decoupling_orth_mat$Q_unstable
    n_unstable <- nrow(Q_u)
  } else {
    Q_u <- matrix(0, nrow = 1, ncol = n)
    n_unstable <- 0
  }

  # Singular Value Decompositions X = u d Conj(t(v)) ----

  # 1) Unstalbe part of Pi (endogenous forecast errors)
  QuPi <- Q_u %*% my_Pi
  rank_QuPi <- base::qr(QuPi)$rank
  list_svd_unstable_endog <- svd(QuPi, nv = n_endog_fc_err)
  # Needed for correct "degree of indeterminacy" (intersection with Vs_1 below)
  # Needed for determining the dimension of  the kernel of Sims existence equation
  if (n_endog_fc_err > rank_QuPi){
    Vu_2 <- list_svd_unstable_endog$v[, (rank_QuPi+1):n_endog_fc_err, drop = FALSE]
  } else {
    Vu_2 <- NULL
  }
  # Needed for evaluating existence according to Sims (2001), page 13, formulae (46) and (48)
  QuPi_svd_u <- list_svd_unstable_endog$u


  # 2) Stable part of Pi (endogenous forecast errors)
  QsPi <- Q_s %*% my_Pi
  rank_QsPi <- base::qr(QsPi)$rank
  list_svd_stable_end <- svd(QsPi, nv = rank_QsPi) # deleted argument
  # Needed for "degree of indeterminacy"
  Vs_1 <- list_svd_stable_end$v
  if (length(list_svd_stable_end$d) > 1){
    Ds <- diag(list_svd_stable_end$d)
  } else {
    Ds <- list_svd_stable_end$d
  }

  # 3) Unstable part of Psi (exogenous variables)
  # Needed for evaluating existence according to Sims (2001), page 13, formulae (46) and (48)
  QuPsi <- Q_u %*% my_Psi
  rank_QuPsi <- base::qr(QuPsi, tol = tol)$rank
  list_svd_stable_exog <- svd(QuPsi, nv = rank_QuPsi)
  QuPsi_svd_u <- list_svd_stable_exog$u

  # Print some info if required ----
  if (verbose == TRUE){
    cat("Matrices in Sims' canonical form: \n")
    print(list_lin_params)

    cat("Orthogonal matrices for stable and unstable subspaces as well as the generalized EV: \n")
    print(list_decoupling_orth_mat)

    cat("Matrices for determining existence and uniqueness: \n\n")
    cat("QuPsi: Exogenous influence on unstable part: \n")
    print(QuPsi)
    cat("QuPi: Endogenous influence on unstable part: \n")
    print(QuPi)
    cat("QsPi: Endogenous influence on stable part: \n")
    print(QsPi)
  }

  # Evaluations ----
  # There is a bug in base::qr()$rank and Matrix::rankMatrix() which returns 1 for input 1e-100.
  # The following lines of code make up for this bug while it is repaired (bug report sent to package maintainers on 22 June 17)

  if(is.null(Vu_2)){
    intersect_spaces <- NULL
  } else {
    intersect_spaces <- t(Conj(Vs_1)) %*% Vu_2
    if (max(dim(intersect_spaces)) == 1 && intersect_spaces < tol) {
      intersect_spaces_dim <- NULL
    }

  }


  return(dplyr::data_frame(
    existence_boolean = ifelse(n_unstable > 0,
                               base::norm( (diag(rep(1, n_unstable)) - QuPi_svd_u %*% t(QuPi_svd_u) ) %*% QuPsi_svd_u, type = "m") < tol,
                               TRUE),
    existence_dim_kernel = n_endog_fc_err - rank_QuPi,
    uniqueness_boolean = ifelse(!(is.null(Vu_2) || is.null(intersect_spaces)),
                                base::norm( intersect_spaces, type = "m" ) < tol,
                                TRUE),
    indeterminacy_dim = ifelse(!is.null(intersect_spaces_dim),
                               base::qr(as.matrix(intersect_spaces))$rank,
                               0),
    indeterminacy_smallest_sv = ifelse(!is.null(intersect_spaces),
                                       svd(intersect_spaces)$d %>% utils::tail(n = 1),
                                       -1)
  ))
}
