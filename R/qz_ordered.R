#' Decoupling into Stable and Unstable Systems
#'
#' Uses the \link[QZ]{qz} decomposition of \eqn{\Gamma_0 , \Gamma_1} from Sims' canoncial form to decouple into a stable and an unstable part.
#' The function arranges the equations such that all generalized eigenvalues in the QZ-decomposition which are larger than the given \code{threshold} are in the bottom right block.
#' To be more precise, the ratios of the diagonal elements of the upper-triangular matrices in the QZ-decomposition for which
#' \deqn{ \Omega_{i,i} / \Lambda_{i,i} > threshold }
#' holds, are in the bottom right corner.
#' The functions ensures that the products \code{t(Q)*} \eqn{\Lambda} \code{*t(Z)} and \code{t(Q)*} \eqn{  \Omega } \code{*t(Z)} remains unchanged.
#'
#' Note that the function \link[QZ]{qz} returns the QZ-decomposition such that
#' \deqn{( \Gamma_0, \Gamma_1 ) = ( VSL *  S  * Conj(t(VSR)), VSL * T * Conj(t(VSR)) ),}
#' where \code{VSL} and \code{VSR} are Hermitian matrices, and \code{S} and \code{T} are upper triangular matrices (output of Fortran routines).
#' However, the QZ-decomposition corresponds in Sims' notation ( see page 9, formula (32) in \href{https://doi.org/10.1023/A:1020517101123}{Sims (2001)} ) to
#' \deqn{( \Gamma_0, \Gamma_1 ) = ( Conj(t(Q)) * \Lambda * Conj(t(Z)) , Conj(t(Q)) * \Omega * Conj(t(Z)) ).}
#' The code is intended to keep the notation close to Sims (2001) and adjusts the outputs of the \link[QZ]{qz} function accordingly.
#'
#' Some parts of this function are based on R code by Christopher A. Sims (22 Feb 2004) which in turn is based on earlier matlab code (finished 27 Apr 2000).
#' Retrieved on 14 June 2017 from \url{http://sims.princeton.edu/yftp/gensys/Rfiles/}.
#'
#' @param threshold Threshold to define when a generalized eigenvalue is considered as unstable
#' @param my_Gamma_0 Matrix from Sims' canonical form
#' @param my_Gamma_1 Matrix from Sims' canonical form
#' @param tol Numeric tolerance level. Default value set to the square root of machine precision.
#'
#' @return List object containing
#'   \itemize{
#'     \item \strong{Q_unstable:} The rows of the orthogonal matrix \code{Q} pertaining to the unstable generalized eigenvalues.
#'     \item \strong{Q_stable:} Same for stable generalized eigenvalues.
#'     \item \strong{generalized_ev:} The generalized eigenvalues.
#'   }
#' @importFrom purrr %>%
#' @export

qz_ordered <- function (threshold = 1.01,
                        my_Gamma_0,
                        my_Gamma_1,
                        tol = sqrt(.Machine$double.eps)) {
  # Apply the QZ decomposition and save the output
  # The list objects are such that (Gamma_0, \Gamma_1) = ( VSL %*% S %*% Conj(t(VSR)), VSL %*% T %*% Conj(t(VSR)) )
  # More info under http://www.netlib.org/lapack/complex16/zgges.f

  qz_list <- QZ::qz(my_Gamma_0, my_Gamma_1)
  # message("Gentle Reminder: \n
  #         Bernd, please checke what happens when there are complex generalized eigenvalues. \n
  #         You use the generic \"qz\" function which calls qz.dgges for real input. \n
  #         This generates the real QZ decomposition which is not necessarily upper-triangular. \n
  #         Try using qz.zgges and check for complex numbers. \n")

  # Make it consistent with Sims' notation
  my_qz_list <- list()
  my_qz_list$Lambda <- qz_list$S
  my_qz_list$Omega <-qz_list$T
  my_qz_list$Q <- qz_list$VSL %>% t() %>% Conj()
  my_qz_list$Z <- qz_list$VSR
  # generalized eigenvalues to be ordered
  my_qz_list$ev_pairs <- cbind(qz_list$ALPHA, qz_list$BETA)
  # Some tweaking to go around division by zero if Lambda is singular
  my_qz_list$gev_abs <- ifelse(abs(qz_list$ALPHA) < tol, -1, abs(qz_list$BETA/qz_list$ALPHA) )
  rm(qz_list)

  n  <- nrow(my_qz_list$Lambda) # dimension of system

  # Starting with the last generalized eigenvalue ----
  for (i in  n:1) {
    m <- 0
    # Search for the "first" unstable root, starting in the bottom right corner
    # If it is in the bottom right corner, nothing needs to be done
    # If it is not in the bottom right corner, it needs to be switched with all stable ones which are more bottom right
    for (j in i:1) {
      if (my_qz_list$gev_abs[j] > threshold || my_qz_list$gev_abs[j] < -.1) { # note that only those generalized eigenvalues are negative that have a zero on the left hand side upper-triangular matrix
        m <- j
        break     #found an unstable root.  Now check for stable ones below.
      }
    }
    if (m == 0) {
      break #quit sort because only stable roots left
    } else {
      # if the unstable roots found at index m is not as much bottom right as it should be, switch blocks until it is as much bottom right as it should be
      if (m < i) {
        for (k in m:(i - 1)) {
          my_qz_list <- qz_switch(my_qz_list, k)
        }
      }
    }
  }

  # Different cases for unstable roots
  # 1) no unstable roots -> length of "first_unstable_index" is zero
  # 2) at least one unstable root -> length of "first_unstable_index" equal to number of unstable roots -> choose first one
  first_unstable_index <- which(my_qz_list$gev_abs > threshold)
  if (length(first_unstable_index) > 0){
    first_unstable_index <- first_unstable_index[1]
  }

  # checked that it should be in accordance with Sims' notation 15 Jun 2017
  if(length(first_unstable_index) == 0){
    # only stable roots
    Q_unstable <- NULL
    Q_stable <- my_qz_list$Q
  } else if (first_unstable_index == 1){
    # only unstable roots
    Q_unstable <- my_qz_list$Q[, first_unstable_index:n, drop = FALSE]
    Q_stable <- NULL
  } else {
    # base case: stable and unstable generalized EV
    Q_unstable <- my_qz_list$Q[first_unstable_index:n, , drop = FALSE]
    Q_stable <- my_qz_list$Q[1:(first_unstable_index-1), , drop = FALSE]
  }

  return(list(
    Q_unstable = Q_unstable,
    Q_stable = Q_stable,
    generalized_ev = my_qz_list$gev_abs
  ))
}


################################################################################

#' Switch diagonal elements of upper-triangular matrices
#'
#' This function switches the diagonal elements of the upper-triangular matrices \eqn{\Lambda} and \eqn{\Omega} in the QZ-decomposition.
#' It is a helper function for making the QZ-decomposition unique in the case where all generalized eigenvalues are different from each other.
#'
#' Based on the code by Christopher A. Sims, retrieved on 18 June 2017 from \url{http://sims.princeton.edu/yftp/gensys/Rfiles/}.
#' See also \href{https://doi.org/10.1007/978-94-015-8196-7_11}{Kagstrom (1993)}
#'
#' @param i Index of the diagonal element to be switched with the following one, i.e. \code{i+1}
#' @param my_qz_list List object containing the orthogonal and upper-triangular matrices returned by the QZ-decomposition in R.
#' @param tol Numeric tolerance level. Default set to square root of machine precision.
#'
#' @return List object \code{qz_list} with elements
#'   \itemize{
#'     \item \strong{Lambda:} Upper-triangular matrix in Sims' notation.
#'     \item \strong{Omega:} Upper-triangular matrix in Sims' notation.
#'     \item \strong{Z:} Orthogonal matrix in Sims' notation.
#'     \item \strong{Q:} Orthogonal matrix in Sims' notation.
#'     \item \strong{gev_abs:} The absolute values of generalized eigenvalues.
#'     \item \strong{ev_pairs:} The pairs of diagonal elements of \eqn{\Lambda} and \eqn{\Omega}.
#'   }
#' @importFrom purrr %>%
#' @export

qz_switch <- function(my_qz_list, i = 1, tol = sqrt(.Machine$double.eps)){

  lambda_11 <- my_qz_list$Lambda[i, i]
  lambda_12 <- my_qz_list$Lambda[i, (i + 1)]
  lambda_22 <- my_qz_list$Lambda[i + 1, i + 1]

  omega_11 <- my_qz_list$Omega[i, i]
  omega_12 <- my_qz_list$Omega[i, (i + 1)]
  omega_22 <- my_qz_list$Omega[i + 1, i + 1]

  w2 <- rep(0, 2)

  # Handle cases with coincident zeros (matrix pencil has determinant identically zero)
  if (abs(lambda_22) < tol & abs(omega_22) < tol){
    # % (2,2) elements of Lambda and Omega are zero
    if (abs(lambda_11) < tol) {
      # (1,1) element of Lambda zero too -> do nothing
      return(my_qz_list)
    } else {
      # (1,1) element of non-zero -> orthogonal transformation (such that 0 in (1,1) element?)

      # Orthogonal right-multiplication of Lambda and Omega
      orth_right <- c(lambda_12, -lambda_11) # normal vector of first row of Lambda
      orth_right <- orth_right / sqrt(sum(Conj(orth_right) * orth_right))
      orth_right <- matrix(c(orth_right, Conj(orth_right[2]), -Conj(orth_right[1])), nrow = 2) # byrow = FALSE is default

      # Orthogonal left-multiplication of Lambda and Omega
      orth_left <- diag(2)
    }
  } else {
    if (abs(lambda_11) < tol && abs(omega_11) < tol){
      if (abs(lambda_22) < tol){
        ## Lambda has a zero at (2,2) element -> do nothing
        return(my_qz_list)
      } else {
        # Element (2,2) of Lambda is non-zero -> create orthogonal matrix

        # Orthogonal right-multiplication of Lambda and Omega
        orth_right <- diag(2)

        # Orthogonal left-multiplication of Lambda and Omega
        orth_left <- c(lambda_22, -lambda_12) # normal vector of second column of Lambda
        orth_left <- orth_left / sqrt(sum(orth_left * Conj(orth_left)))
        orth_left <- matrix(c(Conj(orth_left[2]), -Conj(orth_left[1]), orth_left), byrow = TRUE, nrow = 2, ncol = 2)
      }
    } else {
      # usual case

      # Orthogonal right-multiplication of Lambda and Omega
      orth_right <- c(lambda_22 * omega_12 - omega_22 * lambda_12, Conj(lambda_22 * omega_11 - omega_22 * lambda_11))

      # Orthogonal left-multiplication of Lambda and Omega
      orth_left <- c(Conj(lambda_12 * omega_11 - omega_12 * lambda_11), Conj(lambda_22 * omega_11 - omega_22 * lambda_11))

      # Helpers for normalization of wz and xy
      n <- 1/sqrt(sum(orth_right * Conj(orth_right)))
      m <- 1/sqrt(sum(orth_left * Conj(orth_left)))

      if (Re(m) < tol){
        ## all elements of a and b proportional
        return(my_qz_list)
      }
      orth_right <- orth_right * n
      orth_left <- orth_left * m
      orth_right <- matrix(c(orth_right, Conj(orth_right[2]), -Conj(orth_right[1])),
                byrow = TRUE,
                nrow = 2)
      orth_left <-  matrix( c(orth_left, Conj(orth_left[2]), -Conj(orth_left[1])),
                byrow = TRUE,
                nrow = 2)
    }
  }
  my_qz_list$Lambda[i:(i + 1), ] <- orth_left %*% my_qz_list$Lambda[i:(i + 1), , drop = FALSE]
  my_qz_list$Omega[i:(i + 1), ] <- orth_left %*% my_qz_list$Omega[i:(i + 1), , drop = FALSE]
  my_qz_list$Lambda[, i:(i + 1)] <- my_qz_list$Lambda[, i:(i + 1), drop = FALSE] %*% orth_right
  my_qz_list$Omega[, i:(i + 1)] <- my_qz_list$Omega[, i:(i + 1), drop = FALSE] %*% orth_right
  my_qz_list$Z[, i:(i + 1)] <- my_qz_list$Z[, i:(i + 1), drop = FALSE] %*% orth_right
  my_qz_list$Q[i:(i + 1), ] <- orth_left %*% my_qz_list$Q[i:(i + 1), , drop = FALSE]
  my_qz_list$gev_abs[c(i, i+1)] <- my_qz_list$gev_abs[c(i+1, i)]
  my_qz_list$ev_pairs[c(i, i+1), ] <- my_qz_list$ev_pairs[c(i+1, i), ]

  return(my_qz_list)
}
