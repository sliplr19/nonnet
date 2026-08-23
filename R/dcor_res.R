#' Distance Correlation Test After Residualization
#'
#' Residualizes two random variables with respect to a set of covariates using
#' generalized additive models, then tests the association between the
#' residualized variables using a permutation-based distance correlation test.
#'
#'
#' @param dat A data frame
#' @param A Character string specifying the name of the first variable to be used in the distance correlation.
#' @param C Character string specifying the name of the second variable to be used in the distance correlation.
#' @param B_variables Character vector containing the names of the covariates
#'   used as predictors when residualizing `A` and `C`.
#' @param nperm Integer. Number of permutations used in the distance
#'   correlation test.
#'
#' @return A numeric value giving the permutation p-value from
#'   `energy::dcor.test()`.
#'
#' @export

dcor_res <- function(dat, A, C, B_variables, nperm){
  dat <- data.frame(dat)
  colnames(dat)[colnames(dat) == A] <- "A"
  colnames(dat)[colnames(dat) == C] <- "C"


  Bs <- paste0(B_variables, collapse = "+")

  Aform <- stats::as.formula(paste("A ~", Bs))
  Cform <- stats::as.formula(paste("C ~", Bs))

  A_fit <- mgcv::gam(Aform, data = dat, method = "REML")
  C_fit <- mgcv::gam(Cform, data = dat, method = "REML")

  A_resid <- stats::resid(A_fit)
  C_resid <- stats::resid(C_fit)

  A_resid_lm <- stats::resid(stats::lm(A_resid ~ C_resid))

  datnew <- data.frame("A" = A_resid_lm, "C" = C_resid)

  distancepvAC = energy::dcor.test(datnew$A, datnew$C, R = nperm)$p.value

  return(distancepvAC)
}
