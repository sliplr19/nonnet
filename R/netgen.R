#' Generate Simulated Nonlinear Network Data
#'
#' Simulates a dataset containing a central variable, peripheral variables,
#' and a target variable under different functional relationships.
#'
#' @param n Integer. Number of observations to generate.
#' @param mean_cent Mean of the central variable `A`.
#' @param sd_cent Standard deviation of the central variable `A`.
#' @param mean_peri_non Numeric vector. Means of the peripheral nonlinear variables.
#' @param mean_peri_lin Numeric vector. Means of the peripheral linear variables.
#' @param sd_peri_non Numeric vector. Standard deviations of the peripheral nonlinear variables.
#' @param sd_peri_lin Numeric vector. Standard deviations of the peripheral linear variables.
#' @param mean_target Mean of the random error term for the target variable `C`.
#' @param sd_target Standard deviation of the random error term for the target variable `C`.
#' @param beta_peri_lin Numeric vector. Regression coefficients for effects of the peripheral linear variables on `C`.
#' @param beta_cent_lin Numeric vector. Regression coefficients for effects of `A` on the peripheral linear variables.
#' @param beta_lin Coefficient for the linear effect of `A` on `C`.
#' @param beta_non Coefficient for the nonlinear effect of `A` on `C`.
#' @param beta_cent_non Numeric vector. Regression coefficients for effects of `A` on the peripheral nonlinear variables.
#' @param beta_peri_non Numeric vector. Regression coefficients for effects of the peripheral nonlinear variables on `C`.
#' @param inter_cent_lin Integer vector. Indices identifying peripheral variables involved in interactions with `A`.
#' @param beta_con Regression coefficient for the confounding effect.
#' @param func Character. Functional form relating `A` to `C`.
#' @param con_func Character. Functional form relating `A` to the nonlinear peripheral variables.
#'
#' @return A data frame containing `A`, `C`, linear peripheral variables
#'   `B1`, `B2`, ..., and nonlinear peripheral variables `D1`, `D2`, ....
#'
#' @examples
#' set.seed(123)
#'
#' dat <- netgen(
#'   n = 100,
#'   mean_cent = 0,
#'   sd_cent = 1,
#'   mean_peri_non = rep(0, 4),
#'   mean_peri_lin = rep(0, 3),
#'   sd_peri_non = rep(1, 4),
#'   sd_peri_lin = rep(1, 3),
#'   mean_target = 0,
#'   sd_target = 1,
#'   beta_peri_lin = c(0.3, 0.4, 0.5),
#'   beta_cent_lin = rep(0.5, 3),
#'   beta_lin = 0.5,
#'   beta_non = 0.75,
#'   beta_cent_non = c(0.5, 0.6, 0.7, 0.8),
#'   beta_peri_non = c(0.3, 0.4, 0.5, 0.6),
#'   inter_cent_lin = c(1, 2, 3),
#'   beta_con = 0,
#'   func = "quad",
#'   con_func = "quad"
#' )
#'
#' head(dat)
#'
#' @export


netgen <- function(n, mean_cent, sd_cent, mean_peri_non, mean_peri_lin, sd_peri_non,
                   sd_peri_lin, mean_target, sd_target,

                   beta_peri_lin, beta_cent_lin, beta_lin, beta_non, beta_cent_non,
                   beta_peri_non, inter_cent_lin, beta_con,
                   func,
                   con_func){
  A <- stats::rnorm(n, mean_cent, sd_cent)
  if(func == "quad"){
    if(con_func == "linear"){
      peri_non <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_non)))
      for(i in 1:length(mean_peri_non)){
        peri_non[,i] <- data.frame(beta_cent_non[i]*A + stats::rnorm(n, mean_peri_non[i], sd_peri_non[i]))
      }
      peri_lin <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_lin)))
      for(j in 1:length(mean_peri_lin)){
        peri_lin[,j] <- beta_cent_lin[j]*A + stats::rnorm(n, mean_peri_lin[j], sd_peri_lin[j])
      }
      C <- beta_non*A^2 + beta_lin*A + stats::rnorm(n, mean_target, sd_target)
      for(k in 1:ncol(peri_lin)){
        C = C + beta_peri_lin[k]*peri_lin[k]
      }
      for(l in 1:ncol(peri_non)){
        C= C + beta_peri_non[l]*peri_non[l]
      }
    } else if(con_func == "quad"){
      peri_non <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_non)))
      for(i in 1:length(mean_peri_non)){
        peri_non[,i] <- data.frame(beta_cent_non[i]*A^2 + stats::rnorm(n, mean_peri_non[i], sd_peri_non[i]))
      }
      peri_lin <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_lin)))
      for(j in 1:length(mean_peri_lin)){
        peri_lin[,j] <- beta_cent_lin[j]*A + stats::rnorm(n, mean_peri_lin[j], sd_peri_lin[j])
      }
      C <- beta_non*A^2 + beta_lin*A + stats::rnorm(n, mean_target, sd_target)
      for(k in 1:ncol(peri_lin)){
        C = C + beta_peri_lin[k]*peri_lin[k]
      }
      for(l in 1:ncol(peri_non)){
        C= C + beta_peri_non[l]*peri_non[l]
      } } else if(con_func == "inter"){
        peri_lin <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_lin)))
        for(j in 1:length(mean_peri_lin)){
          peri_lin[,j] <- beta_cent_lin[j]*A + stats::rnorm(n, mean_peri_lin[j], sd_peri_lin[j])
        }
        peri_non <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_non)))
        for(i in 1:length(mean_peri_non)){
            if(i %in% inter_cent_lin){
            peri_non[,i] <- beta_cent_non[i] * A * peri_lin[,i] +
              stats::rnorm(n, mean_peri_non[i], sd_peri_non[i])
          } else {
            peri_non[,i] <- beta_cent_non[i] * A +
              stats::rnorm(n, mean_peri_non[i], sd_peri_non[i])
          }
        }
        C <- beta_non*A^2 + beta_lin*A + stats::rnorm(n, mean_target, sd_target)
        for(m in 1:length(inter_cent_lin)){
          C = C + beta_non* A * peri_lin[,inter_cent_lin[m]]
        }
        for(l in 1:length(inter_cent_lin)){
          C = C + beta_peri_non[l]*peri_non[,l]
        }
      } else if(con_func == "log"){
        peri_non <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_non)))
        for(i in 1:length(mean_peri_non)){
          peri_non[,i] <- data.frame(beta_cent_non[i]*log(abs(A)) + stats::rnorm(n, mean_peri_non[i], sd_peri_non[i]))
        }
        peri_lin <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_lin)))
        for(j in 1:length(mean_peri_lin)){
          peri_lin[,j] <- beta_cent_lin[j]*A + stats::rnorm(n, mean_peri_lin[j], sd_peri_lin[j])
        }
        C <- beta_non*A^2 + beta_lin*A + stats::rnorm(n, mean_target, sd_target)
        for(k in 1:ncol(peri_lin)){
          C = C + beta_peri_lin[k]*peri_lin[k]
        }
        for(l in 1:ncol(peri_non)){
          C= C + beta_peri_non[l]*peri_non[l]
        }
      }
  }
  else if(func == "inter"){
    if(con_func == "linear"){
      peri_non <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_non)))
      for(i in 1:length(mean_peri_non)){
        peri_non[,i] <- data.frame(beta_cent_non[i]*A + stats::rnorm(n, mean_peri_non[i], sd_peri_non[i]))
      }
      peri_lin <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_lin)))
      for(j in 1:length(mean_peri_lin)){
        peri_lin[,j] <- beta_cent_lin[j]*A + stats::rnorm(n, mean_peri_lin[j], sd_peri_lin[j])
      }
      C <- beta_lin*A + + stats::rnorm(n, mean_target, sd_target)
      for(m in 1:length(inter_cent_lin)){
        C = C + beta_non[inter_cent_lin[m]]*A*peri_lin[inter_cent_lin[m]]
      }
      for(l in 1:ncol(peri_non)){
        C= C + beta_peri_non[l]*peri_non[l]
      }
    } else if(con_func == "quad"){
      peri_non <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_non)))
      for(i in 1:length(mean_peri_non)){
        peri_non[,i] <- data.frame(beta_cent_non[i]*A^2 + stats::rnorm(n, mean_peri_non[i], sd_peri_non[i]))
      }
      peri_lin <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_lin)))
      for(j in 1:length(mean_peri_lin)){
        peri_lin[,j] <- beta_cent_lin[j]*A + stats::rnorm(n, mean_peri_lin[j], sd_peri_lin[j])
      }
      C <- beta_lin*A + + stats::rnorm(n, mean_target, sd_target)
      for(m in 1:length(inter_cent_lin)){
        C = C + beta_non[inter_cent_lin[m]]*A*peri_lin[inter_cent_lin[m]]
      }
      for(l in 1:ncol(peri_non)){
        C= C + beta_peri_non[l]*peri_non[l]
      }} else if(con_func == "inter"){
        peri_lin <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_lin)))
        for(j in 1:length(mean_peri_lin)){
          peri_lin[,j] <- beta_cent_lin[j]*A + stats::rnorm(n, mean_peri_lin[j], sd_peri_lin[j])
        }
        peri_non <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_non)))
        for(i in 1:length(mean_peri_non)){
          if(i %in% inter_cent_lin){
            peri_non[,i] <- beta_cent_non[i] * A * peri_lin[,i] +
              stats::rnorm(n, mean_peri_non[i], sd_peri_non[i])
          } else {
            peri_non[,i] <- beta_cent_non[i] * A +
              stats::rnorm(n, mean_peri_non[i], sd_peri_non[i])
          }
        }
        C <- beta_lin*A + stats::rnorm(n, mean_target, sd_target)
        for(m in 1:length(inter_cent_lin)){
          C = C + beta_non*A*peri_lin[,inter_cent_lin[m]]
        }
        for(l in 1:length(inter_cent_lin)){
          C = C + beta_peri_non[l]*peri_non[,l]
        }
      } else if(con_func == "log"){
        peri_non <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_non)))
        for(i in 1:length(mean_peri_non)){
          peri_non[,i] <- data.frame(beta_cent_non[i]*log(abs(A)) + stats::rnorm(n, mean_peri_non[i], sd_peri_non[i]))
        }
        peri_lin <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_lin)))
        for(j in 1:length(mean_peri_lin)){
          peri_lin[,j] <- beta_cent_lin[j]*A + stats::rnorm(n, mean_peri_lin[j], sd_peri_lin[j])
        }
        C <- beta_lin*A + + stats::rnorm(n, mean_target, sd_target)
        for(m in 1:length(inter_cent_lin)){
          C = C + beta_non[inter_cent_lin[m]]*A*peri_lin[inter_cent_lin[m]]
        }
        for(l in 1:ncol(peri_non)){
          C= C + beta_peri_non[l]*peri_non[l]
        }}

  }
  else if(func == "linear"){
    if(con_func == "linear"){
      peri_non <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_non)))
      for(i in 1:length(mean_peri_non)){
        peri_non[,i] <- data.frame(beta_cent_non[i]*A + stats::rnorm(n, mean_peri_non[i], sd_peri_non[i]))
      }
      peri_lin <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_lin)))
      for(j in 1:length(mean_peri_lin)){
        peri_lin[,j] <- beta_cent_lin[j]*A + stats::rnorm(n, mean_peri_lin[j], sd_peri_lin[j])
      }
      C <- beta_non*A + beta_lin*A + stats::rnorm(n, mean_target, sd_target)
      for(k in 1:ncol(peri_lin)){
        C = C + beta_peri_lin[k]*peri_lin[k]
      }
      for(l in 1:ncol(peri_non)){
        C= C + beta_peri_non[l]*peri_non[l]
      }} else if(con_func == "quad"){
        peri_non <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_non)))
        for(i in 1:length(mean_peri_non)){
          peri_non[,i] <- data.frame(beta_cent_non[i]*A^2 + stats::rnorm(n, mean_peri_non[i], sd_peri_non[i]))
        }
        peri_lin <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_lin)))
        for(j in 1:length(mean_peri_lin)){
          peri_lin[,j] <- beta_cent_lin[j]*A + stats::rnorm(n, mean_peri_lin[j], sd_peri_lin[j])
        }
        C <- beta_non*A + beta_lin*A + stats::rnorm(n, mean_target, sd_target)
        for(k in 1:ncol(peri_lin)){
          C = C + beta_peri_lin[k]*peri_lin[k]
        }
        for(l in 1:ncol(peri_non)){
          C= C + beta_peri_non[l]*peri_non[l]
        }} else if(con_func == "inter"){
          peri_lin <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_lin)))
          for(j in 1:length(mean_peri_lin)){
            peri_lin[,j] <- beta_cent_lin[j]*A + stats::rnorm(n, mean_peri_lin[j], sd_peri_lin[j])
          }
          peri_non <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_non)))
          for(i in 1:length(mean_peri_non)){
            if(i %in% inter_cent_lin){
              peri_non[,i] <- beta_cent_non[i] * A * peri_lin[,i] +
                stats::rnorm(n, mean_peri_non[i], sd_peri_non[i])
            } else {
              peri_non[,i] <- beta_cent_non[i] * A +
                stats::rnorm(n, mean_peri_non[i], sd_peri_non[i])
            }
          }
          C <- beta_non*A + beta_lin*A + stats::rnorm(n, mean_target, sd_target)
          for(m in 1:length(inter_cent_lin)){
            C = C + beta_non*A*peri_lin[,inter_cent_lin[m]]
          }
          for(l in 1:length(inter_cent_lin)){
            C = C + beta_peri_non[l]*peri_non[,l]
          }
        } else if(con_func == "log"){
          peri_non <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_non)))
          for(i in 1:length(mean_peri_non)){
            peri_non[,i] <- data.frame(beta_cent_non[i]*log(abs(A)) + stats::rnorm(n, mean_peri_non[i], sd_peri_non[i]))
          }
          peri_lin <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_lin)))
          for(j in 1:length(mean_peri_lin)){
            peri_lin[,j] <- beta_cent_lin[j]*A + stats::rnorm(n, mean_peri_lin[j], sd_peri_lin[j])
          }
          C <- beta_non*A + beta_lin*A + stats::rnorm(n, mean_target, sd_target)
          for(k in 1:ncol(peri_lin)){
            C = C + beta_peri_lin[k]*peri_lin[k]
          }
          for(l in 1:ncol(peri_non)){
            C= C + beta_peri_non[l]*peri_non[l]
          }}
  }
  else if(func == "log"){
    if(con_func == "linear"){
      peri_non <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_non)))
      for(i in 1:length(mean_peri_non)){
        peri_non[,i] <- data.frame(beta_cent_non[i]*A + stats::rnorm(n, mean_peri_non[i], sd_peri_non[i]))
      }
      peri_lin <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_lin)))
      for(j in 1:length(mean_peri_lin)){
        peri_lin[,j] <- beta_cent_lin[j]*A + stats::rnorm(n, mean_peri_lin[j], sd_peri_lin[j])
      }
      C <- beta_non*log(abs(A)) + beta_lin*A + stats::rnorm(n, mean_target, sd_target)
      for(k in 1:ncol(peri_lin)){
        C = C + beta_peri_lin[k]*peri_lin[k]
      }
      for(l in 1:ncol(peri_non)){
        C= C + beta_peri_non[l]*peri_non[l]
      }} else if(con_func == "quad"){
        peri_non <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_non)))
        for(i in 1:length(mean_peri_non)){
          peri_non[,i] <- data.frame(beta_cent_non[i]*A^2 + stats::rnorm(n, mean_peri_non[i], sd_peri_non[i]))
        }
        peri_lin <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_lin)))
        for(j in 1:length(mean_peri_lin)){
          peri_lin[,j] <- beta_cent_lin[j]*A + stats::rnorm(n, mean_peri_lin[j], sd_peri_lin[j])
        }
        C <- beta_non*log(abs(A)) + beta_lin*A + stats::rnorm(n, mean_target, sd_target)
        for(k in 1:ncol(peri_lin)){
          C = C + beta_peri_lin[k]*peri_lin[k]
        }
        for(l in 1:ncol(peri_non)){
          C= C + beta_peri_non[l]*peri_non[l]
        }} else if(con_func == "inter"){
          peri_lin <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_lin)))
          for(j in 1:length(mean_peri_lin)){
            peri_lin[,j] <- beta_cent_lin[j]*A + stats::rnorm(n, mean_peri_lin[j], sd_peri_lin[j])
          }
          peri_non <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_non)))
          peri_non <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_non)))
          for(i in 1:length(mean_peri_non)){
            if(i %in% inter_cent_lin){
              peri_non[,i] <- beta_cent_non[i] * A * peri_lin[,i] +
                stats::rnorm(n, mean_peri_non[i], sd_peri_non[i])
            } else {
              peri_non[,i] <- beta_cent_non[i] * A +
                stats::rnorm(n, mean_peri_non[i], sd_peri_non[i])
            }
          }
          C <- beta_non*log(abs(A)) + beta_lin*A + stats::rnorm(n, mean_target, sd_target)
          for(m in 1:length(inter_cent_lin)){
            C = C + beta_non*A*peri_lin[,inter_cent_lin[m]]
          }
          for(l in 1:length(inter_cent_lin)){
            C = C + beta_peri_non[l]*peri_non[,l]
          }
        } else if(con_func == "log"){
          peri_non <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_non)))
          for(i in 1:length(mean_peri_non)){
            peri_non[,i] <- data.frame(beta_cent_non[i]*log(abs(A)) + stats::rnorm(n, mean_peri_non[i], sd_peri_non[i]))
          }
          peri_lin <- data.frame(matrix(NA, nrow = n, ncol = length(mean_peri_lin)))
          for(j in 1:length(mean_peri_lin)){
            peri_lin[,j] <- beta_cent_lin[j]*A + stats::rnorm(n, mean_peri_lin[j], sd_peri_lin[j])
          }
          C <- beta_non*log(abs(A)) + beta_lin*A + stats::rnorm(n, mean_target, sd_target)
          for(k in 1:ncol(peri_lin)){
            C = C + beta_peri_lin[k]*peri_lin[k]
          }
          for(l in 1:ncol(peri_non)){
            C= C + beta_peri_non[l]*peri_non[l]
          }}
  }
  colnames(peri_lin) <- paste0("B", 1:length(mean_peri_lin))
  colnames(peri_non) <- paste0("D", 1:length(mean_peri_non))
  dat <- cbind(A,C,peri_lin,peri_non)
  colnames(dat)[1] <- "A"
  colnames(dat)[2] <- "C"
  return(dat)
}
