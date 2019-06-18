#
#' Generate smooth graphs for the APF and FDR estimates
#'
#' \code{apf_plot} returns a graph with Average Power Function (APF),
#' Bayes False Discovery Rate (FDR) and APF vs. FDR. In addition, when
#' \code{tab = TRUE}, a table containing APF, FDR, tau and gamma values
#'  for a selected subset of APF and FDR is printed.
#'
#' @param APF_lst The output from the \code{apf_fdr} function.
#' @param tab If \code{TRUE}, a table with relevant values of APF, FDR,
#' tau and gamma is printed.
#' @param APF_inf Sets the minimum value of APF to appear in the table
#' when \code{tab = TRUE}.
#' @param FDR_sup Sets the maximum value of Bayes FDR to appear in the
#' table when \code{tab = TRUE}.
#'
#' @return Smooth graphs for APF vs Gamma (left), FDR vs Gamma (centre)
#' and APF vs FDR (right). Regions where FDR \eqn{\le} \code{FDR_sup} and
#' APF \eqn{\ge} \code{APF_inf} (if presents) are highlighted in yellow
#' and printed in a table (if \code{tab = TRUE}) together with the relative
#' values of \emph{gamma} and \emph{tau}.
#'
#' @examples
#' data("Ex2")
#' apf_plot(Ex2)
#' # Ex2 is an example of output obtained
#' # from the apf_fdr() function.
#'
#' @references Quatto, P, Margaritella, N, et al. Brain networks
#' construction using Bayes FDR and average power function. \emph{Stat Methods Med Res}.
#' Published online May 14th, 2019; DOI:10.1177/0962280219844288.
#'
#' @export
#
apf_plot <- function(APF_lst, tab = TRUE, APF_inf = 0.5,
                     FDR_sup = 0.05){
  graphics::par(mfrow = c(1, 3))
  graphics::plot(APF_lst$gammaval, stats::smooth.spline(APF_lst$APF_gamma,
       spar = 0.8)$y, type = "l", lwd = 2,
       main = "Average Power Function (APF)", ylab = "APF", xlab = "Gamma")
  graphics::plot(APF_lst$gammaval, stats::smooth.spline(APF_lst$FDR_gamma,
       spar = 0.8)$y, type = "l", lwd = 2, main = "Bayes FDR", ylab = "FDR",
       xlab = "Gamma")
  graphics::plot(stats::smooth.spline(APF_lst$FDR_gamma, spar = 0.8)$y,
       stats::smooth.spline(APF_lst$APF_gamma, spar = 0.8)$y, type = "n",
       lwd = 2, main = "APF vs. Bayes FDR", xlab = "Bayes FDR", ylab = "APF")
  graphics::rect(0, 0.5, 0.05, 1, col = "yellow", border = NA, density = 4)
  graphics::lines(stats::smooth.spline(APF_lst$FDR_gamma, spar = 0.8)$y,
       stats::smooth.spline(APF_lst$APF_gamma, spar = 0.8)$y, lwd = 2)
  graphics::abline(h = 0.5, lwd = 2, lty = 2, col = "blue")
  graphics::abline(v = c(0.05, 0.1), lwd = 2, lty = 2, col = c("red", "orange"))
  #
  if (tab == TRUE) {
   sel  <- which(stats::smooth.spline(APF_lst$APF_gamma, spar = 0.8)$y
                 >= APF_inf & stats::smooth.spline(APF_lst$FDR_gamma,
                                                 spar = 0.8)$y <= FDR_sup)
   if ( length(sel) > 0) {
     gamma_sel <- APF_lst$gammaval[sel]
     tau_sel <- APF_lst$tau_gamma[sel]
     APF_sel <- stats::smooth.spline(APF_lst$APF_gamma,  spar = 0.8)$y[sel]
     FDR_sel <- stats::smooth.spline(APF_lst$FDR_gamma, spar = 0.8)$y[sel]
     results_tab <- as.data.frame(round(cbind(gamma_sel, tau_sel, FDR_sel,
                                       APF_sel), 3))
     colnames(results_tab) <- c("Gamma", "tau", "95th FDR", "5th APF")
     return(results_tab)
   }
  }
}
