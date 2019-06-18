#
#' Implementation of APF and FDR robust estimation
#'
#' \code{apf_fdr} returns robust estimates of the Average Power Function (APF)
#' and Bayes False Discovery Rate (FDR) for each value of the threshold Gamma
#' on the p-value and Tau on the correlation coefficient.
#'
#' @param data Either a vector, matrix or dataframe.
#' @param type Set \code{"datf"} if \code{data} is a matrix or dataframe containing
#'  the raw data (columns); \code{"pvl"} for a vector of p-values.
#' @param corr The type of correlation to use when \code{type = "datf"}. It can be
#'  set to either \code{"spearman"} or \code{"pearson"}.
#' @param lobs When \code{type = "pvl"}, it indicates the number of datapoints used to
#'  compute the correlations.
#' @param seed A seed (natural number) for the resampling.
#' @param gamm The threshold gamma on the p-values to explore (typically \eqn{\le} 0.05 or
#'  0.1). A min, max and step length value need to be set.
#'
#' @return A list with four elements. A vector \code{APF_gamma} containing the robust
#' estimates of APF (5th quantiles) for all the gamma values set in \code{gamm}. A vector
#' \code{FDR_gamma} with the robust estimates of Bayes FDR (95th quantiles). A
#' vector \code{tau_gamma} with the correlation coefficients \emph{tau} for each gamma
#' value explored and another vector with the relative values gamma (\code{gammaval}).
#'
#' @examples
#' \donttest{data("Ex1")
#' APF_lst <- apf_fdr(Ex1,"pvl",lobs=100)}
#' # The example uses the dataset Ex1 (in the APFr package) which is
#' # a vector of 100 p-values. The number of datapoints used to
#' # compute each p-value in this example is set to 100. As a result,
#' # a list of 4 objects is returned.
#' \dontshow{
#' prova <- runif(3)
#' apf_fdr(prova,"pvl",lobs=15,gamm=c(0.01,0.05,0.02))}
#'
#' @references Quatto, P, Margaritella, N, et al. Brain networks
#' construction using Bayes FDR and average power function. \emph{Stat Methods Med Res}.
#' Published online May 14th, 2019; DOI:10.1177/0962280219844288.
#'
#' @export
#
apf_fdr <- function(data, type="datf", corr="spearman", lobs=0,
                     seed=111, gamm=c(0.0001, 0.1, 0.002)) {
#
#--ERRORS:
  if (type == "datf") {
   if (is.data.frame(data) == FALSE && is.matrix(data) == FALSE) {
    stop("data are not in matrix nor data frame form")
   }
   if (dim(data)[1] < 3) {
    stop("not enough observations to compute tests")
     }
    }
  if (type == "pvl" && is.vector(data) == FALSE) {
   stop("data are not in vector form")
  }
  if (gamm[1] < 0 || gamm[2] > 1 || gamm[1] > gamm[2] ||
      gamm[3] < 0 || gamm[3] > 1) {
    stop("gamm values are not in [0,1] or gamm[1] > gamm[2]")
  }
  if (seed < 1 || all.equal(seed, as.integer(seed)) != TRUE){
    stop("seed should be a natural positive number")
  }
#-------
  if (type == "datf") {
   # Perform multiple correlation tests
   lobs <- dim(data)[1]
   corre <- stats::cor(as.matrix(data), method = corr,
                       use = "all.obs")
   vect_cor <- corre[upper.tri(corre)]
   vect_t <- vect_cor * sqrt( (lobs - 2) / (1 - (vect_cor ^ 2)))
   samplep <- (2 - 2 * stats::pt(abs(vect_t), lobs - 2))
   # number of tests
   nt <- (ncol(data)) * (ncol(data) - 1) / 2
   #
   rm(vect_cor, vect_t)
   } else {
      samplep <-  data
      nt  <-  length(samplep)
   }
  if (lobs == 0) {
    stop("no. datapoints (lobs) must be specified when type=pvl")
  }
#NEW
  if (length(which(samplep <= gamm[1])) ==
      length(samplep)){
    stop("all p-values <= gamm[1]")
  }
#
  message("Creating Bayes FDR...")
#1) Resample 1000 times + original data -> lst_p -------------------------------
#
  set.seed(seed)
  resamples <- lapply(1:1000, function(i) sample(samplep, replace = T))
  hh <- length(resamples) + 1
  vuoto <- rep(NA, length(samplep))
  lst_p <- as.list(replicate(hh, vuoto, simplify = F))#
  lst_p[[1]] <- samplep
  for (k in 2:length(lst_p)) {
    lst_p[[k]] <- resamples[[k - 1]]
  }
    rm(k, hh, vuoto)
#2) Count # p-values below LAMBDA (lambd) which varies between 0.01 and 0.98
#NEW
  if (length(which(samplep <= 0.0001)) == length(samplep)){
    lambd <- seq(from = 0.0000001, to = 0.0001, by = 0.000001)
  } else {
    lambd <- seq(from = 0.0001, to = 0.9901, by = 0.01)
  }
  R_lambd <- rep(NA, length(lambd))
  lst_R_lambd <- as.list(replicate(length(lst_p), R_lambd, simplify = F))
  for (q in 1:length(lst_p)) {
   for (i in 1:length(lambd)) {
     lst_R_lambd[[q]][i] <-
       max(sum(as.numeric(lst_p[[q]] <= lambd[i]), na.rm = TRUE), 1)
     }
    #NEW
    if (length(which(lst_R_lambd[[q]] == length(samplep))) > 0 &&
        min(which(lst_R_lambd[[q]] == length(samplep))) == 1){
      stop("all p-values <= 1e-07")
      }
  }
  rm(R_lambd, i, q)

#3) Avoid pi_0=0 by carrying over the last R_lambd<length(samplep)
  for (q in 1:length(lst_p)){
    if (length(which(lst_R_lambd[[q]] == length(samplep))) > 0){
      lst_R_lambd[[q]][which(lst_R_lambd[[q]] == length(samplep))] <-
        lst_R_lambd[[q]][min(which(lst_R_lambd[[q]] == length(samplep))) - 1]
      }
  }
#
#5) From lst_R_lambd to lst_W_lambd
  lst_W_lambd <- lapply(lst_R_lambd, function(x) nt - x)
#
#6) Estimate pi_0(lambd)
  lst_pg0_lambd <- lapply(lst_W_lambd, function(x) x /
                            ( (1 - lambd) * nt))
#
#7) R_gamm: no. p-values <= threshold gammval
# (gammaval can be changed a posteriori to focus on a specific
# subset of the gamma range)
#
  gammaval <- seq(from = gamm[1], to = gamm[2], by = gamm[3])
  R_gamm <- rep(NA, length(gammaval) )
  lst_R_gamm <- as.list(replicate(length(lst_p), R_gamm, simplify = F))
  for (q in 1:length(lst_p) ) {
    for (i in 1:length(gammaval) ) {
      lst_R_gamm[[q]][i] <- max(sum(as.numeric(lst_p[[q]] <= gammaval[i]),
                                    na.rm = TRUE), 1)
    }
  }
  rm(R_gamm, i, q)
#
# Estimate P(p-values<= gammaval).
  lst_Pr_P_le_gamm <- lapply(lst_R_gamm, function(x) x / nt)
#
# 8) FDR-----------------------------------------------------------------------
  # create a list with FDR matrices of size lambd X gammaval
  mtx <- matrix(0, nrow = length(lambd), ncol = length(gammaval))
  lst_fdr <- as.list(replicate(length(lst_p), mtx, simplify = F))
  for (q in 1:length(lst_p) ) {
    for (i in 1:length(lambd) ) {
      for (j in 1:length(gammaval) ) {
        lst_fdr[[q]][i, j] <- min( ( ( (lst_pg0_lambd[[q]][i]) * gammaval[j]) /
                                   lst_Pr_P_le_gamm[[q]][j]), 1)
      }
    }
  }
  rm(i, j, q, mtx)
# 9) bootstrap plugin FDR------------------------------------------------------
# compare lst_fdr[[i]], i=2,...,1001 with lst_fdr[[1]] to obtain MSE
# Find FDR that minimises MSE for every gammaval
#
# MSE matrix
  hh <- length(resamples) + 1
  vect <- rep(0, length(resamples))
  MSE <- matrix(0, length(lambd), length(gammaval))
  for (i in 1:length(lambd)) {
    for (j in 1:length(gammaval)) {
      vect <- sapply(lst_fdr[2 : hh], "[[", i, j)
      MSE[i, j] <- (1 / length(resamples)) *
        sum( (vect - lst_fdr[[1]][i, j]) ^ 2)
    }
  }
  rm(hh, vect, i, j)
#
#
  fdr_g <- rep(0, length(lst_p))
  lst_fdr_gamm1 <- as.list(replicate(length(fdr_g), rep(0, length(gammaval)),
                                     simplify = F))
  lst_1 <- as.list(replicate(length(gammaval),
                             matrix(0, nrow = length(lambd), ncol = 3),
                             simplify = F))
  lst_MSE_fdr <- as.list(replicate(length(lst_p), lst_1, simplify = F))
  lst_MSE_fdr_ord <- rep(lst_MSE_fdr, 1)
  cont_lambd <- c(1 : length(lambd))
  cont_lambd_mat <- matrix(cont_lambd, nrow = length(lambd),
                           ncol = length(gammaval))
#
  for (j in 1 : length(lst_p)) {
    for (i in 1 : length(gammaval)) {
      lst_MSE_fdr[[j]][[i]] <- as.matrix(cbind(MSE[, i], lst_fdr[[j]][, i],
                                               cont_lambd_mat[, i]))
    }
  }
  for (j in 1 : length(lst_p)) {
    for (i in 1 : length(gammaval)) {
      lst_MSE_fdr_ord[[j]][[i]] <-
        lst_MSE_fdr[[j]][[i]][order(lst_MSE_fdr[[j]][[i]][, 1]), ]
    }
  }
  for (j in 1 : length(lst_p)) {
    lst_fdr_gamm1[[j]] <- sapply(lst_MSE_fdr_ord[[j]], "[[", 1, 2)
  }
#
  lst_fdr_gamm2 <- as.list(replicate(length(gammaval), fdr_g, simplify = F))
#
  for (i in 1 : length(gammaval)) {
    lst_fdr_gamm2[[i]] <- sapply(lst_fdr_gamm1, "[[", i)
  }
#
# For every lst_fdr_gamm2[[i]] there are 1001 resamples of FDR.
# Extract the 95th quantile.
  FDR_gamma <- rep(0, length(gammaval))
  for (i in 1 : length(gammaval)) {
    FDR_gamma[i] <- stats::quantile(sort(lst_fdr_gamm2[[i]]), .95)
  }
#
  rm(lst_MSE_fdr, lst_MSE_fdr_ord, lst_fdr_gamm1, lst_fdr_gamm2, lst_1, fdr_g)
#
#### FDR_gamma contains 95th FDR percentile for every gammaval###
#
# 10) Average Power Function (APF) --------------------------------------------
#
  message("Creating APF...")
# APF numerator
  lst_APF1 <- as.list(replicate(length(lst_p),
                                matrix(0, nrow = length(lambd),
                                       ncol = length(gammaval)), simplify = F))
  for (q in 1 : length(lst_p)) {
    for (j in 1 : length(gammaval)) {
      for (i in 1 : length(lambd)) {
        lst_APF1[[q]][i, j] <- (lst_R_gamm[[q]][j] -
                                  (gammaval[j] * lst_pg0_lambd[[q]][i] * nt))
      }
    }
  }
#
# APF denominator
  lst_APF2 <- lapply(lst_pg0_lambd, function(x) (1 - x) * nt)
#
  lst_APF <- as.list(replicate(length(lst_p),
                               matrix(0, nrow = length(lambd),
                                      ncol = length(gammaval)), simplify = F))
  for (q in 1 : length(lst_p)) {
    for (i in 1 : length(lambd)) {
      for (j in 1: length(gammaval)) {
        lst_APF[[q]][i, j] <- ( (lst_APF1[[q]][i, j]) / (lst_APF2[[q]][i]))
      }
    }
  }
  for (q in 1 : length(lst_p)) {
    for (i in 1 : length(lambd)) {
      for (j in 1: length(gammaval)) {
        lst_APF[[q]][i, j] <- min(lst_APF[[q]][i, j], 1)
      }
    }
  }
  rm(lst_APF1, lst_APF2)
#
# 11) bootstrap plugin APF------------------------------------------------------
#
# MSE matrix
  hh <- length(resamples) + 1
  vect <- rep(0, length(resamples))
  MSE2 <- matrix(0, length(lambd), length(gammaval))
  for (i in 1:length(lambd)) {
    for (j in 1:length(gammaval)) {
      vect <- sapply(lst_APF[2 : hh], "[[", i, j)
      MSE2[i, j] <- (1 / length(resamples)) *
        sum( (vect - lst_APF[[1]][i, j]) ^ 2)
    }
  }
  rm(hh, vect, i, j)
#
#
  APF_g <- rep(0, length(lst_p))
  lst_APF_gamm1 <- as.list(replicate(length(APF_g),
                                     rep(0, length(gammaval)), simplify = F))
  lst_1 <- as.list(replicate(length(gammaval),
                             matrix(0, nrow = length(lambd), ncol = 3),
                             simplify = F))
  lst_MSE2APF <- as.list(replicate(length(lst_p), lst_1, simplify = F))
  lst_MSE2APF_ord <- rep(lst_MSE2APF, 1)
#
  for (j in 1 : length(lst_p)) {
    for (i in 1 : length(gammaval)) {
      lst_MSE2APF[[j]][[i]] <- as.matrix(cbind(MSE[, i], lst_APF[[j]][, i],
                                               cont_lambd_mat[, i]))

    }
  }
  for (j in 1 : length(lst_p)) {
    for (i in 1 : length(gammaval)) {
      lst_MSE2APF_ord[[j]][[i]] <-
        lst_MSE2APF[[j]][[i]][order(lst_MSE2APF[[j]][[i]][, 1]), ]
    }
  }
  for (j in 1 : length(lst_p)) {
    lst_APF_gamm1[[j]] <- sapply(lst_MSE2APF_ord[[j]], "[[", 1, 2)
  }
#
  lst_APF_gamm2 <- as.list(replicate(length(gammaval), APF_g, simplify = F))
#
  for (i in 1 : length(gammaval)) {
  lst_APF_gamm2[[i]] <- sapply(lst_APF_gamm1, "[[", i)
  }
#
# For each lst_APF_gamm2[[i]] there are 1001 resamples of APF
# Extract the 5th percentile
  APF_gamma <- rep(0, length(gammaval))
  for (i in 1 : length(gammaval)) {
    APF_gamma[i] <- stats::quantile(sort(lst_APF_gamm2[[i]]), .05)
  }
  rm(APF_g, lst_MSE2APF, lst_MSE2APF_ord, lst_APF_gamm1, lst_APF_gamm2, lst_1)
#
#### APF_gamma contains 5th APF percentile for every gammaval###
#
# 12) compute tau from gamma values.
  tau_gamma <- abs(stats::qt(gammaval / 2, lobs - 2) /
                   sqrt(lobs - 2 + (stats::qt(gammaval / 2, lobs - 2)) ^ 2))
#
# 13) create the final list
  APF_lst <- list(APF_gamma, FDR_gamma, tau_gamma, gammaval)
  names(APF_lst) <- c("APF_gamma", "FDR_gamma", "tau_gamma", "gammaval")
  message("Done!")
  return(APF_lst)
}
