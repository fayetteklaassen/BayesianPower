# Bayes_sampsize ====
#' Determine the required sample size for a Bayesian hypothesis test
#'
#' @param h1 A constraint matrix defining H1.
#' @param h2 A constraint matrix defining H2.
#' @param m1 A vector of expected population means under H1 (standardized).
#' @param m2 A vector of expected populations means under H2 (standardized).
#' \code{m2} must be of same length as \code{m1}
#' @param sd1 A vector of standard deviations under H1. Must be a single number (equal
#' standard deviation under all populations), or a vector of the same length as \code{m1}
#' @param sd2 A vector of standard deviations under H2. Must be a single number (equal
#' standard deviation under all populations), or a vector of the same length as \code{m2}
#' @param scale A number specifying the prior scale
#' @param type A character. The type of error to be controlled
#' options are: \code{"1", "2", "de", "aoi", "med.1", "med.2"}
#' @param cutoff A number. The cutoff criterion for type.
#' If \code{type} is \code{"1", "2", "de", "aoi"}, \code{cutoff} must be between 0 and 1
#' If \code{type} is \code{"med.1" or "med.2"}, \code{cutoff} must be larger than 1
#' @param bound1 A number. The boundary above which BF12 favors H1
#' @param bound2 A number. The boundary below which BF12 favors H2
#' @param datasets A number. The number of datasets to compute the error probabilities
#' @param nsamp A number. The number of prior or posterior samples to determine the
#' fit and complexity
#' @param minss A number. The minimum sample size to consider
#' @param maxss A number. The maximum sample size to consider
#' @param seed A number. The random seed to be set
#' @return The sample size for which the chosen type of error probability
#' is at the set cutoff, and the according error probabilities and median Bayes factors
#' @examples
#' # Short computation example NOT SUFFICIENT SAMPLES
#' h1 <- matrix(c(1,-1), nrow= 1, byrow= TRUE)
#' h2 <- 'c'
#' m1 <- c(.4, 0)
#' m2 <- c(0, .1)
#' bayes_sampsize(h1, h2, m1, m2, sd1 = 1, sd2 = 1, scale = 1000,
#' type = "de", cutoff = .125, nsamp = 50, datasets = 50,
#' minss = 40, maxss = 70)
#' @export
bayes_sampsize <- function(h1, h2, m1, m2, sd1 = 1, sd2 = 1, scale = 1000,
                        type = 1, cutoff, bound1 = 1, bound2 = 1 / bound1,
                        datasets = 1000, nsamp = 1000,
                        minss = 2, maxss = 1000, seed = 31) {
  # Errors ====
  if(!is.numeric(cutoff)) stop("expected numeric value")
  if(!is.numeric(minss) || !is.numeric(maxss)) stop("minss and maxss must be integer")

  types <- c("1", "2", "de", "aoi", "med.1", "med.2")
  if(all(types != type)) stop("Incorrect type specified")

  if(type == "med.1" || type == "med.2") {
    if(cutoff < 1) stop("Cutoff must be larger than 1 for controlling median BF")
  } else {
    if(cutoff > 1) stop( "Cutoff must be smaller than 1 for controlling error probabilities")
  }

  if((maxss-minss) < 20 ) stop("Difference between minss and maxss must be at least 20")

  # Function ====
  times <- c(ceiling(log(maxss) / log(minss + 1)) + 1)
  timesT <- 0
  {
    pb <- utils::txtProgressBar(min = 0, max = times)

    set.seed(seed)
    n <- round(maxss / 2)
    oldn <- n + 10

    ngroup = length(m1)

    while ((abs(oldn - n) > 1)) {
      errors <- bayes_power(n = n, h1 = h1, h2 = h2, m1 = m1, m2 = m2,
                            sd1 = sd1, sd2 = sd2,
                            ngroup = ngroup, scale = scale,
                            bound1 = bound1, bound2 = bound2,
                            datasets = datasets, nsamp = nsamp)
      if (errors[[paste(type)]] < cutoff) {
        oldn <- n
        maxss <- n
        n <- round((n + minss) / 2)
      }
      if (errors[[paste(type)]] > cutoff) {
        oldn <- n
        minss <- n
        n <- round((n + maxss) / 2)
      }
      if (errors[[paste(type)]] == cutoff) break

      timesT <- timesT + 1
        utils::setTxtProgressBar(pb, timesT)
    }
    close(pb)
  }
  return(c("Sample size" = as.integer(n), "Type 1 error" = errors[[1]],
           "Type 2 error" = errors[[2]], "Decision error" = errors[[3]],
           "Indecision error" = errors[[4]],
           "Median BF12 under H1" = errors[[5]],
           "Median BF21 under H2"= 1 / errors[[6]],
         "10quantile under H1" = errors[[7]],
         "90quantile under H1" = errors[[8]],
         "10quantile under H2" = errors[[9]],
         "90quantile under H2" = errors[[10]]))
}

# Bayes_power ====
#' Determine the 'power' for a Bayesian hypothesis test
#'
#' @param n A number. The sample size
#' @param h1 A constraint matrix defining H1
#' @param h2 A constraint matrix defining H2
#' @param m1 A vector of expected population means under H1
#' @param m2 A vector of expected populations means under H2
#' \code{m2} must be of same length as \code{m1}
#' @param sd1 A vector of standard deviations under H1. Must be a single number (equal
#' standard deviation under all populations), or a vector of the same length as \code{m1}
#' @param sd2 A vector of standard deviations under H2. Must be a single number (equal
#' standard deviation under all populations), or a vector of the same length as \code{m2}
#' @param ngroup A number or \code{NULL}. The number of groups
#' If \code{NULL} the number of groups is determined from the length of \code{m1}
#' @param scale A number specifying the prior scale
#' @param bound1 A number. The boundary above which BF12 favors H1
#' @param bound2 A number. The boundary below which BF12 favors H2
#' @param datasets A number. The number of datasets to compute the error probabilities
#' @param nsamp A number. The number of prior or posterior samples to determine the
#' fit and complexity
#' @param seed A number. The random seed to be set
#' @return The Type 1, Type 2, Decision error and Area of Indecision probability and
#' the median BF12s under H1 and H2
#' @examples
#' # Short example WITH SMALL AMOUNT OF SAMPLES
#' h1 <- matrix(c(1,-1,0,0,1,-1), nrow= 2, byrow= TRUE)
#' h2 <- "c"
#' m1 <- c(.4,.2,0)
#' m2 <- c(.2,0,.1)
#' bayes_power(40, h1, h2, m1, m2, datasets = 50, nsamp = 50)
#' @export
bayes_power <- function(n, h1, h2, m1, m2,
                        sd1 = 1, sd2 = 1,
                        ngroup = NULL, scale = 1000,
                        bound1 = 1, bound2 = 1/bound1,
                        datasets = 1000, nsamp = 1000,
                        seed = NULL){

  # Errors ====
  if(!is.numeric(bound1) || !is.numeric(bound2) ||
     !is.numeric(m1) || !is.numeric(m2) || !is.numeric(h1) ||
     !is.numeric(datasets) || !is.numeric(n) ||
     !is.numeric(nsamp) || !is.numeric(sd1) || !is.numeric(sd2) ||
     !is.numeric(scale)) {
    stop("expected numeric value")
  }
  if(round(n) != n) stop("n must be integer")
  if(length(dim(h1)) != 2) stop("h1 must be a matrix")
  if(any(round(h1) != h1)) stop("h1 can only contain integers")
  if(any(h1 > 1) || any(h1 < -1)) stop( "h1 must consist of 1 -1 and 0")

  if(is.numeric(h2)){
    if(length(dim(h2)) != 2) stop("h2 must be a matrix")
    if(any(round(h2) != h2)) stop("h2 can only contain integers")
    if(any(h2 > 1) || any(h2 < -1)) stop( "h2 must consist of 1 -1 and 0")
  } else {
    if(h2 != "c" && h2 != "u") stop("Use 'u' or 'c' for unconstrained or complement.")
  }

  if(length(m1) != length(m2)) stop("m1 and m2 must be the same length")
  if(length(m1) != ncol(h1)) stop("h1 and m1 do not match")

  if(length(sd1) != 1){
    if(length(sd1 != length(m1))) stop("sd1 must be a single number or a vector of same length as m1")
  }
  if(length(sd2) != 1){
    if(length(sd2 != length(m2))) stop("sd2 must be a single number or a vector of same length as m2")
  }

  if(bound1 < bound2) stop("bound1 must be larger than bound2")

   # Function ====
  if(!is.null(seed)) {
    if(!is.integer(seed)) stop("seed must be integer")
    set.seed(seed)
  }
  if(is.null(ngroup)) {
    ngroup <- length(m1)
  } else {
    if(!is.integer(ngroup)) stop("ngroup must be integer")
  }

  BF1 <- samp_bf(datasets, n, ngroup, means = m1, sds = sd1, h1, h2, scale, nsamp)
  BF2 <- samp_bf(datasets, n, ngroup, means = m2, sds = sd2, h1, h2, scale, nsamp)
  errors <- bayes_error(BF1, BF2, bound1, bound2)
  return(errors)
}
