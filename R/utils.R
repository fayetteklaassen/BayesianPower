# Documentation eval_const() ====
#' Evaluate a constraint matrix for a set of prior/posterior samples
#'
#' @param hyp A constraint matrix defining a hypothesis.
#' @param samples A matrix. Prior or posterior samples, the number of columns
#' corresponds to the number of groups, the number of rows the number of samples
#' @return A number between 0 and 1. The proportion of samples in which the
#' constraints are met.
eval_const <- function(hyp, samples) {
  # Errors ====
  if(length(dim(hyp)) != 2) {
    stop("Hypothesis must be a matrix")
  }
  if(length(dim(samples)) != 2) {
    stop("Samples must be in a matrix form")
  }
  if(ncol(samples) != ncol(hyp)) {
    stop("Not the right dimensions")
  }
  # Function ====
  nconst <- nrow(hyp)

  out <- mean(apply(samples, 1, function(i) {
    rowCheck <- apply(hyp, 1, function(j) {
      sum(i * j) > 0
    }
    )
    sum(rowCheck) == nconst
  }
  )
  )
  return(out)
}

# Documentation samp_dist() ====
#' Sample from prior or posterior distribution
#'
#' @param nsamp A number. The number of prior or posterior samples to determine the
#' fit and complexity
#' @param means A vector. The prior or posterior means for each group
#' @param sds A number or a vector. The standard deviations for each group
#' If a number is used, the same prior or posterior standard deviation is
#' used for each group.
#' @return A matrix of \code{nsamp} rows and as many columns as the
#' length of \code{means}.
samp_dist <- function(nsamp, means, sds) {
  # Function ====
  ngroup <- length(means)

  samples <- matrix(stats::rnorm(nsamp * ngroup,
                          mean = means,
                          sd = sds),
                    ncol = ngroup,
                    byrow = TRUE)
  return(samples)
}

# Documentation calc_fc() ====
#' Compute the complexity or fit for two hypotheses.
#' @param hyp A constraint matrix defining H1.
#' @param hyp2 A constraint matrix defining H2 OR a character \code{'u'}
#' or \code{'c'} specifying an unconstrained or complement hypothesis
#' @param means A vector of posterior or prior means
#' @param sds A vector or posterior or prior standard deviation
#' @param nsamp A number. The number of prior or posterior samples to determine the
#' fit and complexity
#' @return A vector.
#' The proportion of posterior samples in agreement with H1 and with H2
calc_fc <- function(hyp, hyp2, means, sds, nsamp = 1000) {
  # Errors ====
  if(is.numeric(hyp2)){
    if(length(dim(hyp2)) != 2) stop("Hypothesis must be a matrix")
  } else {
    if(hyp2 != "c" && hyp2 != "u") stop("Use 'u' or 'c' for unconstrained or complement.")
  }
  if(length(dim(hyp)) != 2) {
    stop("Hypothesis must be a matrix")
  }
  # Function ====
  ngroup <- length(means)
  samples <- samp_dist(nsamp, means = means, sds = sds)
  fc1 <- eval_const(hyp, samples)
  if (fc1 == 0) fc1 <- 1/nsamp
  fc2 <- if (is.character(hyp2)) {
    if (hyp2 == 'u') 1
    if (hyp2 == 'c') 1 - fc1
  } else {
    eval_const(hyp2, samples)
  }
  if (fc2 == 0) fc2 <- 1 / nsamp
  out <- c(fc1, fc2)
  return(out)
}

# Documentation calc_bf()====
#' Compute a Bayes factor
#'
#' @param data A matrix. The dataset for which the BF must be computed
#' @param h1 A constraint matrix defining H1.
#' @param h2 A constraint matrix defining H2.
#' @param scale A number specifying the prior scale.
#' @param nsamp A number. The number of prior or posterior samples to determine the
#' @return BF12, that is, the evidence for H1 relative to H2
calc_bf <- function(data, h1, h2, scale, nsamp = 1000) {
  # Function ====
  postmeans <- colMeans(data)
  postsds <- apply(data, 1, stats::sd) / sqrt(nrow(data))
  ngroup <- length(postmeans)
  priormeans <- rep(0, ngroup)
  priorsds <- postsds*scale
  comp <- calc_fc(h1, h2, priormeans, priorsds, nsamp)
  fit <- calc_fc(h1, h2, postmeans, postsds, nsamp)
  bf <- (fit[1] / comp[1]) / (fit[2] / comp[2])
  return(bf)
}

# Documentation samp_bf() ====
#' Sample multiple datasets and compute the Bayes factor in each
#'
#' @param datasets A number. The number of datasets to simulate for each
#' sample size \code{n}
#' @param n A number. The group sample size to be used in data simulation
#' @param ngroup A number. The number of groups.
#' @param means A vector of expected population means.
#' @param sds A vector of expected population standard deviations
#' Note, when standardized, this is a vector of 1s
#' @param h1 A constraint matrix defining H1.
#' @param h2 A constraint matrix defining H2.
#' @param scale A number specifying the prior scale.
#' @param nsamp A number. The number of samples for the fit and complexity
#' See \code{?BayesianPower::calc_fc}
#' @return A vector of Bayes factors BF12 for each of the simulated datasets
samp_bf <- function(datasets, n, ngroup, means, sds, h1, h2, scale, nsamp) {
  # Function ====
  out <- sapply(1:datasets, function(i){
    data <- matrix(stats::rnorm(n * ngroup, mean = means, sd = sds),
                   ncol = ngroup,
                   byrow = TRUE)
    calc_bf(data = data, h1, h2, scale, nsamp)
  }
  )
  return(out)
}

# Documentation bayes_error() ====
#' Determine the unconditional error probabilities for a set of simulated
#' Bayes factors.
#'
#' @param BFs1 A vector. Simulated BF12 under H1 for a given n
#' @param BFs2 A vector. Simulated BF12 under H2 for a given n
#' @param bound1 A number. The boundary above which BF12 favors H1
#' @param bound2 A number. The boundary below which BF12 favors H2
#' @return A named vector. The Type 1, Type 2, Decision error and Area of Indecision probabilities
#' and the median Bayes factors under H1 and H2
bayes_error <- function(BFs1, BFs2, bound1 = 1, bound2 = 1/bound1) {
  # Function ====
  type1 <- mean(BFs1 < bound1)
  type2 <- mean(BFs2 > bound2)
  de <- (type1 + type2)/2
  aoi <- if (bound1 != bound2) {
    (mean((BFs1 > bound1) & (BFs1 < bound2)) +
       mean((BFs2 > bound1) & (BFs2 < bound2))) / 2
  } else {
    0
  }
  med.2.inv <- 1 / stats::median(BFs2)
  med.1 <- stats::median(BFs1)
  return(c("1" = type1, "2" = type2, "de" = de, "aoi" = aoi,
           "med.1" = med.1, "med.2 inverse" = med.2.inv))
}
