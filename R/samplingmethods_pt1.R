
#------------------------------------

#' Inverse CDF Sampling
#'
#' This function generates random samples from a given continuous probability distribution
#' using the inverse cumulative distribution function (CDF) method.
#'
#' @param n The number of samples to generate.
#' @param inv A function representing the inverse CDF of the target distribution.
#' @return A numeric vector of sampled values.
#' @examples
#' # Define the inverse CDF for an Exponential Distribution
#' inv_cdf_exp <- function(u) { -log(1 - u) }
#'
#' # Generate 1000 samples using Inverse CDF Sampling
#' samples <- inv_cdf(inv_cdf_exp, n = 1000)
#'
#' # Plot the comparison between the exponential PDF and sampled values
#' plot_sampling(samples, function(x) dexp(x))  # Compare with exponential PDF
#'
#' @export
inv_cdf <- function(inv, n = 100000) {
  u <- runif(n)  # Generate uniform samples in (0,1)
  inv(u)     # Apply inverse CDF to transform uniform samples
}

#------------------------------------

#' Box-Muller Normal Sampling
#'
#' This function generates random samples from a standard normal distribution
#' using the Box-Muller transform.
#'
#' @param n The number of samples to generate (must be even for full pairs).
#' @return A numeric vector of sampled values from N(0,1).
#' @examples
#' # Generate 1000 samples using the Box-Muller transform
#' samples <- box_muller(1000)
#'
#' # Plot the comparison between the normal PDF and sampled values
#' plot_sampling(samples, function(x) dnorm(x))  # Compare with normal PDF
#'
#' @export
box_muller <- function(n = 100000) {
  if (n %% 2 == 1) stop("n must be even for full pairs.")

  u1 <- runif(n / 2)
  u2 <- runif(n / 2)

  r <- sqrt(-2 * log(u1))
  theta <- 2 * pi * u2

  z1 <- r * cos(theta)
  z2 <- r * sin(theta)

  return(c(z1, z2))  # Return the combined vector of normal samples
}



#------------------------------------

#' General-Purpose Rejection Sampling
#'
#' @description Samples from arbitrary univariate target distributions using rejection sampling.
#'
#' @param target Target density function (univariate: `function(x)`)
#' @param proposal_sampler Function generating proposal samples (`function(n)`)
#' @param proposal_density Proposal density function (`function(x)`)
#' @param n Number of samples to generate (default: 1)
#' @param M Optional bound for target/proposal ratio (will calculate if NULL)
#' @param lower Lower bound for optimization interval (default: -10)
#' @param upper Upper bound for optimization interval (default: 10)
#'
#' @return A list containing:
#' \itemize{
#'   \item `samples`: Vector of accepted samples
#'   \item `theoretical_efficiency`: Theoretical acceptance probability (1/M)
#'   \item `empirical_efficiency`: Empirical acceptance rate (n_samples/total_attempts)
#'   \item `M`: Used bound value
#' }
#'
#'
#' @examples
#' # Example 1: Sampling from Gamma(3,2) using Normal(1.5,1) proposal
#' target <- function(x) dgamma(x, shape = 3, rate = 2)
#' proposal_density <- function(x) dnorm(x, mean = 1.5, sd = 1)
#' proposal_sampler <- function(n) rnorm(n, mean = 1.5, sd = 1)
#'
#' result <- rejection_sampler(
#'   target = target,
#'   proposal_sampler = proposal_sampler,
#'   proposal_density = proposal_density,
#'   n = 1000,
#'   lower = 0,
#'   upper = 10
#' )
#'
#' # View results
#' hist(result$samples, breaks = 30, main = "Gamma(3,2) Samples")
#' cat("Theoretical efficiency:", result$theoretical_efficiency, "\n")
#' cat("Empirical efficiency:", result$empirical_efficiency, "\n")
#' cat("Used M value:", result$M, "\n")
#'
#' # Example 2: Sampling from Beta(2,5) using Uniform(0,1) proposal
#' target <- function(x) dbeta(x, 2, 5)
#' proposal_density <- function(x) dunif(x, 0, 1)
#' proposal_sampler <- function(n) runif(n, 0, 1)
#'
#' result <- rejection_sampler(
#'   target = target,
#'   proposal_sampler = proposal_sampler,
#'   proposal_density = proposal_density,
#'   n = 1000,
#'   lower = 0,
#'   upper = 1
#' )
#'
#' @export

rejection_sampler<-function(target, proposal_sampler, proposal_density, n, M = NULL, lower = -10, upper = 10){
  sampled =numeric(n)
  accepted = 0  #to calculate empirical efficiency and as a placeholder for our conditionals
  tries = 0 #to calculate empirical efficiency and as a placeholder for our conditionals
  optimize_not = is.null(M)


  ratio_function<-function(x){
    target(x)/proposal_density(x)
  }
  if (optimize_not){
    M = optimise(f = ratio_function, interval = c(upper, lower),maximum = TRUE)$objective
  }
  while(accepted < n){
    u = runif(1)
    x_sample = proposal_sampler(1)
    tries = tries+1

    accept_not = (u<=target(x_sample)/(M*proposal_density(x_sample)))

    if (accept_not){
      sampled[accepted + 1] <- x_sample
      accepted = accepted+1
    }
 }

  empirical_efficiency = n/tries
  theoretical_efficiency = 1/M

  return(list("sampled" = sampled, "theoretical_efficiency" = theoretical_efficiency, "empirical_efficiency" = empirical_efficiency))

}



#------------------------------------
#' Metropolis Algorithm
#'
#' This function implements the Metropolis algorithm for sampling from a given target distribution.
#' Currently, Metropolis-Hastings is not implemented, but the function structure allows future extension.
#'
#' @param target A function representing the target distribution (unnormalized density).
#' @param initial A function that generates the initial sample.
#' @param proposal A function that proposes a new sample given the current sample.
#' @param n Integer, number of iterations for the algorithm. Default is 100000.
#'
#' @return A list containing the sampled values.
#'
#' @examples
#' # Define target distribution (unnormalized density)
#' target <- function(x) { dnorm(x, mean = 0, sd = 1) }
#'
#' # Define proposal function (Normal random walk)
#' proposal <- function(x) { rnorm(1, mean = x, sd = 1) }
#'
#' # Define initial value function
#' initial <- function() { 0 }
#'
#' # Run Metropolis algorithm
#' samples <- metropolis(target, initial, proposal, n = 10000, hastings = FALSE)
#'
#' # Plot the histogram of samples
#' hist(unlist(samples), probability = TRUE, main = "Metropolis Samples",
#'      xlab = "x", col = "lightblue", border = "black")
#'
#' @export
metropolis <- function(target, initial, proposal, n = 100000, hastings = FALSE){
  samples <- list(initial())

  for (i in 1:n) {
    current <- samples[[length(samples)]]
    proposed <- proposal(current)

    if (hastings == FALSE){
      if (runif(1) < target(proposed) / target(current)) {
        samples <- append(samples, list(proposed))
      } else {
        samples <- append(samples, list(current))
      }

    }

    if (hastings == TRUE){
      NULL
    }
  }

  return(samples)

}

#------------------------------------
#' Flexible Gibbs Sampling with Optional Rejection Steps
#'
#' @description General-purpose Gibbs sampler supporting both direct sampling and rejection sampling steps.
#'
#' @param conditional_samplers List of conditional samplers (functions or specifications)
#' @param init Initial parameter vector
#' @param n_iter Total iterations (default: 1000)
#' @param burn_in Burn-in period (default: 100)
#' @param thinning Thinning interval (default: 1)
#'
#' @return Matrix where rows are samples and columns are parameters
#'
#' @export
#'
#' @examples
#' # Bivariate normal example
#' sampler <- list(
#'   function(y) rnorm(1, 0.5*y, sqrt(0.75)),
#'   function(x) rnorm(1, 0.5*x, sqrt(0.75))
#' )
#' samples <- gibbs_sampler(sampler, init = c(0,0))
gibbs_sampler <- function(conditional_samplers, init, n_iter = 1000, burn_in = 100, thinning = 1) {

  # Input validation
  if (length(conditional_samplers) != length(init)) {
    stop("Length of conditional_samplers must match init vector length")
  }

  n_params <- length(conditional_samplers)
  samples <- matrix(NA, nrow = n_iter, ncol = n_params)
  samples[1,] <- init

  for (i in 2:n_iter) {
    current <- samples[i-1,]

    for (j in 1:n_params) {
      cond_sampler <- conditional_samplers[[j]]

      if (is.function(cond_sampler)) {
        # Direct sampling case
        current[j] <- cond_sampler(current[-j])
      } else if (is.list(cond_sampler) && cond_sampler$method == "rejection") {
        # Rejection sampling case
        current[j] <- rejection_sampler(
          target = function(x) cond_sampler$target(x, current[-j]),
          proposal = cond_sampler$proposal,
          proposal_gen = cond_sampler$proposal_gen,
          bounds = cond_sampler$bounds
        )$samples[1]
      } else {
        stop("Invalid sampler specification for parameter ", j)
      }
    }

    samples[i,] <- current
  }

  # Apply burn-in and thinning
  keep <- seq(burn_in + 1, n_iter, by = thinning)
  return(samples[keep,])
}




