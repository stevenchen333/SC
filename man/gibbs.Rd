% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samplingmethods_pt1.R
\name{gibbs}
\alias{gibbs}
\title{Gibbs Sampling MCMC Algorithm}
\usage{
gibbs(conditional_samplers, n, init = NULL, burn_in = 0, thinning = 1)
}
\arguments{
\item{conditional_samplers}{A list of functions where each function takes the current
state vector and returns a sample from its conditional distribution}

\item{n}{Number of samples to return (after burn-in and thinning)}

\item{init}{Optional initial values (defaults to vector of 0s)}

\item{burn_in}{Number of initial samples to discard (default: 0)}

\item{thinning}{Keep only every k-th sample (default: 1, no thinning)}
}
\value{
A matrix where each row represents one sample from the multivariate distribution,
with columns corresponding to the variables in the order of conditional_samplers
}
\description{
Implements the Gibbs sampling algorithm, a Markov Chain Monte Carlo (MCMC) method
for sampling from multivariate probability distributions when direct sampling is difficult.
}
\examples{
# Example 1: Bivariate Normal Distribution
samplers <- list(
  function(state) rnorm(1, mean = 0.5 * state[2], sd = sqrt(1 - 0.5^2)),  # X|Y ~ N(0.5Y, 0.75)
  function(state) rnorm(1, mean = 0.5 * state[1], sd = sqrt(1 - 0.5^2))   # Y|X ~ N(0.5X, 0.75)
)
samples <- gibbs(samplers, n = 1000, init = c(0, 0), burn_in = 200)
plot(samples, main = "Bivariate Normal Samples")

# Example 2: Using Rejection Sampling for One Conditional
target_y <- function(y, x) {
  ifelse(y > 0 & y < 10, exp(-abs(x) * y), 0)
}
conditional_y <- function(state) {
  x <- state[1]
  rejection_sampler(
    target = function(y) target_y(y, x),
    proposal_sampler = function(n) runif(n, 0, 10),
    proposal_density = function(y) dunif(y, 0, 10),
    n = 1,
    lower = 0,
    upper = 10
  )$sampled
}
samplers <- list(
  function(state) rnorm(1, mean = -state[2], sd = sqrt(1/2)),  # X|Y ~ N(-Y, 1/2)
  conditional_y                                               # Y|X via rejection
)
samples <- gibbs(samplers, n = 500, burn_in = 100, thinning = 2)

# Example 3: Three-Variable System
samplers <- list(
  function(state) rnorm(1, mean = (state[2] + state[3])/2),  # X|Y,Z
  function(state) rexp(1, rate = abs(state[1]) + 1),        # Y|X,Z
  function(state) rpois(1, lambda = 2 + abs(state[1]))      # Z|X,Y
)
samples <- gibbs(samplers, n = 1000)
pairs(samples, labels = c("X", "Y", "Z"))

}
