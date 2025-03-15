# SC
This R-Package is used for statistical computing (i.e., sampling methods, bayesian inference, Optimizations and EM)

Provides a comprehensive suite of tools for Monte Carlo methods, Markov Chain Monte Carlo (MCMC) sampling, and Bayesian inference. Implements various sampling techniques such as Gibbs Sampling, the Metropolis algorithm, and the Metropolis-Hastings algorithm. Supports random number generation, expectation computation, and resampling methods. Also includes optimization methods such as Expectation-Maximization (EM) and applications of Gaussian Processes for probabilistic modeling.


To use, follow the following instruction:

install.packages("devtools")
devtools::install_github(stevenchen333/SC)


DO NOT use require() or library() to call for my package it will break (currently a bug wll be fixed soon)

use SC::<the function that you want> (e.g. SC:: inverse_cdf_sampling())
