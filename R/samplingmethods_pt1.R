
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

#' Rejection Sampler for Univariate and Multivariate Distributions
#'
#' Samples from target probability densities using rejection sampling.
#' Works for both univariate and multivariate distributions.
#'
#' @param target Target density function (must return scalar for input vector)
#' @param proposal_sampler Function that generates proposal samples
#' @param proposal_density Proposal density function (must match sampler)
#' @param n Number of desired samples
#' @param M Optional upper bound for target/(M*proposal). If NULL, computes automatically.
#' @param lower Lower bound(s) for optimization when finding M
#' @param upper Upper bound(s) for optimization when finding M
#'
#' @return List with three components:
#' \itemize{
#'   \item sampled - Matrix of samples (n rows × d columns)
#'   \item theoretical_efficiency - 1/M (optimal acceptance rate)
#'   \item empirical_efficiency - Actual acceptance rate (n/total_tries)
#' }
#'
#' @examples
#' # =============================================
#' # Example 1: Univariate - Beta(2,5) Distribution
#' # =============================================
#'
#' # Define components for Beta(2,5) sampling
#' target_beta <- function(x) {
#'   ifelse(x >= 0 & x <= 1, dbeta(x, 2, 5), 0)
#' }
#'
#' set.seed(123)
#' beta_result <- rejection_sampler(
#'   target = target_beta,
#'   proposal_sampler = function(n) runif(n, 0, 1),
#'   proposal_density = function(x) dunif(x, 0, 1),
#'   n = 1000,
#'   lower = 0,
#'   upper = 1
#' )
#'
#' # Plot results
#' hist(beta_result$sampled, freq = FALSE, col = "lightblue",
#'      main = "Beta(2,5) Samples via Rejection Sampling")
#' curve(dbeta(x, 2, 5), add = TRUE, col = "red", lwd = 2)
#'
#' # =============================================
#' # Example 2: Multivariate - Bivariate Normal
#' # =============================================
#'
#' # Define components for standard bivariate normal
#' target_normal <- function(x) {
#'   exp(-0.5 * sum(x^2)) / (2 * pi)
#' }
#'
#' set.seed(456)
#' mvnormal_result <- rejection_sampler(
#'   target = target_normal,
#'   proposal_sampler = function(n) matrix(runif(n*2, -5, 5), ncol = 2),
#'   proposal_density = function(x) prod(dunif(x, -5, 5)),
#'   n = 1000,
#'   lower = c(-5, -5),
#'   upper = c(5, 5)
#' )
#'
#' # Plot results
#' plot(mvnormal_result$sampled, pch = 20, col = rgb(0,0,1,0.3),
#'      main = "Bivariate Normal Samples",
#'      xlab = "X1", ylab = "X2")
#' points(0, 0, pch = "X", col = "red", cex = 2)
#'
#' @export
rejection_sampler <- function(target, proposal_sampler, proposal_density, n, M = NULL) {
  # Get dimension from proposal sampler
  test_sample <- proposal_sampler(1)
  d <- length(test_sample)  # Determine dimension

  # Set lower and upper bounds AFTER determining d
  lower <- rep(-100, d)
  upper <- rep(100, d)

  accepted <- 0
  sampled <- matrix(NA, nrow = n, ncol = d)
  tries <- 0
  optimize_not <- is.null(M)

  # Function to find maximum ratio
  ratio_function <- function(x) {
    target_val <- target(x)
    proposal_val <- proposal_density(x)
    if (proposal_val <= 0 || !is.finite(target_val) || !is.finite(proposal_val)) {
      return(0)  # Fallback for invalid cases
    }
    target_val / proposal_val
  }


  # Find M if not provided
  if (optimize_not) {
    opt_result <- optim(
      par = rep(0, d),
      fn = function(x) -ratio_function(x),  # Minimize negative ratio
      lower = lower,
      upper = upper,
      method = "L-BFGS-B"
    )

    if (opt_result$convergence != 0) {
      stop("Optimization failed: ", opt_result$message)
    }

    M <- -opt_result$value  # Convert back to maximum
  }

  # Sampling loop
  while (accepted < n) {
    x_sample <- proposal_sampler(1)  # Get proposal (should be a vector)
    u <- runif(1)
    tries <- tries + 1

    # Calculate acceptance probability
    accept_prob <- target(x_sample) / (M * proposal_density(x_sample))

    if (u <= accept_prob) {
      accepted <- accepted + 1
      sampled[accepted, ] <- x_sample
    }
  }

  # Calculate efficiencies
  empirical_efficiency <- n / tries
  theoretical_efficiency <- 1 / M

  return(list(
    sampled = sampled,
    theoretical_efficiency = theoretical_efficiency,
    empirical_efficiency = empirical_efficiency
  ))
}


#------------------------------------
#' Metropolis Algorithm for Multivariate Distributions
#'
#' Samples from a target probability density using the Metropolis algorithm
#' (special case of Metropolis-Hastings with symmetric proposals).
#'
#' @param target Target density function (must accept and return numeric)
#' @param proposal_sampler Function that generates symmetric proposals (must return same dimension as init)
#' @param init Initial values (defaults to vector of 0s)
#' @param n_iter Total number of iterations
#' @param burn_in Number of initial samples to discard (default: 0)
#' @param thinning Keep only every k-th sample (default: 1)
#' @return Matrix of samples (n rows × d columns)
#'
#' @examples
#' # Bivariate normal target
#' target <- function(x) exp(-0.5 * sum(x^2)) / (2*pi)
#'
#' # Symmetric normal proposal (sd = 0.5 for both dimensions)
#' proposal_sampler <- function(n, current) {
#'   matrix(rnorm(n*2, mean = 0, sd = 0.5), ncol = 2)
#' }
#'
#' samples <- metropolis(
#'   target = target,
#'   proposal_sampler = proposal_sampler,
#'   init = c(0, 0),
#'   n_iter = 5000,
#'   burn_in = 1000,
#'   thinning = 2
#' )
#'
#' plot(samples, main = "Bivariate Normal Samples", pch = 20)
#' @export
metropolis <- function(target, proposal_sampler, init = NULL, n_iter, burn_in = 0, thinning = 1) {

  # Input validation
  if (n_iter <= burn_in) stop("n_iter must be greater than burn_in")

  # Get dimension from init or proposal
  if (!is.null(init)) {
    d <- length(init)
  } else {
    test_sample <- proposal_sampler(1, current = NULL)
    d <- length(test_sample)
    init <- numeric(d)
  }

  # Initialize storage
  samples <- matrix(NA, nrow = n_iter, ncol = d)
  samples[1, ] <- init
  accepted <- 0

  # Main MCMC loop
  for (i in 2:n_iter) {
    current <- samples[i-1, ]
    proposal <- current + proposal_sampler(1, current)  # Symmetric proposal

    # Calculate acceptance probability
    target_current <- target(current)
    target_proposal <- target(proposal)

    # Handle possible invalid densities
    if (target_current <= 0 || is.na(target_current)) {
      warning("Invalid current target density at iteration ", i)
      samples[i, ] <- current
      next
    }

    prob_accept <- min(1, target_proposal / target_current)

    # Accept/reject
    if (runif(1) <= prob_accept) {
      samples[i, ] <- proposal
      accepted <- accepted + 1
    } else {
      samples[i, ] <- current
    }
  }

  # Apply burn-in and thinning
  if (burn_in > 0) samples <- samples[-(1:burn_in), ]
  if (thinning > 1) {
    keep <- seq(1, length(samples), by = thinning)
    samples <- samples[keep ]
  }

  # Calculate and print acceptance rate
  acceptance_rate <- accepted / (n_iter - 1)
  message("Acceptance rate: ", round(acceptance_rate * 100, 1), "%")

  return(samples)
}

#------------------------------------
#' Gibbs Sampling MCMC Algorithm
#'
#' Implements the Gibbs sampling algorithm, a Markov Chain Monte Carlo (MCMC) method
#' for sampling from multivariate probability distributions when direct sampling is difficult.
#'
#' @param conditional_samplers A list of functions where each function takes the current
#'        state vector and returns a sample from its conditional distribution
#' @param n Number of samples to return (after burn-in and thinning)
#' @param init Optional initial values (defaults to vector of 0s)
#' @param burn_in Number of initial samples to discard (default: 0)
#' @param thinning Keep only every k-th sample (default: 1, no thinning)
#'
#' @return A matrix where each row represents one sample from the multivariate distribution,
#'         with columns corresponding to the variables in the order of conditional_samplers
#'
#' @examples
#' # Example 1: Bivariate Normal Distribution
#' samplers <- list(
#'   function(state) rnorm(1, mean = 0.5 * state[2], sd = sqrt(1 - 0.5^2)),  # X|Y ~ N(0.5Y, 0.75)
#'   function(state) rnorm(1, mean = 0.5 * state[1], sd = sqrt(1 - 0.5^2))   # Y|X ~ N(0.5X, 0.75)
#' )
#' samples <- gibbs(samplers, n = 1000, init = c(0, 0), burn_in = 200)
#' plot(samples, main = "Bivariate Normal Samples")
#'
#' # Example 2: Using Rejection Sampling for One Conditional
#' target_y <- function(y, x) {
#'   ifelse(y > 0 & y < 10, exp(-abs(x) * y), 0)
#' }
#' conditional_y <- function(state) {
#'   x <- state[1]
#'   rejection_sampler(
#'     target = function(y) target_y(y, x),
#'     proposal_sampler = function(n) runif(n, 0, 10),
#'     proposal_density = function(y) dunif(y, 0, 10),
#'     n = 1,
#'     lower = 0,
#'     upper = 10
#'   )$sampled
#' }
#' samplers <- list(
#'   function(state) rnorm(1, mean = -state[2], sd = sqrt(1/2)),  # X|Y ~ N(-Y, 1/2)
#'   conditional_y                                               # Y|X via rejection
#' )
#' samples <- gibbs(samplers, n = 500, burn_in = 100, thinning = 2)
#'
#' # Example 3: Three-Variable System
#' samplers <- list(
#'   function(state) rnorm(1, mean = (state[2] + state[3])/2),  # X|Y,Z
#'   function(state) rexp(1, rate = abs(state[1]) + 1),        # Y|X,Z
#'   function(state) rpois(1, lambda = 2 + abs(state[1]))      # Z|X,Y
#' )
#' samples <- gibbs(samplers, n = 1000)
#' pairs(samples, labels = c("X", "Y", "Z"))
#'
#' @export
gibbs<- function(conditional_samplers, n , init = NULL, burn_in = 0, thinning = 1){

  if(!is.list(conditional_samplers)){
    stop("conditional samplers must be a list of functions")
  }

  m <- length(conditional_samplers)

  if (is.null(init)){
    init<- numeric(m)
  }

  if (length(init) != m){
    stop("length of initial values must be the same as number of conditional samplers")
  }

  n = n*thinning +burn_in
  sampled = matrix(0, n, m)
  sampled[1,] = init

  # we start from second row sonce we have initial value in first row
  for (i in 2:n) {
    # Create a copy of the previous sample
    current <- sampled[i-1, ]
    # Update each component using its conditional sampler
    for (j in 1:m) {
      # The sampler function should take the current state as input
      current[j] <- conditional_samplers[[j]](current)
    }

    sampled[i, ] <- current
  }


  #apply burn_in
  if (burn_in > 1){
    sampled<- tail(sampled, (nrow(sampled)-burn_in))
  }
  # apply thinning(it is upto rows of sampled after burn in)
  if (thinning > 1){
    index = seq(from = 1, to = nrow(sampled), by = thinning)
    sampled = sampled[index,]
  }


  return(sampled)

}

#------------------------------------

#' Visualize 2D Target PDF and Sampled Data
#'
#' Creates a contour plot of a 2D target probability density function (PDF)
#' overlaid with a scatter plot of sampled data points.
#'
#' @param target_pdf A function that takes `x` and `y` as inputs and returns the PDF value at `(x,y)`.
#' @param sampled_data A `data.frame` or `matrix` with 2 columns (x, y) of sampled points.
#' @param xlim Numeric vector of length 2 (`c(min, max)`) for the x-axis range.
#' @param ylim Numeric vector of length 2 (`c(min, max)`) for the y-axis range.
#' @param n_grid Number of grid points along each axis for evaluating the PDF (default: `100`).
#' @param contour_bins Number of contour levels (default: `15`).
#' @param point_alpha Transparency of sampled points (0 = fully transparent, 1 = opaque, default: `0.5`).
#' @param point_size Size of sampled points (default: `1`).
#' @param fill_contours If `TRUE`, fills contours with a gradient (default: `FALSE`).
#' @param show_points If `FALSE`, hides the sampled points (default: `TRUE`).
#'
#' @return A `ggplot2` plot object.
#'
#' @examples
#' # Example 1: Bivariate normal
#' target_normal <- function(x, y) {
#'   exp(-0.5 * (x^2 + y^2)) / (2 * pi)
#' }
#' sampled_data <- MASS::mvrnorm(1000, mu = c(0, 0), Sigma = diag(2))
#' plot_2d_pdf(target_normal, sampled_data, xlim = c(-3, 3), ylim = c(-3, 3))
#'
#' # Example 2: Custom PDF (like yours)
#' target_custom <- function(x, y) {
#'   ifelse(x > -1 & y > 0, y * exp(-y * (x + 1)), 0)
#' }
#' sampled_data <- data.frame(x = rexp(1000) - 1, y = rexp(1000))  # Not true samples!
#' plot_2d_pdf(target_custom, sampled_data, xlim = c(-1, 5), ylim = c(0, 5))
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
plot_2d_pdf <- function(
    target_pdf,
    sampled_data,
    xlim,
    ylim,
    n_grid = 100,
    contour_bins = 15,
    point_alpha = 0.5,
    point_size = 1,
    fill_contours = TRUE,
    show_points = TRUE
) {
  # Create grid
  x_seq <- seq(xlim[1], xlim[2], length.out = n_grid)
  y_seq <- seq(ylim[1], ylim[2], length.out = n_grid)
  grid <- expand.grid(x = x_seq, y = y_seq)

  # Evaluate PDF on grid (vectorized)
  grid$density <- target_pdf(grid$x, grid$y)

  # Remove NA/Inf values
  grid <- grid[is.finite(grid$density) & grid$density > 0, ]

  # Convert sampled data to data.frame
  if (!is.data.frame(sampled_data)) {
    sampled_data <- as.data.frame(sampled_data)
  }
  colnames(sampled_data) <- c("x", "y")

  # Base plot
  gg <- ggplot2::ggplot() +
    ggplot2::labs(
      title = "2D Target PDF with Sampled Data",
      x = "x",
      y = "y"
    ) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
    ggplot2::theme_minimal()

  # Add filled contours if requested
  if (fill_contours) {
    gg <- gg +
      ggplot2::geom_raster(
        data = grid,
        aes(x = x, y = y, fill = density),
        alpha = 0.7
      ) +
      ggplot2::scale_fill_viridis_c(name = "Density")
  }

  # Add contour lines
  gg <- gg +
    ggplot2::geom_contour(
      data = grid,
      aes(x = x, y = y, z = density),
      bins = contour_bins,
      color = "darkblue",
      linewidth = 0.5
    )

  # Add sampled points if requested
  if (show_points) {
    gg <- gg +
      ggplot2::geom_point(
        data = sampled_data,
        aes(x = x, y = y),
        alpha = point_alpha,
        size = point_size,
        color = "red"
      )
  }

  return(gg)
}
