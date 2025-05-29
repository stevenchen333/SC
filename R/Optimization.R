# example of function
# f <- function(u) {
#   x <- u[1]; y <- u[2]; z <- u[3]
#   return(x^3 + y^2 + z^2 +12*x*y + 2*z)
# }
gradient_descent<-function(fun,init = c(0,0),n_iter = 10000, tol = 1e-6, lr = 0.05){
  path_hist = vector("list", n_iter)

  x = init
  iter = 1
  path_hist[[iter]] = list(x = x, y = fun(x))

  while  (iter < n_iter){
    grad_x = pracma::grad(fun, x)
    grad_norm = norm(grad_x, type = "2")

    if (grad_norm<tol){

      break
    }
    x = x - lr*grad_x
    iter = iter + 1
    path_hist[[iter]] = list(x = x, y = fun(x))


  }
  path_hist = path_hist[1:iter]
  return(list(
    optimum = x,
    final_objective = fun(x),
    number_iteration = iter,
    path = path_hist,
    converged = grad_norm < tol
  ))
}


newton_method<- function(fun,init = c(0,0),n_iter = 10000, tol = 1e-6){
  path_hist = vector("list", n_iter)

  x = init
  iter = 1
  path_hist[[iter]] = list(x = x, y = fun(x))

  while  (iter < n_iter){
    grad_x = pracma::grad(fun, x)
    grad_norm = norm(grad_x,type = "2")
    hessian_x = pracma::hessian(fun, x)
    if (grad_norm<tol){

      break
    }
    if (rcond(hessian_x) < 1e-15) {
      hessian_x = hessian_x + diag(1e-6, nrow(hessian_x))
    }

    delta = solve(hessian_x, grad_x)
    x = x - delta
    iter = iter + 1
    path_hist[[iter]] = list(x = x, y = fun(x))


  }
  path_hist = path_hist[1:iter]
  return(list(
    optimum = x,
    final_objective = fun(x),
    number_iteration = iter,
    path = path_hist,
    converged = grad_norm < tol
  ))
}

# function example for bfgs:
# f <- function(x) sum((x - c(1,2))^2)
quasi_newton_bfgs <- function(fun, init = c(0,0), n_iter = 10000, tol = 1e-6, lr = NULL) {
  # Environment to store path history during optimization
  history_env <- new.env()
  history_env$path <- list()
  history_env$iter <- 0

  # Wrapped function to record history
  wrapped_fun <- function(x) {
    history_env$iter <- history_env$iter + 1
    val <- fun(x)
    history_env$path[[history_env$iter]] <- list(x = x, y = val)
    val
  }

  # Run optim with BFGS
  res <- optim(
    par = init,
    fn = wrapped_fun,
    method = "BFGS",
    control = list(maxit = n_iter, reltol = tol, trace = 0)
  )

  # Extract path history (may be longer than actual iterations due to internal calls)
  path_hist <- history_env$path

  return(list(
    optimum = res$par,
    final_objective = res$value,
    number_iteration = res$counts["function"],
    path = path_hist,
    converged = res$convergence == 0
  ))
}


stoch_grad <- function(fun, df, init = c(0, 0), n_iter = 10000, tol = 1e-6,
                       lr_fun = NULL, type = c("basic", "batch"), num_batch = NULL) {

  # Validate and match SGD type
  type <- match.arg(type)

  # Default constant learning rate
  if (is.null(lr_fun)) {
    lr_fun <- function(t) 0.01
  }

  # Initialize
  path_hist <- vector("list", n_iter)
  x <- init
  iter <- 1
  path_hist[[iter]] <- list(x = x, y = fun(x, df[1, , drop = FALSE]))  # dummy first sample

  converged <- FALSE  # default

  if (type == "basic") {
    while (iter < n_iter) {
      index <- sample(1:nrow(df), 1)
      sample_point <- df[index, , drop = FALSE]

      grad_x <- pracma::grad(fun, x, sample_point)
      grad_norm <- norm(grad_x, type = "2")

      if (grad_norm < tol) {
        converged <- TRUE
        break
      }

      x <- x - lr_fun(iter) * grad_x
      iter <- iter + 1
      path_hist[[iter]] <- list(x = x, y = fun(x, sample_point))
    }
  }

  if (type == "batch") {
    if (is.null(num_batch)) {
      stop("`num_batch` must be specified for batch SGD.")
    }
    num_batch <- round(num_batch)

    while (iter < n_iter) {
      indices <- sample(1:nrow(df), num_batch)
      sample_points <- df[indices, , drop = FALSE]

      grad_sum <- 0
      for (i in 1:num_batch) {
        sample_i <- sample_points[i, , drop = FALSE]
        grad_sum <- grad_sum + pracma::grad(fun, x, sample_i)
      }

      direction <- grad_sum / num_batch
      grad_norm <- norm(direction, type = "2")

      if (grad_norm < tol) {
        converged <- TRUE
        break
      }

      x <- x - lr_fun(iter) * direction
      iter <- iter + 1
      path_hist[[iter]] <- list(x = x, y = fun(x, sample_points[1, , drop = FALSE]))  # just log one sample
    }
  }

  # Trim path history
  path_hist <- path_hist[1:iter]

  return(list(
    optimum = x,
    final_objective = fun(x, df[1, , drop = FALSE]),  # dummy for output
    number_iteration = iter,
    path = path_hist,
    converged = converged
  ))
}

coordinate_descent <- function(fun, init, n_iter = 1000, tol = 1e-6, lower = -Inf, upper = Inf) {
  x <- init
  d <- length(x)
  path_hist <- vector("list", n_iter)
  path_hist[[1]] <- list(x = x, y = fun(x))

  converged <- FALSE
  iter <- 1

  while (iter < n_iter) {
    x_old <- x

    for (j in 1:d) {
      # Define 1D function fixing all other coordinates
      f1d <- function(val) {
        x_temp <- x
        x_temp[j] <- val
        fun(x_temp)
      }

      # Optimize along coordinate j using Brent's method
      res <- optim(par = x[j], fn = f1d, method = "Brent", lower = lower, upper = upper)
      x[j] <- res$par
    }

    # Store path
    iter <- iter + 1
    path_hist[[iter]] <- list(x = x, y = fun(x))

    # Check convergence
    if (norm(x - x_old, type = "2") < tol) {
      converged <- TRUE
      break
    }
  }

  path_hist <- path_hist[1:iter]

  return(list(
    optimum = x,
    final_objective = fun(x),
    number_iteration = iter,
    path = path_hist,
    converged = converged
  ))
}

# TODO: Finish KKT solver
kkt_solver<-function(obj, equality_const, inequality_const, init =c(0,0), n_iter = 10000, tol = 1e-6){

  # initialize
  x = init


  }




# TODO: EM Algorithm




