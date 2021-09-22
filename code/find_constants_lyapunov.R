# functions to find constants maximizing/minimizing certain functions
default_constants <- function(alpha) return(-1)
minimize_variance_constants <- function(alpha) return(-1 / var(alpha))
maximize_variance_constants <- function(alpha) return(-var(alpha))
maximize_maximum_constants <- function(alpha) return(-max(alpha))
minimize_minimum_constants <- function(alpha) return(-1/min(alpha))

# Input:
# QP dynamical system
# \dot{y_i} = y_i (\lambda_i + \sum_{j = 1}^{m} M_{ij} \Prod_{k = 1}^{n} y_k^{B_{jk}})
# defined by 
# lambda := growth/death rates; length n
# M := interactions; size n x m 
# B := exponents; size m x n

find_constants_GLV_embedding <- function(
  lambda, # growth rates of original system [required]
  M, # matrix of interactions [required]
  B, # matrix of exponents [required]
  ystar = NULL, # stable equilibrium if NULL, find equilibrium numerically
  Y = NULL, # points to test if NULL, sample values uniformly in a n-dimensional box
  THRESH = 10^(-16), # consider extinct below this threshold when integrating dynamics
  EQUIL = 10^(-14), # consider dy/dt|y* = 0 if abs(dy/dt|y*) < EQUIL when searching numerically for equilibrium
  n_points = 100000, # number of points to consider when searching numerically for Lyapunov function
  k = 10, # take points in the n-dimensional box defined by the edges (0, k y_i*)
  n_rounds = 3, # number of iterations for numerical search
  num_steps = 5000, # number of iterations for each optim algorithm
  show_progress = TRUE, # show trace of numerical search
  function_constants_name = "default_constants" # choose from above
){
  function_constants <- match.fun(function_constants_name)
  
  # dimensions
  n <- length(lambda)
  m <- ncol(M)

  # if we are not provided with an equilibrium, find it numerically
  # (assuming there's only one and is globally stable)
  if(is.null(ystar)){
    cat("Finding equilibrium by integrating ODE... ")
    library(deSolve)
    QP <- function(t, y, pars){
      y[y < THRESH] <- 0
      M <- pars$M
      B <- pars$B
      lambda <- pars$lambda
      # compute z_i = \Prod_{k = 1}^{n} y_k^{B_{jk}}
      z <- exp(B %*% log(y))
      dydt <- y * (lambda + as.vector(M %*% z))
      return(list(dydt))
    }
    ystar <- runif(n)
    pr <- list(M = M, B = B, lambda = lambda)
    # keep integrating ODE until convergence
    success <- FALSE
    while(!success){
      out <- ode(y = ystar, times = seq(0, 100, by = 0.01), func = QP, 
               parms = pr, 
               method = "ode45")
      ystar <- out[nrow(out), -1]
      if (all(abs(QP(0, ystar, pr)[[1]]) < EQUIL)) success <- TRUE
    }
    cat("done\n")
  }

  # Transformed variables and parameters
  A <- B %*% M
  r <- as.vector(B %*% lambda)
  # if points are not provided, sample in a box
  if (is.null(Y)){
    # sample y
    Y <- diag(ystar) %*% matrix(runif(n * n_points, 0, k), n, n_points)
  }
  # build a matrix of equil values
  Ystar <- matrix(ystar, n, n_points)
  # Transform variables
  Z <- exp(B %*% log(Y))
  Zstar <- exp(B %*% log(Ystar))
  # Deviations from equilibrium
  DZ <- Z - Zstar

  # Now search for constants for candidate Lyapunov function
  find_constants <- function(alpha, DZ, A){
    alpha <- abs(alpha) # nonnegative constants
    alpha <- alpha / mean(alpha) # set mean to 1 to avoid a = rep(0, m)
    G <- diag(alpha) %*% A
    G <- (G + t(G)) / 2
    dVdt <- colSums(DZ * (G %*% DZ))
    if (any(dVdt > 0)){
      return(sum(dVdt[dVdt > 0]))
    }
    # use specific function
    return(function_constants(alpha))
  }

  cat("Searching numerically for constants of Lyapunov function... ")
  tmp <- list(par = runif(m))
  for (i in n_rounds){
    tmp <- optim(par = tmp$par, fn = find_constants, method = "Nelder-Mead", 
                 A = A, DZ = DZ, control = list(maxit = num_steps, trace = show_progress))
    tmp <- optim(par = tmp$par, fn = find_constants, method = "BFGS", 
                 A = A, DZ = DZ, control = list(maxit = num_steps, trace = show_progress))
  }
  cat("done\n")
  
  if (tmp$value < 0){
    cat("Success!\n")
    alpha <- tmp$par
    alpha <- abs(alpha) 
    alpha <- alpha / mean(alpha)
    print(alpha)
    return(list(lambda = lambda, M = M, B = B, 
                r = r, A = A, DZ = DZ,
                alpha = alpha))
  } else {
    cat("Search failed!\n")
    return(list(lambda = lambda, M = M, B = B, 
                r = r, A = A, DZ = DZ))
  }
}
