# Input:
# QP dynamical system
# \dot{y_i} = y_i (\lambda_i + \sum_{j = 1}^{m} M_{ij} \Prod_{k = 1}^{n} y_k^{B_{jk}})
# defined by 
# lambda := growth/death rates; length n
# M := interactions; size n x m 
# B := exponents; size m x n
# k := distance from equilibrium (i.e., sample y_i uniformly between 0 and k y_i*)
# num_points := number of points to sample

source("example.R")

# dimensions
n <- length(lambda)
m <- ncol(M)

# if we are not provided with an equilibrium, find it numerically
# (assuming there's only one and is globally stable)
if(is.null(ystar)){
  print("Finding equilibrium by integrating ODE")
  library(deSolve)
  THRESH <- 10^(-16)
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
    if (all(abs(QP(0, ystar, pr)[[1]]) < THRESH * 100)) success <- TRUE
  }
  print("done")
}

# Transformed parameters
A <- B %*% M
r <- as.vector(B %*% lambda)
if (is.null(Y)){
  # sample y
  Y <- diag(ystar) %*% matrix(runif(n * np, 0, k), n, np)
}
# build a matrix of equil values
Ystar <- matrix(ystar, n, np)
# Transform variables
Z <- exp(B %*% log(Y))
Zstar <- exp(B %*% log(Ystar))
# Deviation from equilibrium
DZ <- Z - Zstar

# Now search for constants for candidate Lyapunov function
find_constants <- function(alpha, DZ, A){
  alpha <- abs(alpha) # nonnegative constants
  alpha <- alpha / mean(alpha) # set mean to 1 to avoid a = rep(0, m)
  G <- diag(alpha) %*% A
  G <- (G + t(G)) / 2
  # compute DZ^t G DZ
  # 
  dVdt <- colSums(DZ * (G %*% DZ))
  if (any(dVdt > 0)){
    return(sum(dVdt[dVdt > 0]))
  }
  # use function instead!
  return(-1)
}

nrounds <- 3
tmp <- list(par = runif(m))
for (i in nrounds){
  tmp <- optim(par = tmp$par, fn = find_constants, method = "Nelder-Mead", 
               A = A, DZ = DZ, control = list(maxit = 5000, trace = TRUE))
  tmp <- optim(par = tmp$par, fn = find_constants, method = "BFGS", 
               A = A, DZ = DZ, control = list(maxit = 5000, trace = TRUE))
}
if (tmp$value < 0){
  print("Success!")
  alpha <- tmp$par
  alpha <- abs(alpha) 
  alpha <- alpha / mean(alpha)
  print(alpha)
}