source("../find_constants_lyapunov.R")

# find params with 3 equilibria (and not 4)
set.seed(1)
success <- FALSE
while(!success){
  r <- runif(1)
  a <- runif(1)
  g <- runif(1)
  d <- runif(1)
  m <- runif(1)
  if ((-4 * d^2 * g + m^2) > 0){
    otherz1s <- (m - sqrt(-4 * d^2 * g + m^2)/(2 * d))
    otherz2s <- (r - a * otherz1s) * (g + otherz1s^2)
    z1s <- (m + sqrt(-4 * d^2 * g + m^2)/(2 * d))
    z2s <- (r - a * z1s) * (g + z1s^2)
    if ((otherz1s < 0) & (z1s > 0) & (z2s > 0) )success <- TRUE
  }
}

k <- 5
np <- 10^5

y1 <- runif(np) * k * z1s
y2 <- runif(np) * k * z2s
y1s <- z1s
y2s <- z2s

# polynomial system
Y <- rbind(y1, 
           y2)

ystar <- c(y1s, y2s)

lambda <- c(g * r, - d * g)

B <- matrix(c(
  1,0,
  0,1,
  2,0,
  3,0
), 4, 2, byrow = TRUE)

M <- matrix(c(
  -a * g, -1, r,-a,
  m, 0, -d, 0
),2, 4, byrow = TRUE)


tmp <- find_constants_GLV_embedding(
  lambda = lambda, 
  M = M, 
  B = B, 
  ystar = ystar,
  Y = Y, 
  THRESH = 10^(-16), 
  EQUIL = 10^(-14), 
  n_points = 10000, # number of points to consider when searching numerically for Lyapunov function
  k = 5, # take points in the n-dimensional box defined by the edges (0, k y_i*)
  n_rounds = 5, # number of iterations for numerical search
  num_steps = 50000, # number of iterations for each optim algorithm
  show_progress = TRUE, # show trace of numerical search
  function_constants_name = "default_constants", # choose from above
  weight_pattern = c(1,1,0,0) # weights for the constants (use to set some to zero)
)

