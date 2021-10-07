source("../find_constants_lyapunov.R")
d <- 0.3
k <- 20
np <- 10^6
y1s <- (1-d)/2
y2s <- (1+d)/2
xs <- (1 - d^2)

y1 <- runif(np) * k * y1s
y2 <- runif(np) * k * y2s
x <- runif(np) * k * xs

# polynomial system
Y <- rbind(x, y1 + y2, ((1 + d) * y1 + (1-d) * y2) / (y1 + y2))
ystar <- c(xs, y1s + y2s, ((1 + d) * y1s + (1-d) * y2s) / (y1s + y2s))

lambda <- c(1,0,-2)

B <- matrix(c(
  1, 0, 0,
  0, 1, 0,
  0, 0, 1,
  0, 0, -1,
  1, 0, -1
), 5, 3, byrow = TRUE)

M <- matrix(c(
  0, -1, 0, 0, 0,
  1, 0, -1, 0, 0,
  -1, 0, 1, 1 - d^2, 1
),3, 5, byrow = TRUE)

find_constants_GLV_embedding(
  lambda = lambda, 
  M = M, 
  B = B, 
  ystar = ystar,
  Y = Y, 
  THRESH = 10^(-16), 
  EQUIL = 10^(-14), 
  n_points = 10000, # number of points to consider when searching numerically for Lyapunov function
  k = 10, # take points in the n-dimensional box defined by the edges (0, k y_i*)
  n_rounds = 10, # number of iterations for numerical search
  num_steps = 50000, # number of iterations for each optim algorithm
  show_progress = TRUE, # show trace of numerical search
  function_constants_name = "default_constants", # choose from above
  weight_pattern = c(1,1,1,0,0) # weights for the constants (use to set some to zero)
)
