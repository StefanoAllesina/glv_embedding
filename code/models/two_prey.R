source("../find_constants_lyapunov.R")
e <- 0.5
k <- 5
np <- 10^5
x1s <- 1/(2 * (1 + e))
x2s <- 1/(2 * (1 - e))
ys <- 1/(1 - e^2)

x1 <- runif(np) * k * x1s
x2 <- runif(np) * k * x2s
y <- runif(np) * k * ys

# polynomial system
Y <- rbind(x1 + x2, y, ((1 + e) * x1 + (1-e) * x2) / (x1 + x2))
ystar <- c(x1s + x2s, ys, ((1 + e) * x1s + (1-e) * x2s) / (x1s + x2s))

lambda <- c(1,-1,-1)

B <- matrix(c(
  0, 1, 0,
  0, 1, 1, 
  1, 1, 0, 
  0, 0, -1, 
  0, 1, -1
), 5, 3, byrow = TRUE)

M <- matrix(c(
  0, -1, 0, 0, 0,
  0, 0, 1, 0, 0,
  -2, 1, 0, 1, 1 - e^2
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
