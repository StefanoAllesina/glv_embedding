# # QP formulation
lambda <- c(60, 72, -52)
# 
M <- matrix(c(
  -65, -5, -3, -20, -5, -25, 
  14, 0, -8, -4, -18, -22, 
  14, 59, 0, 0, 0, 0
), 3, 6, byrow = TRUE)

B <- matrix(c(
  1, 0, 0,
  0, 1, 0,
  0, 0, 1,
  2, 0, 0,
  0, 2, 0,
  1, 1, 0
), 6, 3, byrow = TRUE)
# 
source("../find_constants_lyapunov.R")
#example finding equilibrium numerically

# BUG: for some reasons this take a long/infinite time

# result <- find_constants_GLV_embedding(lambda = lambda,
#                                        M = M,
#                                        B = B,
#                                        function_constants_name = "minimize_minimum_constants")
# example providing equilibrium
result <- find_constants_GLV_embedding(lambda = lambda,
                                       M = M,
                                       B = B,
                                       ystar = c(2311/4464, 1693/2232, 1661795/1107072),
                                       #function_constants_name = "maximize_maximum_constants",
                                       function_constants_name = "maximize_variance_constants",
                                       show_progress = TRUE,
                                       weight_pattern = c(1,1,1,0,0,0), n_rounds = 10)
