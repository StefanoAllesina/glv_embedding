# # QP formulation
lambda <- c(0, 0, 0)
# 
M <- matrix(c(
  0, -1, 2, -2, 1, 0,
  2, 0, -1, 1, 0, -2,
  -1, 2, 0, 0, -2, 1
), 3, 6, byrow = TRUE)
# 
B <- matrix(c(
  2, 0, 0,
  0, 2, 0,
  0, 0, 2,
  1, 1, 0,
  1, 0, 1,
  0, 1, 1
), 6, 3, byrow = TRUE)
# 
source("../../find_constants_lyapunov.R")
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
                                       ystar = rep(1/3,3),
                                       function_constants_name = "minimize_minimum_constants",
                                       #function_constants_name = "maximize_maximum_constants",
                                       show_progress = FALSE)
DZ <- result$DZ
A <- result$A
G <- diag(c(1,1,1,2,2,2))
M <- G %*% A
M <- M + t(M)
dVdt <- sapply(1:ncol(DZ), function(i) DZ[,i] %*% M %*% DZ[,i])
