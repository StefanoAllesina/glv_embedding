# SIR model with input:
# Eq. 4.1 of 
# Shuai, Z. and van den Driessche, P. (2013). 
# Global stability of infectious disease models using Lyapunov functions. 
# SIAM Journal on Applied Mathematics, 73(4), 1513-1532.

# with y1 = S, y2 = I, and y3 = R of the original equation

## System of equations
# \dot{y1} = y1 (-\delta + \mu / y1 - \beta y2)
# \dot{y2} = y2 (-(\delta + \gamma + \alpha) + \beta y1)
# \dot{y3} = y3 (-delta + \gamma y2 / y3)
# 
## Equilibrium
# y1^\star = \frac{\alpha + \gamma + \delta}{\beta}
# y2^\star = \frac{\mu}{\alpha + \gamma + \delta} - \frac{delta}{\beta}
# y3^\star = \frac{\gamma \mu}{(\alpha + \gamma + \delta) \delta} - \frac{gamma}{\beta}
# Note:
# y3^\star / y2^\star = \frac{\gamma}{\delta}

# Take random parameters yielding feasible equilibrium
set.seed(11)

success <- FALSE
while(!success){
  mu <- runif(1)
  beta <- runif(1)
  gamma <- runif(1)
  delta <- runif(1)
  alpha <- runif(1)
  y1s <- (alpha + gamma + delta) / (beta)
  y2s <- (mu)/(alpha + gamma + delta) - (delta)/(beta)
  y3s <- y2s * gamma / delta
  if (all(c(y1s, y2s, y3s) > 0)) success <- TRUE
}

# # QP formulation
lambda <- c(-delta, -(alpha + gamma + delta), -delta)
# 
M <- matrix(c(
   0, -beta, mu, 0,
   beta, 0, 0, 0,
   0, 0, 0, gamma
), 3, 4, byrow = TRUE)
# 
B <- matrix(c(
   1,0,0,
   0,1,0,
   -1,0,0,
   0,1,-1
), 4, 3, byrow = TRUE)
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
                                       ystar = c(y1s, y2s, y3s),
                                       function_constants_name = "minimize_minimum_constants",
                                       #function_constants_name = "maximize_maximum_constants",
                                       show_progress = FALSE)
