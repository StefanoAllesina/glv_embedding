# Leslie-Gower model:
# Eq. 1 of 

# A. Korobeinikov
# A Lyapunov function for Leslie–Gower predator–prey models
# Appl. Math. Lett., 14 (6) (2001), pp. 697-699

# with x = H and y = P of the original equation

## System of equations
# \dot{x} = x (r_1 - a_1 y)
# \dot{y} = y (r_2 - a_2 \frac{y}{x})

## Equilibrium
# x^\star = \frac{r_1 a_2}{a_1 r_2}
# y^\star = \frac{r_1}{a_1}

# Known Lyapunov function:
# V =  \frac{x^\star}{x} + \log \frac{x}{x^\star} + 
#      \frac{a_1 x^\star}{a_2} \left(\frac{y^\star}{y} + \log \frac{x}{x^\star})

# Take random parameters
r1 <- runif(1)
r2 <- runif(1)
a1 <- runif(1)
a2 <- runif(1)
# QP formulation
lambda <- c(r1, r2)

M <- matrix(c(
  0, -a1, 0, 
  0, 0, -a2
), 2, 3, byrow = TRUE)

B <- matrix(c(
  1,0,
  0,1,
  -1,1
), 3, 2, byrow = TRUE)

source("../find_constants_lyapunov.R")
set.seed(5)
# example finding equilibrium numerically
result <- find_constants_GLV_embedding(lambda = lambda, 
                                       M = M, 
                                       B = B, 
                                       function_constants_name = "minimize_minimum_constants")
# example providing equilibrium 
result <- find_constants_GLV_embedding(lambda = lambda, 
                                       M = M, 
                                       B = B, 
                                       ystar = c(r1 * a1 / (r2 * a2), r1 / a1),
                                       function_constants_name = "minimize_minimum_constants", 
                                       show_progress = FALSE)