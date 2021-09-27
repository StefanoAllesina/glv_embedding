# Leslie-Gower model:
# Eq. 2 of 
# A. Korobeinikov
# A Lyapunov function for Leslie–Gower predator–prey models
# Appl. Math. Lett., 14 (6) (2001), pp. 697-699

# with y1 = H and y2 = P of the original equation

## System of equations
# \dot{y1} = y1 (r_1 - y1 - a_1 y2)
# \dot{y2} = y2 (r_2 - a_2 \frac{y2}{y1})

## Equilibrium
# y1^\star = \frac{r_1 a_2}{a_1 r_2 + a2}
# y2^\star = \frac{r_1 r2}{a_1 r2 + a2}

# Known Lyapunov function:
# V =  \frac{y1^\star}{y1} + \log \frac{y1}{y1^\star} + 
#      \frac{a_1 y1^\star}{a_2} \left(\frac{y2^\star}{y2} + \log \frac{y2}{y2^\star})

# Take random parameters
r1 <- runif(1)
r2 <- runif(1)
a1 <- 1#runif(1)
a2 <- 1#runif(1)
# QP formulation
lambda <- c(r1, r2)

M <- matrix(c(
  -1, -a1, 0, 
  0, 0, -a2
), 2, 3, byrow = TRUE)

B <- matrix(c(
  1,0,
  0,1,
  -1,1
), 3, 2, byrow = TRUE)

source("../find_constants_lyapunov.R")
#set.seed(5)
# example finding equilibrium numerically
# result <- find_constants_GLV_embedding(lambda = lambda, 
#                                        M = M, 
#                                        B = B, 
#                                        function_constants_name = "minimize_minimum_constants")
# example providing equilibrium 
result <- find_constants_GLV_embedding(lambda = lambda, 
                                       M = M, 
                                       B = B, 
                                       ystar = c((r1 * a2)/(a1 * r2 + a2),(r1 * r2)/(a1 * r2 + a2)),
                                       function_constants_name = "minimize_minimum_constants", 
                                       #function_constants_name = "maximize_maximum_constants", 
                                       show_progress = FALSE)

DZ <- result$DZ
A <- result$A
test <- function(w1,w2,w3){
  G <- diag(c(w1,w2,w3)) %*% A
  G <- G + t(G)
  dVdt <- sapply(1:ncol(DZ), function(i) DZ[,i] %*% G %*% DZ[,i])
  return(range(dVdt))
}

