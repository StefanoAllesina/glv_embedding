source("../find_constants_lyapunov.R")
r1 <- runif(1)
r2 <- runif(1)
a1 <- runif(1)
a2 <- runif(1)
b1 <- runif(1)
b2 <- runif(1)
g1 <- runif(1)
g2 <- runif(1)

k <- 5
np <- 10^5

z1s <- (1/(2  * a1  * (b2+a2 * g2+r2))) * 
  (b1 * b2+b2 * r1+a2 * g2 * r1+b1 * r2+r1 * r2-a1 * g1 * (a2 * g2+r2)+
     sqrt(4 * a1 *  g1 * (b2+a2 * g2+r2) * (a2 * g2 * r1+(b1+r1) * r2)+
            (b1 * (b2+r2)-a1 * g1 * (a2 * g2+r2)+r1 * (b2+a2 * g2+r2))^2))
z2s <- (1/(2 * a2 * (b1+a1 * g1+r1))) * 
  (-a1 * a2 * g1 * g2+b2 * r1-a2 * g2 * r1+a1 * g1* r2+r1* r2+b1* (b2+r2)+
     sqrt(4 * a1* g1 * (b2+a2*  g2+r2)* (a2* g2* r1+(b1+r1)* r2)+(b1 *(b2+r2)-a1 * g1* (a2* g2+r2)+r1* (b2+a2* g2+r2))^2))

z1 <- runif(np) * k * z1s
z2 <- runif(np) * k * z2s
x1 <- z1
x2 <- z2
x3 <- 1 / (g1 + x1)
x4 <- 1 / (g2 + x2)
x1s <- z1s
x2s <- z2s
x3s <- 1 / (g1 + x1s)
x4s <- 1 / (g2 + x2s)

# polynomial system
Y <- rbind(x1, 
           x2,
           x1 * x3,
           x2 * x4,
           x1^2 * x3,
           x2^2 * x4,
           x1*x2*x3*x4)
Y <- rbind(x1, 
           x2,
           x3,
           x4)

ystar <- c(x1s, x2s, x3s, x4s)

lambda <- c(r1, r2, 0, 0)

B <- matrix(c(
  1, 0, 0, 0, 
  0, 1, 0, 0, 
  1, 0, 1, 0, 
  0, 1, 0, 1, 
  2, 0, 1, 0,
  0, 2, 0, 1, 
  1, 1, 1, 1
), 7, 4, byrow = TRUE)

M <- matrix(c(
  -a1, 0, 0, b1, 0, 0, 0,
  0, -a2, b2, 0, 0, 0, 0,
  0, 0, -r1, 0, a1, 0, -b1,
  0, 0, 0, -r2, 0, a2, -b2
),4, 7, byrow = TRUE)


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
  weight_pattern = c(0,0,0,0,0,1,1) # weights for the constants (use to set some to zero)
)
# 1,2
# 1,4
# 1,6
# 1,7
# 2,3
# 2,5
# 2,7
# 3,4
# 3,6
# 3,7
# 4,5
# 4,7
# 5,6
# 5,7
# 6,7





# A <- tmp$A
# DZ <- tmp$DZ
# alpha <- c(g2 + y2, 1/a2, 0,0,0,0,0)
# weight_pattern = c(1,1,0,0,0,0,0)
# find_constants <- function(alpha, DZ, A){
#   alpha <- abs(alpha) * weight_pattern # nonnegative constants
#   alpha <- alpha / mean(alpha) # set mean to 1 to avoid a = rep(0, m)
#   G <- diag(alpha) %*% A
#   G <- (G + t(G)) / 2
#   dVdt <- colSums(DZ * (G %*% DZ))
#   if (any(dVdt > 0)){
#     return(sum(dVdt[dVdt > 0]))
#   }
#   # use specific function
#   return(function_constants(alpha))
# }
