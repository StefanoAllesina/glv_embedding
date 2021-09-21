## Leslie-Gower model
# \dot{x} = x (1 - x - z)
# \dot{z} = z (delta - beta z / x)
# In QP form, we have:
# y = {x, y}
# lambda = (1, delta)
# M = {{-1, -1, 0}, {0, 0, -beta}}
# B = {{1,0},{0,1},{-1,1}}
# take random parameters
delta <- runif(1)
beta <- runif(1)

lambda <- c(1, delta)
M <- matrix(c(
  -1, -1, 0, 
  0, 0, -beta
), 2, 3, byrow = TRUE)

B <- matrix(c(
  1,0,
  0,1,
  -1,1
), 3, 2, byrow = TRUE)


k <- 5
np <- 1000
ystar <- NULL # find equilibrium numerically
Y <- NULL # sample values uniformly in the rectangular space 