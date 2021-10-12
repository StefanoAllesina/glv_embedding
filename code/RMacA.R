# goal: find stable parameterization for Rosenzweig-MacArthur model with two competing resources and one predator
rm(list = ls())
library(deSolve)
THRESH <- 10^-6

# Eqs:
rmca <- function(t, x, pars = NULL){
  x[x < THRESH] <- 0
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  dx1 <- x1 * (r1 - A[1,1] * x1 - A[1,2] * x2) - b1 * x1 * x3 / (h + x1 + x2)
  dx2 <- x2 * (r2 - A[2,1] * x1 - A[2,2] * x2) - b2 * x2 * x3 / (h + x1 + x2)
  dx3 <- -d * x3 + e * (b1 * x1 + b2 * x2) * x3 / (h + x1 + x2)
  return(list(c(dx1, dx2, dx3)))
}

possiblevals <- (1:20)

success <- FALSE
my_seed <- 0
while(!success){
  my_seed <- my_seed + 1
  set.seed(my_seed)
  print(my_seed)
  r1 <- sample(possiblevals, 1)
  r2 <- sample(possiblevals, 1)
  b1 <- sample(possiblevals, 1)
  b2 <- sample(possiblevals, 1)
  d <- sample(possiblevals, 1)
  h <- sample(possiblevals, 1)
  e <- sample(possiblevals, 1)
  A <- matrix(sample(possiblevals, 4, replace = TRUE), 2, 2)
  x0 <- runif(3)
  tmp <- ode(y = x0, times = seq(0, 100, by = 1), parms = NULL, method = "ode45", func = rmca)
  finaldens <- tmp[nrow(tmp), -1]
  if (all(finaldens > 0.05)){
    if (sum(unlist(rmca(0, finaldens, NULL))^2) < 10^-10){
      plot(tmp)
      success <- TRUE
    }
  }
}
