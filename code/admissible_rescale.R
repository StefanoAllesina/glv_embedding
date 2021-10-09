check_TRLD <- function(pars, M, B){
  # weights 
  n <- ncol(M)
  w <- abs(pars[1:n])
  w <- w / mean(w)
  # scales 
  ts <- pars[n + 1:ncol(B)]
  #ts <- ts / max(ts)
  # build new M
  B <- B + matrix(ts, n, ncol(B), byrow = TRUE)
  A <- B %*% M
  G <- diag(w) %*% A
  G <- G + t(G)
  eG <- eigen(G, only.values = TRUE, symmetric = TRUE)$values
  if (all(eG <= 10^-16)){
    return(-1)
  }
  return(sum(eG[eG > 0]))
}

search_admissible_via_rescale <- function(s, M, B, nrounds = 5, mxit = 5000){
  M <- cbind(M, s)
  B <- rbind(B, 0)
  n <- ncol(M)
  tmp <- list(par = rnorm(ncol(M) + ncol(B)))
  for (i in 1:nrounds){
    tmp <- optim(par = tmp$par, fn = check_TRLD, method = "Nelder-Mead",
                 control = list(maxit = mxit, trace = TRUE), M = M, B = B)
    tmp <- optim(par = tmp$par, fn = check_TRLD, method = "BFGS",
                 control = list(maxit = mxit, trace = TRUE), M = M, B = B)
  }
  n <- ncol(M)
  pars <- tmp$par
  w <- abs(pars[1:n])
  w <- w / mean(w)
  # scales 
  ts <- pars[n + 1:ncol(B)]
  #ts <- ts / max(ts)
  # build new M
  B <- B + matrix(ts, n, ncol(B), byrow = TRUE)
  A <- B %*% M
  G <- diag(w) %*% A
  G <- G + t(G)
  eG <- eigen(G, only.values = TRUE, symmetric = TRUE)$values
  print(eG)
  if (tmp$value < 0){
    print("Success")
  } else {
    print("Failed")
  }
  return(list(ts = ts, w = w))
}

# Two prey
e <- 0.5
# polynomial system
s <- c(1,-1,-1)

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

# Two pred
d <- 0.3

s <- c(1,0,-2)

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

# find params with 3 equilibria (and not 4)
#set.seed(1)
success <- FALSE
while(!success){
  r <- runif(1)
  a <- runif(1)
  g <- runif(1)
  d <- runif(1)
  m <- runif(1)
  if ((-4 * d^2 * g + m^2) > 0){
    otherz1s <- (m - sqrt(-4 * d^2 * g + m^2)/(2 * d))
    otherz2s <- (r - a * otherz1s) * (g + otherz1s^2)
    z1s <- (m + sqrt(-4 * d^2 * g + m^2)/(2 * d))
    z2s <- (r - a * z1s) * (g + z1s^2)
    if ((otherz1s < 0) & (z1s > 0) & (z2s > 0) )success <- TRUE
  }
}

s <- c(g * r, - d * g)

B <- matrix(c(
  1,0,
  0,1,
  2,0,
  3,0
), 4, 2, byrow = TRUE)

M <- matrix(c(
  -a * g, -1, r,-a,
  m, 0, -d, 0
),2, 4, byrow = TRUE)



print(search_admissible_via_rescale(s, M, B))