check_LD <- function(pars, A, twoeigen = FALSE){
  pars <- abs(pars)
  pars <- pars / mean(pars)
  G <- diag(pars) %*% A
  G <- G + t(G)
  if (twoeigen){
    max(Re(eigen(G, only.values = TRUE, symmetric = TRUE)$values[-nrow(A)]))
  }
  return(max(Re(eigen(G, only.values = TRUE, symmetric = TRUE)$values)))
}

check_TRLD <- function(pars, A, twoeigen = FALSE){
  # weights 
  n <- ncol(A)
  w <- abs(pars[1:n])
  w <- w / mean(w)
  # scales 
  ts <- pars[n + 1:(n-1)]
  # build new A
  B <- diag(rep(1,n-1))
  B <- rbind(B, 0)
  B <- B + matrix(ts, n, n-1, byrow = TRUE)
  M <- B %*% A
  G <- diag(w) %*% M
  G <- G + t(G)
  eG <- eigen(G, only.values = TRUE, symmetric = TRUE)$values
  if (all(eG <= 10^-16)){
    return(-1)
  }
  return(sum(eG[eG > 0]))
}


A <- matrix(c(
  -1/2, 0, -1/3,
  -1/4, -1/3, 0,
  0, 0, -1/4
), 3, 3, byrow = TRUE)

global_success <- FALSE
i <- 188
while(!global_success){
  i <- i + 1
  set.seed(i)
  
  # step 1: choose A23, A32 such that
  # A is stable
  # A is not Lyapunov-Diagonally Stable
  success <- FALSE
  print("Find matrix")
  while(!success){
    A[2,3] <- round(rnorm(1), 2)
    A[3,2] <- round(rnorm(1), 2)
    maxReL1 <- max(Re(eigen(A, only.values = TRUE, symmetric = FALSE)$values))
    if (maxReL1 < 0){
      tmp <- optim(par = rep(1, 3), fn = check_LD, method = "Nelder-Mead", control = list(maxit = 10000), A = A)
      if (tmp$value > 0) success <- TRUE
    }
  }
  
  # step 2: choose r such that the equilibrium is feasible
  success <- FALSE
  print("Find equil")
  while(!success){
    r <- sample(1:10,3)/5
    xstar <- solve(A, -r)
    if (all(xstar > 0.05)) success <- TRUE
  }
  
  print("Find TR")
  # step 3: a time-reparametrization makes it possible to prove global stability
  Atilde <- cbind(A, r)
  Atilde <- rbind(Atilde, 0)
  # new equil is an eigenvector
  xstilde <- Re(eigen(Atilde)$vectors[,4])
  xstilde <- xstilde / xstilde[4]
  # the last eigenvalue of Mtilde will always be zero
  tmp <- optim(par = rnorm(7), fn = check_TRLD, method = "Nelder-Mead", control = list(maxit = 10000, trace = FALSE), A = Atilde[1:3,])
  tmp <- optim(par = tmp$par, fn = check_TRLD, method = "BFGS", control = list(maxit = 10000, trace = FALSE), A = Atilde[1:3,])
  tmp <- optim(par = tmp$par, fn = check_TRLD, method = "Nelder-Mead", control = list(maxit = 10000, trace = FALSE), A = Atilde[1:3,])
  if (tmp$value < 0) {
    print(paste(i, "success"))
    global_success <- TRUE
  } else {
    print(paste(i, "failure"))
  }
  
}
Att <- Atilde[1:3,]
n <- ncol(Att)
pars <- tmp$par
w <- abs(pars[1:n])
w <- w / mean(w)
# scales 
ts <- pars[n + 1:(n-1)]
# build new A
B <- diag(rep(1,n-1))
B <- rbind(B, 0)
B <- B + matrix(ts, n, n-1, byrow = TRUE)
M <- B %*% Att
G <- diag(w) %*% M
G <- G + t(G)
eG <- eigen(G, only.values = TRUE, symmetric = TRUE)$values
print(eG)
print(A)
print(r)
