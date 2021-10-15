# 1) find admissible A
set.seed(7)
is_admissible <- function(par, A){
  par <- abs(par)
  par <- par / mean(par)
  par <- round(par, 2)
  B <- diag(par) %*% A
  B <- B + t(B)
  eB <- eigen(B, symmetric = TRUE, only.values = TRUE)$values
  if (max(eB) > 0) return(max(eB))
  return(-1)
}

success <- FALSE
n <- 2
while(!success){
  A <- matrix(round(sample(-10:10,n,n)), n, n)
  diag(A) <- diag(A) - 2
  tmp <- optim(par = rnorm(n), fn = is_admissible, method = "Nelder-Mead", A = A)
  found_admissible <- FALSE
  if (tmp$value < 0){
    w <- abs(tmp$par)
    w <- w / mean(w)
    w <- round(w, 2)
    found_admissible <- TRUE
    print("found admissible")
  }
  if (found_admissible){
    # now choose r that make the equilibrium feasible
    xstar <- sample(1:10, n) / 10
    r <- - as.vector(A %*% xstar)
    # now cycle through time rescales and check whether 
    # the resulting matrix is not admissible
    B <- diag(rep(1,n))
    B <- rbind(B, 0)
    M <- cbind(A, r)
    for (i in 1:n){
      print(i)
      Bprime <- B - matrix(B[i,], nrow(B), ncol(B),byrow = TRUE)
      Aprime <- Bprime %*% M
      newM <- Aprime[-i,]
      newr <- newM[,i]
      newA <- newM[,-i]
      if (abs(det(newA)) > 10^-7 ){
        print("good det")
        tmp <- optim(par = rnorm(n), fn = is_admissible, method = "Nelder-Mead", A = newA)
        if (tmp$value > 0){
          if (all(diag(newA) <= 0)){
            print("success!")
            success <- TRUE
            break
          }
        }
      }
    }
    if (!success)  print("no luck")
  }
}

print(newA)

