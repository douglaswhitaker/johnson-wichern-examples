matStats <- function(X){
  X <- as.matrix(X)
  n <- nrow(X)
  x.bar <- t(matrix(1,ncol=n) %*% X)/n
  D <- X - matrix(1,nrow=n) %*% t(x.bar)
  S <- (n-1)^(-1) * t(D)%*%D
  Sinv <- solve(S)
  return(list(x.bar=x.bar,D=D,S=S,Sinv=Sinv,n=n))
}
matStats2 <- function(X){
  X <- as.matrix(X)
  n <- nrow(X)
  x.bar <- t(matrix(1,ncol=n) %*% X)/n
  D <- X - matrix(1,nrow=n) %*% t(x.bar)
  S <- (n-1)^(-1) * t(D)%*%D
  #Sinv <- solve(S)
  return(list(x.bar=x.bar,D=D,S=S,
              #Sinv=Sinv,
              n=n))
}