# Examples for Chapter 2
rm(list=ls())

# Example 2.1

x <- matrix(c(1,3,2),ncol=1)
y <- matrix(c(-2,1,-1),ncol=1)

vect <- c(1,3,2)
t(vect)
t(t(vect))

as.matrix(vect)

3*x
x+y

( Lx <- sqrt( t(x) %*% x ) )
(Ly <- sqrt(t(y)%*%y))

(cosTheta <- (t(x)%*%y)/(Lx*Ly))
(acos(cosTheta))

r2d <- function(r){
  return(r*180/pi)
  }

r2d(acos(cosTheta))

L3x <- sqrt( t(3*x) %*% (3*x) )
L3x == 3*Lx

# Example 2.3
rm(list=ls())
A <- matrix(c(3,-1,2,
              1,5,4),byrow=TRUE,ncol=3)
t(A)

# Example 2.4
rm(list=ls())
A <- matrix(c(0,3,1,
              1,-1,1),byrow=TRUE,ncol=3)
B <- matrix(c(1,-2,-3,
              2,5,1),byrow=TRUE,ncol=3)
4 * A

A+B

# Example 2.5
rm(list=ls())
A <- matrix(c(3,-1,2,
              1,5,4),byrow=TRUE,ncol=3)
B <- matrix(c(-2,7,9),ncol=1)
C <- matrix(c(2,0, # note about name
              1,-1),byrow=TRUE,ncol=2)

A %*% B
C %*% A
A %*% C # error

# Example 2.7
rm(list=ls())
A <- matrix(c(3,5,5,-2),byrow=TRUE,ncol=2)

A==t(A)
isSymmetric(A)

# Example 2.8
rm(list=ls())
A <- matrix(c(3,2,4,1),byrow=TRUE,ncol=2)

(Ainv <- solve(A))

library(MASS)
fractions(Ainv)

# Example 2.10
rm(list=ls())
A <- matrix(c(13,-4,2,
              -4,13,-2,
              2,-2,10),byrow=TRUE,ncol=3)

eigen(A) # Note different order of eigenvalues and different signs for eigen

lambda <- eigen(A)$values # Because two are the same, eigenvectors for these are not unique
e <- eigen(A)$vectors 

lambda[1]*e[,1]%*%t(e[,1])+
  lambda[2]*e[,2]%*%t(e[,2])+
  lambda[3]*e[,3]%*%t(e[,3])
A == lambda[1]*e[,1]%*%t(e[,1])+lambda[2]*e[,2]%*%t(e[,2])+lambda[3]*e[,3]%*%t(e[,3])
A == round(lambda[1]*e[,1]%*%t(e[,1])+lambda[2]*e[,2]%*%t(e[,2])+lambda[3]*e[,3]%*%t(e[,3]))

# Example 2.12
rm(list=ls())
x1 <- c(-1,0,1)
px1 <- c(.3,.3,.4)
x2 <- c(0,1)
px2 <- c(.8,.2)

(Ex1 <- sum(x1*px1))
(Ex2 <- sum(x2*px2))
(EX <- matrix(c(Ex1,Ex2),ncol=1))

# Example 2.13
px12 <- matrix(c(.24,.06,
                 .16,.14,
                 .40,.00),ncol=2,byrow=TRUE)

(s11 <- sum((x1-Ex1)^2 * px1))
(s22 <- sum((x2-Ex2)^2 * px2))

s12 <- 0
for (i in 1:3){
  for (j in 1:2){
    s12 <- s12 + (x1[i]-Ex1)*(x2[j]-Ex2)*px12[i,j]
  }
}
(s21 <- s12) 

# Could now make matrices if we wanted

# Example 2.14
rm(list=ls())
S <- matrix(c(4,1,2,
              1,9,-3,
              2,-3,25),byrow=TRUE,ncol=3)
(V12 <- diag(sqrt(diag(S))))
library(MASS)
(V12.inv <- fractions(diag(sqrt(1/diag(S)))))

(rho <- V12.inv %*% S %*% V12.inv)
fractions(rho)
