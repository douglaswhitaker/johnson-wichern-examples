# Examples for Chapter 3
rm(list=ls())

library(rgl)
library(matlib)
library(MASS)

# Example 3.1
X <- matrix(c(4,-1,3,
              1,3,5),
            ncol=2,byrow=FALSE)

plot(X)
smean <- apply(X, FUN = mean, 2)
points(x=smean[1],y=smean[2],pch=2)

# Example 3.2
open3d()
vectors3d(X[,1],col="blue",lwd=2)
vectors3d(X[,2],col="red",lwd=2)
axes3d()

# Example 3.3
v1 <- matrix(rep(1,nrow(X)),ncol=1)
mean1 <- as.vector(smean[1]*v1)
mean2 <- as.vector(smean[2]*v1)

vectors3d(mean1,col="black",lwd=4)
vectors3d(mean2,col="black",lwd=1)

dev1 <- as.vector((X[,1]-mean1))
dev2 <- as.vector((X[,2]-mean2))
vectors3d(mean1+dev1,origin=mean1,col="blue",lwd=1)
vectors3d(mean2+dev2,origin=mean2,col="red",lwd=1)

# Example 3.4
open3d()
vectors3d(dev1,col="blue",lwd=2)
vectors3d(dev2,col="red",lwd=2)
axes3d()

Ld1 <- sqrt(t(dev1)%*%dev1)
Ld2 <- sqrt(t(dev2)%*%dev2)
cosTheta <- (t(dev1)%*%dev2)/(Ld1*Ld2) 

r2d <- function(r){
  return(r*180/pi)
}

r2d(acos(cosTheta)) # Matches the graph

(L2d1 <- t(dev1)%*%t(t(dev1)))
(L2d2 <- t(dev2)%*%t(t(dev2)))

n <- nrow(X)
(s11 <- L2d1/n)
(s22 <- L2d2/n)
(s12 <- t(dev1)%*%t(t(dev2))/n)
(r12 <- s12/(sqrt(s11)*sqrt(s22))) # cos(theta)

r2d(acos(r12)) # Hey, look at that

Sn <- matrix(c(s11,s12,
               s12,s22),
             ncol=2)
print(fractions(Sn))
R <- matrix(c(1,r12,
              r12,1),
            ncol=2) 
print(R)
  
# Example 3.7
rm(list=ls())

S <- matrix(c(252.04,-68.43,
              -68.43,123.67),byrow = TRUE, ncol=2)

det(S) # Generalized variance

# Example 3.9
X <- matrix(c(1,4,4,
              2,1,0,
              5,6,4),byrow=FALSE,ncol=3)
n <- nrow(X)

(x.bar <- t(matrix(1,ncol=n) %*% X)/n)

(D <- X - matrix(1,nrow=n) %*% t(x.bar))

(S <- (n-1)^(-1) * t(D)%*%D)

D[,1]+2*D[,2]

det(S)

open3d()
vectors3d(D[,1],col="blue",lwd=2)
vectors3d(D[,2],col="red",lwd=2)
vectors3d(D[,3],col="purple",lwd=2)
axes3d()

# Example 3.10 
rm(list=ls())

X <- matrix(c(1,4,2,5,3,
              9,12,10,8,11,
              10,16,12,13,14),byrow=FALSE,ncol=3)
n <- nrow(X)
(x.bar <- t(matrix(1,ncol=n) %*% X)/n)
(D <- X - matrix(1,nrow=n) %*% t(x.bar))
(S <- (n-1)^(-1) * t(D)%*%D)
det(S)

# Example 3.11
rm(list=ls())

S <- matrix(c(4,3,1,
              3,9,2,
              1,2,1),byrow=TRUE,ncol=3)
V12 <- diag(sqrt(diag(S)))
V12.inv <- diag(sqrt(1/diag(S)))
R <- V12.inv %*% S %*% V12.inv
fractions(R)

det(S)
det(R)
prod(diag(S))*det(R)

# Example 3.12
rm(list=ls())

S1 <- matrix(c(252.04,-68.43,
              -68.43,123.67),byrow = TRUE, ncol=2)
tr(S1)

S2 <- matrix(c(3,-3/2,0,
               -3/2,1,1/2,
               0,1/2,1),byrow=TRUE,ncol=3)
tr(S2)