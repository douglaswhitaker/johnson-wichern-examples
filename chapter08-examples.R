# Example 8.1
rm(list=ls())

# Note that this is the covariance matrix from Exercise 4.3 about independence of MV Normal random variables 
S <- matrix(c(1,-2,0,
              -2,5,0,
              0,0,2),ncol=3,byrow=TRUE)

(l <- eigen(S)$values)
(e <- eigen(S)$vectors)

# Coefficients for principle components (linear combinations)
#                1st PC      2nd PC   3rd PC
#       X1 coef -0.3826834    0       0.9238795
#       X2 coef  0.9238795    0       0.3826834
#       X3 coef  0.0000000    1       0.0000000

# Total variance
(tv <- sum(diag(S)))

# Proportion of variance explained by each PC:
(l/tv)

# Correlation between Y1 and X1
e[1,1]*sqrt(l[1])/sqrt(S[1,1])
# Correlation between Y1 and X2
e[2,1]*sqrt(l[1])/sqrt(S[2,2])
# Correlation between Y2 and X1
e[1,2]*sqrt(l[2])/sqrt(S[1,1])
# Correlation between Y2 and X2
e[2,2]*sqrt(l[2])/sqrt(S[2,2])
# Correlation between Y2 and X3
e[3,2]*sqrt(l[2])/sqrt(S[3,3])

# Example 8.2
rm(list=ls())

S <- matrix(c(1,4,4,100),ncol=2,byrow=TRUE)
(R <- cov2cor(S))

# PCA using Covariance Matrix
(l.s <- eigen(S)$values)
(e.s <- eigen(S)$vectors)
# Proportion of variability explained
(l.s/sum(l.s))

# PCA using Correlation Matrix
(l.r <- eigen(R)$values)
(e.r <- eigen(R)$vectors)
# Proportion of variability explained
p <- ncol(R)
(l.r/p)

# Example 8.3
rm(list=ls())

dat <- read.table("data/T8-5.DAT")
names(dat) <- c("TotalPop","ProfDegree","Employed","GovEmp","MedianHome")

X <- as.matrix(dat)
n <- nrow(X)
p <- ncol(X)
(x.bar <- t(matrix(1,ncol=n) %*% X)/n)
D <- X - matrix(1,nrow=n) %*% t(x.bar)
(S <- (n-1)^(-1) * t(D)%*%D)
Sinv <- solve(S)

(l <- eigen(S)$values)
(e <- eigen(S)$vectors)

l/sum(l)
sum((l/sum(l))[1:2]) # 92% of variability explained by two principle components

e[,1]*sqrt(l[1])/sqrt(diag(S))
e[,2]*sqrt(l[2])/sqrt(diag(S))

# Example 8.4
rm(list=ls())

dat <- read.table("data/T6-9.dat")
X1 <- dat[which(dat[,4]=="male"),1:3]
X <- sapply(X1,log)
n <- nrow(X)
p <- ncol(X)
(x.bar <- t(matrix(1,ncol=n) %*% X)/n)
D <- X - matrix(1,nrow=n) %*% t(x.bar)
(S <- (n-1)^(-1) * t(D)%*%D)


(l <- eigen(S)$values)
(e <- eigen(S)$vectors)

l/sum(l)

e[,1]*sqrt(l[1])/sqrt(diag(S))
e[,2]*sqrt(l[2])/sqrt(diag(S))
e[,3]*sqrt(l[3])/sqrt(diag(S))

plot(1:length(l),l,type="l")

# Using built-in functions
prcomp(X,scale=FALSE) # scale=FALSE is for using covariance matrix
# note that there is another function called princomp() - it is included for compatibility reasons. prcomp() is preferred.
screeplot(prcomp(X),type="lines")

# Example 8.5 
rm(list=ls())

dat <- read.table("data/T8-4.DAT")
names(dat) <- c("JPMorgan","Citibank","WellsFargo","RoyalDutchShell","ExxonMobil")
X <- as.matrix(dat)
n <- nrow(X)
p <- ncol(X)
(x.bar <- t(matrix(1,ncol=n) %*% X)/n)
D <- X - matrix(1,nrow=n) %*% t(x.bar)
(S <- (n-1)^(-1) * t(D)%*%D)
(R <- cov2cor(S))


(l <- eigen(R)$values)
(e <- eigen(R)$vectors)

prcomp(dat,scale=TRUE) # scale=TRUE is for using correlation matrix

l/p
sum(l[1:2])/p

# Example 8.7
rm(list=ls())

dat <- read.table("data/T6-9.dat")
X1 <- dat[which(dat[,4]=="male"),1:3]
X <- sapply(X1,log)
n <- nrow(X)
p <- ncol(X)
(x.bar <- t(matrix(1,ncol=n) %*% X)/n)
D <- X - matrix(1,nrow=n) %*% t(x.bar)
(S <- (n-1)^(-1) * t(D)%*%D)


(l <- eigen(S)$values)
(e <- eigen(S)$vectors)

l/sum(l)

e[,1]*sqrt(l[1])/sqrt(diag(S))
e[,2]*sqrt(l[2])/sqrt(diag(S))
e[,3]*sqrt(l[3])/sqrt(diag(S))

qqnorm(X %*% e[,2])
plot(X %*% e[,2], X %*% e[,1])

# Using built-in functions
pr.fit <- prcomp(X,scale=FALSE)
qqnorm(pr.fit$x[,2])
plot(pr.fit$x[,2],pr.fit$x[,1])

# Example 8.8
rm(list=ls())

dat <- read.table("data/T8-4.DAT")
names(dat) <- c("JPMorgan","Citibank","WellsFargo","RoyalDutchShell","ExxonMobil")
X <- as.matrix(dat)
n <- nrow(X)
p <- ncol(X)
(x.bar <- t(matrix(1,ncol=n) %*% X)/n)
D <- X - matrix(1,nrow=n) %*% t(x.bar)
(S <- (n-1)^(-1) * t(D)%*%D)

# Note this example is based on Sigma, not Rho
(l <- eigen(S)$values)
(e <- eigen(S)$vectors)

alpha <- 0.05
(lb <- l[1]/(1+qnorm(1-alpha/2)*sqrt(2/n)))
(ub <- l[1]/(1-qnorm(1-alpha/2)*sqrt(2/n)))

# Example 8.6
rm(list=ls())

x.bar <- matrix(c(39.88, 45.08, 48.11, 49.95),ncol=1)
R <- matrix(c(1.000, .7501, .6329, .6363,
              .7501, 1.000, .6925, .7386,
              .6329, .6925, 1.000, .6625,
              .6363, .7386, .6625, 1.000), ncol=4, byrow=TRUE)

(l <- eigen(R)$values)
(e <- eigen(R)$vectors)

(r.bar <- mean(R[lower.tri(R)]))
(1+(4-1)*r.bar)
l[1]
l[2:4]

# Example 8.9 
rm(list=ls())

R <- matrix(c(1.000, .7501, .6329, .6363,
              .7501, 1.000, .6925, .7386,
              .6329, .6925, 1.000, .6625,
              .6363, .7386, .6625, 1.000), ncol=4, byrow=TRUE)
n <- 150
p <- 4
alpha <- 0.05

rk <- c()
(rk[1] <- (1/(p-1))*(sum(R[,1])-1))
(rk[2] <- (1/(p-1))*(sum(R[,2])-1))
(rk[3] <- (1/(p-1))*(sum(R[,3])-1))
(rk[4] <- (1/(p-1))*(sum(R[,4])-1))
(r.bar <- mean(R[lower.tri(R)]))
(part1 <- sum((R[lower.tri(R)]-r.bar)^2))
(part2 <- sum((rk-r.bar)^2))
(gammahat <- ((p-1)^2 * (1-(1-r.bar)^2))/(p-(p-2)*(1-r.bar)^2))
(T2 <- (n-1)/(1-r.bar)^2 * (part1 - gammahat*part2))
(T2.crit <- qchisq(1-alpha,(p+1)*(p-2)/2))

# Example 8.10
rm(list=ls())

dat <- read.table("data/T5-8.DAT")
names(dat) <- c("Legal","Extraordinary","Holdover","COA","Meeting")

X <- as.matrix(dat)
n <- nrow(X)
p <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D

(l <- eigen(S)$values)
(e <- eigen(S)$vectors)

alpha <- 0.05

c2 <- qchisq(1-alpha,df=2)

library(plotrix)
angle <- 0 # principal components uncorrelated
plot(0,pch='',ylab='',xlab='',xlim=c(-6000,5000),ylim=c(-4000,3000))
lengths <- c(sqrt(c2*eigen(S)$values[1]),
             sqrt(c2*eigen(S)$values[2]))
# Centered at 0
draw.ellipse(x=0,y=0,a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)
# Mean-centered 
points(t(e[,1])%*%t(X-(t(matrix(1,ncol=n)) %*% t(X.mean))),
       t(e[,2])%*%t(X-(t(matrix(1,ncol=n)) %*% t(X.mean))))

# Example 8.11
T2j <- (t(e[,3])%*%t(X-(t(matrix(1,ncol=n)) %*% t(X.mean))))^2/l[3] + 
        (t(e[,4])%*%t(X-(t(matrix(1,ncol=n)) %*% t(X.mean))))^2/l[4] + 
        (t(e[,5])%*%t(X-(t(matrix(1,ncol=n)) %*% t(X.mean))))^2/l[5]
plot(y=T2j,x=1:length(T2j),type="l")
alpha <- 0.05
p <- 5
UCL <- qchisq(1-alpha,p-2)
abline(h=UCL)

# Example 8.12 
# observation 11 out of control

X <- as.matrix(dat[-11,])
n <- nrow(X)
p <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D

(l <- eigen(S)$values)
(e <- eigen(S)$vectors)


alpha <- 0.01

c2 <- qchisq(1-alpha,df=2)

library(plotrix)
angle <- 0 # principal components uncorrelated
plot(0,pch='',ylab='',xlab='',xlim=c(-6000,5000),ylim=c(-4000,3000))
lengths <- c(sqrt(c2*eigen(S)$values[1]),
             sqrt(c2*eigen(S)$values[2]))
# Centered at 0
draw.ellipse(x=0,y=0,a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)
# Mean-centered 
points(t(e[,1])%*%t(X-(t(matrix(1,ncol=n)) %*% t(X.mean))),
       t(e[,2])%*%t(X-(t(matrix(1,ncol=n)) %*% t(X.mean))))