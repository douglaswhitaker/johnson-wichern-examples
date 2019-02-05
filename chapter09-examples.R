# Example 9.1
rm(list=ls())

Sigma <- matrix(c(19,30,2,12,
                  30,57,5,23,
                  2,5,38,47,
                  12,23,47,68),ncol=4,byrow=TRUE)
L <- matrix(c(4,1,
              7,2,
              -1,6,
              1,8),ncol=2,byrow=TRUE)
Psi <- diag(c(2,4,1,3))
(L%*%t(L)+Psi)
Sigma==L%*%t(L)+Psi

# Communality of X1
sum((L[1,])^2)

# Variance decomposition 
Sigma[1,1] == sum((L[1,])^2) + Psi[1,1]

# Example 9.4 (based on Example 8.5)
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

l/p
sum(l[1:2])/p


# One-Factor Solution
(L1 <- t(t(sqrt(l[1])*e[,1])))
(Psi1 <- diag(1-(L1^2)[,1]))
R - L1%*%t(L1) - Psi1 # Residual matrix
#   Communalities
rowSums(L1^2)
# Specific Variances
1-rowSums(L1^2)
# Proportion of variability explained
sum(l[1])/p

# Two-Factor Solution
(L2 <- cbind(t(t(sqrt(l[1])*e[,1])),t(t(sqrt(l[2])*e[,2]))))
(Psi2 <- diag(1-rowSums(L2^2)))
R - L2%*%t(L2) - Psi2 # Residual matrix
# Communalities
rowSums(L2^2)
# Specific Variances
1-rowSums(L2^2)
# Proportion of variability explained
sum(l[1:2])/p


# Example 9.10
factanal(dat,2,rotation="none") # uses MLE
factanal(dat,2,rotation="varimax")


# Example Big 5 http://openpsychometrics.org/_rawdata/BIG5.zip 
rm(list=ls())

dat <- read.table(file="data/BIG5/data.csv",sep="\t",header=TRUE)
factanal(dat[,8:57],5,rotation="none")$loadings
fit <- factanal(dat[,8:57],5,rotation="varimax")
print(fit$loadings,cutoff=0.3)
