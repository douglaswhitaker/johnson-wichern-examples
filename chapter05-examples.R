# Example 5.1
rm(list=ls())

X <- matrix(c(6,10,8,
              9,6,3),byrow=FALSE,ncol=2)
mu0 <- matrix(c(9,5),ncol=1)
n <- nrow(X)
p <- ncol(X)
(X.mean <- t(matrix(1,ncol=n) %*% X)/n)
(D <- X - matrix(1,nrow=n) %*% t(X.mean))
(S <- (n-1)^(-1) * t(D)%*%D)
Sinv <- solve(S)

T2 <- n * t(X.mean-mu0) %*% Sinv %*% (X.mean-mu0)

# A priori, T^2 is distributed ((n-1)*p)/(n-p)*qf(1-alpha,p,n-p)

(3-1)*2/(3-2)*qf(1-0.05,2,3-2)

# Example 5.2
rm(list=ls())

alpha <- 0.10
mu0 <- matrix(c(4,50,10),ncol=1)

dat <- read.table("data/T5-1.DAT")
names(dat) <- c("Sweat_Rate","Sodium","Potassium")

qqnorm(dat$Sweat_Rate)
qqnorm(dat$Sodium)
qqnorm(dat$Potassium)

X <- as.matrix(dat)
n <- nrow(X)
p <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D
Sinv <- solve(S)

(T2 <- n * t(X.mean-mu0) %*% Sinv %*% (X.mean-mu0)) # observed
(T.crit <- ((n-1)*p)/(n-p)*qf(1-alpha,p,n-p)) # critical

# Example 5.3
rm(list=ls())

alpha <- 0.05

dat.closed <- read.table("data/T4-1.DAT") # Door closed
dat.open <- read.table("data/T4-5.DAT") # Door Open
dat <- data.frame(closed=dat.closed[,1]^.25, 
                  open=dat.open[,1]^.25)

X <- as.matrix(dat)
n <- nrow(X)
p <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D
Sinv <- solve(S)

inRegion <- function(mu, X.mean, Sinv, n, alpha){
  critical.value <- (p*(n-1))/(n-p)*qf(1-alpha,p,n-p)
  return(n*t(X.mean-mu)%*%Sinv%*%(X.mean-mu) <= critical.value)
}

inRegion(mu=matrix(c(0.562,0.589),
                   ncol=1),
         X.mean=X.mean,
         Sinv=Sinv,
         n=n,
         alpha=alpha)

library(plotrix)
angle <- atan(eigen(S)$vectors[2,1]/eigen(S)$vectors[1,1]) # sohcahtoa
plot(0,pch='',ylab='',xlab='',xlim=c(0.5,0.65),ylim=c(0.55,0.65))
axis1 <- sqrt(eigen(S)$values[1])*
  sqrt(
    (p*(n-1))/
      (n*(n-p))* # not the same as critical.value above
      qf(1-alpha,p,n-p)) 
axis2 <- sqrt(eigen(S)$values[2])*
  sqrt((p*(n-1))/(n*(n-p))*qf(1-alpha,p,n-p))
lengths <- c(axis1,axis2)
draw.ellipse(x=X.mean[1,1],y=X.mean[2,1],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)

sqrt(eigen(S)$values[1])/sqrt(eigen(S)$values[2])

# Example 5.4
a1 <- matrix(c(1,0),ncol=1)
a2 <- matrix(c(0,1),ncol=1)

l1 <- X.mean[1,1] - sqrt((p*(n-1))/(n*(n-p))*qf(1-alpha,p,n-p)*t(a1)%*%S%*%a1)
u1 <- X.mean[1,1] + sqrt((p*(n-1))/(n*(n-p))*qf(1-alpha,p,n-p)*t(a1)%*%S%*%a1)

l2 <- X.mean[2,1] - sqrt((p*(n-1))/(n*(n-p))*qf(1-alpha,p,n-p)*t(a2)%*%S%*%a2)
u2 <- X.mean[2,1] + sqrt((p*(n-1))/(n*(n-p))*qf(1-alpha,p,n-p)*t(a2)%*%S%*%a2)

l1;u1
l2;u2

abline(v=l1,lty=3);abline(v=u1,lty=3)
abline(h=l2,lty=3);abline(h=u2,lty=3)

# Example 5.5
rm(list=ls())

alpha <- 0.05

dat <- read.table(file="T5-2.DAT")
names(dat) <- c("Social","Verbal","Science")

X <- as.matrix(dat)
n <- nrow(X)
p <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D
Sinv <- solve(S)

a1 <- matrix(c(1,0,0),ncol=1)
a2 <- matrix(c(0,1,0),ncol=1)
a3 <- matrix(c(0,0,1),ncol=1)

l1 <- X.mean[1,1] - sqrt((p*(n-1))/(n*(n-p))*qf(1-alpha,p,n-p)*t(a1)%*%S%*%a1)
u1 <- X.mean[1,1] + sqrt((p*(n-1))/(n*(n-p))*qf(1-alpha,p,n-p)*t(a1)%*%S%*%a1)

l2 <- X.mean[2,1] - sqrt((p*(n-1))/(n*(n-p))*qf(1-alpha,p,n-p)*t(a2)%*%S%*%a2)
u2 <- X.mean[2,1] + sqrt((p*(n-1))/(n*(n-p))*qf(1-alpha,p,n-p)*t(a2)%*%S%*%a2)

l3 <- X.mean[3,1] - sqrt((p*(n-1))/(n*(n-p))*qf(1-alpha,p,n-p)*t(a3)%*%S%*%a3)
u3 <- X.mean[3,1] + sqrt((p*(n-1))/(n*(n-p))*qf(1-alpha,p,n-p)*t(a3)%*%S%*%a3)

l1;u1
l2;u2
l3;u3

a23 <- matrix(c(0,1,-1),ncol=1)
l23 <- (X.mean[2,1] - X.mean[3,1]) - sqrt((p*(n-1))/(n*(n-p))*qf(1-alpha,p,n-p)*t(a23)%*%S%*%a23)
u23 <- (X.mean[2,1] - X.mean[3,1]) + sqrt((p*(n-1))/(n*(n-p))*qf(1-alpha,p,n-p)*t(a23)%*%S%*%a23)
l23;u23

S23 <- S[-1,-1]

library(plotrix)
angle <- atan(eigen(S23)$vectors[2,1]/eigen(S23)$vectors[1,1]) # sohcahtoa
plot(0,pch='',ylab='',xlab='',xlim=c(50,59),ylim=c(23,27))
axis1 <- sqrt(eigen(S23)$values[1])*sqrt((p*(n-1))/(n*(n-p))*qf(1-alpha,p,n-p)) 
axis2 <- sqrt(eigen(S23)$values[2])*sqrt((p*(n-1))/(n*(n-p))*qf(1-alpha,p,n-p))
lengths <- c(axis1,axis2)
draw.ellipse(x=X.mean[2,1],y=X.mean[3,1],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)

# Example 5.6
rm(list=ls())
# Run all of the Example 5.3 and 5.4 code first

b <- 2
t.bonf <- qt(1-alpha/(2*b),df=n-1)

(l1b <- X.mean[1,1] - t.bonf * sqrt(S[1,1]/n))
(u1b <- X.mean[1,1] + t.bonf * sqrt(S[1,1]/n))
(l2b <- X.mean[2,1] - t.bonf * sqrt(S[2,2]/n))
(u2b <- X.mean[2,1] + t.bonf * sqrt(S[2,2]/n))


abline(v=l1b,lty=2);abline(v=u1b,lty=2)
abline(h=l2b,lty=2);abline(h=u2b,lty=2)

# Example 5.7
rm(list=ls())

means <- c(28.1,26.6,35.4,34.2,23.6,22.0,22.7)
sds <- c(5.76,5.85,3.82,5.12,3.76,3.93,4.03)

alpha <- 0.10
n <- 96

confval <- sqrt(qchisq(1-alpha,df=7))

for (i in 1:7){
  print(i)
  print(means[i]-confval*sds[i]/sqrt(n))
  print(means[i]+confval*sds[i]/sqrt(n))
}

# Example 5.8
rm(list=ls())

dat <- read.table("data/T5-8.DAT")
names(dat) <- c("Legal","Extraordinary","Holdover","COA","Meeting")
n <- nrow(dat)

center <- mean(dat$Legal)
UL <- center + 3*sd(dat$Legal)
LL <- center - 3*sd(dat$Legal)

plot(x=1:n,dat$Legal,type="l",ylim=c(1500,5500),xlab="Observation Number",ylab="Individual Value")
abline(h=center);abline(h=UL);abline(h=LL)

# Example 5.9
X <- as.matrix(dat[,1:2])
n <- nrow(X)
p <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D

alpha <- 0.01
c2 <- qchisq(1-alpha,df=2)

library(plotrix)
angle <- atan(eigen(S)$vectors[2,1]/eigen(S)$vectors[1,1]) # sohcahtoa
plot(0,pch='',ylab='',xlab='',xlim=c(1000,6000),ylim=c(-2000,6000))
points(X)
lengths <- c(sqrt(c2*eigen(S)$values[1]),
             sqrt(c2*eigen(S)$values[2]))
draw.ellipse(x=X.mean[1,1],y=X.mean[2,1],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)

center <- mean(dat$Extraordinary)
UL <- center + 3*sd(dat$Extraordinary)
LL <- center - 3*sd(dat$Extraordinary)

plot(x=1:n,dat$Legal,type="l",ylim=c(-3000,6000),xlab="Observation Number",ylab="Individual Value")
abline(h=center);abline(h=UL);abline(h=LL)

# Example 5.10
Sinv <- solve(S)
d <- c()
for (i in 1:n){
  d[i] <- (X[i,]-t(X.mean))%*%Sinv%*%t(X[i,]-t(X.mean))
}

UCL <- qchisq(1-alpha,df=2)

plot(x=1:length(d),y=d,ylim=c(0,12),xlab="Period",ylab="T2")
abline(h=UCL,lty=3)
identify(x=1:length(d),y=d)

# Example 5.11
rm(list=ls())

dat <- read.table("data/T5-9.dat")
names(dat) <- c("Voltage","Current","Feed_Speed","Gas_Flow")

qqnorm(dat[,1]);qqnorm(dat[,2]);qqnorm(dat[,3]);qqnorm(dat[,4])
qqnorm(log(dat$Gas_Flow))

dat <- data.frame(dat[,1:3],log_Gas_Flow=log(dat$Gas_Flow))

X <- as.matrix(dat)
n <- nrow(X)
p <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D
Sinv <- solve(S)

# T^2 chart
d <- c()
for (i in 1:n){
  d[i] <- (X[i,]-t(X.mean))%*%Sinv%*%t(X[i,]-t(X.mean))
}

UCL95 <- qchisq(1-0.05,df=p)
UCL99 <- qchisq(1-0.01,df=p)

plot(x=1:length(d),y=d,ylim=c(0,14),xlab="Case",ylab="T2")
abline(h=UCL95,lty=2)
abline(h=UCL99,lty=1)

# 99% quality control ellipse
X14 <- X[,c(1,4)]
X14.mean <- t(matrix(1,ncol=n) %*% X14)/n
D14 <- X14 - matrix(1,nrow=n) %*% t(X14.mean)
S14 <- (n-1)^(-1) * t(D14)%*%D14

c2 <- qchisq(1-0.01,df=2)
library(plotrix)
angle <- atan(eigen(S14)$vectors[2,1]/eigen(S14)$vectors[1,1]) # sohcahtoa
plot(0,pch='',ylab='',xlab='',ylim=c(3.85,4.05),xlim=c(20.5,24.0))
points(X14)
lengths <- c(sqrt(c2*eigen(S14)$values[1]),
             sqrt(c2*eigen(S14)$values[2]))
draw.ellipse(x=X14.mean[1,1],y=X14.mean[2,1],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)

# Univariate Xbarbar chart
center <- mean(dat$log_Gas_Flow)
UL <- center + 3*sd(dat$log_Gas_Flow)
LL <- center - 3*sd(dat$log_Gas_Flow)

plot(x=1:n,dat$log_Gas_Flow,type="l",ylim=c(3.85,4.05),xlab="Observation Number",ylab="Individual Value")
abline(h=center);abline(h=UL);abline(h=LL)


# Example 5.12
rm(list=ls())

# Run code for Examples 5.8-5.10 here
dat.stable <- dat[-11,]


X <- as.matrix(dat.stable[,1:2])
n <- nrow(X)
p <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D
Sinv <- solve(S)


alpha <- 0.05
c2 <- (2*(n^2-1))/(n*(n-2))*qf(1-alpha,2,n-2)

library(plotrix)
angle <- atan(eigen(S)$vectors[2,1]/eigen(S)$vectors[1,1]) # sohcahtoa
plot(0,pch='',ylab='',xlab='',xlim=c(1500,5500),ylim=c(-500,3000))
points(X)
lengths <- c(sqrt(c2*eigen(S)$values[1]),
             sqrt(c2*eigen(S)$values[2]))
draw.ellipse(x=X.mean[1,1],y=X.mean[2,1],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)

# Example 5.13
rm(list=ls())

X <- matrix(c(NA,7,5,NA,
              0,2,1,NA,
              3,6,2,5),byrow=FALSE,ncol=3)

# Note: the code below is designed to be illustrative, not practical
# In practice, a function should be made that performs the steps automatically and 
# stops when some stop condition is met. 
# (e.g. estimated mu and Sigma matrices change by less than epsilon on an iteration)


# Estimation Step 1

mu1 <- mean(na.omit(X[,1]))
mu2 <- mean(na.omit(X[,2]))
mu3 <- mean(na.omit(X[,3]))

X0 <- X
X0[is.na(X[,1]),1] <- mu1
X0[is.na(X[,2]),2] <- mu2
X0[is.na(X[,3]),3] <- mu3

mut <- matrix(c(mu1,mu2,mu3),ncol=1)

X2 <- X0
n <- nrow(X2)
p <- ncol(X2)
X2.mean <- t(matrix(1,ncol=n) %*% X2)/n
D2 <- X2 - matrix(1,nrow=n) %*% t(X2.mean)
Snt <- (n)^(-1) * t(D2)%*%D2
 
# Prediction Step 1

# Missing first column entry in row 1 of X
x11 <- mut[1,1] + 
  Snt[c(1),c(2:3)] %*% 
  solve(Snt[2:3,2:3]) %*%
  (X[1,2:3]-mut[2:3,1])

x2.11 <- Snt[1,1] - 
  Snt[1,2:3] %*% 
  solve(Snt[2:3,2:3]) %*%
  Snt[2:3,1] + x11^2

x11x12 <- c(x11) * X[1,2:3]

# Missing two components in row 4 of X

x41 <- mut[1:2,1] + 
  Snt[1:2,3] %*% 
  solve(Snt[3,3]) %*%
  (X[4,3]-mut[3,1])

x2.41 <- Snt[1:2,1:2] - 
  Snt[1:2,3] %*% 
  solve(Snt[3,3]) %*%
  Snt[3,1:2] + 
  (x41 %*% t(x41))

x41xxx <- x41 * c(X[4,3])

X2[1,1] <- x11
X2[4,1:2] <- c(x41)

T1 <- t(t(colSums(X2)))

(T2.11 <- x2.11 + X[2,1]^2 + X[3,1]^2 + x2.41[1,1])
(T2.22 <- X[1,2]^2 + X[2,2]^2 + X[3,2]^2 + x2.41[2,2])
(T2.33 <- sum(X[,3]^2))

(T2.21 <- x11x12[1] + X[2,1]*X[2,2] + X[3,1]*X[3,2] + x2.41[1,2])
(T2.31 <- x11x12[2] + X[2,1]*X[2,3] + X[3,1]*X[3,3] + x41xxx[1])
(T2.32 <- X[1,2]*X[1,3] + X[2,2]*X[2,3] + X[3,2]*X[3,3] + x41xxx[2])

T2 <- matrix(c(T2.11, T2.21, T2.31,
               T2.21, T2.22, T2.32,
               T2.31, T2.32, T2.33),byrow=TRUE,ncol=3)

# Estimation Step 2, using formula (5-40)
(mut <- (1/n)*T1)
(Snt <- (1/n)*T2-mut%*%t(mut))




##### Copy and pasted code several times below here, no new code, just running until mut and Snt stop changing much
# Prediction Step 2 (copied code from above)

# Missing first column entry in row 1 of X
x11 <- mut[1,1] + 
  Snt[c(1),c(2:3)] %*% 
  solve(Snt[2:3,2:3]) %*%
  (X[1,2:3]-mut[2:3,1])

x2.11 <- Snt[1,1] - 
  Snt[1,2:3] %*% 
  solve(Snt[2:3,2:3]) %*%
  Snt[2:3,1] + x11^2

x11x12 <- c(x11) * X[1,2:3]

# Missing two components in row 4 of X

x41 <- mut[1:2,1] + 
  Snt[1:2,3] %*% 
  solve(Snt[3,3]) %*%
  (X[4,3]-mut[3,1])

x2.41 <- Snt[1:2,1:2] - 
  Snt[1:2,3] %*% 
  solve(Snt[3,3]) %*%
  Snt[3,1:2] + 
  (x41 %*% t(x41))

x41xxx <- x41 * c(X[4,3])

X2[1,1] <- x11
X2[4,1:2] <- c(x41)

T1 <- t(t(colSums(X2)))

(T2.11 <- x2.11 + X[2,1]^2 + X[3,1]^2 + x2.41[1,1])
(T2.22 <- X[1,2]^2 + X[2,2]^2 + X[3,2]^2 + x2.41[2,2])
(T2.33 <- sum(X[,3]^2))

(T2.21 <- x11x12[1] + X[2,1]*X[2,2] + X[3,1]*X[3,2] + x2.41[1,2])
(T2.31 <- x11x12[2] + X[2,1]*X[2,3] + X[3,1]*X[3,3] + x41xxx[1])
(T2.32 <- X[1,2]*X[1,3] + X[2,2]*X[2,3] + X[3,2]*X[3,3] + x41xxx[2])

T2 <- matrix(c(T2.11, T2.21, T2.31,
               T2.21, T2.22, T2.32,
               T2.31, T2.32, T2.33),byrow=TRUE,ncol=3)

# Estimation Step 3, using formula (5-40)
(mut <- (1/n)*T1)
(Snt <- (1/n)*T2-mut%*%t(mut))

# Prediction Step 3 (copied code from above)

# Missing first column entry in row 1 of X
x11 <- mut[1,1] + 
  Snt[c(1),c(2:3)] %*% 
  solve(Snt[2:3,2:3]) %*%
  (X[1,2:3]-mut[2:3,1])

x2.11 <- Snt[1,1] - 
  Snt[1,2:3] %*% 
  solve(Snt[2:3,2:3]) %*%
  Snt[2:3,1] + x11^2

x11x12 <- c(x11) * X[1,2:3]

# Missing two components in row 4 of X

x41 <- mut[1:2,1] + 
  Snt[1:2,3] %*% 
  solve(Snt[3,3]) %*%
  (X[4,3]-mut[3,1])

x2.41 <- Snt[1:2,1:2] - 
  Snt[1:2,3] %*% 
  solve(Snt[3,3]) %*%
  Snt[3,1:2] + 
  (x41 %*% t(x41))

x41xxx <- x41 * c(X[4,3])

X2[1,1] <- x11
X2[4,1:2] <- c(x41)

T1 <- t(t(colSums(X2)))

(T2.11 <- x2.11 + X[2,1]^2 + X[3,1]^2 + x2.41[1,1])
(T2.22 <- X[1,2]^2 + X[2,2]^2 + X[3,2]^2 + x2.41[2,2])
(T2.33 <- sum(X[,3]^2))

(T2.21 <- x11x12[1] + X[2,1]*X[2,2] + X[3,1]*X[3,2] + x2.41[1,2])
(T2.31 <- x11x12[2] + X[2,1]*X[2,3] + X[3,1]*X[3,3] + x41xxx[1])
(T2.32 <- X[1,2]*X[1,3] + X[2,2]*X[2,3] + X[3,2]*X[3,3] + x41xxx[2])

T2 <- matrix(c(T2.11, T2.21, T2.31,
               T2.21, T2.22, T2.32,
               T2.31, T2.32, T2.33),byrow=TRUE,ncol=3)

# Estimation Step 4, using formula (5-40)
(mut <- (1/n)*T1)
(Snt <- (1/n)*T2-mut%*%t(mut))


###########
library(mvnmle)
mlest(X)
mlest(X,iterlim=1000)
