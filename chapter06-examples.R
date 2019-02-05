
# Example 6.1
rm(list=ls())

dat <- read.table("data/T6-1.dat")
names(dat) <- c("BOD.com","SS.com","BOD.sta","SS.sta")

alpha <- 0.05
n <- nrow(dat)
p <- 2
mu0 <- matrix(c(0,0),ncol=1)

dj1 <- dat$BOD.com - dat$BOD.sta
dj2 <- dat$SS.com - dat$SS.sta
d <- matrix(c(dj1,dj2),ncol=2,nrow=n)
d.bar <- matrix(c(mean(dj1),mean(dj2)),ncol=1)
Sd <- var(d)

(T2 <- n*t(d.bar)%*%solve(Sd)%*%d.bar)

(T2.crit <- p*(n-1)/(n-p)*qf(1-alpha,p,n-p))

ifelse(T2>T2.crit,"Reject H0","Fail to Reject H0")

d.bar[1,1] - sqrt(p*(n-1)/(n-p)*qf(1-alpha,p,n-p))*sqrt(Sd[1,1]/n)
d.bar[1,1] + sqrt(p*(n-1)/(n-p)*qf(1-alpha,p,n-p))*sqrt(Sd[1,1]/n)

d.bar[2,1] - sqrt(p*(n-1)/(n-p)*qf(1-alpha,p,n-p))*sqrt(Sd[2,2]/n)
d.bar[2,1] + sqrt(p*(n-1)/(n-p)*qf(1-alpha,p,n-p))*sqrt(Sd[2,2]/n)


# Example 6.2
rm(list=ls())
alpha <- 0.05

dat <- read.table("data/T6-2.dat")
C <- matrix(c(-1,-1,1,1,
              1,-1,1,-1,
              1,-1,-1,1),byrow=TRUE,ncol=4)
X <- as.matrix(dat)
n <- nrow(X)
q <- ncol(X)
X.mean <- t(matrix(1,ncol=n) %*% X)/n
D <- X - matrix(1,nrow=n) %*% t(X.mean)
S <- (n-1)^(-1) * t(D)%*%D
Sinv <- solve(S)

C%*%X.mean
C%*%S%*%t(C)

(T2 <- n*t(C%*%X.mean)%*%solve(C%*%S%*%t(C))%*%(C%*%X.mean))

(Critical.T2 <- ((n-1)*(q-1)/(n-q+1))*qf(1-alpha,q-1,n-q+1))

# Confidence Intervals
# Interaction
C[3,]%*%X.mean - sqrt(((n-1)*(q-1)/(n-q+1))*qf(1-alpha,q-1,n-q+1))*sqrt(t(C[3,])%*%S%*%C[3,]/n)
C[3,]%*%X.mean + sqrt(((n-1)*(q-1)/(n-q+1))*qf(1-alpha,q-1,n-q+1))*sqrt(t(C[3,])%*%S%*%C[3,]/n)

# Halothane Influence
C[1,]%*%X.mean - sqrt(((n-1)*(q-1)/(n-q+1))*qf(1-alpha,q-1,n-q+1))*sqrt(t(C[1,])%*%S%*%C[1,]/n)
C[1,]%*%X.mean + sqrt(((n-1)*(q-1)/(n-q+1))*qf(1-alpha,q-1,n-q+1))*sqrt(t(C[1,])%*%S%*%C[1,]/n)

# CO2 pressure influence
C[2,]%*%X.mean - sqrt(((n-1)*(q-1)/(n-q+1))*qf(1-alpha,q-1,n-q+1))*sqrt(t(C[2,])%*%S%*%C[2,]/n)
C[2,]%*%X.mean + sqrt(((n-1)*(q-1)/(n-q+1))*qf(1-alpha,q-1,n-q+1))*sqrt(t(C[2,])%*%S%*%C[2,]/n)

# Example 6.3
rm(list=ls())

alpha <- 0.05
n1 <- 50
n2 <- 50
p <- 2
xbar1 <- matrix(c(8.3,4.1),ncol=1)
xbar2 <- matrix(c(10.2,3.9),ncol=1)
S1 <- matrix(c(2,1,1,6),ncol=2)
S2 <- matrix(c(2,1,1,4),ncol=2)

(Spooled <- (n1-1)/(n1+n2-2) * S1 + (n2-1)/(n1+n2-2) * S2)
(xbard <- xbar1 - xbar2)


library(plotrix)
angle <- atan(eigen(Spooled)$vectors[2,1]/eigen(Spooled)$vectors[1,1]) # sohcahtoa
plot(0,pch='',ylab='',xlab='',xlim=c(-3,1),ylim=c(-1,2))
c2 <- ((n1+n2-2)*p)/(n1+n2-p-1) * qf(1-alpha,p,n1+n2-p-1)
axis1 <- sqrt(eigen(Spooled)$values[1])*sqrt((1/n1+1/n2)*c2) 
axis2 <- sqrt(eigen(Spooled)$values[2])*sqrt((1/n1+1/n2)*c2) 

lengths <- c(axis1,axis2)
draw.ellipse(x=xbard[1,1],y=xbard[2,1],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)

# Example 6.4
rm(list=ls())
alpha <- 0.05
xbar1 <- matrix(c(204.4,556.6),ncol=1)
xbar2 <- matrix(c(130.0,355.0),ncol=1)
S1 <- matrix(c(13825.3,23823.4,23823.4,73107.4),ncol=2)
S2 <- matrix(c(8632.0,19616.7,19616.7,55964.5),ncol=2)
n1 <- 45
n2 <- 55
p <- 2
(Spooled <- (n1-1)/(n1+n2-2) * S1 + (n2-1)/(n1+n2-2) * S2)
(xbard <- xbar1 - xbar2)
(c2 <- ((n1+n2-2)*p)/(n1+n2-p-1) * qf(1-alpha,p,n1+n2-p-1))

xbard[1,1] - sqrt(c2)*sqrt((1/n1+1/n2)*Spooled[1,1])
xbard[1,1] + sqrt(c2)*sqrt((1/n1+1/n2)*Spooled[1,1])

xbard[2,1] - sqrt(c2)*sqrt((1/n1+1/n2)*Spooled[2,2])
xbard[2,1] + sqrt(c2)*sqrt((1/n1+1/n2)*Spooled[2,2])

library(plotrix)
angle <- atan(eigen(Spooled)$vectors[2,1]/eigen(Spooled)$vectors[1,1]) # sohcahtoa
plot(0,pch='',ylab='',xlab='',xlim=c(0,300),ylim=c(0,350))
c2 <- ((n1+n2-2)*p)/(n1+n2-p-1) * qf(1-alpha,p,n1+n2-p-1)
axis1 <- sqrt(eigen(Spooled)$values[1])*sqrt((1/n1+1/n2)*c2) 
axis2 <- sqrt(eigen(Spooled)$values[2])*sqrt((1/n1+1/n2)*c2) 

lengths <- c(axis1,axis2)
draw.ellipse(x=xbard[1,1],y=xbard[2,1],a=lengths[1],b=lengths[2],angle=angle,deg=FALSE)

# Example 6.5
1/n1 * S1 + 1/n2 * S2

a1 <- matrix(c(1,0),ncol=1)
a2 <- matrix(c(0,1),ncol=1)

xbard[1,1] - sqrt(qchisq(1-alpha,p))*sqrt(t(a1)%*%(1/n1 * S1 + 1/n2 * S2)%*%a1)
xbard[1,1] + sqrt(qchisq(1-alpha,p))*sqrt(t(a1)%*%(1/n1 * S1 + 1/n2 * S2)%*%a1)

xbard[2,1] - sqrt(qchisq(1-alpha,p))*sqrt(t(a2)%*%(1/n1 * S1 + 1/n2 * S2)%*%a2)
xbard[2,1] + sqrt(qchisq(1-alpha,p))*sqrt(t(a2)%*%(1/n1 * S1 + 1/n2 * S2)%*%a2)

(T2 <- t(xbard)%*%solve(1/n1 * S1 + 1/n2 * S2)%*%xbard)
(T2.critical <- qchisq(1-alpha,p))

(ahat <- solve(1/n1 * S1 + 1/n2 * S2)%*%xbard)

# Example 6.6
library(psych) # for trace function tr() (could also use something like sum(diag(X)) or sum(eigen(X)$values))
denom1 <- 1/n1 * (tr((1/n1 * S1 %*% solve(1/n1 * S1 + 1/n2 * S2)) %*% 
                       (1/n1 * S1 %*% solve(1/n1 * S1 + 1/n2 * S2)))+ (tr(1/n1 * S1 %*% solve(1/n1 * S1 + 1/n2 * S2)))^2 )
denom2 <- 1/n2 * (tr((1/n2 * S2 %*% solve(1/n1 * S1 + 1/n2 * S2)) %*% 
                       (1/n2 * S2 %*% solve(1/n1 * S1 + 1/n2 * S2)))+ (tr(1/n2 * S2 %*% solve(1/n1 * S1 + 1/n2 * S2)))^2 )

(v <- (p+p^2)/(denom1+denom2))

(T2)
(F.crit <- (v*p)/(v-p+1)*qf(1-alpha,p,v-p+1))

# Example 6.10
rm(list=ls())
# Examples 6.10-6.12 are limited by the summary statistics rather than raw data - calculations are bit off due to rounding
alpha <- 0.01

n1 <- 271
n2 <- 138
n3 <- 107
n <- n1 + n2 + n3
p <- 4
g <- 3

xbar1 <- matrix(c(2.066,0.480,0.082,0.360),ncol=1)
xbar2 <- matrix(c(2.167,0.596,0.124,0.418),ncol=1)
xbar3 <- matrix(c(2.273,0.521,0.125,0.383),ncol=1)

S1 <- matrix(c(.291,-.001,.002,.010,
               -.001,.011,.000,.003,
               .002,.000,.001,.000,
               .010,.003,.000,.010),ncol=4)
S2 <- matrix(c(.561,.011,.001,.037,
               .011,.025,.004,.007,
               .001,.004,.005,.002,
               .037,.007,.002,.019),ncol=4)
S3 <- matrix(c(.261,.030,.003,.018,
               .030,.017,-.000,.006,
               .003,-.000,.004,.001,
               .018,.006,.001,.013),ncol=4)

W <- (n1-1)*S1 + (n2-1)*S2 + (n3-1)*S3
xbar <- (1/(n1+n2+n3))*(n1*xbar1 + n2*xbar2 + n3*xbar3)
B <- n1 * (xbar1-xbar)%*%t(xbar1-xbar) + n2 * (xbar2-xbar)%*%t(xbar2-xbar) + n3 * (xbar3-xbar)%*%t(xbar3-xbar)

(Lambda <- det(W)/det(B+W))

(test.statistic <- (n-p-2)/p * (1-sqrt(Lambda))/sqrt(Lambda))
(critical.value <- qf(1-alpha,2*p,2*(n-p-2)))

(test.statistic.LS <- -(n-1-(p+g)/2)*log(Lambda))
(critical.value.LS <- qchisq(1-alpha,p*(g-1)))

# Example 6.11
tau1 <- (xbar1-xbar)
tau2 <- (xbar2-xbar)
tau3 <- (xbar3-xbar)
alpha <- 0.05

(tau1[3,1] - tau3[3,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n3)*W[3,3]/(n-g))
(tau1[3,1] - tau3[3,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n3)*W[3,3]/(n-g))


(tau1[3,1] - tau2[3,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n3)*W[3,3]/(n-g))
(tau1[3,1] - tau2[3,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n3)*W[3,3]/(n-g))

(tau2[3,1] - tau3[3,1]) - qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n3)*W[3,3]/(n-g))
(tau2[3,1] - tau3[3,1]) + qt(1-alpha/(p*g*2),n-g)*sqrt((1/n1 + 1/n3)*W[3,3]/(n-g))

# Example 6.12 

(Spooled <- (n1-1)/(n1+n2+n3-g) * S1 + (n2-1)/(n1+n2+n3-g) * S2 + (n3-1)/(n1+n2+n3-g) * S3)

u <- (1/(n1-1) + 1/(n2-1) + 1/(n3-1) - 1/(n1 + n2 + n3 - g)) * ((2*p^2 + 3*p - 1)/(6*(p+1)*(g-1)))
(M <- (n1 + n2 + n3 - g)*log(det(Spooled)) - ( (n1-1)*log(det(S1)) + (n2-1)*log(det(S2)) + (n3-1)*log(det(S3)) ) )

C <- (1-u)*M
v <- (1/2)*p*(p+1)*(g-1)
C.crit <- qchisq(1-alpha,v)

# Example 6.13
rm(list=ls())

dat <- read.table("data/T6-4.dat")
names(dat) <- c("ExtrusionRate","Additive","TearRes","Gloss","Opacity")

# For computerized approaches - make appropriate columns factors, but this can be really annoying so don't do it in dat
dat2 <- dat
dat2[,1] <- as.factor(dat2[,1])
dat2[,2] <- as.factor(dat2[,2])

# Part (a)
alpha <- 0.05

# Using built-in functions 
fit <- manova(cbind(TearRes,Gloss,Opacity) ~ ExtrusionRate * Additive, data=dat2) # if you use dat instead of dat2, residual df=8 instead of 6
(fit.summary <- summary(fit)) 
fit.summary$SS
anova(fit)

# Example 6.14 
rm(list=ls())

x.bar1 <- matrix(c(6.833,7.033,3.967,4.700),ncol=1)
x.bar2 <- matrix(c(6.633,7.000,4.000,4.533),ncol=1)
Spooled <- matrix(c(.606,.262,.066,.161,
                    .262,.637,.173,.143,
                    .066,.173,.810,.029,
                    .161,.143,.029,.306),ncol=4,byrow=TRUE)
C <- matrix(c(-1,1,0,0,
              0,-1,1,0,
              0,0,-1,1),ncol=4,byrow=TRUE)
n1 <- n2 <- 30
p <- 4
alpha <- 0.05

# Test for Parallel Profiles for Two Normal Populations
(T2 <- (t(x.bar1-x.bar2)%*%t(C))%*%solve((1/n1+1/n2)*C%*%Spooled%*%t(C))%*%(C%*%(x.bar1-x.bar2)))
(c2 <- ((n1+n2-2)*(p-1))/(n1+n2-p) * qf(1-alpha,p-1,n1+n2-p))
# Fail to Reject, conduct next test

# Test for Coincident Profiles
one <- matrix(rep(1,p),ncol=1)
(T2 <- (t(one)%*%(x.bar1-x.bar2))%*%
          solve((1/n1+1/n2)*
                  t(one)%*%Spooled%*%one)%*%
    (t(one)%*%(x.bar1-x.bar2)))
(c2 <- qf(1-alpha,1,n1+n2-2))
# Fail to Reject

