# Examples

library(MASS) # Venables & Ripley, (2002), Modern Applied Statistics with S

# Example 4.1
set.seed(1)
mu <- c(0,0)
Sigmaa <- matrix(c(1,0,0,1),ncol=2)
Sigmab <- matrix(c(1,.75,.75,1),ncol=2)

dat1 <- mvrnorm(5000, mu = mu, Sigma = Sigmaa)
dat2 <- mvrnorm(5000, mu = mu, Sigma = Sigmab)

kde1 <- kde2d(dat1[,1],dat1[,2],n=100)
kde2 <- kde2d(dat2[,1],dat2[,2],n=100)

persp(kde1,phi=45,theta=30,shade=0.3,border=NA)
filled.contour(kde1)
persp(kde2,phi=45,theta=30,shade=0.3,border=NA)
filled.contour(kde2)

library(rgl)
colors1 <- heat.colors(length(kde1$z))[rank(kde1$z)]
persp3d(kde1, col=colors1)
colors2 <- heat.colors(length(kde2$z))[rank(kde2$z)]
persp3d(kde2, col=colors2)


# Example 4.2
e1 <- matrix(c(1/sqrt(2),1/sqrt(2)),ncol=1)
e2 <- matrix(c(1/sqrt(2),-1/sqrt(2)),ncol=1)
eigenvectors <- cbind(e1,e2)

angle <- atan(eigenvectors[2,1]/eigenvectors[1,1]) # sohcahtoa

library(plotrix)
mu <- c(0,0)

s11 <- 1
s12 <- 0.75

plot(0,pch='',ylab='',xlab='',xlim=c(-5,5),ylim=c(-5,5))
c50 <- qchisq(0.50,df=2)
lengths50 <- c(c50*sqrt(s11+s12),c50*sqrt(s11-s12))
draw.ellipse(x=mu[1],y=mu[2],a=lengths50[1],b=lengths50[2],angle=angle,deg=FALSE)

c90 <- qchisq(0.90,df=2)
lengths90 <- c(c90*sqrt(s11+s12),c90*sqrt(s11-s12))
draw.ellipse(x=mu[1],y=mu[2],a=lengths90[1],b=lengths90[2],angle=angle,deg=FALSE)

# Example 4.9
xo <- c(-1.00,-.10,.16,.41,.62,.80,1.26,1.54,1.71,2.30)
n <- length(xo)
problevels <- ((1:n)-0.5)/n
q <- qnorm(problevels)
plot(q,xo,pch=19)


# Example 4.10
dat <- read.table("data/T4-1.DAT")
names(dat) <- c("Radiation")

xo <- sort(dat$Radiation)
n <- length(xo)
problevels <- ((1:n)-0.5)/n
q <- qnorm(problevels)
plot(q,xo,pch=19)

qqnorm(dat$Radiation) # built-in functions


# Example 4.11
xo <- c(-1.00,-.10,.16,.41,.62,.80,1.26,1.54,1.71,2.30)
n <- length(xo)
problevels <- ((1:n)-0.5)/n
q <- qnorm(problevels)
q.bar <- mean(q)
x.bar <- mean(xo)
critval <- 0.9351 # from table 4.2 p. 181
# reject normality if rQ is below this critical value
(rq1 <- sum((xo-x.bar)*(q-q.bar))/(sqrt(sum((xo-x.bar)^2))*sqrt(sum((q-q.bar)^2))))


# Example 4.12
x1 <- matrix(c(108.28,152.36,95.04,65.45,62.97,263.99,265.19,285.06,92.01,165.68),ncol=1)
x2 <- matrix(c(17.05,16.59,10.91,14.14,9.52,25.33,18.54,15.73,8.10,11.13),ncol=1)
X <- cbind(x1,x2)

n <- nrow(X)
p <- ncol(X)
(x.bar <- t(matrix(1,ncol=n) %*% X)/n)
(D <- X - matrix(1,nrow=n) %*% t(x.bar))
(S <- (n-1)^(-1) * t(D)%*%D)
Sinv <- solve(S)

d <- c()
for (i in 1:n){
  d[i] <- (X[i,]-t(x.bar))%*%Sinv%*%t(X[i,]-t(x.bar))
}

print(d)
(cutoff <- qchisq(0.50,2))
which(d<cutoff)
sum(d<cutoff)/n # 40% of data fall within 50% contour




# Example 4.13
chisq.quantiles <- qchisq(((1:n)-0.5)/n,df=2)
plot(y=sort(d),x=chisq.quantiles,pch=19)

# Example 4.14 
dat <- read.table("data/T4-3.DAT")
names(dat) <- c("x1","x2","x3","x4","d2")
n <- nrow(dat)
chisq.quantiles <- qchisq(((1:n)-0.5)/n,df=4)
plot(y=sort(dat$d2),x=chisq.quantiles,pch=19)

# Example 4.15
library(lattice)
library(hexbin)
splom(dat[,-5])
trellis.focus('panel',1,1)
clicked.points <- panel.link.splom(pch=19,col='green')

# Example 4.16
dat <- read.table("data/T4-1.DAT")
names(dat) <- c("Radiation")
n <- nrow(dat)

l.lambda <- function(lambda,n,x){
  if(lambda==0){
    x.lambda <- log(x)
  }
  else {
    x.lambda <- (x^lambda - 1)/lambda
  }
  return(
    (-n/2)*log(
      sum(
        (x.lambda-mean(x.lambda))^2)
      )+
      (lambda-1)*sum(log(x)))
}

lambdas <- seq(from=-1,to=1.50,by=.01)
l.lambdas <- c()
for (i in 1:length(lambdas)){
  l.lambdas[i] <- l.lambda(lambdas[i],n,dat$Radiation)
}
plot(lambdas,l.lambdas)
lambdas[which(l.lambdas==max(l.lambdas))] # this is the transformation we should use

qqnorm(dat$Radiation)
qqnorm((dat$Radiation^0.3 - 1)/0.3)

# Using optimization function nlm

l.lambda2 <- function(lambda,n,x){
  return(-l.lambda(lambda=lambda,n=n,x=x)) # nlm only minimizes, so we turn this into a minimization problem
}

nlm(l.lambda2,p=-2,n=n,x=dat$Radiation)$estimate

# Example 4.17
dat2 <- read.table(file="T4-5.DAT")
dat1 <- dat; rm(dat)
dat <- data.frame(closed=dat1$Radiation,open=dat2$V1)

l.lambda.biv <- function(lambdas,n,x){
  if(lambdas[1]==0){
    x1.lambda <- log(x[,1])
  }
  else {
    x1.lambda <- (x[,1]^lambdas[1] - 1)/lambdas[1]
  }
  
  if(lambdas[2]==0){
    x2.lambda <- log(x[,2])
  }
  else {
    x2.lambda <- (x[,2]^lambdas[2] - 1)/lambdas[2]
  }
  
  Slambda <- cov(as.matrix(cbind(matrix(x1.lambda,ncol=1),matrix(x2.lambda,ncol=1))))

  
  return(-(n/2)*log(det(Slambda))+(lambdas[1]-1)*sum(log(x[,1]))+(lambdas[2]-1)*sum(log(x[,2])))
}

lambdas1 <- lambdas2 <- seq(from=0,to=.5,by=0.02) # choice of step is important, 0.01 results in slightly different values
l.lambdas.biv <- matrix(NA,nrow=length(lambdas1),ncol=length(lambdas2))
for (i in 1:length(lambdas1)){
  for (j in 1:length(lambdas2)){
    l.lambdas.biv[i,j] <- l.lambda.biv(lambdas=c(lambdas1[i],lambdas2[j]),n=n,x=dat)
  }
}

filled.contour(lambdas1,lambdas2,l.lambdas.biv)
maximizing.index <- which(l.lambdas.biv==max(l.lambdas.biv), arr.ind=TRUE)
lambdas1[maximizing.index[1]]
lambdas2[maximizing.index[2]]
