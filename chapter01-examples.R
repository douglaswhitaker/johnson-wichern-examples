# Examples for Chapter 1

# Example 1.1
X <- matrix( c(42,52,48,58,
               4,5,4,3), byrow = FALSE, ncol = 2 )
print(X)

# could also have done
X <- matrix(c(42,4, # this could all be one one line
              52,5,
              48,4,
              58,3),byrow = TRUE, ncol = 2)

x1 <- c(42,52,48,58)
x2 <- c(4,5,4,3)
Xdf <- data.frame(Variable1=x1,Variable2=x2)

print (X) 

# Example 1.2
# A bit about accessing matrix elements
X[,1] # all rows (blank), column 1
X[3,] # row 3, all columns (blank)
X[2,1] # x21 
X[1:3,2] # rows 1 through 3, column 2

# Calculating means
(xbar1 <- sum(X[,1])/length(X[,1]))
(xbar2 <- sum(X[,2])/length(X[,2]))

x.bar <- matrix(c(xbar1,xbar2),ncol=1) # for column or row matrix, byrow doesn't matter

(xbar1 <- mean(X[,1]))
(xbar2 <- mean(X[,2]))

mean(Xdf$Variable1)
mean(Xdf$Variable2)
colMeans(Xdf)

# Computing summary statistics element-wise
n <- length(x1) # all variables have the same length
s12 <- (1/n)*sum((x1-mean(x1))*(x2-mean(x2))) # same as s21
print(s12)
s11 <- (1/n)*sum((x1-mean(x1))*(x1-mean(x1)))
s22 <- (1/n)*sum((x2-mean(x2))*(x2-mean(x2)))
r12 <- s12/(sqrt(s11)*sqrt(s22)) # same as r21
print(r12)

# Creating matrices, because they are symmetric byrow = FALSE works fine
x.bar <- matrix(c(mean(x1),
                  mean(x2)),
                ncol=1)
print(X)
Sn <- matrix(c(s11,s12,
               s12,s22),
             ncol=2)
print(Sn)
R <- matrix(c(1,r12,
              r12,1),
            ncol=2) 
print(R)

cov(X) # different, because of 1/n vs. 1/(n-1) denominator
cor(X) # same as R


# Could also use loops
Sn2 <- matrix(rep(NA,2*2),ncol=2)
for(i in 1:2){
  for(j in 1:2){
    Sn2[i,j] <- (1/4)*sum((X[,i]-x.bar[i])*(X[,j]-x.bar[j]))
  }
}
print(Sn==Sn2)

# Example (like 1.3, from Exercise 1.27)
dat <- read.table("data/T1-11.dat")
head(dat)
names(dat) <- c("Size","Visitors")
head(dat)

plot(x=dat$Size,y=dat$Visitors)

plot(x=dat$Size,y=dat$Visitors,
     xlim=c(0,2500),
     ylim=c(0,10),
     main="Attendance and Size of National Parks",
     xlab="Size (acres)",
     ylab="Visitors (millions)",
     pch=19)


# install.packages("ggExtra")
library(ggplot2)
library(ggExtra)

# Note that ggplot2 has a very different syntax
# dotplot would be of limited use as sample size increases
ggMarginal(
  ggplot(dat,aes(Size,Visitors)) + 
    geom_point() + 
    theme_gray() + 
    ggtitle("Attendance and Size of National Parks") + 
    xlab("Size (acres)") + 
    ylab("Visitors (millions)"), 
  type="histogram", fill="steelblue",col="darkblue")

# Example 1.5
rm(list=ls())

dat <- read.table("data/T1-2.DAT")
names(dat) <- c("Density","MachineD","CrossD")
pairs(~Density+MachineD+CrossD,data=dat,main="Scatterplot Matrix")
pairs(~Density+MachineD+CrossD,data=dat,main="Scatterplot Matrix", pch=19) # note the ... in the ? description

# Lots of ways to do this!
library(car) # 1 of over 10,000 R packages...
scatterplot.matrix(~Density+MachineD+CrossD,data=dat,main="Scatterplot Matrix")

# Example 1.6
rm(list=ls())
dat <- read.table("data/T1-3.dat")
names(dat) <- c("Mass","SVL","HLS")

# Lots of ways 
library(scatterplot3d)
scatterplot3d(dat)
scatterplot3d(scale(dat))

library(rgl) # Interactive
plot3d(dat$Mass,dat$SVL,dat$HLS)

library(car) # Builds on rgl
scatter3d(dat$Mass,dat$SVL,dat$HLS)

# Example 1.7

sex <- c(1,0,1,1,0,1,0,1,0,1,0,1,0,
         0,0,0,1,0,0,0,1,1,0,1,1)
sex[which(sex==0)] <- 19 
scatterplot3d(dat,pch=sex)
plot3d(dat$Mass,dat$SVL,dat$HLS,pch=sex) # this library doesn't support pch... different tools for different purposes

# Example 1.8
rm(list=ls())
dat <- read.table("data/T1-2.DAT")
names(dat) <- c("Density","MachineD","CrossD")
pairs(~Density+MachineD+CrossD,data=dat,main="Scatterplot Matrix")

library(lattice)
library(hexbin)
splom(dat)
trellis.focus('panel',1,1)
clicked.points <- panel.link.splom(pch=19,col='green')

# Example 1.9
rm(list=ls())
dat <- read.table("data/T4-3.DAT")
library(rgl)
plot3d(dat$V2,dat$V3,dat$V4)
library(scatterplot3d)
scatterplot3d(dat$V2,dat$V3,dat$V4)
scatterplot3d(dat$V3,dat$V4,dat$V2)
scatterplot3d(dat$V4,dat$V2,dat$V3)

# Example 1.10
rm(list=ls())
dat <- read.table("data/T1-4.dat") # the filenames generally aren't case-sensitive, at least on Windows
names(dat) <- c("Wt2","Wt3","Wt4","Wt5","Lngth2","Lngth3","Lngth4","Lngth5") # no spaces in column names... using "." instead is okay (that's what they are automatically changed to when importing at least)

library(reshape)
dat2 <- melt(dat[,1:4])
dat2 <- cbind(bear=rep(1:(nrow(dat2)/4),4),dat2)

library(lattice)
xyplot(value ~ variable, group=bear,type="l",data=dat2)

dat3 <- cbind(dat2,year=rep(2:5,each=length(unique(dat2$bear))))
par(mfrow=c(2,4))
for (i in unique(dat2$bear)){
  plot(value ~ year, data=dat3[which(dat3$bear==i),], type="o",main=paste("Bear",i,"Weight"))
}
par(mfrow=c(1,1))

# For length
dat4 <- melt(dat[,5:8])
dat4 <- cbind(bear=rep(1:(nrow(dat4)/4),4),dat4)
dat4 <- cbind(dat4,year=rep(2:5,each=length(unique(dat4$bear))))
par(mfrow=c(2,4))
for (i in unique(dat2$bear)){
  plot(value ~ year, data=dat4[which(dat4$bear==i),], type="o",main=paste("Bear",i,"Length"))
}
par(mfrow=c(1,1))

# Example 1.11
rm(list=ls())
dat <- read.table("data/T12-4.DAT")

library(fmsb)
dat.tmp <- dat[,1:8]
dat.tmp <- rbind(floor(sapply(dat[,1:8],min)),dat.tmp)
dat.tmp <- rbind(ceiling(sapply(dat[,1:8],max)),dat.tmp)
radarchart(dat.tmp,maxmin=TRUE)
par(mfrow=c(2,3))
for (i in 1:5){
  radarchart(rbind(dat.tmp[1:2,],dat[i+2,1:8]),maxmin=FALSE)
}
par(mfrow=c(1,1))

# Example 1.12
library(aplpack)
faces(dat[,1:8])
