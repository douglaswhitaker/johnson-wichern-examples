# Example 4.2

library(plotrix) # needed for draw.ellipse()

# Eigenvectors, given in the example but you could get them any way you like
e1 <- matrix(c(1/sqrt(2),1/sqrt(2)),ncol=1)
e2 <- matrix(c(1/sqrt(2),-1/sqrt(2)),ncol=1)
eigenvectors <- cbind(e1,e2)

# Eigenvalues, given in the example but you could get them any way you like
lambda1 <- s11+s12 
lambda2 <- s11-s12

# angle of rotation for the ellipse
angle <- atan(eigenvectors[2,1]/eigenvectors[1,1]) # sohcahtoa, in radians (see deg=FALSE below)

mu <- c(0,0) # center of ellipse

# Given variance and covariance values
s11 <- 1
s12 <- 0.75

# Create an empty plot - draw.ellipse adds to an existing plot but won't create a plot itself
plot(0,pch='',ylab='',xlab='',xlim=c(-5,5),ylim=c(-5,5))

c50 <- qchisq(0.50,df=2) # c-value, determined by the specifics of the problem

lengths50 <- c(c50*sqrt(lambda1),
               c50*sqrt(lambda2)) # lengths of the axes

draw.ellipse(x=mu[1],y=mu[2], # Center
             a=lengths50[1],b=lengths50[2], # Axes
             angle=angle,deg=FALSE) # Rotation

# Another ellipse
c90 <- qchisq(0.90,df=2)
lengths90 <- c(c90*sqrt(s11+s12),c90*sqrt(s11-s12))
draw.ellipse(x=mu[1],y=mu[2],a=lengths90[1],b=lengths90[2],angle=angle,deg=FALSE)