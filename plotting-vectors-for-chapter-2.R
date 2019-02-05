rm(list=ls())

library(rgl)
library(matlib)
library(MASS)

# Example 3.2 
X <- matrix(c(4,-1,3,
              1,3,5),
            ncol=2,byrow=FALSE)

# This method works well but there is no arrow/cone to show direction
open3d()
vectors3d(X[,1],col="blue",lwd=2)
vectors3d(X[,2],col="red",lwd=2)
axes3d() 


# Using extra code to add cones to the end to show direction
source("ch2-hw-helper-functions.R")

plot.vec(vec=X[,1],col="blue")
plot.vec(vec=X[,2],col="red")
axes3d()

snapshot3d("example-3.2-vectors.png") # save screenshot