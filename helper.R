# Helper functions for Chapter 2
projAonB <- function(A,B){
  return(fractions(det((t(A)%*%B)/(t(B)%*%B))*B))
}

l.vec <- function(v){
  return(sqrt(t(v)%*%v))
}

r2d <- function(r){return(r*180/pi)}

plot.vec <- function(vec,col="black"){
  v <- rbind(c(0,0,0),vec)
  segments3d(v,col=col)
  cone3d(base=v[2,]-(v[1,]+v[2,]/6),
         rad=0.5,tip=v[2,],col=col,front="lines",back="lines")
}

# from rgl demo file lollipop3d.R
cone3d <- function(base,tip,rad,n=30,...) {
  degvec <- seq(0,2*pi,length=n)
  ax <- tip-base
  ## what do if ax[1]==0?
  if (ax[1]!=0) {
    p1 <- c(-ax[2]/ax[1],1,0)
    p1 <- p1/sqrt(sum(p1^2))
    if (p1[1]!=0) {
      p2 <- c(-p1[2]/p1[1],1,0)
      p2[3] <- -sum(p2*ax)
      p2 <- p2/sqrt(sum(p2^2))
    } else {
      p2 <- c(0,0,1)
    }
  } else if (ax[2]!=0) {
    p1 <- c(0,-ax[3]/ax[2],1)
    p1 <- p1/sqrt(sum(p1^2))
    if (p1[1]!=0) {
      p2 <- c(0,-p1[3]/p1[2],1)
      p2[3] <- -sum(p2*ax)
      p2 <- p2/sqrt(sum(p2^2))
    } else {
      p2 <- c(1,0,0)
    }
  } else {
    p1 <- c(0,1,0); p2 <- c(1,0,0)
  }
  ecoord2 <- function(theta) {
    base+rad*(cos(theta)*p1+sin(theta)*p2)
  }
  for (i in 1:(n-1)) {
    li <- ecoord2(degvec[i])
    lj <- ecoord2(degvec[i+1])
    triangles3d(c(li[1],lj[1],tip[1]),c(li[2],lj[2],tip[2]),c(li[3],lj[3],tip[3]),...)
  }   
}