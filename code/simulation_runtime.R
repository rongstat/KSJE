######## simulation for KSBE: denoising (numerical results)

library(Rfast)
library(RSpectra)
library(BiocNeighbors)

###load data: smelly face

n=18000
eyes.r = 0.1*sqrt(runif(n/4, 0, 1))
eyes.theta = runif(n/4,0, 2*pi)
eyes.sample = cbind(eyes.r*sin(eyes.theta),eyes.r*cos(eyes.theta))
eyes.sample[1:(n/8),]= eyes.sample[1:(n/8),] + c(0.25,0.25)
eyes.sample[-(1:(n/8)),1]= eyes.sample[-(1:(n/8)),1] - 0.25
eyes.sample[-(1:(n/8)),2]= eyes.sample[-(1:(n/8)),2] + 0.25
face.r = sqrt(runif(n/2, 0.9^2, 1))
face.theta = runif(n/2,0, 2*pi)
face.sample = cbind(face.r*sin(face.theta),face.r*cos(face.theta))
mouth.r = sqrt(runif(n/4, 0.45^2, 0.55^2))
mouth.theta = runif(n/4,0, pi)
mouth.sample = cbind(mouth.r*cos(mouth.theta), -mouth.r*sin(mouth.theta))
data=rbind(eyes.sample, face.sample, mouth.sample)
data=data*2
r=2

###load data: Setting 2 - smelly face - no mouth

n=18000
eyes.r = 0.1*sqrt(runif(n/4, 0, 1))
eyes.theta = runif(n/4,0, 2*pi)
eyes.sample = cbind(eyes.r*sin(eyes.theta),eyes.r*cos(eyes.theta))
eyes.sample[1:(n/8),]= eyes.sample[1:(n/8),] + c(0.25,0.25)
eyes.sample[-(1:(n/8)),1]= eyes.sample[-(1:(n/8)),1] - 0.25
eyes.sample[-(1:(n/8)),2]= eyes.sample[-(1:(n/8)),2] + 0.25
face.r = sqrt(runif(n/2, 0.9^2, 1))
face.theta = runif(n/2,0, 2*pi)
face.sample = cbind(face.r*sin(face.theta),face.r*cos(face.theta))
data2=rbind(eyes.sample, face.sample)
data2=data2*2
r=2

### knn dist

knn_dist <- function(X,Y){
  #knn graphs
  gx = findKmknn(X, k = 50)$index
  gy = findKmknn(Y, k = 50)$index
  
  out=c()
  for(i in 1:dim(gx)[1]){
    out[i]=length(intersect(gx[i,],gy[i,]))/length(union(gx[i,],gy[i,]))
    #out[i]= cor(as.matrix(dist(X[c(i,gx[i,]),]))[1,-1], as.matrix(dist(Y[c(i,gy[i,]),]))[1,-1])
  }
  return(mean(out))
}


###

sig1 = 1
time=matrix(ncol=8,nrow=100)
error=matrix(ncol=8,nrow=100)
error.out=matrix(ncol=8,nrow=4)
n=1500
p=500
pp.v=seq(800,12000,length.out=8)
for(pp in 1:50){
  n=pp.v[pp]
  #q=1
    
 
    #setting 4: similar but distinct structure
    data.X1=data[sample(dim(data)[1],n),]
    data.X2=data[sample(dim(data2)[1],n),]
    data.Z=matrix(rnorm(n*p,0,sig1), ncol=p)
    data.Y1 = data.Z
    theta=n^(2/5)
    data.Y1[,1:r] = data.Y1[,1:r]+as.matrix(data.X1*theta)
    
    sig2 = theta*2 #augmented structures
    data.Z=matrix(rnorm(n*p,0,sig1*3), ncol=p)
    data.Y2 = data.Z
    data.Y2[,1:r] = data.Y2[,1:r]+as.matrix(data.X2*theta)
    data.Y2[,-(1:r)] = data.Y2[,-(1:r)]+matrix(rnorm(n*(p-r),0,sig2),ncol=(p-r))
    
    
    
    
    #various kernel matrices
    data.Y1=t(data.Y1)
    data.Y2=t(data.Y2)
    

    #LBDM
    start_time <- Sys.time()
    dist.mat = Dist(t(cbind(data.Y1, data.Y2)))
    cross.dist = as.matrix(dist.mat)[1:dim(data.Y1)[2], (dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])]
    K.mat = exp(-cross.dist^2/quantile(dist.mat,0.5)^2)  
    alpha=1
    A = K.mat
    D1 = rowSums(K.mat)
    D2 = colSums(K.mat)
    A = t(t(D1^(-1/2) * A) * (D2^(-1/2)))
    lbdm.svd = svds(A, k=r+1)
    X.combined = rbind(lbdm.svd$u[,1:r+1], lbdm.svd$v[,1:r+1])
    V = c(D1,D2)^(-1/2) * X.combined %*% diag(lbdm.svd$d[1:r+1]^alpha)
    out.lbdm = V[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),]
    end_time <- Sys.time()
    #plot(out.lbdm)
    time[pp,6]=as.numeric(difftime(end_time, start_time, units = "secs"))
    
    #roseland
    start_time <- Sys.time()
    dist.mat = Dist(t(cbind(data.Y1, data.Y2)))
    cross.dist = as.matrix(dist.mat)[1:dim(data.Y1)[2], (dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])]
    K.mat = exp(-cross.dist^2/quantile(dist.mat,0.5)^2)  
    D.r = rowSums(t(K.mat) %*% (K.mat))
    r.svd.r = svds(((D.r)^(-1/2)) * t(K.mat), k=r+1)
    out.rl=(D.r)^(-1/2) * r.svd.r$u[,2:(r+1)] %*% diag(r.svd.r$d[2:(r+1)])^2
    end_time <- Sys.time()
    #plot(out.rl)
    time[pp,7]=as.numeric(difftime(end_time, start_time, units = "secs"))
    
    #Prop
    start_time <- Sys.time()
    dist.mat = Dist(t(cbind(data.Y1, data.Y2)))
    cross.dist = as.matrix(dist.mat)[1:dim(data.Y1)[2], (dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])]
    K.mat = exp(-cross.dist^2/quantile(dist.mat,0.5)^2)  
    prop.svd = svds(K.mat, k=r+2)
    out.prop = prop.svd$v[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
    end_time <- Sys.time()
    #plot(out.prop,xlab="",ylab="",col=info2)
    time[pp,8]=as.numeric(difftime(end_time, start_time, units = "secs"))
    
  print(pp)
}

