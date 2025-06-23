############# negative control for duo-landmark

library(Rfast)
library(RSpectra)
library(BiocNeighbors)
library(pcaPP)
library(tdaunif)

ali.scr=c()
for(nn in 1:100){
data.Y1 <- sample_klein_tube(3000, sd = .00)
#pairs(data.Y1, asp = 1, pch = 19, cex = .5)
dim(data.Y1)


n=2000
data.Y2=cbind( runif(n,-1,1), rep(0,n), rep(0,n), rep(0,n))



#various kernel matrices
data.Y1=t(data.Y1)
data.Y2=t(data.Y2)
dist.mat = Dist(t(cbind(data.Y1, data.Y2)))
cross.dist = as.matrix(dist.mat)[1:dim(data.Y1)[2], (dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])]
K.mat = exp(-cross.dist^2/quantile(dist.mat,0.5)^2)  

K.mat.all = exp(-as.matrix(dist.mat)^2/quantile(dist.mat,0.5)^2)


#joint PCA
par(mfrow=c(1,3))

#joint kPCA
r=4
kpca.svd = eigs(K.mat.all, k=20)
X.combined = kpca.svd$vectors[,1:r+1] %*% diag(kpca.svd$values[1:r+1])^(1/2)
out.kpca2 = X.combined[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),]
out.kpca1 = X.combined[(1):(dim(data.Y1)[2]),]

cl = c(rep(0,dim(K.mat)[1]), rep(1,dim(K.mat)[2]))

#plot(X.combined[,1:2], col=factor(cl),xlab="KEF1",ylab="KEF2" )
#plot(X.combined[,2:3], col=factor(cl),xlab="KEF2",ylab="KEF3" )
#plot(X.combined[,3:4], col=factor(cl),xlab="KEF3",ylab="KEF4" )

library(bluster)
#boxplot(neighborPurity(X.combined,cl)[,1])
ali.scr[nn]=median(neighborPurity(X.combined,cl)[,1])

#prop
prop.svd = svds(K.mat, k=20)
out.prop2 = prop.svd$v[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
out.prop1 = prop.svd$u[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
X.combined = rbind(out.prop1,out.prop2)


#plot(X.combined[,1:2], col=factor(cl),xlab="KEF1",ylab="KEF2" )
#plot(X.combined[,2:3], col=factor(cl),xlab="KEF2",ylab="KEF3" )
#plot(X.combined[,3:4], col=factor(cl),xlab="KEF3",ylab="KEF4" )


}

summary(ali.scr)



######## signal & noise

ali.scr=c()
for(nn in 1:100){
  n=3000
  u = runif(n, 0, 2*pi)
  v = runif(n, 0, 2*pi)
  data = cbind( (2+0.8*cos(u))*cos(v),
                (2+0.8*cos(u))*sin(v),
                0.8*sin(u))
  data.Y1=data
  
  
  n=2000
  data.Y2=cbind( rnorm(n,0,1), rnorm(n,0,1), rnorm(n,0,1))/2
  
  #various kernel matrices
  data.Y1=t(data.Y1)
  data.Y2=t(data.Y2)
  dist.mat = Dist(t(cbind(data.Y1, data.Y2)))
  cross.dist = as.matrix(dist.mat)[1:dim(data.Y1)[2], (dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])]
  K.mat = exp(-cross.dist^2/quantile(dist.mat,0.5)^2)  
  
  K.mat.all = exp(-as.matrix(dist.mat)^2/quantile(dist.mat,0.5)^2)
  
  
  #joint PCA
  par(mfrow=c(1,3))
  
  #joint kPCA
  r=4
  kpca.svd = eigs(K.mat.all, k=20)
  X.combined = kpca.svd$vectors[,1:r+1] %*% diag(kpca.svd$values[1:r+1])^(1/2)
  out.kpca2 = X.combined[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),]
  out.kpca1 = X.combined[(1):(dim(data.Y1)[2]),]
  
  cl = c(rep(0,dim(K.mat)[1]), rep(1,dim(K.mat)[2]))
  
  #plot(X.combined[,1:2], col=factor(cl),xlab="KEF1",ylab="KEF2" )
  #plot(X.combined[,2:3], col=factor(cl),xlab="KEF2",ylab="KEF3" )
  #plot(X.combined[,3:4], col=factor(cl),xlab="KEF3",ylab="KEF4" )
  
  library(bluster)
  #boxplot(neighborPurity(X.combined,cl)[,1])
  ali.scr[nn]=median(neighborPurity(X.combined,cl)[,1])
  
  #prop
  prop.svd = svds(K.mat, k=20)
  out.prop2 = prop.svd$v[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
  out.prop1 = prop.svd$u[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
  X.combined = rbind(out.prop1,out.prop2)
  
  
  #plot(X.combined[,1:2], col=factor(cl),xlab="KEF1",ylab="KEF2" )
  #plot(X.combined[,2:3], col=factor(cl),xlab="KEF2",ylab="KEF3" )
  #plot(X.combined[,3:4], col=factor(cl),xlab="KEF3",ylab="KEF4" )
  

}

summary(ali.scr)


######## noise & noise
 
ali.scr=c()
for(nn in 1:100){
  n=3000
  data.Y1=cbind( rnorm(n,0,1), rnorm(n,0,1), rnorm(n,0,1))/2
  
  
  n=2000
  data.Y2=cbind( rnorm(n,0,1), rnorm(n,0,1), rnorm(n,0,1))/2
  
  #various kernel matrices
  data.Y1=t(data.Y1)
  data.Y2=t(data.Y2)
  dist.mat = Dist(t(cbind(data.Y1, data.Y2)))
  cross.dist = as.matrix(dist.mat)[1:dim(data.Y1)[2], (dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])]
  K.mat = exp(-cross.dist^2/quantile(dist.mat,0.5)^2)  
  
  K.mat.all = exp(-as.matrix(dist.mat)^2/quantile(dist.mat,0.5)^2)
  
  
  #joint PCA
  par(mfrow=c(1,3))
  
  #joint kPCA
  r=4
  kpca.svd = eigs(K.mat.all, k=20)
  X.combined = kpca.svd$vectors[,1:r+1] %*% diag(kpca.svd$values[1:r+1])^(1/2)
  out.kpca2 = X.combined[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),]
  out.kpca1 = X.combined[(1):(dim(data.Y1)[2]),]
  
  cl = c(rep(0,dim(K.mat)[1]), rep(1,dim(K.mat)[2]))
  
  #plot(X.combined[,1:2], col=factor(cl),xlab="KEF1",ylab="KEF2" )
  #plot(X.combined[,2:3], col=factor(cl),xlab="KEF2",ylab="KEF3" )
  #plot(X.combined[,3:4], col=factor(cl),xlab="KEF3",ylab="KEF4" )
  
  library(bluster)
  #boxplot(neighborPurity(X.combined,cl)[,1])
  ali.scr[nn]=median(neighborPurity(X.combined,cl)[,1])
  
  #prop
  prop.svd = svds(K.mat, k=20)
  out.prop2 = prop.svd$v[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
  out.prop1 = prop.svd$u[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
  X.combined = rbind(out.prop1,out.prop2)
  
  
  #plot(X.combined[,1:2], col=factor(cl),xlab="KEF1",ylab="KEF2" )
  #plot(X.combined[,2:3], col=factor(cl),xlab="KEF2",ylab="KEF3" )
  #plot(X.combined[,3:4], col=factor(cl),xlab="KEF3",ylab="KEF4" )
  
  
}

boxplot(ali.scr)

hist(ali.scr, xlab="median KNN purity")
