######## simulation for KSBE: denoising (numerical results)

library(Rfast)
library(RSpectra)
library(BiocNeighbors)
library(pcaPP)
###load data: torus

n=5000
u = runif(n, 0, 2*pi)
v = runif(n, 0, 2*pi)
data = cbind( (2+0.8*cos(u))*cos(v),
              (2+0.8*cos(u))*sin(v),
              0.8*sin(u))
data1=data
colfunc <- colorRampPalette(c("blue","forestgreen", "yellow"))
scatterplot3d(x=data[,1], y=data[,2], z=data[,3], angle=35, color = colfunc(n)[rank(-data[,1])],pch=20,
              xlab=NA, ylab=NA, zlab=NA)

u = runif(n, 0, pi)
v = runif(n, 0, 2*pi)
data = cbind( (2+0.8*cos(u))*cos(v),
              (2+0.8*cos(u))*sin(v),
              0.8*sin(u))
data2=data
colfunc <- colorRampPalette(c("blue","forestgreen", "yellow"))
scatterplot3d(x=data[,1], y=data[,2], z=data[,3], angle=35, color = colfunc(n)[rank(-data[,1])],pch=20,
              xlab=NA, ylab=NA, zlab=NA)

r=3
### knn dist


knn_dist <- function(X,Y){
  #knn graphs
  gx = findKNN(X, k = 50)$index
  gy = findKNN(Y, k = 50)$index
  
  out=c()
  for(i in 1:dim(gx)[1]){
    out[i]=length(intersect(gx[i,],gy[i,]))/length(union(gx[i,],gy[i,]))
    #out[i]= cor(as.matrix(dist(X[c(i,gx[i,]),]))[1,-1], as.matrix(dist(Y[c(i,gy[i,]),]))[1,-1])
  }
  return(mean(out))
}


### varying sample size

sig1 = 0.2
pp.v=seq(400,1000, length.out= 20)
error=matrix(ncol=5,nrow=100)
error.out=matrix(ncol=5,nrow=length(pp.v))
ali.scr=matrix(ncol=length(pp.v),nrow=100)
#p=floor(pp.v[pp])
p=800



for(pp in 1:length(pp.v)){
  n=floor(pp.v[pp])
  theta=n^(1/2)*0.2
  for(NN in 1:100){
    
    # #setting 1: basic
    data.X1=data1[sample(dim(data1)[1],n),]
    data.X2=data1[sample(dim(data1)[1],n),]
    data.Z=matrix(rnorm(n*p,0,0.4), ncol=p)
    data.Y1 = data.Z
    data.Y1[,1:r] = data.Y1[,1:r]+as.matrix(data.X1*theta)
    
    data.Z=matrix(rnorm(n*p,0,1), ncol=p)
    data.Y2 = data.Z
    data.Y2[,1:r] = data.Y2[,1:r]+as.matrix(data.X2*theta)
    data.Y2[,(r+1):(r+20)] = data.Y2[,(r+1):(r+20)]+matrix(runif(n*20,-8,8),ncol=20) #noise in Y2
    
    
    data.Y1 = scale(data.Y1, scale=FALSE)
    #data.Y10=t(data.Y10)
    data.Y1=t(data.Y1)
    data.Y2 = scale(data.Y2, scale=FALSE)
    #data.Y20=t(data.Y20)
    data.Y2=t(data.Y2)
    dist.mat = Dist(t(cbind(data.Y1, data.Y2)))
    cross.dist = as.matrix(dist.mat)[1:dim(data.Y1)[2], (dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])]
   
  
    
    
    
    #Prop-0.8
    K.mat = exp(-cross.dist^2/quantile(cross.dist,0.8)^2)  
    prop.svd = svds(K.mat, k=r+2)
    out.kpca = prop.svd$v[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
    
    
    #Prop-0.65
    K.mat = exp(-cross.dist^2/quantile(cross.dist,0.65)^2)  
    prop.svd = svds(K.mat, k=r+2)
    out.pca = prop.svd$v[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
   
    
    #Prop-0.35
    K.mat = exp(-cross.dist^2/quantile(cross.dist,0.35)^2)  
    prop.svd = svds(K.mat, k=r+2)
    out.lbdm = prop.svd$v[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
    
    #Prop-0.2
    K.mat = exp(-cross.dist^2/quantile(cross.dist,0.2)^2)  
    prop.svd = svds(K.mat, k=r+2)
    out.rl = prop.svd$v[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
   
    
    #Prop-0.5
    K.mat = exp(-cross.dist^2/quantile(cross.dist,0.5)^2)  
    prop.svd = svds(K.mat, k=r+2)
    out.prop = prop.svd$v[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)


  
      error[NN,]=c(knn_dist(data.X2,out.kpca),knn_dist(data.X2,out.pca),
                   knn_dist(data.X2,out.lbdm),knn_dist(data.X2,out.rl),knn_dist(data.X2,out.prop))

    
    print(NN)
  }
  print(pp)
  error.out[pp,] = colMeans(error,na.rm = T)
}


data_p = data.frame(rand.index = c(error.out),
                    method = factor(rep(c("prop-w-0.8", "prop-w-0.65", "prop-w-0.35", "prop-w-0.2",
                                          "prop-w-0.5"), 
                                        each=length(pp.v))),
                    sigma = rep(pp.v, times = 5))

data_p$method <- factor(data_p$method, levels = c("prop-w-0.2","prop-w-0.35","prop-w-0.5","prop-w-0.65", "prop-w-0.8"))

plt <- ggplot(data=data_p, aes(x=sigma, y=rand.index), group=factor(method)) +
  geom_line(aes(linetype=method,color=method))+
  geom_point(aes(shape=method,color=method) ,size=5) + ylab("concordance") + xlab("n")+
  theme(axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14,face="bold"), #change legend title font size
        legend.text = element_text(size=14,face="bold")) #600*500

plt 


######### check sensitivity to r

for(pp in 1:length(pp.v)){
  n=floor(pp.v[pp])
  theta=n^(1/2)*0.2
  for(NN in 1:100){
    
    # #setting 1: basic
    r=3
    data.X1=data1[sample(dim(data1)[1],n),]
    data.X2=data1[sample(dim(data1)[1],n),]
    data.Z=matrix(rnorm(n*p,0,0.4), ncol=p)
    data.Y1 = data.Z
    data.Y1[,1:r] = data.Y1[,1:r]+as.matrix(data.X1*theta)
    
    data.Z=matrix(rnorm(n*p,0,1), ncol=p)
    data.Y2 = data.Z
    data.Y2[,1:r] = data.Y2[,1:r]+as.matrix(data.X2*theta)
    data.Y2[,(r+1):(r+20)] = data.Y2[,(r+1):(r+20)]+matrix(runif(n*20,-8,8),ncol=20) #noise in Y2
    
    
    data.Y1 = scale(data.Y1, scale=FALSE)
    #data.Y10=t(data.Y10)
    data.Y1=t(data.Y1)
    data.Y2 = scale(data.Y2, scale=FALSE)
    #data.Y20=t(data.Y20)
    data.Y2=t(data.Y2)
    dist.mat = Dist(t(cbind(data.Y1, data.Y2)))
    cross.dist = as.matrix(dist.mat)[1:dim(data.Y1)[2], (dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])]
    
    
    
    K.mat = exp(-cross.dist^2/quantile(cross.dist,0.5)^2)  
    
    #Prop-r-2
    r=2
    prop.svd = svds(K.mat, k=r+2)
    out.kpca = prop.svd$v[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
    
    
    #Prop-r-3
    r=3
    prop.svd = svds(K.mat, k=r+2)
    out.pca = prop.svd$v[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
    
    
    #Prop-r-4
    r=4
    prop.svd = svds(K.mat, k=r+2)
    out.lbdm = prop.svd$v[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
    
    #Prop-r-5
    r=5
    prop.svd = svds(K.mat, k=r+2)
    out.rl = prop.svd$v[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
    
    
    #Prop-r-6
    r=6
    prop.svd = svds(K.mat, k=r+2)
    out.prop = prop.svd$v[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
    
    
    
    error[NN,]=c(knn_dist(data.X2,out.kpca),knn_dist(data.X2,out.pca),
                 knn_dist(data.X2,out.lbdm),knn_dist(data.X2,out.rl),knn_dist(data.X2,out.prop))
    
    
    print(NN)
  }
  print(pp)
  error.out[pp,] = colMeans(error,na.rm = T)
}


data_p = data.frame(rand.index = c(error.out),
                    method = factor(rep(c("prop-r-2", "prop-r-3", "prop-r-4", "prop-r-5",
                                          "prop-r-6"), 
                                        each=length(pp.v))),
                    sigma = rep(pp.v, times = 5))

data_p$method <- factor(data_p$method, levels = c("prop-r-2", "prop-r-3", "prop-r-4", "prop-r-5",
                                                  "prop-r-6"))

plt <- ggplot(data=data_p, aes(x=sigma, y=rand.index), group=factor(method)) +
  geom_line(aes(linetype=method,color=method))+
  geom_point(aes(shape=method,color=method) ,size=5) + ylab("concordance") + xlab("n")+
  theme(axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14,face="bold"), #change legend title font size
        legend.text = element_text(size=14,face="bold")) #600*500

plt 
