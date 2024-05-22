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
  gx = findKmknn(X, k = 50)$index
  gy = findKmknn(Y, k = 50)$index
  
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
time=matrix(ncol=3,nrow=100)
time.out=matrix(ncol=3,nrow=length(pp.v))
error=matrix(ncol=7,nrow=100)
error.out=matrix(ncol=7,nrow=length(pp.v))
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
    K.mat = exp(-cross.dist^2/quantile(cross.dist,0.5)^2)  
    K.mat.rl = exp(-cross.dist^2/quantile(dist.mat,0.05)^2) #common manifold model & discussion
    
    K.mat.all = exp(-as.matrix(dist.mat)^2/quantile(dist.mat,0.5)^2)
    
    K.mat1=exp(-as.matrix(dist.mat)[1:dim(data.Y1)[2],1:dim(data.Y1)[2]]^2/quantile(dist.mat[1:dim(data.Y1)[2],1:dim(data.Y1)[2]],0.5)^2)
    K.mat2=exp(-as.matrix(dist.mat)[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])]^2/quantile(dist.mat[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])],0.5)^2)
    
    
    
    #par(mfrow = c(2, 4), mar = c(4,1,1,1) + 0.1)
    
    #PCA
    sepl.svd = svds(data.Y2,k=r+1)
    out.sepl2 = sepl.svd$v[,1:r] %*% diag(sepl.svd$d[1:r])^(1/2)
    sepl.svd = svds(data.Y1,k=r+1)
    out.sepl1 = sepl.svd$v[,1:r] %*% diag(sepl.svd$d[1:r])^(1/2)
    #plot(out.sepl1)
    
    
    #KPCA
    sep.svd = eigs(K.mat2,k=r+1)
    out.sep2 = sep.svd$vectors[,1:r+1] %*% diag(sep.svd$values[1:r+1])^(1/2)
    sep.svd = eigs(K.mat2,k=r+1)
    out.sep1 = sep.svd$vectors[,1:r+1] %*% diag(sep.svd$values[1:r+1])^(1/2)
    #plot(out.sep2)
    
    
    #joint PCA
    pca.svd = svds(cbind(data.Y1,data.Y2), k=r)
    X.combined =pca.svd$v[,1:r] %*% diag(pca.svd$d[1:r])^(1/2)
    out.pca2 = X.combined[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),]
    out.pca1 = X.combined[(1):(dim(data.Y1)[2]),]
    #plot(out.pca)
    
    
    #joint kPCA
    kpca.svd = eigs(K.mat.all, k=r+1)
    X.combined = kpca.svd$vectors[,1:r+1] %*% diag(kpca.svd$values[1:r+1])^(1/2)
    out.kpca2 = X.combined[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),]
    out.kpca1 = X.combined[(1):(dim(data.Y1)[2]),]
    #plot(out.kpca)
    
    
    #LBDM
    start_time <- Sys.time()
    alpha=1
    A = K.mat
    D1 = rowSums(K.mat)
    D2 = colSums(K.mat)
    A = diag(D1^(-1/2)) %*% A %*% diag(D2^(-1/2))
    lbdm.svd = svds(A, k=r+1)
    X.combined = rbind(lbdm.svd$u[,1:r+1], lbdm.svd$v[,1:r+1])
    V = diag(c(D1,D2)^(-1/2)) %*% X.combined %*% diag(lbdm.svd$d[1:r+1]^alpha)
    out.lbdm = V[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),]
    end_time <- Sys.time()
    #plot(out.lbdm)
    time[NN,1]=as.numeric(difftime(end_time, start_time, units = "secs"))
    
    #roseland
    start_time <- Sys.time()
    D.r = rowSums(t(K.mat.rl) %*% (K.mat.rl))
    r.svd.r = svds(diag((D.r)^(-1/2)) %*% t(K.mat.rl), k=r+1)
    out.rl=diag((D.r)^(-1/2)) %*% r.svd.r$u[,2:(r+1)] %*% diag(r.svd.r$d[2:(r+1)])^2
    end_time <- Sys.time()
    #plot(out.rl)
    time[NN,2]=as.numeric(difftime(end_time, start_time, units = "secs"))
    
    #Prop
    start_time <- Sys.time()
    prop.svd = svds(K.mat, k=r+2)
    out.prop = prop.svd$v[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
    end_time <- Sys.time()
    #plot(out.prop,xlab="",ylab="")
    time[NN,3]=as.numeric(difftime(end_time, start_time, units = "secs"))
    
    
    error[NN,]=c(knn_dist(data.X2,out.sepl2),knn_dist(data.X2,out.sep2),knn_dist(data.X2,out.pca2),knn_dist(data.X2,out.kpca2),
                 knn_dist(data.X2,out.lbdm),knn_dist(data.X2,out.rl),knn_dist(data.X2,out.prop))
    print(NN)
  }
  print(pp)
  time.out[pp,]=colMeans(time,na.rm = T)
  error.out[pp,] = colMeans(error,na.rm = T)
}


data_p = data.frame(rand.index = c(error.out),
                    method = factor(rep(c("pca", "kpca", "j-pca", "j-kpca",
                                          "lbdm","rl","prop"), 
                                        each=length(pp.v))),
                    sigma = rep(pp.v, times = 7))

data_p$method <- factor(data_p$method, levels = c("prop","lbdm","rl","j-pca", "j-kpca", "kpca", "pca"))

plt <- ggplot(data=data_p, aes(x=sigma, y=rand.index), group=factor(method)) +
  geom_line(aes(linetype=method,color=method))+
  geom_point(aes(shape=method,color=method) ,size=5) + ylab("Jaccard index") + xlab("n")+
  theme(axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14,face="bold"), #change legend title font size
        legend.text = element_text(size=14,face="bold")) #600*500

plt 

plot_data <- ggplot_build(plt)
color_codes <- plot_data$data[[1]]$colour
unique_colors <- unique(color_codes)
unique_colors


data_p = data.frame(time = c(time.out),
                    method = factor(rep(c("lbdm","rl","prop"), 
                                        each=length(pp.v))),
                    dim = rep(pp.v, times = 3))
data_p$method <- factor(data_p$method, levels = c("prop","lbdm","rl"))


method_colors <- c("prop" = "#F8766D", "lbdm" = "#C49A00", "rl" = "#53B400", 
                   "j-pca" = "purple", "j-kpca" = "blue", "kpca" = "pink", "pca" = "brown")
method_linetypes <- c("prop" = "solid", "lbdm" = "dashed", "rl" = "dotted", 
                      "j-pca" = "dotdash", "j-kpca" = "longdash", "kpca" = "twodash", "pca" = "blank")


plt2 <- ggplot(data=data_p, aes(x=dim, y=time, group=factor(method))) + 
  geom_line(aes(linetype=method, color=method)) +
  geom_point(aes(shape=method, color=method), size=5) +
  ylab("running time") + xlab("n") +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  theme(axis.text=element_text(size=15, face="bold"),
        axis.title=element_text(size=15, face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=14, face="bold"))

plt2
