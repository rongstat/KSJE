######## simulation for KSBE: denoising (numerical results)

library(Rfast)
library(RSpectra)
library(BiocNeighbors)
library(fossil)
library(ggplot2)
library(cluster)
library(uwot)
library(bluster)
###
r=6
c.prop = c(1,1,1,1,1,1)
c.prop = c.prop/sum(c.prop)#class proportion
rho=5

#set mean vectors
mu=rep(0,r)
mu1=mu
mu1[1]=rho
mu2=mu
mu2[2]=rho
mu3=mu
mu3[3]=rho
mu4=mu
mu4[4]=rho
mu5=mu
mu5[5]=rho
mu6=mu
mu6[6]=rho
sigma=diag(rep(1,r))


error=matrix(ncol=7,nrow=100)
error.out=matrix(ncol=7,nrow=20)
error2=matrix(ncol=7,nrow=100)
ali.scr=matrix(ncol=20,nrow=100)
error.out2=matrix(ncol=7,nrow=20)
s.ind1 = matrix(ncol=7,nrow=100)
s.out1 = matrix(ncol=7,nrow=20)
s.ind2 = matrix(ncol=7,nrow=100)
s.out2 = matrix(ncol=7,nrow=20)
sig.v=seq(6,9, length.out= 20)

for(pp in 1:length(sig.v)){
  theta=3
  n0=600
  n=sum(round(n0*c.prop))
  p=800
  #epsilon = (log(n)/n)^(1/(4*r+5))
  for(NN in 1:100){
    
    
    #setting 1:
    #generate 6 clusters
    data1=rmvnorm(round(n0*c.prop[1]),mu1,sigma)
    data2=rmvnorm(round(n0*c.prop[2]),mu2,sigma)
    data3=rmvnorm(round(n0*c.prop[3]),mu3,sigma)
    data4=rmvnorm(round(n0*c.prop[4]),mu4,sigma)
    data5=rmvnorm(round(n0*c.prop[5]),mu5,sigma)
    data6=rmvnorm(round(n0*c.prop[6]),mu6,sigma)
    # data set X
    data.X1=rbind(data1,data2,data3,data4,data5,data6)
    data.Z1=rmvnorm(n,rep(0,p),diag(rep(0.5,p)))
    data.Y1 = data.Z1
    data.Y1[,1:6] = data.Y1[,1:6]+data.X1*theta
    info1=rep(1:6, times = round(n0*c.prop))

    data1=rmvnorm(round(n0*c.prop[1]),mu1,sigma)
    data2=rmvnorm(round(n0*c.prop[2]),mu2,sigma)
    data3=rmvnorm(round(n0*c.prop[3]),mu3,sigma)
    data4=rmvnorm(round(n0*c.prop[4]),mu4,sigma)
    data5=rmvnorm(round(n0*c.prop[5]),mu5,sigma)
    data6=rmvnorm(round(n0*c.prop[6]),mu6,sigma)
    # data set X
    data.X2=rbind(data1,data2,data3,data4,data5,data6)
    data.Z2=rmvnorm(n,rep(0,p),diag(rep(1,p)))
    data.Y2 = data.Z2
    data.Y2[,1:6] = data.Y2[,1:6]+data.X2*theta
    data.Y2[,(r):(r+20)] = data.Y2[,(r):(r+20)]+matrix(runif(n*21,-sig.v[pp]*3,sig.v[pp]),ncol=21) #noise in Y2
    info2=rep(1:6, times =  round(n0*c.prop))

    # 
    #setting 2: distinct cluster structure
    # #generate 5 clusters
    # data1=rmvnorm(n*0.25,mu3,sigma)
    # data2=rmvnorm(n*0.25,mu4,sigma)
    # data4=rmvnorm(n*0.25,mu5,sigma)
    # data6=rmvnorm(n*0.25,mu6,sigma)
    # # data set X
    # data.X1=rbind(data1,data2,data4,data6)
    # data.Z1=rmvnorm(n,rep(0,p),diag(rep(0.5,p)))
    # data.Y1 = data.Z1
    # data.Y1[,1:6] = data.Y1[,1:6]+data.X1*theta
    # info1 = rep(1:4, times=c(150,150,150,150))
    # #plot(umap(data.Y1),col=info1)
    # 
    # data1=rmvnorm(round(n0*c.prop[1]),mu1,sigma)
    # data2=rmvnorm(round(n0*c.prop[2]),mu2,sigma)
    # data3=rmvnorm(round(n0*c.prop[3]),mu3,sigma)
    # data4=rmvnorm(round(n0*c.prop[4]),mu4,sigma)
    # data5=rmvnorm(round(n0*c.prop[5]),mu5,sigma)
    # data6=rmvnorm(round(n0*c.prop[6]),mu6,sigma)
    # # data set X
    # data.X2=rbind(data1,data2,data3,data4,data5,data6)
    # data.Z2=rmvnorm(n,rep(0,p),diag(rep(8,p)))
    # data.Y2 = data.Z2
    # data.Y2[,1:6] = data.Y2[,1:6]+data.X2*theta
    # data.Y2[,(r):(r+20)] = data.Y2[,(r):(r+20)]+matrix(runif(n*21,-sig.v[pp]*3,sig.v[pp]),ncol=21) #noise in Y2
    # info2=rep(1:6, times =  round(n0*c.prop))



    #various kernel matrices
    data.Y1 = scale(data.Y1, scale=FALSE)
    data.Y1=t(data.Y1)
    data.Y2 = scale(data.Y2, scale=FALSE)
    data.Y2=t(data.Y2)
    dist.mat = Dist(t(cbind(data.Y1, data.Y2)))
    cross.dist = as.matrix(dist.mat)[1:dim(data.Y1)[2], (dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])]
    K.mat = exp(-cross.dist^2/quantile(cross.dist,0.5)^2)  
    K.mat.rl = exp(-cross.dist^2/quantile(dist.mat,0.05)^2) #common manifold model & discussion
    
    K.mat.all = exp(-as.matrix(dist.mat)^2/quantile(dist.mat,0.5)^2)
    
    K.mat1=exp(-as.matrix(dist.mat)[1:dim(data.Y1)[2],1:dim(data.Y1)[2]]^2/quantile(dist.mat[1:dim(data.Y1)[2],1:dim(data.Y1)[2]],0.5)^2)
    K.mat2=exp(-as.matrix(dist.mat)[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])]^2/quantile(dist.mat[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])],0.5)^2)
    
    
    
    #PCA
    sepl.svd = svds(data.Y2,k=r+1)
    out.sepl = sepl.svd$v[,1:r] %*% diag(sepl.svd$d[1:r])^(1/2)
    cls2=kmeans(out.sepl,centers=6)$cluster
    s.ind2[NN,1] = mean(silhouette( info2, dist(out.sepl[1:n,]))[,3])
    sepl.svd = svds(data.Y1,k=r+1)
    out.sepl = sepl.svd$v[,1:r] %*% diag(sepl.svd$d[1:r])^(1/2)
    cls1=kmeans(out.sepl,centers=6)$cluster
    s.ind1[NN,1] = mean(silhouette( info1, dist(out.sepl))[,3])
    error[NN,1] = (rand.index(cls2[1:n],info2)+rand.index(cls1,info1))/2
    error2[NN,1] = rand.index(cls2[1:n],info2)
    
    
    #KPCA
    sep.svd = eigs(K.mat2,k=r+1)
    out.sep = sep.svd$vectors[,1:r+1] %*% diag(sep.svd$values[1:r+1])^(1/2)
    cls2=kmeans(out.sep,centers=6)$cluster
    s.ind2[NN,2] = mean(silhouette(info2, dist(out.sep[1:n,]))[,3])
    #plot(umap(out.sep,spread=2),col=info2)
    sep.svd = eigs(K.mat1,k=r+1)
    out.sep = sep.svd$vectors[,1:r+1] %*% diag(sep.svd$values[1:r+1])^(1/2)
    cls1=kmeans(out.sep,centers=6)$cluster
    s.ind1[NN,2] = mean(silhouette( info1, dist(out.sep))[,3])
    error[NN,2] =(rand.index(cls2[1:n],info2)+rand.index(cls1,info1))/2
    error2[NN,2] = rand.index(cls2[1:n],info2)
    
    #joint PCA
    pca.svd = svds(cbind(data.Y1,data.Y2), k=r)
    X.combined =pca.svd$v[,1:r] %*% diag(pca.svd$d[1:r])^(1/2)
    out.pca = X.combined[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),]
    cls2=kmeans(out.pca,centers=6)$cluster
    s.ind2[NN,3] = mean(silhouette( info2, dist(out.pca[1:n,]))[,3])
    #plot(umap(out.pca,spread=2),col=info2)
    out.pca = X.combined[1:(dim(data.Y1)[2]),]
    cls1=kmeans(out.pca,centers=6)$cluster
    s.ind1[NN,3] = mean(silhouette( info1, dist(out.pca))[,3])
    error[NN,3] = (rand.index(cls2[1:n],info2)+rand.index(cls1,info1))/2
    error2[NN,3] = rand.index(cls2[1:n],info2)
    
    #joint kPCA
    kpca.svd = eigs(K.mat.all, k=r+1)
    X.combined = kpca.svd$vectors[,1:r+1] %*% diag(kpca.svd$values[1:r+1])^(1/2)
    out.kpca = X.combined[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),]
    cls2=kmeans(out.kpca,centers=6)$cluster
    s.ind2[NN,4] = mean(silhouette( info2, dist(out.kpca[1:n,]))[,3])
    #plot(umap(out.kpca,spread=2),col=info2)
    out.kpca = X.combined[1:(dim(data.Y1)[2]),]
    cls1=kmeans(out.kpca,centers=6)$cluster
    s.ind1[NN,4] = mean(silhouette( info1, dist(out.kpca))[,3])
    error[NN,4] = (rand.index(cls2[1:n],info2)+rand.index(cls1,info1))/2
    error2[NN,4] = rand.index(cls2[1:n],info2)
    
    cl = c(rep(0,dim(K.mat)[1]), rep(1,dim(K.mat)[2]))
    #boxplot(neighborPurity(X.combined,cl)[,1])
    #alignability screening step here
    ali.scr[NN,pp] = median(neighborPurity(X.combined,cl)[,1])
    
    #LBDM
    alpha=1
    A = K.mat.rl
    D1 = rowSums(K.mat.rl)
    D2 = colSums(K.mat.rl)
    A = diag(D1^(-1/2)) %*% A %*% diag(D2^(-1/2))
    lbdm.svd = svds(A, k=r+1)
    X.combined = rbind(lbdm.svd$u[,1:r+1], lbdm.svd$v[,1:r+1])
    V = diag(c(D1,D2)^(-1/2)) %*% X.combined %*% diag(lbdm.svd$d[1:r+1]^alpha)
    out.lbdm = V[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),]
    cls2=kmeans(out.lbdm,centers=6)$cluster
    s.ind2[NN,6] = mean(silhouette( info2, dist(out.lbdm[1:n,]))[,3])
    #plot(umap(out.lbdm),col=info2)
    out.lbdm = V[1:(dim(data.Y1)[2]),]
    cls1=kmeans(out.lbdm,centers=6)$cluster
    s.ind1[NN,6] = mean(silhouette( info1, dist(out.lbdm))[,3])
    error[NN,6] = (rand.index(cls2[1:n],info2)+rand.index(cls1,info1))/2
    error2[NN,6] = rand.index(cls2[1:n],info2)
    
    #roseland
    D.r = rowSums(t(K.mat.rl) %*% (K.mat.rl))
    r.svd.r = svds(diag((D.r)^(-1/2)) %*% t(K.mat.rl), k=r+1)
    out.rl=diag((D.r)^(-1/2)) %*% r.svd.r$u[,2:(r+1)] %*% diag(r.svd.r$d[2:(r+1)])^2
    cls2=kmeans(out.rl,centers=6)$cluster
    s.ind2[NN,7] = mean(silhouette( info2, dist(out.rl[1:n,]))[,3])
    #plot(umap(out.rl),col=info2)
    D.r = rowSums((K.mat.rl) %*% t(K.mat.rl))
    r.svd.r = svds(diag((D.r)^(-1/2)) %*% (K.mat.rl), k=r+1)
    out.rl=diag((D.r)^(-1/2)) %*% r.svd.r$u[,2:(r+1)] %*% diag(r.svd.r$d[2:(r+1)])^2
    cls1=kmeans(out.rl,centers=6)$cluster
    s.ind1[NN,7] = mean(silhouette( info1, dist(out.rl))[,3])
    error[NN,7] = (rand.index(cls2[1:n],info2)+rand.index(cls1,info1))/2
    error2[NN,7] = rand.index(cls2[1:n],info2)
    
    #Prop
    
    
    if(ali.scr[NN,pp]!=1){
      prop.svd = svds(K.mat, k=r+2)
      out.prop = prop.svd$v[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
      #plot(umap(out.prop))
      cls2=kmeans(out.prop,centers=6)$cluster
      s.ind2[NN,5] = mean(silhouette( info2, dist(out.prop[1:n,]))[,3])
      #plot(umap(out.prop),col=info2)
      out.prop = prop.svd$u[,1:r+1] %*% diag(prop.svd$d[1:r+1])^(1/2)
      cls1=kmeans(out.prop,centers=6)$cluster
      s.ind1[NN,5] = mean(silhouette( info1, dist(out.prop))[,3])
      error[NN,5] = (rand.index(cls2[1:n],info2)+rand.index(cls1,info1))/2
      error2[NN,5] = rand.index(cls2[1:n],info2)
    }else{
      error[NN,5] = NA
      error2[NN,5] = NA
    }
    print(NN)
  }
  error.out[pp,] = colMeans(error,na.rm=T)
  error.out2[pp,] = colMeans(error2,na.rm=T)
  s.out1[pp,] = colMeans(s.ind1,na.rm=T)
  s.out2[pp,] = colMeans(s.ind2,na.rm=T)
}

#error.out

data = data.frame(rand.index = c(error.out),
                  method = factor(rep(c("pca", "kpca", "j-pca", "j-kpca","prop",
                                        "lbdm","rl"), 
                                      each=length(sig.v))),
                  dim = rep(sig.v, times = 7))

data$method <- factor(data$method, levels = c("prop","lbdm","rl","j-pca", "j-kpca", "kpca", "pca"))

ggplot(data, aes(x=dim, y=rand.index, group=factor(method), colour = method)) +
  geom_line(aes(linetype=method))+xlab("tau")+
  geom_point(aes(shape=method) ,size=5) + ylab("rand index") + 
  theme(axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14,face="bold"), #change legend title font size
        legend.text = element_text(size=14,face="bold")) #600*500


data = data.frame(rand.index = c(error.out2),
                  method = factor(rep(c("pca", "kpca", "j-pca", "j-kpca","prop",
                                        "lbdm","rl"), 
                                      each=length(sig.v))),
                  dim = rep(sig.v, times = 7))

data$method <- factor(data$method, levels = c("prop","lbdm","rl","j-pca", "j-kpca", "kpca", "pca"))

ggplot(data, aes(x=dim, y=rand.index, group=factor(method), colour = method)) +
  geom_line(aes(linetype=method))+xlab("tau")+
  geom_point(aes(shape=method) ,size=5) + ylab("rand index (dataset2)") + 
  theme(axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14,face="bold"), #change legend title font size
        legend.text = element_text(size=14,face="bold")) #600*500


boxplot(ali.scr,  xlab="tau",ylab="KNN purity")
ali.scr=as.data.frame(ali.scr)
colnames(ali.scr)=round(sig.v,2)
