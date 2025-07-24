
library(Rfast)
library(RSpectra)
library(BiocNeighbors)
library(pcaPP)

#main function: duo-landkmark joint embedding

#X - centered and normalized dataset1, n1 x p
#Y - centered and normalized dataset2, n2 x p 
#w - numerical value between 0 and 1, determing the bandwidth selection. Default = 0.5.
#Gamma - a collecting of integers smaller than n1 and n2, being the index set of latent embedding. We recommend choosing Gamma=c(1:(r+1)) where r is the desired latent dimension.

DL.embed <- function(X,Y, w=0.5, Gamma=c(1:3)){
  dist.mat = Dist(rbind(data.Y1, data.Y2))
  cross.dist = as.matrix(dist.mat)[1:dim(data.Y1)[1], (dim(data.Y1)[1]+1):(dim(data.Y1)[1]+dim(data.Y2)[1])]
  
  K.mat = exp(-cross.dist^2/quantile(cross.dist,0.5)^2)  
  prop.svd = svds(K.mat, k=max(Gamma))
  
  out1= sqrt(dim(prop.svd$v[,Gamma])[1])*prop.svd$v[,Gamma] %*% diag(prop.svd$d[1:Gamma])^(1/2)
  out2= sqrt(dim(prop.svd$u[,Gamma])[1])*prop.svd$u[,Gamma] %*% diag(prop.svd$d[1:Gamma])^(1/2)
  
  return(list(embed.data1=out1, embed.data2=out2))
}
