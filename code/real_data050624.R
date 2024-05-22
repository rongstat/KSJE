library(Seurat)
library(SeuratData)
library(uwot)
library(Rfast)
library(RSpectra)
library(fossil)
library(ggplot2)
library(HDF5Array)
library(zellkonverter)
library(cluster)
############################IFNB
# install dataset
InstallData("ifnb")
# load dataset
LoadData("ifnb")
ifnb=UpdateSeuratObject(ifnb)
# split the dataset into a list 
data.list <- SplitObject(ifnb, split.by = "stim")

# normalize and identify variable features for each dataset independently
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = data.list, nfeatures =1000)

table(data.list$CTRL@meta.data$seurat_annotations)
#CD14 Mono  CD4 Naive T CD4 Memory T    CD16 Mono            B        CD8 T  T activated           NK           DC 
#2215          978          859          507          407          352          300          298          258 
#B Activated           Mk          pDC        Eryth 
#185          115           51           23 
table(data.list$STIM@meta.data$seurat_annotations)
#CD14 Mono  CD4 Naive T CD4 Memory T    CD16 Mono            B        CD8 T  T activated           NK           DC 
#2147         1526          903          537          571          462          333          321          214 
#B Activated           Mk          pDC        Eryth 
#203          121           81           32 


data.X1 = as.matrix(data.list$CTRL@assays$RNA@data)
ct.X1 = data.list$CTRL@meta.data$seurat_annotations
data.X2 = as.matrix(data.list$STIM@assays$RNA@data)
ct.X2 = data.list$STIM@meta.data$seurat_annotations

########### feature selection
data.X1 = data.X1[match(features,rownames(data.X1)),]
data.X2 = data.X2[match(features,rownames(data.X2)),]
ct.X1 = ct.X1[which(colSums(data.X1)!=0)]
ct.X2 = ct.X2[which(colSums(data.X2)!=0)]
data.X1 = data.X1[,which(colSums(data.X1)!=0)]
data.X2 = data.X2[,which(colSums(data.X2)!=0)]

dim(data.X1)
dim(data.X2)



############# kernel joint embedding


data.1 = t(scale(t(data.X1), center = T, scale = F))
data.2 = t(scale(t(data.X2), center = T, scale = F))
ct.X1=factor(ct.X1)
ct.X2=factor(ct.X2)


dist.mat = Dist(t(cbind(data.1, data.2)))
dim(as.matrix(dist.mat)[1:dim(data.1)[2], (dim(data.1)[2]+1):(dim(data.1)[2]+dim(data.2)[2])])
cross.dist = as.matrix(dist.mat)[1:dim(data.1)[2], (dim(data.1)[2]+1):(dim(data.1)[2]+dim(data.2)[2])]

K.mat = exp(-cross.dist^2/quantile(cross.dist,0.5)^2)
K.mat.rl = exp(-cross.dist^2/quantile(dist.mat,0.05)^2) #common manifold model & discussion

K.mat.all = exp(-as.matrix(dist.mat)^2/quantile(dist.mat,0.5)^2)

K.mat1=exp(-as.matrix(dist.mat)[1:dim(data.1)[2],1:dim(data.1)[2]]^2/quantile(dist.mat[1:dim(data.1)[2],1:dim(data.1)[2]],0.5)^2)
K.mat2=exp(-as.matrix(dist.mat)[(dim(data.1)[2]+1):(dim(data.1)[2]+dim(data.2)[2]),
                                (dim(data.1)[2]+1):(dim(data.1)[2]+dim(data.2)[2])]^2/
             quantile(dist.mat[(dim(data.1)[2]+1):(dim(data.1)[2]+dim(data.2)[2]),(dim(data.1)[2]+1):(dim(data.1)[2]+
                                                                                                        dim(data.2)[2])],0.5)^2)


alpha=1
A = K.mat
D1 = rowSums(K.mat)
D2 = colSums(K.mat)
A = diag(D1^(-1/2)) %*% A %*% diag(D2^(-1/2))

D.l = rowSums(K.mat.rl %*% t(K.mat.rl))
r.svd.l = svds(diag((D.l)^(-1/2)) %*% K.mat.rl, k=21)
D.r = rowSums(t(K.mat) %*% (K.mat))
r.svd.r = svds(diag((D.r)^(-1/2)) %*% t(K.mat.rl), k=21)


svd.out = svds(cbind(data.1,data.2), k=50)
plot(svd.out$d[1:45]/svd.out$d[2:46])
which(svd.out$d[1:45]/svd.out$d[2:46]>1.02)

r.index.kcca1=c()
r.index.kpca1=c()
r.index.cca1=c()
r.index.pca1=c()
r.index.kcca2=c()
r.index.kpca2=c()
r.index.cca2=c()
r.index.pca2=c()
r.index.sep1=c()
r.index.sep2=c()
r.index.sepl1=c()
r.index.sepl2=c()
r.index.rl1=c()
r.index.rl2=c()
r.index.lbdm1=c()
r.index.lbdm2=c()
k=0
method="average"
for(r in 5:20){ 
  print(r)
  k=k+1
  #### KCCA 
  cca.out = svds(K.mat, k=r+1)
  hc=hclust(dist(cca.out$u[,1:r+1] %*% diag(cca.out$d[1:r+1])^(1/2)), method=method)
  memb = cutree(hc,k=r+1)
  r.index.kcca1[k]=rand.index(memb, as.numeric(ct.X1))
  hc=hclust(dist(cca.out$v[,1:r+1] %*% diag(cca.out$d[1:r+1])^(1/2)), method=method)
  memb = cutree(hc,k=r+1)
  r.index.kcca2[k]=rand.index(memb, as.numeric(ct.X2))
  
  
  #### j-KPCA
  eig.out = eigs(K.mat.all, k=r+1)
  X.combined = eig.out$vectors[,1:r+1] %*% diag(eig.out$values[1:r+1])^(1/2)
  hc=hclust(dist(X.combined[1:dim(data.1)[2],]), method=method)
  memb = cutree(hc,k=r+1)
  r.index.kpca1[k]=rand.index(memb, as.numeric(ct.X1))
  hc=hclust(dist(X.combined[(dim(data.1)[2]+1):(dim(data.1)[2]+dim(data.2)[2]),]), method=method)
  memb = cutree(hc,k=r+1)
  r.index.kpca2[k]=rand.index(memb, as.numeric(ct.X2))
  
  ##CCA
  svd.out = svds(t(data.1) %*% (data.2), k=r, nu = r, nv = r)
  hc=hclust(dist(svd.out$u[,1:r] %*% diag(svd.out$d[1:r])^(1/2)), method=method)
  memb = cutree(hc,k=r+1)
  r.index.cca1[k]=rand.index(memb, as.numeric(ct.X1))
  hc=hclust(dist(svd.out$v[,1:r] %*% diag(svd.out$d[1:r])^(1/2)), method=method)
  memb = cutree(hc,k=r+1)
  r.index.cca2[k]=rand.index(memb, as.numeric(ct.X2))
  
  ##j-PCA
  svd.out = svds(cbind(data.1,data.2), k=r, nu = r, nv = r)
  X.combined =svd.out$v[,1:r] %*% diag(svd.out$d[1:r])^(1/2)
  hc=hclust(dist(X.combined[1:dim(data.1)[2],]), method=method)
  memb = cutree(hc,k=r+1)
  r.index.pca1[k]=rand.index(memb, as.numeric(ct.X1))
  hc=hclust(dist(X.combined[(dim(data.1)[2]+1):(dim(data.1)[2]+dim(data.2)[2]),]), method=method)
  memb = cutree(hc,k=r+1)
  r.index.pca2[k]=rand.index(memb, as.numeric(ct.X2))
  
  ##kpca
  eig.out = eigs(K.mat1, k=r+1)
  X.combined = eig.out$vectors[,1:r+1] %*% diag(eig.out$values[1:r+1])^(1/2) 
  hc=hclust(dist(X.combined), method=method)
  memb = cutree(hc,k=r+1)
  r.index.sep1[k]=rand.index(memb, as.numeric(ct.X1))
  
  eig.out = eigs(K.mat2, k=r+1)
  X.combined = eig.out$vectors[,1:r+1] %*% diag(eig.out$values[1:r+1])^(1/2)
  hc=hclust(dist(X.combined), method=method)
  memb = cutree(hc,k=r+1)
  r.index.sep2[k]=rand.index(memb, as.numeric(ct.X2))
  
  ##pca
  svd.out = svds(data.1, k=r)
  X.combined = svd.out$v[,1:r] %*% diag(svd.out$d[1:r])^(1/2)
  hc=hclust(dist(X.combined), method=method)
  memb = cutree(hc,k=r+1)
  r.index.sepl1[k]=rand.index(memb, as.numeric(ct.X1))
  
  svd.out = svds(data.2, k=r, nu = r, nv = r)
  X.combined = svd.out$v[,1:r] %*% diag(svd.out$d[1:r])^(1/2) 
  hc=hclust(dist(X.combined), method=method)
  memb = cutree(hc,k=r+1)
  r.index.sepl2[k]=rand.index(memb, as.numeric(ct.X2))
  
  ## LBDM
  lbdm.out = svds(A, k=r+1)
  X.combined = rbind(lbdm.out$u[,1:r+1], lbdm.out$v[,1:r+1])
  V = diag(c(D1,D2)^(-1/2)) %*% X.combined %*% diag(lbdm.out$d[1:r+1]^alpha)
  hc=hclust(dist(V[1:dim(data.1)[2],]), method=method)
  memb = cutree(hc,k=r+1)
  r.index.lbdm1[k]=rand.index(memb, as.numeric(ct.X1))
  hc=hclust(dist(V[(dim(data.1)[2]+1):(dim(data.1)[2]+dim(data.2)[2]),]), method=method)
  memb = cutree(hc,k=r+1)
  r.index.lbdm2[k]=rand.index(memb, as.numeric(ct.X2))
  
  #Roseland
  hc=hclust(dist(diag((D.l)^(-1/2)) %*% r.svd.l$u[,2:(r+1)] %*% diag(r.svd.l$d[2:(r+1)])^2), method=method)
  memb = cutree(hc,k=r+1)
  r.index.rl1[k]=rand.index(memb, as.numeric(ct.X1))
  hc=hclust(dist(diag((D.r)^(-1/2)) %*% r.svd.r$u[,2:(r+1)] %*% diag(r.svd.r$d[2:(r+1)])^2), method=method)
  memb = cutree(hc,k=r+1)
  r.index.rl2[k]=rand.index(memb, as.numeric(ct.X2))
  
}


fig.data = data.frame(rand.index = c(r.index.kcca1+r.index.kcca2, r.index.cca1+r.index.cca2, 
                                     r.index.kpca1+r.index.kpca2, r.index.pca1+r.index.pca2, r.index.sep1+r.index.sep2, 
                                     r.index.sepl1+ r.index.sepl2, r.index.lbdm1+r.index.lbdm2,r.index.rl1+r.index.rl2)/2,
                      method = rep(c("prop","seurat", "j-kpca", "j-pca", "kpca","pca","lbdm","rl"), 
                                   each = length(r.index.kcca1)), 
                      r = rep(5:20, 8)) 

ggplot(data=fig.data, aes(x=reorder(method, rand.index, FUN = median), y=rand.index)) + geom_boxplot()+ xlab("")+
  ylab("rand index")+
  theme(axis.text=element_text(size=15,face="bold"),
        axis.text.x = element_text(angle = 35, vjust=0.7),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14,face="bold"), #change legend title font size
        legend.text = element_text(size=14,face="bold")) #5*4


#######################################ATAC

data0=readH5AD("/Users/rom2554/Dropbox/Xiucai-Rong/Project II/Datasets/small_atac_gene_activity.h5ad", reader = "python")
counts = data0@assays@data@listData$counts
dim(counts)
colnames(counts) = rownames(data0@colData)
rownames(counts) = names(data0@rowRanges)

meta.data= data.frame(method = data0@colData$batchname, cell_type=data0@colData$final_cell_label)
row.names(meta.data)=c(colnames(counts))
all.data <- CreateSeuratObject(counts = counts, project = "brain", min.cells = 3,
                               min.features = 10, meta.data = meta.data)
ct.X1 = meta.data$cell_type[which(meta.data$method=="10x Genomics")]
ct.X2 = meta.data$cell_type[which(meta.data$method=="Fang et al.")]
table(ct.X1)
table(ct.X2)
data.X1 = counts[,which(meta.data$method=="10x Genomics")]
data.X2 =  counts[,which(meta.data$method=="Fang et al.")]
all.RNA = cbind(data.X1, data.X2)
dim(all.RNA)
meta.data=data.frame(cell_type = factor(c(ct.X1,ct.X2)),
                     method = factor(c(rep("data1",length(ct.X1)), rep("data2",length(ct.X2)))))
row.names(meta.data)=c(colnames(all.RNA))
all.data <- CreateSeuratObject(counts = all.RNA, project = "align", min.cells = 3,
                               min.features = 10, meta.data = meta.data)
all.data <- SplitObject(all.data, split.by = "method")
all.data <- lapply(X = all.data, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
  x <- NormalizeData(x)
})
features <- SelectIntegrationFeatures(object.list = all.data, nfeatures =1000)
data.X1 = as.matrix(all.data$data1@assays$RNA$data)
data.X2 = as.matrix(all.data$data2@assays$RNA$data)
data.X1 = data.X1[match(features,rownames(data.X1)),]
data.X2 = data.X2[match(features,rownames(data.X2)),]


############# kernel joint embedding


data.1 = t(scale(t(data.X1), center = T, scale = F))
data.2 = t(scale(t(data.X2), center = T, scale = F))
ct.X1=factor(ct.X1)
ct.X2=factor(ct.X2)


dist.mat = Dist(t(cbind(data.1, data.2)))
dim(as.matrix(dist.mat)[1:dim(data.1)[2], (dim(data.1)[2]+1):(dim(data.1)[2]+dim(data.2)[2])])
cross.dist = as.matrix(dist.mat)[1:dim(data.1)[2], (dim(data.1)[2]+1):(dim(data.1)[2]+dim(data.2)[2])]

K.mat = exp(-cross.dist^2/quantile(cross.dist,0.5)^2)
K.mat.rl = exp(-cross.dist^2/quantile(dist.mat,0.05)^2) #common manifold model & discussion

K.mat.all = exp(-as.matrix(dist.mat)^2/quantile(dist.mat,0.5)^2)

K.mat1=exp(-as.matrix(dist.mat)[1:dim(data.1)[2],1:dim(data.1)[2]]^2/quantile(dist.mat[1:dim(data.1)[2],1:dim(data.1)[2]],0.5)^2)
K.mat2=exp(-as.matrix(dist.mat)[(dim(data.1)[2]+1):(dim(data.1)[2]+dim(data.2)[2]),
                                (dim(data.1)[2]+1):(dim(data.1)[2]+dim(data.2)[2])]^2/
             quantile(dist.mat[(dim(data.1)[2]+1):(dim(data.1)[2]+dim(data.2)[2]),(dim(data.1)[2]+1):(dim(data.1)[2]+
                                                                                                        dim(data.2)[2])],0.5)^2)


alpha=1
A = K.mat
D1 = rowSums(K.mat)
D2 = colSums(K.mat)
A = diag(D1^(-1/2)) %*% A %*% diag(D2^(-1/2))

D.l = rowSums(K.mat.rl %*% t(K.mat.rl))
r.svd.l = svds(diag((D.l)^(-1/2)) %*% K.mat.rl, k=21)
D.r = rowSums(t(K.mat) %*% (K.mat))
r.svd.r = svds(diag((D.r)^(-1/2)) %*% t(K.mat.rl), k=21)


svd.out = svds(cbind(data.1,data.2), k=50)
plot(svd.out$d[1:45]/svd.out$d[2:46])
which(svd.out$d[1:45]/svd.out$d[2:46]>1.02)

r.index.kcca1=c()
r.index.kpca1=c()
r.index.cca1=c()
r.index.pca1=c()
r.index.kcca2=c()
r.index.kpca2=c()
r.index.cca2=c()
r.index.pca2=c()
r.index.sep1=c()
r.index.sep2=c()
r.index.sepl1=c()
r.index.sepl2=c()
r.index.rl1=c()
r.index.rl2=c()
r.index.lbdm1=c()
r.index.lbdm2=c()
k=0
method="complete"
for(r in 5:20){ 
  print(r)
  k=k+1
  #### KCCA 
  cca.out = svds(K.mat, k=r+1)
  hc=hclust(dist(cca.out$u[,1:r+1] %*% diag(cca.out$d[1:r+1])^(1/2)), method=method)
  memb = cutree(hc,k=r+1)
  r.index.kcca1[k]=rand.index(memb, as.numeric(ct.X1))
  hc=hclust(dist(cca.out$v[,1:r+1] %*% diag(cca.out$d[1:r+1])^(1/2)), method=method)
  memb = cutree(hc,k=r+1)
  r.index.kcca2[k]=rand.index(memb, as.numeric(ct.X2))
  
  
  #### j-KPCA
  eig.out = eigs(K.mat.all, k=r+1)
  X.combined = eig.out$vectors[,1:r+1] %*% diag(eig.out$values[1:r+1])^(1/2)
  hc=hclust(dist(X.combined[1:dim(data.1)[2],]), method=method)
  memb = cutree(hc,k=r+1)
  r.index.kpca1[k]=rand.index(memb, as.numeric(ct.X1))
  hc=hclust(dist(X.combined[(dim(data.1)[2]+1):(dim(data.1)[2]+dim(data.2)[2]),]), method=method)
  memb = cutree(hc,k=r+1)
  r.index.kpca2[k]=rand.index(memb, as.numeric(ct.X2))
  
  ##CCA
  svd.out = svds(t(data.1) %*% (data.2), k=r, nu = r, nv = r)
  hc=hclust(dist(svd.out$u[,1:r] %*% diag(svd.out$d[1:r])^(1/2)), method=method)
  memb = cutree(hc,k=r+1)
  r.index.cca1[k]=rand.index(memb, as.numeric(ct.X1))
  hc=hclust(dist(svd.out$v[,1:r] %*% diag(svd.out$d[1:r])^(1/2)), method=method)
  memb = cutree(hc,k=r+1)
  r.index.cca2[k]=rand.index(memb, as.numeric(ct.X2))
  
  ##j-PCA
  svd.out = svds(cbind(data.1,data.2), k=r, nu = r, nv = r)
  X.combined =svd.out$v[,1:r] %*% diag(svd.out$d[1:r])^(1/2)
  hc=hclust(dist(X.combined[1:dim(data.1)[2],]), method=method)
  memb = cutree(hc,k=r+1)
  r.index.pca1[k]=rand.index(memb, as.numeric(ct.X1))
  hc=hclust(dist(X.combined[(dim(data.1)[2]+1):(dim(data.1)[2]+dim(data.2)[2]),]), method=method)
  memb = cutree(hc,k=r+1)
  r.index.pca2[k]=rand.index(memb, as.numeric(ct.X2))
  
  ##kpca
  eig.out = eigs(K.mat1, k=r+1)
  X.combined = eig.out$vectors[,1:r+1] %*% diag(eig.out$values[1:r+1])^(1/2) 
  hc=hclust(dist(X.combined), method=method)
  memb = cutree(hc,k=r+1)
  r.index.sep1[k]=rand.index(memb, as.numeric(ct.X1))
  
  eig.out = eigs(K.mat2, k=r+1)
  X.combined = eig.out$vectors[,1:r+1] %*% diag(eig.out$values[1:r+1])^(1/2)
  hc=hclust(dist(X.combined), method=method)
  memb = cutree(hc,k=r+1)
  r.index.sep2[k]=rand.index(memb, as.numeric(ct.X2))
  
  ##pca
  svd.out = svds(data.1, k=r)
  X.combined = svd.out$v[,1:r] %*% diag(svd.out$d[1:r])^(1/2)
  hc=hclust(dist(X.combined), method=method)
  memb = cutree(hc,k=r+1)
  r.index.sepl1[k]=rand.index(memb, as.numeric(ct.X1))
  
  svd.out = svds(data.2, k=r, nu = r, nv = r)
  X.combined = svd.out$v[,1:r] %*% diag(svd.out$d[1:r])^(1/2) 
  hc=hclust(dist(X.combined), method=method)
  memb = cutree(hc,k=r+1)
  r.index.sepl2[k]=rand.index(memb, as.numeric(ct.X2))
  
  ## LBDM
  lbdm.out = svds(A, k=r+1)
  X.combined = rbind(lbdm.out$u[,1:r+1], lbdm.out$v[,1:r+1])
  V = diag(c(D1,D2)^(-1/2)) %*% X.combined %*% diag(lbdm.out$d[1:r+1]^alpha)
  hc=hclust(dist(V[1:dim(data.1)[2],]), method=method)
  memb = cutree(hc,k=r+1)
  r.index.lbdm1[k]=rand.index(memb, as.numeric(ct.X1))
  hc=hclust(dist(V[(dim(data.1)[2]+1):(dim(data.1)[2]+dim(data.2)[2]),]), method=method)
  memb = cutree(hc,k=r+1)
  r.index.lbdm2[k]=rand.index(memb, as.numeric(ct.X2))
  
  #Roseland
  hc=hclust(dist(diag((D.l)^(-1/2)) %*% r.svd.l$u[,2:(r+1)] %*% diag(r.svd.l$d[2:(r+1)])^2), method=method)
  memb = cutree(hc,k=r+1)
  r.index.rl1[k]=rand.index(memb, as.numeric(ct.X1))
  hc=hclust(dist(diag((D.r)^(-1/2)) %*% r.svd.r$u[,2:(r+1)] %*% diag(r.svd.r$d[2:(r+1)])^2), method=method)
  memb = cutree(hc,k=r+1)
  r.index.rl2[k]=rand.index(memb, as.numeric(ct.X2))
  
}


fig.data = data.frame(rand.index = c(r.index.kcca1+r.index.kcca2, r.index.cca1+r.index.cca2, 
                                     r.index.kpca1+r.index.kpca2, r.index.pca1+r.index.pca2, r.index.sep1+r.index.sep2, 
                                     r.index.sepl1+ r.index.sepl2, r.index.lbdm1+r.index.lbdm2,r.index.rl1+r.index.rl2)/2,
                      method = rep(c("prop","seurat", "j-kpca", "j-pca", "kpca","pca","lbdm","rl"), 
                                   each = length(r.index.kcca1)), 
                      r = rep(5:20, 8)) 

ggplot(data=fig.data, aes(x=reorder(method, rand.index, FUN = median), y=rand.index)) + geom_boxplot()+ xlab("")+
  ylab("rand index")+
  theme(axis.text=element_text(size=15,face="bold"),
        axis.text.x = element_text(angle = 35, vjust=0.7),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14,face="bold"), #change legend title font size
        legend.text = element_text(size=14,face="bold")) #5*4


