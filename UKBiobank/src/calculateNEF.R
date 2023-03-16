rm(list = setdiff("main.dir", ls()))
print("Calculating number of effective tests - this takes a few minutes ...")

library(Biobase)
load("./data/poag_ukbb_es_example.RData")

mdata = exprs(poag.ukb)
N = dim(mdata)[1]
cov.mat = matrix(data = NA,nrow = N, ncol = N)

for(i in 1:N){
  
  if(i%%10 == 1)
    print(paste("i =",i))
  
  for(j in 1:N){
    
    if(j%%50 == 1)
      print(paste("j =",j))
    
    if(i<j){
      c.corr = cov(x = mdata[i,], y = mdata[j,], use = "pairwise.complete.obs", method = "pearson")
      cov.mat[i,j] = c.corr
      cov.mat[j,i] = c.corr
    }  
  }
}

diag(cov.mat) = 1

length(which(is.na(cov.mat)))
range(cov.mat)
# running singular value decomposition
pca.res = eigen(cov.mat)
pca.res$values[1:10]
### NEF = 34
ukb.nef = which(cumsum(pca.res$values/sum(pca.res$values))>0.995)[1]

save(ukb.nef, file = "./results/ukb_nef.RData")

print("Done.")