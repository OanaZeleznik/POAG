rm(list = setdiff("main.dir", ls()))

library(Biobase)

print("Calculating the number of effective tests ...")

load("data/mets_eset_with_covars_example.RData")

### all samples
fdata = fData(mets.glauc)
mdata = exprs(mets.glauc)
pdata = pData(mets.glauc)
tapply(pdata$agemeno, pdata$sex, summary)

cov.mat = cov(t(mdata), use = "pairwise.complete.obs")
cov.mat[which(is.na(cov.mat),arr.ind = T)] = 0

# running singular value decomposition
pca.res = eigen(cov.mat)
pca.res$values[1:10]
print("NEF = 132 based on the real data which will be used in this analysis")
print("The example dataset has NEF = ")
print(which(cumsum(pca.res$values/sum(pca.res$values))>0.995)[1])
print("Done")
