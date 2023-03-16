rm(list = setdiff("main.dir", ls()))

library(Biobase)

print("Calculating residuals for categorical variables ...")

load("./data/mets_eset_with_covars_example.RData")

transformMetabolitesToProbitScores = function(eset){
   
  my.pdata = pData(eset)
  my.met.exp = exprs(eset)
  
  for(i in 1:dim(my.met.exp)[1]){
    orig.met = my.met.exp[i,]
    
    # skip NAs in the transformation
    is.not.na.index = which(!is.na(orig.met))
    cur.met = my.met.exp[i,is.not.na.index]
    
    # transform to probit scores
    N = length(cur.met)
    cur.met.r = rank(cur.met)/(N+1)
    probit.met = qnorm(p = cur.met.r)
    
    # set values
    my.met.exp[i, is.not.na.index] = probit.met
  }
  
  result.eset = ExpressionSet(assayData = my.met.exp,
                              phenoData = AnnotatedDataFrame(my.pdata),
                              featureData = AnnotatedDataFrame(fData(eset)))
  
}

createESetBasedOnResiduals = function(my.eset, my.covars){
    
  pdata = pData(my.eset)
  
  levels(pdata$hrace) = c(" : non-hispanic caucasian"," : other"," : other"," : other")
  pData(my.eset) = pdata
  
  mdata = exprs(my.eset)
  mets.residuals.m1 = mdata
  colnames(mets.residuals.m1) = pdata$id
  
  for(i in 1:dim(mdata)[1]){
   
    curr.data.m1 = data.frame(cbind(met = mdata[i,], pdata[,my.covars]))
    curr.data.m1.help = curr.data.m1[complete.cases(curr.data.m1),]
    
    if(dim(curr.data.m1.help)[1]<10){
      mets.residuals.m1[i,] = met = mdata[i,]
      
    }else{
      curr.data.m1 = curr.data.m1.help
      curr.data.m1[,-1] = lapply(X = curr.data.m1[,-1], FUN = function(x) {return(as.factor(as.character(x)))})
      num.levels.m1 = lapply(X = curr.data.m1[,-1], FUN = function(x) {return(length(levels(x)))})
      curr.data.m1 = curr.data.m1[,c(1,which(num.levels.m1>1)+1)]
      
    
      help.m1 = residuals(lm(met~.,data = curr.data.m1))
      mets.residuals.m1[i,names(help)] = help.m1
    }
    
  }   
  
  my.eset.residuals = ExpressionSet(assayData = mets.residuals.m1, phenoData = AnnotatedDataFrame(pdata), featureData = AnnotatedDataFrame(fData(my.eset)))
   
  return(my.eset.residuals)
}

#conditional logistic regression, no matching factors in M1  
mets.glauc.m1 = mets.glauc
mets.glauc.m2 = createESetBasedOnResiduals(my.eset = mets.glauc, 
                                           my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort","smokec")) 

mets.glauc.m3 = createESetBasedOnResiduals(my.eset = mets.glauc, 
                                           my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort","smokec","hgla","hrace")) 
mets.glauc.m4 = mets.glauc.m3

mets.glauc.m5 = createESetBasedOnResiduals(my.eset = mets.glauc, 
                                           my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort","smokec","hgla","hrace","hbph","cholh","sdbh","eversteroid")) 

save(mets.glauc, mets.glauc.m1, mets.glauc.m2,mets.glauc.m3, mets.glauc.m4, mets.glauc.m5, file = "data/mets_eset_with_covars_and_residuals.RData")

load("data/mets_eset_with_covars_by_subtype.RData")
# subtype analyses, unconditional logistic regression, matching factors go into model 1
# unconditional logistic regression, matching factors in M1  
mets.glaucpa.m1 = createESetBasedOnResiduals(my.eset = mets.glaucpa, 
                                             my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort")) 

mets.glaucpa.m2 = createESetBasedOnResiduals(my.eset = mets.glaucpa, 
                                           my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort","smokec")) 

mets.glaucpa.m3 = createESetBasedOnResiduals(my.eset = mets.glaucpa, 
                                           my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort","smokec","hgla","hrace")) 
mets.glaucpa.m4 = mets.glaucpa.m3

mets.glaucpa.m5 = createESetBasedOnResiduals(my.eset = mets.glaucpa, 
                                           my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort","smokec","hgla","hrace","hbph","cholh","sdbh","eversteroid")) 
# unconditional logistic regression, matching factors in M1  
mets.glaucpe.m1 = createESetBasedOnResiduals(my.eset = mets.glaucpe, 
                                             my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort")) 

mets.glaucpe.m2 = createESetBasedOnResiduals(my.eset = mets.glaucpe, 
                                             my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort","smokec")) 

mets.glaucpe.m3 = createESetBasedOnResiduals(my.eset = mets.glaucpe, 
                                             my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort","smokec","hgla","hrace")) 
mets.glaucpe.m4 = mets.glaucpe.m3

mets.glaucpe.m5 = createESetBasedOnResiduals(my.eset = mets.glaucpe, 
                                             my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort","smokec","hgla","hrace","hbph","cholh","sdbh","eversteroid")) 

save(mets.glaucpa, mets.glaucpa.m1, mets.glaucpa.m2, mets.glaucpa.m3, mets.glaucpa.m4, mets.glaucpa.m5, 
     mets.glaucpe, mets.glaucpe.m1, mets.glaucpe.m2, mets.glaucpe.m3, mets.glaucpe.m4, mets.glaucpe.m5, 
     file = "data/mets_eset_with_covars_by_subtype_and_residuals.RData")

################################################################################
# stratified analyses
load("data/mets_eset_with_covars_for_stratified_analyses.RData")
################################################################################

### BMI; unconditional logistic regression, matching factors in M1  
normal.bmi.m1 = createESetBasedOnResiduals(my.eset = normal.bmi, 
                                           my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort")) 
normal.bmi.m5 = createESetBasedOnResiduals(my.eset = normal.bmi, 
                                           my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort","smokec","hgla","hrace","hbph","cholh","sdbh","eversteroid")) 
high.bmi.m1 = createESetBasedOnResiduals(my.eset = high.bmi, 
                                         my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort")) 
high.bmi.m5 = createESetBasedOnResiduals(my.eset = high.bmi, 
                                         my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort","smokec","hgla","hrace","hbph","cholh","sdbh","eversteroid")) 

### age; unconditional logistic regression, matching factors in M1  
age.high.glauc.m1 = createESetBasedOnResiduals(my.eset = age.high.glauc, 
                                           my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort")) 
age.high.glauc.m5 = createESetBasedOnResiduals(my.eset = age.high.glauc, 
                                           my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort","smokec","hgla","hrace","hbph","cholh","sdbh","eversteroid")) 
age.low.glauc.m1 = createESetBasedOnResiduals(my.eset = age.low.glauc, 
                                         my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort")) 
age.low.glauc.m5 = createESetBasedOnResiduals(my.eset = age.low.glauc, 
                                         my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort","smokec","hgla","hrace","hbph","cholh","sdbh","eversteroid")) 

### family history; unconditional logistic regression, matching factors in M1  
fh.glauc.m1 = createESetBasedOnResiduals(my.eset = fh.glauc, 
                                               my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort")) 
fh.glauc.m5 = createESetBasedOnResiduals(my.eset = fh.glauc, 
                                               my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort","smokec","hrace","hbph","cholh","sdbh","eversteroid")) 
nfh.glauc.m1 = createESetBasedOnResiduals(my.eset = nfh.glauc, 
                                         my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort")) 
nfh.glauc.m5 = createESetBasedOnResiduals(my.eset = nfh.glauc, 
                                         my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort","smokec","hrace","hbph","cholh","sdbh","eversteroid")) 

### by sex; conditional logistic regression, no matching factors in M1
# among men, I will not adjust for hrace, sdbh, and eversteroid use, due to very low numbers of, or zero, participants across categories 
# collapsing " : past" and " : current" smoker categories
pdata = pData(men.glauc)
table(pdata$smokec, pdata$caco)
levels(pdata$smokec) = c(" : never",   " : past",    " : past")
pData(men.glauc) = pdata

men.glauc.m1 = men.glauc
men.glauc.m5 = createESetBasedOnResiduals(my.eset = men.glauc, 
                                         my.covars = c("fast", "bldtime", "bldmon", "smokec","hbph","cholh","hgla")) 

women.glauc.m1 = women.glauc
women.glauc.m5 = createESetBasedOnResiduals(my.eset = women.glauc, 
                                          my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort","smokec","hgla","hrace","hbph","cholh","sdbh","eversteroid")) 

save(women.glauc.m1, men.glauc.m1, women.glauc.m5, men.glauc.m5,
     file="data/mets_eset_with_covars_for_stratified_analyses_and_residuals.RData")

### by time blood-DX; conditional logistic regression, no matching factors in M1  
closeDX.glauc.m1 = closeDX.glauc
closeDX.glauc.m5 = createESetBasedOnResiduals(my.eset = closeDX.glauc, 
                                          my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort","smokec","hgla","hrace","hbph","cholh","sdbh","eversteroid")) 

farDX.glauc.m1 = farDX.glauc
farDX.glauc.m5 = createESetBasedOnResiduals(my.eset = farDX.glauc, 
                                            my.covars = c("fast", "menobld", "menodiag","bldtime", "bldmon","cohort","smokec","hgla","hrace","hbph","cholh","sdbh","eversteroid")) 

save(women.glauc.m1, men.glauc.m1, women.glauc.m5, men.glauc.m5,
     age.low.glauc.m1, age.high.glauc.m1, age.low.glauc.m5, age.high.glauc.m5,
     fh.glauc.m1, nfh.glauc.m1, fh.glauc.m5, nfh.glauc.m5,
     normal.bmi.m1, high.bmi.m1, normal.bmi.m5, high.bmi.m5,
     closeDX.glauc.m1, farDX.glauc.m1, closeDX.glauc.m5, farDX.glauc.m5,
     file="data/mets_eset_with_covars_for_stratified_analyses_and_residuals.RData")
print("Done")
