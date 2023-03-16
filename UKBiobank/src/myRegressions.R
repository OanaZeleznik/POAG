library(plyr)
library(stats)
library(Biobase)
library(broom)
library(tibble)
library(survival)
library(NPC)
library(qvalue)
library(multtest)

runRobustUnconditionalLogisticRegression = function(one.exposure.data, regression.formula, 
                                                    debug = F){
  # this function implements a robust regression with tryCatch statements
  # it will be able to store errors/warnings
    
  if(debug) 
    message("Running runRobustUnconditionalLogisticRegression\n")
  
  glm.fit = tryCatch(expr = glm(formula = regression.formula, family = binomial(link = "logit"), 
                                data = one.exposure.data),
                     error = function(e) e,
                     warning = function(w) w)
  
  if(debug)
    message(paste("Result is of class: ",class(glm.fit)[1],"\n"))
  if(debug)
    message("End Running runRobustUnconditionalLogisticRegression\n")
  return(glm.fit)
}

runRobustConditionalLogisticRegression = function(one.exposure.data, regression.formula, 
                                                    debug = F){
  # this function implements a robust regression with tryCatch
  # it will be able to store errors/warnings
  
  if(debug) 
    message("Running runRobustConditionalLogisticRegression\n")
  
  clogit.fit = tryCatch(expr = clogit(formula = regression.formula, data = one.exposure.data),
                     error = function(e) e,
                     warning = function(w) w)
  
  if(debug)
    message(paste("Result is of class: ",class(clogit.fit)[1],"\n"))
  if(debug)
    message("End Running runRobustUnconditionalLogisticRegression\n")
  return(clogit.fit)
}

runRobustUnconditionalLogisticRegressionForIndicatorExposure = function(one.exposure.data,
                                                                        regression.formula, 
                                                                        debug = F){
  # this functionimplements Pete Kraft's concept for testing metabolites in the presence of 
  # missing values
    
  if(debug)
    message("Running runRobustUnconditionalLogisticRegressionForIndicatorExposure\n")
  
  glm.fit.all = runRobustUnconditionalLogisticRegression(one.exposure.data = one.exposure.data,
                                                         regression.formula = regression.formula)
  
  
  # if there were problems with model return
  if(class(glm.fit.all)[1] != "glm"){ 
    if(debug)
      message("Early End due to Error/Warning in 
              runRobustUnconditionalLogisticRegressionForIndicatorExposure\n")
    return(glm.fit.all)
  }
  
  # if there were problems with model return
  if(is.na(glm.fit.all$coefficients["exposure"])){
    if(debug)
      message("Early End due to exposure not in model in 
            runRobustUnconditionalLogisticRegressionForIndicatorExposure\n")
    return(glm.fit.all)
    
  }
     
  
  one.exposure.data = one.exposure.data[,-1*which(colnames(one.exposure.data) %in% 
                                                    c("exposure","ind.var"))] 
  
  glm.fit.small = runRobustUnconditionalLogisticRegression(one.exposure.data = one.exposure.data,
                                                           regression.formula = regression.formula)
  
  # if there were problems with model here return
  if(class(glm.fit.small)[1] != "glm"){
    if(debug)
      message("Early End due to error/warning in smaller model in 
            runRobustUnconditionalLogisticRegressionForIndicatorExposure\n")
    return(glm.fit.small)
  }
  res = anova(glm.fit.all, glm.fit.small, test = "Chisq")
    
  if(debug) 
    message("End Running runRobustUnconditionalLogisticRegressionForIndicatorExposure\n")
  return(list(glm.fit.all = glm.fit.all, anova = res))
  
}

runRobustConditionalLogisticRegressionForIndicatorExposure = function(one.exposure.data,
                                                                        regression.formula, 
                                                                        debug = F){
  # this functionimplements Pete Kraft's concept for testing metabolites in the presence of 
  # missing values
    
  if(debug)
    message("Running runRobustUnconditionalLogisticRegressionForIndicatorExposure\n")
  
  clogit.fit.all = runRobustConditionalLogisticRegression(one.exposure.data = one.exposure.data,
                                                         regression.formula = regression.formula)
  
  
  # if there were problems with model return
  if(class(clogit.fit.all)[1] != "clogit"){ 
    if(debug)
      message("Early End due to Error/Warning in 
              runRobustConditionalLogisticRegressionForIndicatorExposure\n")
    return(clogit.fit.all)
  }
  
  # if there were problems with model return
  if(is.na(clogit.fit.all$coefficients["exposure"])){
    if(debug)
      message("Early End due to exposure not in model in 
              runRobustUnconditionalLogisticRegressionForIndicatorExposure\n")
    return(clogit.fit.all)
    
  }
  
  
  one.exposure.data = one.exposure.data[,-1*which(colnames(one.exposure.data) %in% 
                                                    c("exposure","ind.var"))] 
  
  # remove exposure from formula
  reg.formula.small = paste(as.character(regression.formula[2]),as.character(regression.formula[1]),as.character(regression.formula[3]))
  reg.formula.small = gsub(pattern = " + exposure",replacement = "",x = reg.formula.small, fixed = T)
  reg.formula.small = as.formula(reg.formula.small)
  
  
  clogit.fit.small = runRobustConditionalLogisticRegression(one.exposure.data = one.exposure.data,
                                                           regression.formula = reg.formula.small)
  
  # if there were problems with model here return
  if(class(clogit.fit.small)[1] != "clogit"){
    if(debug)
      message("Early End due to error/warning in smaller model in 
            runRobustConditionalLogisticRegressionForIndicatorExposure\n")
    return(clogit.fit.small)
  }
  res = anova(clogit.fit.all, clogit.fit.small, test = "Chisq")
   
  if(debug) 
    message("End Running runRobustConditionalLogisticRegressionForIndicatorExposure\n")
  return(list(clogit.fit.all = clogit.fit.all, anova = res))
  
}

### start here to review
runUnconditionalLogisticRegression = function(all.data, outcome, covariates,
                                              allow.na.values.in.exposure, debug = F){
  
  #### 
  # this function runs unconditional logistic regressions for each exposure in the dataset
  # exposures, ouctome and covariates are in the columns of the dataset
  # for exposures with missing values, it implemntes Pete Kraft's approach to test associations if 
  # allow.na.values.in.exposure is set to FALSE
  ##### input parameters
  # all.data = dataset with exposures, outcome, covariates in columns
  # outcome = name of column in all.data containing the outcome
  # covariates = names of columns in all.data containing the covariates
  # allow.na.values.in.exposure = if FALSE: implemnte Pete Kraft's approach to test associations
  # allow.na.values.in.exposure = if TRUE: run the regressions and 
  # drop observations with missing values 
  ####
  # returns a list with the results of each regression (the whole glm object)
  ####
  
  # init result to match regression result
  error.result = tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA, 
                        convergence = "error")
  warning.result = tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA,
                          convergence = "warning")
  
  # how many exposures are in the data
  num.exposures = length(setdiff(colnames(all.data),c(outcome, covariates)))
  # index of outcome and covariates
  index.non.exp = which(colnames(all.data) %in% c(outcome,covariates))
  
  # data with exposure vars only
  exp.data = all.data[,-1*index.non.exp]
  
  # data with outcome and covariates
  non.exp.data = all.data[,index.non.exp]
  # make sure that outcome is not a factor
  non.exp.data[,outcome] = as.numeric(as.character(non.exp.data[,outcome]))
  
  
  reg.res = vector("list",num.exposures)
  names(reg.res) = colnames(exp.data)
  
  for(i in 1:num.exposures){
    
    if(i %% 50 == 1)
      message(paste("i = ",i,"\n"))
    reg.data = cbind(non.exp.data,exposure = exp.data[,i])
    reg.formula = as.formula(paste(outcome, "~ ."))
    
    # we implement metabolites with missing values as indicator variables
    if(length(which(is.na(reg.data$exposure)))>0 & allow.na.values.in.exposure == F){
      # add indicator variable: 1 for those with missing values in the exposure and 0 for the others
      reg.data$ind.var = 0
      reg.data$ind.var[which(is.na(reg.data$exposure))] = 1
      
      # set NA in exposure to 0
      reg.data$exposure[which(is.na(reg.data$exposure))] = 0
      
      # multiply exposure by (1-M)
      reg.data$exposure = reg.data$exposure*(1-reg.data$ind.var)
      
      # Model in this case is: g(E[Y]) = a + b M + c (1-M) X + d Z
      # Test whether X is associated with Y, calculate the two d.f. test of H0: b=c=0.
      reg.res[[i]] = runRobustUnconditionalLogisticRegressionForIndicatorExposure(
        one.exposure.data = reg.data, regression.formula = reg.formula)
      
    }
    else{
      reg.res[[i]] = runRobustUnconditionalLogisticRegression(one.exposure.data = reg.data,
                                                              regression.formula = reg.formula)
    }
  }
  return(reg.res)
}

runConditionalLogisticRegression = function(all.data, outcome, covariates,strata, 
                                            allow.na.values.in.exposure, debug = F)
{
  #### 
  # this function runs conditional logistic regressions for each exposure in the dataset
  # exposures, oucome and covariates are in the colums of the dataset
  # for exposures with missings values it implemnte Pete Kraft's approach to test associations if 
  # allow.na.values.in.exposure is set to FALSE
  ##### input parameters
  # all.data = dataset with exposures, outcome, covariates in columns
  # outcome = name of column in all.data containing the outcome
  # covariates = names of columns in all.data containing the covariates
  # allow.na.values.in.exposure = if FALSE: implemnte Pete Kraft's approach to test associations
  # allow.na.values.in.exposure = if TRUE: will not change anythig; just run the regressions and 
  # drop observations with missing values 
  ####
  # returns a list with the results of each regression (the whole glm object)
  #### 
  
  # init result to match regression result
  error.result = tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA, 
                        convergence = "error")
  warning.result = tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA,
                          convergence = "warning")
  
  # how many exposures are in the data
  num.exposures = length(setdiff(colnames(all.data),c(outcome, covariates, strata)))
  # index of outcome and covariates
  index.non.exp = which(colnames(all.data) %in% c(outcome,covariates,strata))
  
  # data with exposure vars only
  exp.data = all.data[,-1*index.non.exp]
  
  # data with outcome and covariates
  non.exp.data = all.data[,index.non.exp]
  # make sure that outcome is not a factor
  non.exp.data[,outcome] = as.numeric(as.character(non.exp.data[,outcome]))
    
  reg.res = vector("list",num.exposures)
  names(reg.res) = colnames(exp.data)
  
  for(i in 1:num.exposures){
  #for(i in c(1,2,93,104)){
    #i = 1
    if(i %% 50 == 1)
      cat(paste("i = ",i,"... "))
    reg.data = cbind(non.exp.data,exposure = exp.data[,i])
        
    if(is.null(covariates))
      reg.formula = as.formula(paste(outcome, "~ strata(",strata,") + exposure",sep=""))
    else
      reg.formula = as.formula(paste(outcome, "~ strata(",strata,") + exposure",
                                     paste("+",covariates,collapse = ""),sep=""))
    
    if(debug)
      print(as.character(reg.formula))
    
    # we implement metabolites with missing values as indicator variables
    if(length(which(is.na(reg.data$exposure)))>0 & allow.na.values.in.exposure == F){
      # add indicator variable: 1 for those with missing values in the exposure and 0 for the others
      if(debug) 
        print("Running CLR on indicator metabolite")
      reg.data$ind.var = 0
      reg.data$ind.var[which(is.na(reg.data$exposure))] = 1
      
      # set NA in exposure to 0
      reg.data$exposure[which(is.na(reg.data$exposure))] = 0
      
      # multiply exposure by (1-M)
      reg.data$exposure = reg.data$exposure*(1-reg.data$ind.var)
      
      # Model in this case is: g(E[Y]) = a + b M + c (1-M) X + d Z
      # Test whether X is associated with Y, calculate the two d.f. test of H0: b=c=0.
      reg.res[[i]] = runRobustConditionalLogisticRegressionForIndicatorExposure(
        one.exposure.data = reg.data, regression.formula = reg.formula, debug = debug)
      
    }
    else{
      if(debug) 
        print("Running CLR on NON-indicator metabolite")
      reg.res[[i]] = runRobustConditionalLogisticRegression(one.exposure.data = reg.data,
                                                              regression.formula = reg.formula,
                                                            debug = debug)
    }
  }
  cat("\n")
  return(reg.res)
  
}



robustConfInt = function(my.glm,parm){
  # caclulates profile CIs, which are more stable than Wald CIs (standard CIs)
    
  if(class(my.glm)[1] == "clogit"){
    return(summary(my.glm)$conf.int[parm,c("lower .95","upper .95")])
  }
  
  robustCI = tryCatch(suppressMessages(confint(my.glm,parm)),
                      error = function(e) e,
                      warning = function(w) w)
  
}

createDatasetForRegression = function(my.eset, covariates, outcome, strata, debug = F){
  # will create a dataset containing all exposures, outcome, covariates and strata
  # if outcome is coded as case/control, it will be turned in 1=case, 0=control
   
  if(debug) print("Start createDatasetForRegression")
  met.vals = exprs(my.eset)
  sample.data = pData(my.eset)
  if(debug) print(colnames(sample.data))
  
  if(length(intersect(c(outcome,covariates,strata),colnames(sample.data)))<
     length(c(outcome,covariates,strata))){
    cat("ERROR: Please make sure that outcome, covariates and/or strata are in the provided ESet\n")
    cat("Currently not in provided ESet:", setdiff(c(outcome,covariates,strata),
                                                   colnames(sample.data)),"\n")
    return(0)
  }
  
  if(debug) print("Init tests ok")
  my.data = data.frame(cbind(sample.data[,c(outcome,covariates, strata)],t(met.vals)))
  if(debug) print(colnames(my.data))
  if(debug) print(table(my.data[,outcome]))
  
  # recode case/control to 1/0 
  if(all.equal(sort(names(table(my.data[,outcome]))), c("case","control")) == T){
    if(debug) 
      message(paste("Recoding outcome variable from", paste(names(table(my.data[,outcome])),
                                                            collapse = "/"),"to "))
    my.data[,outcome] = factor(my.data[,outcome])
    my.data[,outcome] = mapvalues(my.data[,outcome], from = c("case", "control"), to = c("1", "0"))
    if(debug)
      message(paste(names(table(my.data[,outcome])),collapse = "/"))
  }
  
  if(all.equal(sort(names(table(my.data[,outcome]))), c("1","2")) == T){
    if(debug) 
      message(paste("Recoding outcome variable from", paste(names(table(my.data[,outcome])),
                                                            collapse = "/"),"to "))
    my.data[,outcome] = factor(my.data[,outcome])
    my.data[,outcome] = mapvalues(my.data[,outcome], from = c("1", "2"), to = c("1", "0"))
    if(debug)
      message(paste(names(table(my.data[,outcome])),collapse = "/"))
  }
  if(debug) print("All ok in createDatasetForRegression")
  return(my.data)
    
}

createResultTable = function(my.data, term.of.interest, debug = F){
  
  if(debug) print("createResultTable")
  
  class.help = unlist(lapply(my.data, function(x) class(x)[1]))
    
  result.table = data.frame(id = rep(NA,length(my.data)), estimate = NA, std.error = NA, 
                            statistic = NA, p.value = NA, p.value.anova = NA, l.ci = NA, u.ci = NA,
                            convergence = "reason")
  result.table$convergence = as.character(result.table$convergence)
  result.table$id = names(my.data)
  
  for(i in 1:length(my.data)){
    if(debug) cat("i = ",i,"\n")
    
    if(class(my.data[[i]])[1] %in% c("glm","clogit")){
      if(term.of.interest %in% rownames(summary(my.data[[i]])$coefficients)){
        
        help = summary(my.data[[i]])$coefficients[term.of.interest,]
        if("exp(coef)" %in% names(help))
          help = help[-1*which(names(help) == "exp(coef)")]
        result.table[i,c("estimate","std.error", "statistic","p.value" )] = help
          
        robust.ci = robustConfInt(my.data[[i]], parm = "exposure")
        
        if(class(robust.ci)[1] == "numeric"){
          result.table[i,c("l.ci","u.ci")] = robust.ci
          result.table[i,"convergence"] = "converged"
        }else
          result.table[i,"convergence"] = "converged, profile based CIs did not work"
      }else
        result.table[i,"convergence"] = "Exposure not in model"
      
    }
    if(class(my.data[[i]])[1]=="list" & class(my.data[[i]][[1]])[1] %in% c("glm","clogit")){
      if(term.of.interest %in% rownames(summary(my.data[[i]][[1]])$coefficients)){
        
        help = summary(my.data[[i]][[1]])$coefficients[term.of.interest,]
        if("exp(coef)" %in% names(help))
          help = help[-1*which(names(help) == "exp(coef)")]
        
        result.table[i,c("estimate","std.error", "statistic","p.value" )] = help
          
        #result.table[i,"p.value.anova"] = my.data[[i]][["anova"]]$`Pr(>Chi)`[2]
        result.table[i,"p.value"] = my.data[[i]][["anova"]][2,dim(my.data[[i]][["anova"]])[2]]
        
        robust.ci = robustConfInt(my.data[[i]][[1]], parm = "exposure")
        
        if(class(robust.ci)[1] == "numeric"){
          result.table[i,c("l.ci","u.ci")] = robust.ci
          result.table[i,"convergence"] = "converged, run ANOVA"
        }else
          result.table[i,"convergence"] = "converged, run ANOVA, profile based CIs did not converge"
      }else
        result.table[i,"convergence"] = "Exposure not in model"
    } 
    if(class(my.data[[i]])[1]=="simpleWarning" | 
       (class(my.data[[i]])[1]=="list" & class(my.data[[i]][[1]])[1] == "simpleWarning")){
      result.table[i,"convergence"] = my.data[[i]]$message
    }
    if(class(my.data[[i]])[1]=="simpleError" |
       (class(my.data[[i]])[1]=="list" & class(my.data[[i]][[1]])[1] == "simpleError")){
      result.table[i,"convergence"] = my.data[[i]]$message
    }
  }
  
  result.table$p.value.anova = NULL
  if(debug) print("End createResultTable")
  return(result.table) 
}


runRegression = function(regression.type = "unconditional.logistic.regression", 
                         my.eset = NULL, outcome = NULL, covariates = NULL, 
                         strata = NULL, allow.na.values.in.exposure = F, 
                         file.with.result.prefix, store.complete.result = T, store.result = T, debug = F){
  
  # this function works on Expression Sets (ESet) only
  # will run the specified regression on all measured features in the ESet
  # all data have to be transfomed and imputed
  # allow.na.values.in.exposure = F : if a metabolite contains NA values it will be used as an 
  # indicator metabolite (Pete Kraft approach)
  #
  #### INPUT parameter
  # regression type = currently can be: unconditional.logistic.regression
  #                   or conditional.logistic.regression
  # my.eset = ExpressionSet with data
  # outcome = NULL: column name containing the outcome variable in pData of ESet
  # covariates = NULL: column names containing the covariates in pData of ESet
  # strata = NULL: column name containing the strata variable in pData of ESet
  # allow.na.values.in.exposure = F: what to do with missing values. 
  #                                  If False use Pete Kraft's approach
  # file.with.result.prefix: prefix of file with results. 
  #
  #### RETURN values
  # result.data = table with estimates, CI and p-values
  # saves different version of the results
  #
  #### 
  
  if(debug) print("Start")
  all.regression.types =  c("unconditional.logistic.regression", "conditional.logistic.regression")
  
  if(is.null(my.eset)|is.null(outcome)){
    message("ERROR: The parameters my.eset and/or outcome cannot empty\n")
    return(0)
  }
  
  if(!regression.type %in% all.regression.types){
    message(paste("ERROR: Supported regression types are: ", all.regression.types,".\n"))
    return(0)
  }
  
  if(regression.type == "conditional.logistic.regression" & is.null(strata)){
    message("Please specify the strata to use in your conditional logistic regression\n")
    return(0)
  }
  
  if(debug) print("Init tests ok")
  regression.data = createDatasetForRegression(my.eset=my.eset, outcome = outcome, 
                                               covariates = covariates, strata=strata, debug = debug)
  
  cat("Running ", regression.type,"on ", dim(exprs(my.eset))[1],"exposures\n")
  cat("Covariates: ", covariates,"\n")
  cat("Outcome info: \n")
  print(table(regression.data[,outcome]))
  
  print(c("The following variables are treated as factors:", paste(colnames(regression.data)[unlist(lapply(regression.data,is.factor))], collapse = ", ")))
  
  # use switch to be able to add other types of regressions
  result = switch (regression.type,
    unconditional.logistic.regression = 
      runUnconditionalLogisticRegression(all.data = regression.data, outcome = outcome, 
                                         covariates = covariates, 
                                         allow.na.values.in.exposure = allow.na.values.in.exposure, debug = debug),
    conditional.logistic.regression = 
      runConditionalLogisticRegression(all.data = regression.data, outcome = outcome, 
                                       covariates = covariates, strata = strata, 
                                       allow.na.values.in.exposure = allow.na.values.in.exposure, debug = debug),
    
    "No matching regression type specified"
  )
  if(store.complete.result){
    cat("Saving all results from the regressions to: ", paste(file.with.result.prefix,"_all.RData",
                                                              sep = ""),"\n")
    save(result, file = paste(file.with.result.prefix,"_all.RData",sep = ""))
  
  }
  
  
  result.table = createResultTable(my.data = result, term.of.interest = "exposure", debug = debug)
  met.annot = fData(my.eset)
  met.annot$id = rownames(met.annot)
  result.data = cbind(met.annot[,c("METABOLITE","HMDB_ID","is.indicator")],result.table)
  
  if(store.result){
    cat("Saving table with raw main results from the regressions to: ", 
        paste(file.with.result.prefix,"_table_raw.RData",sep=""),"\n")
    save(result.data, file = paste(file.with.result.prefix,"_table_raw.RData",sep=""))
  }
  return(result.data)
}

calculateORCI = function(result.data=NULL, or.scale = 1, file.with.result.prefix=NULL,
                         estimate = "estimate", l.ci = "l.ci", u.ci = "u.ci", exp.ci = F, 
                         add.qval = F, estimate.CI = F){
  
  # this function calculates OR anc CI for OR based on or.scale\
  # it needs a data frame with the estimates, and raw CI
  #
  #### INPUT values
  # result.data = data frame with the estimates, and raw CI
  # or.scale = 1, how to scale OR. 1 is equal to OR per 1 SD increase in exposure values
  # file.with.result.prefix = prefix of file to strore results
  #
  #### RETURN values
  # or.print: data frame with table ready for manuscript
  # stores different versions of the data as RData/csv file
  # _all.RData: complete regression results with all avaliable details as a list
  # _table_raw.RData: a table with main results such as: HMDB IDs, Metabolites, estimates, 
  #                   standard errors, p-values and test statistic
  # _table_OR.RData
  #_manuscript_table.csv
  # 
 
  
  if(is.null(result.data))
    cat("ERROR: result data can't be empty\n")
  
  if(is.null(file.with.result.prefix))
    cat("ERROR: file.with.result.prefix can't be empty\n")
  
  or.data = result.data
  or.data$OR = exp(or.data[,estimate]*or.scale)
  
  if(estimate.CI == T){
    
      if(any(is.na(or.data[,"l.ci"])) | any(is.na(or.data[,"u.ci"]))){
        or.data[,l.ci] = or.data[,"estimate"]-1.96*or.data$std.error
        or.data[,u.ci] = or.data[,"estimate"]+1.96*or.data$std.error
      }
    
    if(exp.ci == F){
      or.data[,l.ci] = exp(or.data[,l.ci])
      or.data[,u.ci] = exp(or.data[,u.ci])
    }
  }
  
  
  if(exp.ci){
    or.data$LCI.OR = exp(or.data[,l.ci]*or.scale)
    or.data$UCI.OR = exp(or.data[,u.ci]*or.scale)
  }else{
    or.data$LCI.OR = exp(log(or.data[,l.ci])*or.scale)
    or.data$UCI.OR = exp(log(or.data[,u.ci])*or.scale)
  }
    
  cat("Saving a table with OR main results from the regressions to: ", 
      paste(file.with.result.prefix,"_table_OR_raw.RData",sep=""),"\n")
  save(or.data, file = paste(file.with.result.prefix,"_table_OR_raw.RData",sep=""))
   
  or.data$OR = round(or.data$OR,2)
  or.data$LCI.OR = round(or.data$LCI.OR,2)
  or.data$UCI.OR = round(or.data$UCI.OR,2)
  or.data$OR.CI = paste(or.data$OR," (",or.data$LCI.OR,"-",or.data$UCI.OR,")",sep = "")
  
  
  if(add.qval == T){
    or.data = or.data[order(or.data$p.value, decreasing = F),]
    help = mt.rawp2adjp(rawp = or.data$p.value, proc = "TSBH", alpha = 0.05)
    or.data$q.value = round(help$adjp[,"TSBH_0.05"],3)
  }

  or.data$p.value = round(or.data$p.value,3)
  
  sel.cols = c("HMDB_ID","METABOLITE","OR.CI","p.value","convergence")
  if("adj.pval" %in% colnames(or.data))
    sel.cols = c(sel.cols, "adj.pval")
  if("q.value" %in% colnames(or.data))
    sel.cols = c(sel.cols, "q.value")
  
  or.print = or.data[,sel.cols]
  or.print = or.print[order(or.print$p.value, decreasing = F),]
  
  if("adj.pval" %in% colnames(or.print))
    or.print = or.print[order(or.print$adj.pval, or.print$p.value, decreasing = F),]
  
  or.print$convergence[which(or.print$convergence == "converged, run ANOVA")] = "converged"
  
  cat("Saving the results table for the manuscript as RData: ", 
      paste(file.with.result.prefix,"_table_OR.RData",sep=""),"\n")
  save(or.print, file = paste(file.with.result.prefix,"_table_OR.RData",sep=""))
  
  cat("Saving the results table for the manuscript as csv file: ", 
      paste(file.with.result.prefix,"_manuscript_table.csv",sep=""),"\n")
  write.table(x = or.print, file = paste(file.with.result.prefix,"_manuscript_table.csv",sep=""),
              col.names = T, row.names = F)
  return(or.print)
}


permuteCACO = function(caco.data, outcome, permute.strata){
   
  pdata = pData(caco.data)
  caco.data.new = pdata
  caco.data.new$caco.new = caco.data.new[,outcome]
  
  caco.data.new = caco.data.new[order(caco.data.new[,permute.strata]),]
  caco.data.new$caco.new = unlist(tapply(X = caco.data.new[,outcome],
                                         INDEX = factor(caco.data.new[,permute.strata]),
                                         FUN = function(y) {z = y[sample(x = length(y))]; return(z)}))
  
  # set to initial order
  caco.data.new = caco.data.new[rownames(pData(caco.data)),] 
  caco.data.new[,outcome] = caco.data.new$caco.new
  caco.data.new$caco.new = NULL
  pData(caco.data) = caco.data.new
    
  return(caco.data)
}


runRegressionWithPermutations = function(regression.type = "unconditional.logistic.regression", 
                                         my.eset = NULL, outcome = NULL, covariates = NULL, 
                                         strata = NULL, allow.na.values.in.exposure = F, 
                                         file.with.result.prefix, num.perm = 10, permute.strata){
   
  cat("Running ", num.perm, "permutations ...\n")
  # perm.file will allways be the same as it stores the whole models which is just too much data
  perm.file = paste(file.with.result.prefix,"_perm", sep = "")
  perm.res = vector("list",num.perm)
  for(i in 1:num.perm){
    cat("============= i = ",i,"\n")
    my.eset = permuteCACO(caco.data = my.eset,outcome = outcome, permute.strata = permute.strata)
    
    perm.res[[i]] = runRegression(regression.type = regression.type, my.eset = my.eset, 
                                  outcome = outcome, covariates = covariates, 
                                  strata = strata, 
                                  allow.na.values.in.exposure = allow.na.values.in.exposure, 
                                  file.with.result.prefix, store.complete.result = F,
                                  store.result = F)
  }
  return(perm.res)
}


mergePermutationBatchesPvalues = function(batch.files, num.batches){
    
  #init result
  cur.file = gsub(pattern = "_B",replacement = "_1",x = batch.files, fixed = T)
  
  help = load(cur.file)
  cur.perm = get(help)
  perm.all = do.call(what = rbind, args = lapply(cur.perm, function(x) return(x$p.value)))
  colnames(perm.all) = cur.perm[[1]]$id
  
  for(i in 2:num.batches){
    cur.file = gsub(pattern = "_B",replacement = paste("_",i,sep = ""),x = batch.files, fixed = T)
    help = load(cur.file)
    cur.perm = get(help)
    
    perm.help = do.call(what = rbind, args = lapply(cur.perm, function(x) return(x$p.value)))
    colnames(perm.help) = cur.perm[[1]]$id
    common.ids = intersect(colnames(perm.help), colnames(perm.all))
    perm.all = rbind(perm.all[,common.ids], perm.help[,common.ids])
  }
  
  to.print = paste("Permutation pvalues dimension:", paste(dim(perm.all)))
  print(to.print)
  return(perm.all)
}

mergePermutationBatchesEstimates = function(batch.files, num.batches){
   
  #init result
  cur.file = gsub(pattern = "B",replacement = "1",x = batch.files, fixed = T)
  help = load(cur.file)
  cur.perm = get(help)
  perm.all = do.call(what = rbind, args = lapply(cur.perm, function(x) return(x$estimate)))
  
  for(i in 2:num.batches){
    cur.file = gsub(pattern = "B",replacement = i,x = batch.files, fixed = T)
    help = load(cur.file)
    cur.perm = get(help)
    perm.help = do.call(what = rbind, args = lapply(cur.perm, function(x) return(x$estimate)))
    perm.all = rbind(perm.all, perm.help)
  }
  
  colnames(perm.all) = cur.perm[[1]]$id
 
  return(perm.all)
}

calcPermPvalGlobalTest = function(perm.data, obs.data){
    
  num.nas = apply(perm.data,2, function(x) {length(which(is.na(x)))})
  table(num.nas)
  
  # change order to match observed data
  perm.data = perm.data[,obs.data$id]
    
  min.pval = apply(X = perm.data, MARGIN = 1, FUN = min)
  print("Quantiles of minimal p-values across permutations:")
  print(round(quantile(x = min.pval, probs = seq(0.05,0.5,0.05)),5))
  print("Summary of minimal p-values in each distribution:")
  print(summary(min.pval))
  
  print("observed minimal p-value:")
  print(min(obs.data$p.value))
  
  # The permutation p-value for test of the overall null 
  # (no metabolite is associated with outcome/ovarian cancer) was estimated as k/(1,001), where k is the 
  # number of permutations where the smallest p-value was smaller than the smallest observed p-value
  overall.perm.pval = length(which(min.pval<min(obs.data$p.value)))/(length(min.pval)+1)
  print("Global test permuation p-value:")
  print(overall.perm.pval)
  return(overall.perm.pval)
}

calcCorrectedPvaluesMinP = function(perm.data, obs.data){
   
  require(NPC)
  perm.help = obs.data
  perm.help = perm.help[order(perm.help$p.value,decreasing = F),]
  
  perm.data = perm.data[,perm.help$id]
  
  all.pv = rbind(perm.help$p.value, perm.data)
  colnames(all.pv) = perm.help$id
  dim(all.pv)
  all.pv[1:5, 1:10]
  
  adj.pv = FWE(PV = all.pv, stepdown = T)
  adj.df = data.frame(adj.pval = round(adj.pv,3), id = rownames(adj.pv))
  
  res.data = merge(x = obs.data, y = adj.df, by.x = "id",by.y = "id") 
  res.data = res.data[order(res.data$adj.pval),]
  
  print("Summary of permutation adjusted p-values:")
  print(summary(res.data$adj.pval))
  return(res.data)
}

evaluatePermutationResults = function(batch.files, num.batches, obs.data, result.file, 
                                      exp.ci = T, or.scale){
  # calculates global test pvalue
  # calculates permuatations adjusted p-values based on Westfall and Young
  # adds adj.pval to results table and writes table to csv file
  # obs.data has to be raw data where p-values have not been rounded 
  
  perm.pvals = mergePermutationBatchesPvalues(batch.files = batch.files,num.batches = num.batches)
  table(apply(perm.pvals, 2, function(x) length(which(is.na(x)))))
  # setting NA to 1 (not converged metabolites pvalue to 1)
  perm.pvals[which(is.na(perm.pvals), arr.ind = T)] = 1
  
  # reorder to match result table
  perm.pvals = perm.pvals[,obs.data$id]
  perm.mets = obs.data$id
  dim(perm.pvals)
  
  # subset to metabolites that converged in the main analysis
  rownames(obs.data) = obs.data$id
  table(obs.data$convergence)
  perm.mets = obs.data$id[which(obs.data$convergence %in% c("converged", "converged, run ANOVA",
                                                            "converged, profile based CIs did not work",
                                                            "converged, run ANOVA, profile based CIs did not converge"))]
  perm.pvals.valid = perm.pvals[,perm.mets]
  dim(perm.pvals.valid)
  table(apply(perm.pvals.valid, 2, function(x) length(which(is.na(x)))))

  overall.perm.pval = calcPermPvalGlobalTest(perm.data = perm.pvals.valid, 
                                             obs.data = obs.data[perm.mets,])
  adj.pvals = calcCorrectedPvaluesMinP(perm.data = perm.pvals.valid,
                                       obs.data = obs.data[perm.mets,])
  
  all.res = adj.pvals
  
  help.pvals = obs.data[which(!obs.data$id %in% perm.mets),]
  if(dim(help.pvals)[1]>=1){
    help.pvals$adj.pval = NA
    all.res = rbind(adj.pvals,help.pvals)
  }
    
  or.res = calculateORCI(result.data = all.res, or.scale = or.scale, file.with.result.prefix = result.file, 
                         estimate = "estimate", l.ci = "l.ci",u.ci = "u.ci", exp.ci = exp.ci)
  
  return(or.res)
}
