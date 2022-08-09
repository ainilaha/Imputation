##################################################################################
################################ENWAS ############################################
##################################################################################

###your parameters###
##Note: particpant IDs should be called "ID", age and sex should be called, "age", "sex"

##1. data: the dataset with the outcome, model covariates, and the post-stratification variable;
##2. ph.list: the covariate list;
##3. cachepath: your cache path;
##4. formula: model covariates;
##5. prop: post-stratification statistics;
##6. post_stratified_cat: post-stratification variable categories; (default variable name: ps_cat)
##7. N_jobs: numbers of parallel jobs you would like to run at the same time;
##8. outcome: your model outcome
##9. sex_specific: sex specific run
library(dplyr)
library(rlang)
library(survival)
library(survey)
library(sandwich)
library(splines)
library(Rsge)
library(purrr)

ph.list = c('num1',"num2","bin")
model_var = c('age',"sex",'pc.0')
data = data_enwas
outcome = "outcome1"
gen_cat = "age_sex"
sex_spec_cat = "age_cat"
ss.out = "No"
ss.out.cat = "M/F"
n_cutoff= 10 

i = 1

ENWAS_general=function(i,
                       ph.list, # list of variables to run through as independent variables
                       data, # dataframe with ID, outcome,age (age centered at 50),exposures, poststratification variable
                       model_var, #character vector, list of variable to include in model
                       gen_cat, #name of post stratification variable
                       outcome,#name of outcome variable
                       sex_spec_cat, # if sex specific, name of post-stratification variable
                       ss.out, # "Yes" or "No", is outcome sex-specific
                       ss.out.cat, # "M/F", "M", "F"
                       n_cutoff # minimum sample size required to run model.
                       ){
  
  keep = TRUE
  
  #placeholder for full model results
  
  full_mod = read.csv(text = c("outcome,phenotype,model_var,Estimate,SE,t_stat,p_value,r2,aic,ll,rank"))
  full_mod[1,] = c(outcome, ph.list[i], rep(NA, ncol(full_mod)-2))

  
  # check if outcome is non-numeric and smallest bucket is less than 100, do not run ENWAS
  if(!is.numeric(data[,outcome])){
    cc_count =  data %>% count(!!sym(outcome))
  if(nrow(cc_count)<2 | (cc_count %>% dplyr::select(n) %>% min(.$n))<100){
    keep = FALSE #set keep = FALSE so that enwas is not run
  }
  }
  
 if(keep==FALSE){ 
  
   # stop analysis for outcomes with total cases <100
   
   result= data.frame(
     outcome = outcome,
     phenotype = ph.list[i],type = NA, N = NA, SD = NA, rsq = NA,aic = NA, ll = NA, rank = NA, cc_pheno1  = NA,  cc_pheno2 = NA, 
     post_stratified_beta=NA,original_se = NA,original_pvalue =NA,
     san_se=NA,san_pvalue=NA,min_weight=NA,mean_weight=NA,max_weight=NA,
     abs_post_beta_SD=NA, sex_specific = NA, sex_cat = NA , error = "Case/Control count <100")
 
   }else{

   # continue on with enwas analysis

   # model form is specified. 
   # currently model form includes a spline term for the age variable. age = current age -50. 
   # knots are at every 10 years from -20 (30) to 20 (70), with boundaries at -30 (20) and 40 (90) 
   # This would need to be adjusted if another age variable is used. 
   # Adjustments should be made on lines 104,115 and 121.
     
  if(ss.out == 'Yes'){ # if outcome is sex-specific
    post_strat = sex_spec_cat
    mod = c("bs(age, knots = seq(-20,20, by = 10), Boundary.knots = c(-30,40))",model_var[which(!model_var %in% c('age','sex'))])
    sex_specific = "Yes"
    sex_cat = ss.out.cat
    X1 = data[,c("ID",ph.list)] # # here it gets dicey - would need to think about structure of users database
      #get.pheno.data(names = c(ph.list[i])) %>% mutate(ID = rownames(.)) %>% dplyr::select(ID, everything()) 
  }else{
    # if not sex specific, pull phenotype and sex, do checks that we have enough of each variable. 
    X=get.pheno.data(names = c(ph.list[i], 'sex'))
    minN = X %>% filter(complete.cases(.)) %>% group_by(sex) %>% count()
  
  if(min(minN$n)<=20|nrow(minN)==1){
    post_strat = sex_spec_cat
    data = data %>% filter(sex==minN$sex[which(minN$n==max(minN$n))])
    mod = c("bs(age, knots = seq(-20,20, by = 10), Boundary.knots = c(-30,40))",model_var[which(!model_var %in% c('age','sex'))])

    sex_specific = "Yes"
    sex_cat = minN$sex[which(minN$n==max(minN$n))]
  }else{
    post_strat = gen_cat
    mod = c("bs(age, knots = seq(-20,20, by = 10), Boundary.knots = c(-30,40))*sex",model_var[which(!model_var %in% c('age','sex'))])
    sex_specific = "No"
    sex_cat = "M/F"
    }
  
  # dataset with predictor
  
  X1= X %>% mutate(ID = rownames(.)) %>% dplyr::select(ID, everything(), -sex)
  }
  
  # define the proportion for post-stratification  
  colnames(data)[which(colnames(data)==post_strat)] = 'ps.cat'
  
  prop = data %>% count(ps.cat) %>% mutate(prop = prop.table(n))
  

  # merge dataset of predictor with the reference dataset
  
  Y = data %>% left_join(X1, by = "ID") %>% filter(complete.cases(.))
  if(nrow(Y)==0){
    result= data.frame(
      outcome = outcome,
      phenotype = ph.list[i],
      type = NA, 
      N = NA, SD = NA, rsq = NA,
      aic = NA, ll = NA, rank = NA, 
      cc_pheno1  = NA,  cc_pheno2 = NA, 
      post_stratified_beta=NA,original_se = NA,original_pvalue =NA,
      san_se=NA,san_pvalue=NA,min_weight=NA,mean_weight=NA,max_weight=NA,
      abs_post_beta_SD=NA, sex_specific, sex_cat , error = "Case/Control count <100")
  }else{
  prop.sample = Y %>% count(ps.cat) %>%mutate(prop = prop.table(n))
  
  colnames(Y)[length(Y)]="Predictor"
  Y$Predictor=as.numeric(Y$Predictor)
  
  ## pull counts of 'controls' vs "cases", basically just take total - number at lowest level, or means
  
  if((!is.numeric(Y[,outcome]))){
    if(length(unique(Y[,outcome]))==2){class_outcome = 'binary'}else{class_outcome = "categorical"}
    cc_pheno1 = nrow(Y) -  (Y %>% count(!!sym(outcome)) %>%slice(1) %>% .$n)
    cc_pheno2 = cc_pheno1/nrow(Y)
  }else{
    class_outcome = "numeric"
    cc_pheno1 = mean(Y[,outcome],na.rm = TRUE)
    cc_pheno2 = sd(Y[,outcome],na.rm = TRUE)
  }
  
  # change factors that are non binary - to be numeric (for now)
  # maybe eventually we think of a better model here. 
  if((!is.numeric(Y[,outcome]) & length(unique(Y[,outcome]))>2)){
    Y[,outcome]=as.numeric(Y[,outcome])
  }   
  
  #SD
  SD=sd(Y$Predictor,na.rm = T)
  
  #N
  N = nrow(Y)
  
  # Cases for non numeric outcomes
  keep = TRUE
  if(!is.numeric(Y[,outcome])){
   cc_count =  Y %>% count(!!sym(outcome))
   if(nrow(cc_count)<2 | (cc_count %>% dplyr::select(n) %>% min(.$n))<100){
     keep = FALSE
   }
  }
  


  if(N<n_cutoff){
    result= data.frame(
      outcome = outcome,
      phenotype = ph.list[i],type = class_outcome,N, SD = NA, rsq = NA, 
      aic = NA, ll = NA, rank = NA, 
      cc_pheno1  = NA,  cc_pheno2 = NA, 
      post_stratified_beta=NA,original_se = NA,original_pvalue =NA,
      san_se=NA,san_pvalue=NA,min_weight=NA,mean_weight=NA,max_weight=NA,
      abs_post_beta_SD=NA, sex_specific, sex_cat , error = paste0("N<",n_cutoff))
  }else{
    if(keep==FALSE){
      # stop analysis for outcomes with outcome/phenotype cases <100
      result= data.frame(
        outcome = outcome, phenotype = ph.list[i],type = class_outcome,N, SD = NA,rsq = NA, aic = NA, ll = NA, rank = NA, cc_pheno1  = NA,  cc_pheno2 = NA,
        post_stratified_beta=NA,original_se = NA,original_pvalue =NA,
        san_se=NA,san_pvalue=NA,min_weight=NA,mean_weight=NA,max_weight=NA,
        abs_post_beta_SD=NA, sex_specific, sex_cat , error = "Case/Control count <100")
    }else{
      
## start running models
  form = reformulate(termlabels = c(mod,"Predictor"), response = outcome)

  #Post-stratification 
  ps.dist <- data.frame(ps.cat = prop$ps.cat, 
                        #Freq = prop$n)
                        Freq = nrow(Y) * prop$prop)
  
  dat.unweighted <- svydesign(ids=~ID, data=Y)
  
  if('try-error' %in% class(try(rake(
    design = dat.unweighted,sample.margins = list(~ps.cat),population.margins = list(ps.dist)),silent = TRUE))){
    result= data.frame(
      outcome = outcome,
      phenotype = ph.list[i],type = class_outcome,N, SD = NA, rsq = NA,aic = NA, ll = NA, rank = NA,cc_pheno1  = NA,  cc_pheno2 = NA,
      post_stratified_beta=NA,original_se = NA,original_pvalue =NA,
      san_se=NA,san_pvalue=NA,min_weight=NA,mean_weight=NA,max_weight=NA,
      abs_post_beta_SD=NA, sex_specific, sex_cat ,error = "Post-stratification error: rake failure")
  }else{
    dat.svy.rake.unstrat = rake(design = dat.unweighted,
                                sample.margins = list(~ps.cat),
                                population.margins = list(ps.dist))
    
    mean_weight = mean(weights(dat.svy.rake.unstrat))
    max_weight = max(weights(dat.svy.rake.unstrat))
    min_weight = min(weights(dat.svy.rake.unstrat))
    upper_limit = 5*mean_weight
    dat.svy.unstrat.trim <- trimWeights(dat.svy.rake.unstrat, upper=upper_limit, strict=TRUE) 
    
    ## if not numeric - run logistic regression
    if (is.numeric(Y[,outcome])==FALSE) { 
      
    f1 =try(glm(form, data=Y,family=binomial(), weights = 1/dat.svy.unstrat.trim[["prob"]]), silent = TRUE)
    
      
    ## now continue with f1
    
    if(!('try-error' %in% class(f1))){

      post_stratified_beta=try(coef(summary(f1))["Predictor","Estimate"], silent = TRUE)
      if(class(post_stratified_beta)=='try-error'){
        post_stratified_beta = NA
      }
      if(!is.na(post_stratified_beta)){
        original_se=coef(summary(f1))["Predictor","Std. Error"]
        original_pvalue=coef(summary(f1))["Predictor","Pr(>|z|)"]
        r2 = 1 - f1$deviance/f1$null.deviance
        aic = f1$aic
        ll = as.numeric(logLik(f1))
        rank = f1$rank
        full_mod = data.frame(
          outcome = outcome, 
          phenotype = ph.list[i],
          as.data.frame(summary(f1)$coefficients) %>% mutate(model_var = row.names(.)), 
          r2 = r2,
          aic = aic,
          ll = ll,
          rank = rank
        ) %>% select(outcome, phenotype, model_var,everything())
        colnames(full_mod) = c("outcome", "phenotype",   "model_var",  "Estimate",   "SE", "t_stat",    "p_value",   "r2",'aic','ll','rank')

      }else{
        f1 = NA
        full_mod = read.csv(text = c("outcome,phenotype,model_var,Estimate,SE,t_stat,p_value,r2,aic,ll,rank"))
        full_mod[1,] = c(outcome, ph.list[i], rep(NA, ncol(full_mod)-2))
       
        }
      }else{
        f1 = NA
        full_mod = read.csv(text = c("outcome,phenotype,model_var,Estimate,SE,t_stat,p_value,r2,aic,ll,rank"))
        full_mod[1,] = c(outcome, ph.list[i], rep(NA, ncol(full_mod)-2))
      }}else {
        # if outcome is numeric
          f1 = try(glm(form, family = gaussian(), data=Y, weights = 1/dat.svy.unstrat.trim[["prob"]])  , silent = TRUE)
          
          if(!('try-error' %in% class(f1))){
           
            post_stratified_beta=try(coef(summary(f1))["Predictor","Estimate"], silent = TRUE)
            if(class(post_stratified_beta)=='try-error'){
              post_stratified_beta = NA
            }
            
          if(!is.na(post_stratified_beta)){
          original_se=coef(summary(f1))["Predictor","Std. Error"]
          original_pvalue=coef(summary(f1))["Predictor","Pr(>|t|)"]
          r2 = 1 - f1$deviance/f1$null.deviance
          aic = f1$aic
          ll = as.numeric(logLik(f1))
          rank = f1$rank
          full_mod = data.frame(
            outcome = outcome, 
            phenotype = ph.list[i],
            as.data.frame(summary(f1)$coefficients) %>% mutate(model_var = row.names(.)), 
            r2 = r2,
            aic = aic, 
            ll = ll, 
            rank = rank
          ) %>% select(outcome, phenotype, model_var,everything())
          colnames(full_mod) = c("outcome", "phenotype",   "model_var",  "Estimate",   "SE", "t_stat",    "p_value",   "r2",'aic','ll','rank')
          
          }else{
            f1 = NA
            full_mod = read.csv(text = c("outcome,phenotype,model_var,Estimate,SE,t_stat,p_value,r2,aic,ll,rank"))
            full_mod[1,] = c(outcome, ph.list[i], rep(NA, ncol(full_mod)-2))

          }
        }else{
          f1 = NA
          full_mod = read.csv(text = c("outcome,phenotype,model_var,Estimate,SE,t_stat,p_value,r2,aic,ll,rank"))
          full_mod[1,] = c(outcome, ph.list[i], rep(NA, ncol(full_mod)-2))

          }
      }
  
  
  if(TRUE %in% is.na(f1)){
    result= data.frame(
      outcome = outcome, 
      phenotype = ph.list[i],type = class_outcome,N, SD = NA,rsq = NA,aic = NA, ll = NA, rank = NA,cc_pheno1  = NA,  cc_pheno2 = NA,
      post_stratified_beta=NA,original_se = NA,original_pvalue =NA,
      san_se=NA,san_pvalue=NA,min_weight=NA,mean_weight=NA,max_weight=NA,
      abs_post_beta_SD=NA, sex_specific, sex_cat ,error = "Model failed to converge")
    }else{
  #Sandwich_var
  sandwich <- as.data.frame(diag(vcovHC(f1, type = "HC")))
  colnames(sandwich)="sandwich_var"
  san_se=sqrt(sandwich$sandwich_var[which(rownames(sandwich)=="Predictor")])
  stat <- post_stratified_beta/san_se
  san_pvalue <- pchisq(q = stat^2, df = 1, lower.tail=FALSE)
        

  
  #Effect size per SD
  abs_post_beta_SD=abs(post_stratified_beta)*SD
  
  result=data.frame(outcome = outcome, phenotype = ph.list[i],type = class_outcome, N, SD, rsq = r2, aic = aic, ll = ll, rank = rank,
                    cc_pheno1  = cc_pheno1,  cc_pheno2 = cc_pheno2, post_stratified_beta,original_se,original_pvalue,
                    san_se, san_pvalue, min_weight, mean_weight, max_weight,
                    abs_post_beta_SD, sex_specific, sex_cat,error = "None")
    }
    }
  }
  }
  }
 }
  list(result = result,full_mod = full_mod)

}



## Make fake ENWAS data
N = 100000
set.seed(234729374)
data_enwas = 
  data.frame(
    ID = 1:N,
    outcome1 = rnorm(n = N,mean = 10, sd = 3),
    outcome2 = rbinom(n = N, size = 1,prob = 0.2),
    age = rnorm(n = N,mean = 60, sd = 10),
    sex = sample(c('M','F'),size = N,replace = TRUE),
    pc.0 = rnorm(n = N, mean = 0, sd = 1),
    num1 = rnorm(n = N, mean = 0, sd = 1),
    num2 = rnorm(n = N, mean = 0, sd = 1),
    bin1 = rbinom(n = N, size = 1,prob = 0.2)
  ) %>% filter(age>20, age<90) %>% 
  mutate(
    age_cat = cut(age,breaks = seq(20, 90,by = 10),include.lowest = TRUE),
    age_sex = paste0(sex,"-",age_cat),
    age_50 = age-50
  )








i=1
ph.list=exposure_vars # list of variables to run through as independent variables
data # dataframe with ID, outcome,age (age centered at 50),exposures, poststratification variable
model_var=c("BMXWT","BMXHT") #character vector, list of variable to include in model********，and just variable?
gen_cat="SDMVSTRA" #name of post stratification variable
outcome="BMXWAIST" #name of outcome variable
sex_spec_cat="RIAGENDR" # if sex specific, name of post-stratification variable
ss.out="Yes" # "Yes" or "No", is outcome sex-specific
ss.out.cat="M/F" # "M/F", "M", "F"
n_cutoff=10 # minimum sample size required to run model.
