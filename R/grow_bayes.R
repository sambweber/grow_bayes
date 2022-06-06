require('R2jags')
require('future')
require('tidyverse')
require('coda')

# ------------------------------------------------------------------------------------
# Growth model - used internally by grow_bayes to build the growth model
# ------------------------------------------------------------------------------------

# Returns the specified growth model equation with group indexing for fitting models with different parameters for subgroups.

growth_model = function(model = c('VB','logistic','Gompertz','Richards','cessation','VBlogK')){
  
  model <- match.arg(model)
  f = switch(model,
             Gompertz  = "Linf[group[i]] * pow(L0[group[i]]/Linf[group[i]],exp(-exp(1) * k[group[i]] * age[i]))",
             VB =,logistic =, Richards  = "Linf[group[i]] * pow(1+(pow((L0[group[i]]/Linf[group[i]]),1-delta[group[i]])-1) * exp(-k[group[i]] * age[i] / pow(delta[group[i]],delta[group[i]]/(1-delta[group[i]]))),1/(1-delta[group[i]]))",
             cessation = "L0[group[i]] + rmax[group[i]] * (((log(exp(-k[group[i]] * t50[group[i]]) + 1) - log(exp(k[group[i]]*(age[i]-t50[group[i]]))+1))/k[group[i]]) + age[i])",
             VBlogK    = "Linf[group[i]] * (1 - exp(-k2[group[i]]*(age[i] - t0[group[i]])) * pow((1+exp(-beta[group[i]]*(age[i]-t0[group[i]]-alpha[group[i]])))/(1+exp(beta[group[i]]*alpha[group[i]])),(k1[group[i]]-k2[group[i]])/beta[group[i]]))"
             )
  if(model=='VB')       f = gsub('delta\\[group\\[i\\]\\]','2/3',f)
  if(model=='logistic') f = gsub('delta\\[group\\[i\\]\\]','2',f)
  f
}

# ----------------------------------------------------------------------------------
# Error model - used internally by grow_bayes to build the growth model
# ----------------------------------------------------------------------------------

# Returns the error model relating expecting length to observed length

error_model = function(model = c('normal','lognormal')){
  model = match.arg(model)
  switch(model,
         normal = 'length[i] ~ dnorm(Lt[i], tau.y)',
         lognormal = 'Log_Lt[i] = log(max(Lt[i],1))\n       
                     length[i] ~ dlnorm(Log_Lt[i], tau.y)'
  )
}


# ------------------------------------------------------------------------------------
# Priors  - used internally by grow_bayes to build the growth model
# ------------------------------------------------------------------------------------

# Returns default priors for the parameters of the specific growth function in BUGS language

gm_priors = function(model=c('VB','logistic','Gompertz','Richards','cessation','VBlogK')){
  
  model <- match.arg(model)
  
  switch(model,
         
         VB =, logistic =, Gompertz =, Richards  = 
           
           {cat("
     Linf[j] ~ dnorm(0,0.0001) I(L0[j],)
     k[j] ~ dgamma(0.001,0.001)
     L0[j] ~ dnorm(0,0.0001) I(0,)
     ")
             if(model=='Richards'){
               cat("delta[j] ~ dgamma(0.001,0.001) # Must be >0
      ")}
           },
         
         cessation = 
           
           cat(" 
      L0[j] ~ dnorm(0,0.0001) I(0,)
      k[j]  ~ dnorm(0,0.0001) I(0,)
      rmax[j] ~ dnorm(0,0.0001)
      t50[j]  ~ dnorm(0,0.0001) I(0,)
      "),
         
         VBlogK =
           
           # Prior distributions based on those provided in Dortel REF for Indian Ocean  
           cat("
     Linf[j] ~ dunif(100,10000)
     k1[j]  ~ dgamma(2.78,4.74) 
     k2[j]  = k1[j] + kappa[j]
     kappa[j] ~ dunif(0,4)
     t0[j]  ~ dunif(-2,0)
     alpha[j]  ~ dgamma(4,1.38)
     beta[j]   ~ dunif(0,30)
     ")
         
)
}

# -------------------------------------------------------------------------------------
# grow_bayes: main fitting function
# -------------------------------------------------------------------------------------

#' This is the principle function used for fitting growth curves. It can fit either single or multiple 
#' growth models, with the option of fitting in parallel to speed up processing times. If a single
#' model is specified the function will return the output of the JAGS run. If multiple models are 
#' specified the function will return a nested tibble containing a list of JAGS model ouputs and will 
#' rank them based on Deviance Information Criterion or Root Mean Square Error.

#' @params size,age Numeric vectors containing sizes and associated ages
#' @params group Optional: a vector of group labels to be modelled using separate parameters. If not supplied
#' a single growth curve will be fit to the whole dataset.
#' @params model The growth model(s) to fit.
#' @params errors The error model to use
#' @params mod.dir Optional: Path to a directory where the BUGS model file will be written for inspection. If not specified
#' a temporary directory will be created.
#' @params n.iter,n.burnin,n.thin The total number of MCMC iteractions to run, the number of burn in iterations and the sample
#' thinning rate, respectively.
#' @params inits A list or function generating initial values for the 3 MCMC chains.
#' @params ncores If ncores > 1 and number of models specified is also >1, models will be run in parallel.

grow_bayes = function(size,age,group=NULL,model = c('VB','logistic','Gompertz','Richards','cessation','VBlogK'),
                      errors = c('normal','lognormal'),mod.dir=NULL,n.iter=25000,n.thin=20,n.burnin=5000,inits=NULL,ncores=1){
  
  if(!is.null(group)){ 
    group = as.factor(group)
    levs  = levels(group)
    group = as.numeric(group)
  } else group = rep(1,length(size))
  
  data = data.frame(size,age,group) # this will stop the function if vectors are different lengths
  keep = complete.cases(data)
  data = data[keep,]
  remove= sum(!keep)
  if(remove>0) cat('Removing',remove,'rows with missing data')
  
  jags.data = list(length = data$size,age = data$age,group = data$group,N=length(data$size),G = length(unique(data$group)))
  
  if(is.null(mod.dir)) mod.dir = tempdir()
  mod.dir = gsub('\\/$', '', mod.dir)
  
  if(length(model)>1){
    
    if(ncores>1){
    plan(tweak(multisession,workers = ncores))
    cat('Fitting',length(model),'models in parallel.....')
    }
      
    fit = tibble(model = model) %>% 
      mutate(jags.fit = furrr::future_map(model, ~fit_gm(jags.data,.x,errors,mod.dir,n.iter,n.thin,n.burnin,inits))) %>%
      mutate(DIC = map_dbl(jags.fit,DIC)) %>%
      mutate(RMSE = map_dbl(jags.fit,RMSE)) %>%
      mutate(MPSRF = map_dbl(jags.fit,MPSRF)) %>%
      arrange(DIC) %>% mutate(model = forcats::fct_reorder(model,DIC))
    
    plan(sequential)
    
  } else {
    fit = fit_gm(jags.data,model,errors,mod.dir,n.iter,n.thin,n.burnin,inits)
  }
  
  return(fit)
}

# ----------------------------------------------------------------------------------------------------------------------------------
# Fit growth model - used internally by grow_bayes
# ----------------------------------------------------------------------------------------------------------------------------------

fit_gm = function(data,model,errors,mod.dir,n.iter=25000,n.thin=20,n.burnin=5000,inits=NULL){                  
  
  f = growth_model(model)
  
  mod.file = paste0(mod.dir,"/grow_bayes_",model,".txt")
  if(file.exists(mod.file)) unlink(mod.file)
  
  # Write the model file
  sink(file = mod.file, append = TRUE)
  
  cat(
    "model{
   for(i in 1:N){
   Lt[i] <- ",f,
    "
   ",
    error_model(errors),
    "
   }
 
   # SD
   tau.y ~ dgamma(0.001,0.001)
   sigma.y <- 1/sqrt(tau.y)

   # Priors
   for(j in 1:G){
   "
  )
  
  gm_priors(model)
  
  cat("
  }
}")
  
  sink()
  
  # Fit the model
  params = all.vars(formula(paste('~',f)))
  params = subset(params,!params %in% c('i','group','age'))
  
  if(!is.null(inits)) inits = lapply(1:3,function(x) inits(data$G))
  
  out = jags(data, inits=inits, params, model.file = mod.file, n.iter=n.iter,n.thin=n.thin,n.burnin=n.burnin,working.directory = mod.dir)
  
  structure(out,class = c("grow_bayes","rjags"))
}

# ----------------------------------------------------------------------------------------------------------------------------------
# Predict method
# ----------------------------------------------------------------------------------------------------------------------------------

# A predict method for objects of class `grow_bayes`. If newdata is not specified, size will be predicted at 1000 regularly spaced intervals
# between the minimum and maximum of age in the original data
                                     
predict.grow_bayes = function(model,newdata=NULL,plot=T){
  
  data = model$model$data()
  data = data.frame(age=data$age,length=data$length,group=data$group)
  
  if(is.null(newdata)){
    newdata = expand.grid(
      age = seq(min(data$age),max(data$age),length.out = 1000),
      group = unique(data$group)
    )
  } else if(is.null(newdata$group)){
    newdata = expand.grid(age=newdata$age,group = unique(data$group))
  }
  
  # Get the mean parameter estimates  
  summary = model$BUGSoutput$summary
  params  = model$parameters.to.save
  
  # This assigns new vector objects containing the means of each estimated parameter
  # whatever they may be
  for(p in params){
    assign(p,summary[which(grepl(p, rownames(summary))),'50%'])
  }
  
  # Retrieve and evaluate the formula from the object
  pow = function(x,p) x^p
  
  form = model$model$model()[3]
  form = gsub('\\[i\\]','',form)
  form = gsub(".*= ","",form)
  age = newdata$age; group = newdata$group
  newdata$length2 = eval(parse(text=form))
  newdata$group = factor(newdata$group)
  
  # Find the credible/confidence intervals for the mean using posterior samples
  post = model$BUGSoutput$sims.matrix
  post = post[sample(1:nrow(post),min(1000,nrow(post))),]
  
  # Now sapply the evaluation across the columns
  samples = sapply(1:nrow(post),function(j){
    for(p in params){
      assign(p,post[j,which(grepl(p, colnames(post)))])
    }
    eval(parse(text=form))
  })
  
  newdata$lower = apply(samples, 1, quantile, probs=0.025)
  newdata$upper = apply(samples, 1, quantile, probs=0.975)
  newdata$length = apply(samples, 1, mean) # Applying the mean estimates from
  # the model actually doesn't give you the best confidence intervals because of 
  # the way growth is calculated. Better to use the median/mean from the samples
  
  # Plot the results
  if(length(unique(newdata$group))>1){
    pl= ggplot(newdata,aes(x=age,y=length,colour=group,fill=group)) 
  } else {pl = ggplot(newdata,aes(x=age,y=length))}
  
  data$group = factor(data$group)
  pl = pl + geom_point(data=data) + geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.1,colour='transparent') +
    geom_line()
  
  if(plot) print(pl)
  return(newdata)
}

# ------------------------------------------------------------
# summary method
# ------------------------------------------------------------                                   

# Summary method for objects returned by `grow_bayes`
                                     
summary.grow_bayes =  function(jags.model){
  
  as.data.frame(jags.model$BUGSoutput$summary) %>% 
    tibble::rownames_to_column('parameter') %>% subset(parameter!='deviance') %>% 
    dplyr::select(parameter,mean,median=`50%`,lcl=`2.5%`,ucl=`97.5%`)
  
}

# ------------------------------------------------------------
# Deviance Information Criterion Method
# ------------------------------------------------------------

# Extracts DIC from a fitted JAGS model.                                     
                                     
DIC = function(model){
  
  if(!inherits(model,'rjags')) stop("model should be an object of class rjags")
  model$BUGSoutput$DIC
  
}

# --------------------------------------------------------------------
# Extract multiple potential scale reduction factor as a convergence diagnostic
# --------------------------------------------------------------------

MPSRF = function(model){
  
  if(!inherits(model,'rjags')) stop("model should be an object of class rjags")
  coda::gelman.diag(as.mcmc(model))$mpsrf

}

# ------------------------------------------------------------
# Extract root mean square error/sum of squares from grow_bayes object
# ------------------------------------------------------------

RMSE = function(model,type = c('RMSE','SS')){
  
  if(!inherits(model,'grow_bayes')) stop("model should be an object of class grow_bayes")
  type = match.arg(type)
  data = model$model$data()
  data = data.frame(age=data$age,length=data$length,group=data$group)
  
  resids = data$length - predict(model,newdata = data,plot=F)$length
  switch(type,RMSE = sqrt(mean(resids^2)), SS = sum(resids^2))
}
                                     

# -------------------------------------------------------------------------------------------------                                     
# Growth model functions
# -------------------------------------------------------------------------------------------------                                    

# To make functions readable in BUGS format
pow = function(x,p) x^p 
                                     
Richards = function(age,Linf,k,a,b){
  Linf*pow(1 - a * exp(-k*age),b)
}
                                     
Richards.t0 = function(age,Linf,k,t0,b){
  Linf*pow(1 - 1/b * exp(-k*(age-t0)),b)
}
                                     
VBlogK    = function(age,Linf,k1,k2,t0,alpha,beta){                                  
  Linf*(1 - exp(-k2*(age - t0)) * pow((1+exp(-beta*(age-t0-alpha)))/(1+exp(beta*alpha)),(k1-k2)/beta))
}
                                     
 Richards.logK    = function(age,Linf,k1,k2,t0,alpha,beta,d){                                  
  Linf*(1 - exp(-k2*(age - t0)) * pow(pow((1+1/d*exp(-beta*(age-t0-alpha))),d)/(1+exp(beta*alpha)),(k1-k2)/beta))
}
  
# -------------------------------------------------------------------------------------------------                                     
# Converting growth model functions into BUGS language to run
# -------------------------------------------------------------------------------------------------                                          

fun2jags = function(fun){
  
  args = formalArgs(fun)
  new  = paste0(args,ifelse(args=='age','[i]','[group[i]]'))
  args = paste0('\\b',args,'\\b')
  
  .f = deparse(get(fun))[-1] %>% subset(!.%in% c('{','}')) %>% stringi::stri_trim() %>% paste(collapse=' ')
  stringi::stri_replace_all_regex(.f, pattern = args, replacement = new, vectorize_all = FALSE) 
  
}
                                                                      

# -------------------------------------------------------------------------------------------------                                     
# An implementation of the Farley et al. 2020 decimal age correction algorithm
# -------------------------------------------------------------------------------------------------                                      

# x is a dataframe where each row is an individual and columns are measurements from the primordium to 
# the start of each opaque zone. The first column is assumed to be the primordium and should be zero
# for each individual. The last column is the measurement to the otolith edge and should have a value for each indidual.
# Intermediate columns may have NAs where measurements are missing for younger individuals.
                                     
decimal_age = function(x,method = 'Farley'){
 
 incs = t(apply(x[,-ncol(x)],1,diff))
 mu = colMeans(incs,na.rm=T)
 
 get_edge = function(i){
   d = diff(na.omit(i))
   n = length(d)
   age2first = 0.7486*d[1]^1.8709
   mir = ifelse(n>1,min(1,d[n]/mu[n]),NA)
   cbind(age2first,max(0,n-2),mir)  
 }
 
 t(apply(x,1,get_edge)) %>%
 cbind(rowSums(.,na.rm=T)) %>%
 data.frame() %>% 
 setNames(c('first','zones','margin','decimal_age'))
 
}
                                                                          
# ------------------------------------------------------------------------------------------
# Function to calculate the Index of Average Percentage Error
# ------------------------------------------------------------------------------------------
                                     
IAPE = function(x){
  
  x   = x[rowSums(!is.na(x))>1,]
  mu  = rowMeans(x,na.rm=T)
  xij = abs(x - mu)/mu
  Xj  = rowMeans(xij,na.rm=T)
  mean(Xj,na.rm=T) * 100
  
}
