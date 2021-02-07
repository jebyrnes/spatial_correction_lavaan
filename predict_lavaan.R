#Code by Dr. Jarrett E.K. Byrnes; Last Modified 12/06/2016
#created in regards to Issue #44 on lavaan
#https://github.com/yrosseel/lavaan/issues/44

predict_lavaan <- function(fit, newdata = NULL){
  stopifnot(inherits(fit, "lavaan"))
  
  #Make sure we can use this
  if(!inspect(fit, "meanstructure")) stop("Need to supply meanstructure = TRUE in fit\n")
  if(is.null(newdata)){
    newdata <- data.frame(inspect(fit, "data"))
    names(newdata) <- lavNames(fit)
  }
  
  if(length(lavNames(fit, type="lv"))!=0) stop("Does not currently work with latent variables\n")
  
  #check for new data
  if(sum(!(lavNames(fit, type="ov.x") %in% names(newdata)))>0) stop("Not all exogenous variables supplied!")
  
  #Add some new columns to newdata
  newdata$Intercept <- 1
  newdata[lavNames(fit, "ov.nox")] <- 0
  
  
  mod_df <- data.frame(lhs = fit@ParTable$lhs,
                       op = fit@ParTable$op,
                       rhs = fit@ParTable$rhs,
                       exo = fit@ParTable$exo,
                       est = fit@ParTable$est,
                       se = fit@ParTable$se,
                       stringsAsFactors=FALSE)
  
  #Drop covariances
  mod_df <- mod_df[-which(mod_df$op=="~~"),]
  mod_df[which(mod_df$op=="~1"),]$rhs <- "Intercept"
  
  #get rid of exogenous on lhs
  mod_df <- mod_df[-which(mod_df$exo==1),]
  
  #Order by lhs
  mod_df <- mod_df[sort(mod_df$lhs, index.return=TRUE)$ix,]
  
  #let us know which variables on the rhs are exogenous
  mod_df$ord <- 0
  mod_df[which(!(mod_df$rhs %in% mod_df$lhs)),]$ord <- 1
  
  #Make a "order"
  ord_current <- 1
  while(sum(mod_df$ord==0)>0){
    for(r in unique(mod_df$lhs)){
      val <-  sum(mod_df[which(mod_df$lhs==r),]$ord==0)
      if(val==0) {
        mod_df[which(mod_df$lhs==r),]$ord <- ord_current
        
        if(sum(mod_df$rhs==r)>0)
          mod_df[which(mod_df$rhs==r),]$ord <- ord_current+1
      }
    }
    ord_current <- ord_current +1
  }
  
  #correct for ragged ordering
  for(r in unique(mod_df$lhs)){
    mod_df[which(mod_df$lhs==r),]$ord <- max(mod_df[which(mod_df$lhs==r),]$ord)
  }
  
  #sort by order 
  mod_df <- mod_df[sort(mod_df$ord, index.return=TRUE)$ix,]
  
  #now do the fitting in order
  fit_df <- data.frame(base = rep(1, nrow(newdata)))
  
  for(r in unique(mod_df$lhs)){
    subdf <- subset(mod_df, mod_df$lhs==r)
    #make a formula
    rhs <- paste0(subdf$rhs, collapse=" + ")
    form <- as.formula(paste0(r, " ~ ", rhs))
    
    #use formula to get right part of the data in right format
    mod_mat <- model.matrix(form, newdata)[,-1]
    new_val = mod_mat %*% subdf$est
    
    fit_df[[r]] <- new_val
    newdata[[r]] <- new_val
  }
  
  return(fit_df[,-1])
  
}


fitted_lavaan <- function(fit){
  predict_lavaan(fit)
}

residuals_lavaan <- function(fit){
  fitted_vals <- fitted_lavaan(fit)
  
  rawdata <- data.frame(inspect(fit, "data"))
  names(rawdata) <- lavNames(fit)
  
  res <- data.frame(base = rep(1, nrow(rawdata)))
  for(vals in names(fitted_vals)){
    res[[vals]] <- rawdata[[vals]] - fitted_vals[[vals]] 
  }
  
  return(res[,-1])
}


############################################
########### EXAMPLE CODE
############################################
# 
#iris_mod <- "
# Sepal.Width ~ Sepal.Length
# Petal.Width ~ Sepal.Width + Petal.Length
# "
 
# iris_fit <- sem(iris_mod, data=iris, meanstructure = TRUE)
 
# fitted_lavaan(iris_fit)
 
# res <- residuals_lavaan(iris_fit)
 
# par(mfrow=c(1,2))
# for(i in 1:2) {
#   qqnorm(res[,i], main=names(res)[i])
#   qqline(res[,i])
# }
# par(mfrow=c(1,1))   
 
# newdata <- data.frame(Sepal.Length=1:10, Petal.Length=11:20)
# pred <- predict_lavaan(iris_fit, newdata=newdata)
 
# plot(Petal.Width ~ Sepal.Length, data=cbind(newdata, pred))
 
# library(mvnormtest)
# mshapiro.test(t(res))
 