## LARS/LASSO stage 1 
## Without imputation (complete-case analysis)
## Written to rerun ARIES analyses 
## Yiwen Zhu
## Oct 23, 2018

## Add mobility hypothesis
## Jiaxuan Jessie Liu
## Apr 5, 2019

## Add FWL confounding adjustment
## Jiaxuan Jessie Liu
## March 16, 2020

## Add CI for SI 
## Adapted from: LARS-noimpute-function-20200520.R (Yiwen)
## Jiaxuan Jessie Liu
## July 22, 2020

## df: analytic sample with complete data
## outcome.vec: a vector of outcome names (e.g. all CpGs)
## adver: short name of adversity
## covars: a list of covariates
## hypos: additional hypotheses to be tested (aside from the individual time point)
## recency.vec: recency vector - weighting by age in months 
## inf.method: inference method, covTest, sI (selective inference)
## exposures: default to individual time points; otherwise specify a list of collapsed time points


## NOTE: Before running this, please check carefully that the variable type in the df 
## are set to the correct type(categorical/numeric)
## otherwise this won't run properly!! 


loadlibs <- function(){
  library(lars)
  library(dplyr)
  library(glmnet)
  library(intervals)
  library(selectiveInference)
}
suppressPackageStartupMessages(loadlibs())

rm(loadlibs)


### The complete function written by Yiwen, modified by Jessie: add  mobility hypothesis, only applicable for period-collapsed exposures
select.LARS.complete.addmob <- function(df, outcome.vec, adver, covars, hypos, recency.vec = NULL, 
                                inf.method = "covTest", 
                                exposures = "default", 
                                FWL = FALSE) {
  
  DNAm <- df[, c("ID", outcome.vec)]
  df <- df[, !grepl("^cg", colnames(df)) & !grepl("^ch", colnames(df))]
  
  
  ## data prep ---- 
  
  # define exposures 
  if (exposures == "default") { 
    exposures <- grep(paste0("^", adver), colnames(df), value = T)
  }
  
  # subset to those with complete data on this exposure and all covariates
  df <- df %>% filter(rowSums(is.na(.[,c(exposures,covars)])) == 0)
  DNAm <- DNAm[DNAm$ID %in% df$ID, ]
  
  # define sample size:
  n <- nrow(df) 
  
  # add ever/never exposed, accumulation, recency, and mobility
  exp <- do.call("cbind", lapply(df[, exposures], function(x) as.numeric(as.character(x))))
  df <- df %>% 
    mutate(accumulation = rowSums(exp)) %>% 
    mutate(ever = ifelse(accumulation > 0, 1, 0))
  
  if ("recency" %in% hypos){
    df$recency <-  rowSums(exp %*% diag(recency.vec))
  } 
  
  if ("mobility" %in% hypos){## this would depends on the exposure variables we input. If we already input collapsed periods, we don't need to create period measures. 
    x1=paste0(adver,"_very_early")
    x2=paste0(adver,"_early")
    x3=paste0(adver,"_middle")
    df$mobility_D12 <- (1-df[,x2])*df[,x1]
    df$mobility_U12 <- (1-df[,x1])*df[,x2]
    df$mobility_D23 <- (1-df[,x3])*df[,x2]
    df$mobility_U23 <- (1-df[,x2])*df[,x3]
    hypos <- c(hypos [hypos != "mobility"],"mobility_D12","mobility_U12","mobility_D23","mobility_U23")
  } 
  
  ## select covariates -- 
  ## what are the categorical covariates
  cat.covars <- do.call("c", lapply(1:length(covars), function(x) {
    cov1 <- covars[x]
    if (is.factor(df[, cov1])) {
      return(cov1)
    } else {
      return(NULL)
    }
  }))
  
  ## number of hypotheses 
  n_hypo <- length(exposures) + length(hypos)
  
  ## dummy code the categorical covars
  cat.covariates <- model.matrix(as.formula(paste0("~", paste(cat.covars, collapse = "+"))), df)[,-1]
  covariates <- as.matrix(cbind(cat.covariates, as.matrix(df[,setdiff(covars, cat.covars)])))
  X_hypos <- do.call("cbind", lapply(df[, c(exposures, hypos)], function(x) as.numeric(as.character(x))))
  X_residual <- lm(X_hypos ~ covariates)$resid
  
  ## Normalize X_hypos for use in max-|t| calculations
  col_mean <- apply(X_residual, 2, mean)
  X_centered <- X_residual - rep(col_mean, rep(n, n_hypo)) #subtract mean
  col_sss <- apply(X_centered, 2, function(x) sqrt(sum(x^2)))
  X_normed <- X_centered / rep(col_sss, rep(n, n_hypo)) #divide by sqrt sum squares
  Xt <- t(X_normed)
  XtX <- Xt %*% X_normed
  
  
  ## create LASSO function to return results
  
  ## covariance test ---- 
  ## y: single locus DNAm 
  run.covTest <- function(i){
    y <- DNAm[, outcome.vec[i]]
    y <- as.numeric(as.character(y))
    
    ## Frish-Waugh-Lovell 
    if (FWL) {
      y <- lm(y~covariates)$resid
    }
    
    if (i %% 1000 == 0){
      print(i)
    }
    
    lasso <- lars(X_residual, y)
    # last.action <- length(lasso$actions)
    
    ## Covariance test results
    # d <- covTest3(lasso, X_residual, y, maxp=last.action)$results
    
    d <- covTest3(lasso, X_residual, y, maxp=1)$results
    d <- as.data.frame(d)
    row.names(d) <- NULL
    d$Variable_Name <- colnames(X_residual)[abs(d$Predictor_Number)]
    d$Improvement_R2 <- lasso$R2[2]
    ## Remove steps where hypo dropped:
    d <- d[!is.na(d$Drop_in_covariance), ]
    
    # ## Make vector of results, up to max # hypos added:
    # k <- nrow(d)
    # d2 <- as.data.frame(t(d[, c(4, 5, 3)]))
    # r <- as.character(unlist(unname(d2)))
    # r2 <- c(r, rep(NA, c((15*3) - length(r))))
    
  }
  
  ## selective inference ---- 
  # create a similar function using selectiveInference: fixedLassoInference
  run.sI <- function(i){
    
    y <- DNAm[, outcome.vec[i]]
    y <- as.numeric(as.character(y))
    
    ## Frish-Waugh-Lovell 
    if (FWL) {
      y <- lm(y~covariates)$resid
    }
    
    if (i %% 1000 == 0){
      print(i)
    }
    
    # prevent p-value calculations from rounding
    options(digits=10)
    
    # Apply LARS
    lasso <- lars(X_residual, y)
    
    # scale X_residual
    sumsq <- lasso$normx
    X_normed <- scale(X_residual, scale = sumsq)
    
    #Covariance test results
    # first significant one
    fli <- fixedLassoInf(X_normed, y, lasso$beta[2, ], lasso$lambda[2], tol.beta = 1e-22) # lowered the tolerance for determining if a coefficient is zero
    scale <- sumsq[fli$vars]
    fli$coef0 <- fli$coef0/scale
    fli$sd <- fli$sd/scale
    fli$ci <- fli$ci/cbind(scale,scale)
    fli$vars <- names(fli$vars)
    fli$probe <- outcome.vec[i]
    fli.res <- as.data.frame(cbind(fli$probe,fli$vars, fli$coef0, fli$sd, fli$pv, fli$ci))
    colnames(fli.res) <- c("Probe","Variable", "Beta", "SE", "Pvalue", "CI.lo", "CI.up")
    fli.res$R2 <- lasso$R2[[2]]
    
    return(fli.res)
    
  }
  
  if (inf.method == "covTest") { 
    results <- as.data.frame(t(sapply(1:length(outcome.vec), run.covTest)))
  }
  
  if (inf.method == "sI") { 
    if (length(outcome.vec)==1) {results <- run.sI(1)} 
    else {results <- data.frame(do.call("rbind", lapply(1:length(outcome.vec), run.sI)))}
  }
  
  return(results)
}



###-----------------------------------------------
## load modified covTest, without rounding p-values:
covTest3 <- function (fitobj, x, y, sigma.est = "full", status = NULL, maxp = min(nrow(x), 
                                                                                  ncol(x))) 
{
  s4 = substring(fitobj$call, 1, 4)[1]
  s7 = substring(fitobj$call, 1, 7)[1]
  s8 = substring(fitobj$call, 1, 8)[1]
  if (s4 == "lars" & s7 != "lars.en" & s8 != "lars.glm") {
    calltype = "lars"
    family = "gaussian"
  }
  if (s7 == "lars.en") {
    calltype = "lars.en"
    family = "gaussian"
  }
  if (s8 == "lars.glm") {
    calltype = "lars.glm"
    family = fitobj$family
  }
  if (family == "cox") 
    stop("Cox model not yet implemented")
  if (calltype == "lars") {
    if (fitobj$type != "LASSO") {
      stop("Call to Lars must use type='LASSO'")
    }
  }
  if (calltype == "lars") {
    type = "lar"
  }
  if (calltype == "lars.en") {
    type = "lars.en"
  }
  if (calltype == "lars.glm") {
    type = "lars.glm"
  }
  if (calltype == "lars.glm" & sigma.est == "full") {
    sigma.est = 1
    cat("glm model; sigma set to 1", fill = TRUE)
  }
  n = nrow(x)
  p = ncol(x)
  my = mean(y)
  lambda.min.ratio = ifelse(nrow(x) < ncol(x), 0.1, 1e-04)
  jlist = unlist(fitobj$act)
  if (type == "lar") 
    lamlist = c(fitobj$lambda, 0)
  if (type == "lars.en") 
    lamlist = c(fitobj$lambda, 0)
  if (type == "lars.glm") 
    lamlist = c(fitobj$lambda, 0)
  maxp.call = maxp
  maxp = length(jlist)
  maxp = min(maxp, which(lamlist == 0))
  maxp = min(maxp, maxp.call)
  jlist = jlist[1:maxp]
  cov0 = cov = sig = rep(NA, maxp)
  yy = y - my
  if (family == "binomial") {
    glmobj = glmnet(x, y, family = "binomial", standardize = fitobj$standardize, 
                    lambda.min.ratio = lambda.min.ratio)
  }
  if (family == "cox") {
    glmobj = glmnet(x, Surv(y, status), family = "cox", standardize = fitobj$standardize, 
                    lambda.min.ratio = lambda.min.ratio)
  }
  if (family == "cox") {
    junk = calcz02(x, y, status)
    sc = junk$sceta
    inf = junk$infeta
    miinf = misqrt((t(inf) + inf)/2)
    yy = t(sc) %*% miinf
  }
  for (j in 1:maxp) {
    if (jlist[j] > 0) {
      lambda = lamlist[j + 1]
      if (type == "lar") 
        yhat = predict(fitobj, x, s = lambda, type = "fit", 
                       mode = "lam")$fit
      if (type == "lars.en") 
        yhat = (1 + fitobj$lambda2) * predict.lars.en(fitobj, 
                                                      x, lambda)
      if (type == "lars.glm" & family == "binomial") {
        yhat = as.vector(predict(glmobj, x, type = "link", 
                                 s = lambda/n))
      }
      if (type == "lars.glm" & family == "cox") {
        yhat = as.vector(predict(glmobj, x, type = "link", 
                                 s = lambda/n))
      }
      cov[j] = sum(yy * yhat)
      if (j == 1) {
        cov0[j] = 0
      }
      if (j > 1) {
        tt0 = which(fitobj$beta[j, ] != 0)
        if (type == "lar") {
          aa = update(fitobj, x = x[, tt0, drop = F])
          yhat0 = predict(aa, x[, tt0], type = "fit", 
                          s = lambda, mode = "lam")$fit
        }
        if (type == "lars.en") {
          aa = update(fitobj, x = x[, tt0, drop = F])
          yhat0 = (1 + fitobj$lambda2) * predict.lars.en(aa, 
                                                         x[, tt0], lambda)
        }
        if (type == "lars.glm") {
          if (family == "binomial") {
            if (length(tt0) == 1) {
              tt0 = c(tt0, tt0)
            }
            glmobj0 = glmnet(x[, tt0, drop = F], y, family = "binomial", 
                             standardize = fitobj$standardize, lambda.min.ratio = lambda.min.ratio)
            yhat0 = as.vector(predict(glmobj0, x[, tt0, 
                                                 drop = F], type = "link", s = lambda/n))
          }
          if (family == "cox") {
            if (length(tt0) == 1) {
              tt0 = c(tt0, tt0)
            }
            glmobj0 = glmnet(x[, tt0, drop = F], Surv(y, 
                                                      status), family = "cox", standardize = fitobj$standardize, 
                             lambda.min.ratio = lambda.min.ratio)
            yhat0 = as.vector(predict(glmobj0, x[, tt0, 
                                                 drop = F], type = "link", s = lambda/n))
          }
        }
        cov0[j] = sum(yy * yhat0)
      }
    }
  }
  if (is.numeric((sigma.est))) {
    sigma = sigma.est
    sigma.type = "known"
    null.dist = "Exp(1)"
    if (sigma.est <= 0) {
      stop("sigma.est must be positive")
    }
  }
  if (sigma.est == "full") {
    if (nrow(x) < ncol(x) + 1) 
      stop("Number of observations must exceed number of variables,\nwhen sigma.est is `full; you need to specify a numeric value for sigma.est")
    sigma.type = "unknown"
    aaa = lsfit(x, y)
    sigma = sqrt(sum(aaa$res^2)/(n - p))
    np = n - p
    null.dist = paste("F(2,", as.character(np), ")", sep = "")
  }
  tt = ((cov - cov0)/sigma^2)
  if (sigma.type == "known") {
    out = cbind(jlist, tt, 1 - pexp(tt, 1))
    dimnames(out)[[2]] = c("Predictor_Number", "Drop_in_covariance", 
                           "P-value")
  }
  if (sigma.type == "unknown") {
    out = cbind(jlist, tt, 1 - pf(tt, 2, n - p))
    dimnames(out)[[2]] = c("Predictor_Number", "Drop_in_covariance", 
                           "P-value")
  }
  dimnames(out)[[1]] = rep("", nrow(out))
  return(list(results = out, 
              sigma = round(sigma,  4), 
              null.dist = null.dist))
}