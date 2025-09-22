library(Formula)
library(VGAM) 
library(tidyverse)
library(hnp) 
library(pscl) 
library(MASS) 
library(reshape2) 
library(reshape)

# This script assumes that the WeibullT functions have already been loaded,
# especially from 'funcoesWeibullT.R'

# Fitting the truncated Weibull model
WeibullT_reg.fit <- function(formula,output,link.mu,link.phi,data){
  #-----------------------------------------------------------------#
  # FORMULA
  #-----------------------------------------------------------------#
  form = Formula(formula)
  model = model.frame(form,data)
  y<-model.response(model) # response variable
  X<-model.matrix(form,data=model, rhs = 1) # design matrix for the mu parameter
  Z<-model.matrix(form,data=model, rhs = 2) # design matrix for the phi parameter
  n <- length(y)
  p = ncol(X)
  q = ncol(Z)
  #-----------------------------------------------------------------#
  # Link function for the mode
  #-----------------------------------------------------------------#
  if(any(link.mu == c("logit","probit","cloglog","identity"))) {
    stats1 <- make.link(link.mu)
  }else stop(paste(link.mu, "link not available, available links are",
                   "\"log\" and \"identity\""))
  
  link1 <- structure(list(link = link.mu, 
                          linkfun = stats1$linkfun,
                          linkinv = stats1$linkinv, 
                          mu.eta = stats1$mu.eta, 
                          diflink = function(t) 1/(stats1$mu.eta(stats1$linkfun(t)))
  ))
  #-----------------------------------------------------------------#
  # Link function for the PHI
  #-----------------------------------------------------------------#
  if(any(link.phi == c("log", "identity", "1/mu^2"))) {
    stats2 <- make.link(link.phi)
  }else stop(paste(link.mu, "link not available, available links are",
                   "\"logit\", \"probit\". \"cauchit\" and \"cloglog\"" ))
  
  link2 <- structure(list(link = link.phi, 
                          linkfun = stats2$linkfun,
                          linkinv = stats2$linkinv, 
                          mu.eta = stats2$mu.eta, 
                          diflink = function(t) 1/(stats2$mu.eta(stats1$linkfun(t)))
  ))
  #-----------------------------------------------------------------#
  # Derivative of the incomplete gamma function
  #-----------------------------------------------------------------#
  derivada_gama <- function(theta, X1, Z1, a, c){
    beta <- theta[1:p]
    nu <- theta[(p+1):length(theta)]
    
    eta.m <- X1%*%as.matrix(beta)
    eta.phi <- Z1%*%as.matrix(nu)
    
    if(link.mu == "identity"){
      m <-link1$linkinv(eta.m)
    }else{
      m <- 2*link1$linkinv(eta.m)-1
    }
    phi <- link2$linkinv(eta.phi)
    
    h<- list()
    g<- numeric()
    delta = ((2/(m+1))^phi)*((phi-1)/phi)
    for(i in 1:length(m)){
      h[[i]] = function(t){
        (t^(a-1) )* (log(t)^c) * exp(-t)
      }
      g[i] = integrate(h[[i]], lower=0, upper=delta[i])$value
    }
    return(g)
  }
  #-----------------------------------------------------------------#
  # Likelihood function
  #-----------------------------------------------------------------#
  fn2 <- function(theta){
    beta <- theta[1:p]
    nu <- theta[(p+1):length(theta)]
    
    eta.m <- X%*%as.matrix(beta)
    eta.phi <- Z%*%as.matrix(nu)
    
    if(link.mu == "identity"){
      m <-link1$linkinv(eta.m)
    }else{
      m <- 2*link1$linkinv(eta.m)-1
    }
    
    phi <- link2$linkinv(eta.phi)
    
    log.likelihood <- suppressWarnings(sum(
      log(phi/(m+1)) + log((phi-1)/phi) + (phi-1)*log((y+1)/(m+1)) - (((y+1)/(m+1))^phi) * ((phi-1)/phi) -
        log(1 - exp(- ((2/(m+1))^phi) * ((phi-1)/phi)))
    ))
    return(-log.likelihood)
  }
  #-----------------------------------------------------------------#
  # Fisher information
  #-----------------------------------------------------------------#
  Fisher <- function(theta, data){
    beta <- theta[1:p]
    nu<- theta[(p+1):length(theta)]
    #
    eta.m <- X%*%as.matrix(beta) # eta_1i = g1(m_i*) = X^T * beta
    if(link.mu == "identity"){
      m <-link1$linkinv(eta.m)
    }else{
      m <- 2*link1$linkinv(eta.m)-1
    }
    #m
    eta.phi <- Z%*%as.matrix(nu) # eta_2i = g2(phi_i) = Z^{T} * nu
    phi <- link2$linkinv(eta.phi) # phi_i = g2^(-1) (eta.phi)
    delta <- ((2/(m+1))^phi) * ((phi-1)/phi)
    
    #-----------------------------------------------------------------#
    mT <- diag(as.vector( 2 * (1/link1$diflink((m+1)/2) ) )) # (dm_i/d eta_1i = 2/ g_1^{'} (m*)) 
    mC <- diag(as.vector( link2$mu.eta(eta.phi)   )  ) # (dphi_i/eta_2i) = 1/g_{2}^{'} (phi)
    
    #-----------------------------------------------------------------#
    # Important results
    #-----------------------------------------------------------------#
    EY =  ((phi/(phi-1))/(1- exp(-delta)) ) * derivada_gama(theta, X, Z, 2, 0)
    #EY =  ((phi/(phi-1))/(1- exp(-delta)) ) * pgamma.deriv.unscaled(delta, 2)[,1]
    EYlog = (   (1/(phi-1))/ (1 - exp(-delta)))*( derivada_gama(theta, X, Z, 2, 1) + log(phi/(phi-1)) * derivada_gama(theta, X, Z, 2, 0) )
    #EYlog = (   (1/(phi-1))/ (1 - exp(-delta)))*( pgamma.deriv.unscaled(delta, 2)[,2] + log(phi/(phi-1)) * pgamma.deriv.unscaled(delta, 2)[,1] )
    EYlog2 = (  (  (1/phi)*(1/(phi-1)) )/(1- exp(-delta))  ) * ( 2 * log(phi/(phi-1)) * derivada_gama(theta,X, Z, 2, 1) + (  log(phi/(phi-1))  )^(2) * derivada_gama(theta,X,Z,2,0) + derivada_gama(theta,X, Z,2,2) )
    #EYlog2 = (  (  (1/phi)*(1/(phi-1)) )/(1- exp(-delta))  ) * ( 2 * log(phi/(phi-1)) * pgamma.deriv.unscaled(delta, 2)[,2] + (  log(phi/(phi-1))  )^(2) * pgamma.deriv.unscaled(delta, 2)[,1] + pgamma.deriv.unscaled(delta, 2)[,3] )
    #-----------------------------------------------------------------#
    # Second-order derivative with respect to m
    h1 <- (   exp(-delta) * ((2/(m+1))^(2*phi)) * (( (phi-1)/(m+1) )^2)-
                exp(-delta) * ((2/(m+1))^phi) * (phi/(m+1)) * ((phi-1)/(m+1)) +
                exp(-2*delta) * ((2/(m+1))^phi) * (phi/(m+1)) * ((phi-1)/(m+1)) -
                exp(-delta) * ((2/(m+1))^phi) * ((phi-1)/((m+1)^2)) +
                exp(-2*delta) * ((2/(m+1))^phi) * ((phi-1)/((m+1)^2))
    )/ ( (1-exp(-delta))^2 )
    # Second-order derivative with respect to phi
    h3 <- -((
      -exp(-delta) * ((2/(m+1))^(2*phi)) * (((phi-1)/phi)^2) * ((log( 2/(m+1)))^2) -
        2 * exp(-delta) * ((2/(m+1))^(2*phi)) * ((phi-1)/phi) * (1/(phi^2)) * log(2/(m+1)) +
        exp(-delta) * ((2/(m+1))^phi) * ((phi-1)/phi) * ((log(2/(m+1)) )^2) -
        exp(-2*delta) * ((2/(m+1))^phi) * ((phi-1)/phi) * ((log(2/(m+1)) )^2) +
        2 * exp(-delta) * ((2/(m+1))^phi) * (1/(phi^2)) * log(2/(m+1)) -
        2 * exp(-2*delta) * ((2/(m+1))^phi) * (1/(phi^2)) * log(2/(m+1)) -
        exp(-delta) * ((2/(m+1))^(2*phi)) * ((1/(phi^2))^2) -
        exp(-delta) * ((2/(m+1))^phi) * (2/(phi^3)) +
        exp(-2*delta) * ((2/(m+1))^phi) * (2/(phi^3))
    )/( ( 1-exp(-delta))^2 ) )
    # Mixed partial derivative
    h2 <- ( -exp(-delta) * ((2/(m+1))^(2*phi)) * log(2/(m+1)) * ((phi-1)/phi) * ((phi-1)/(m+1)) -
              exp(-delta) * ((2/(m+1))^(2*phi)) * (1/(phi^2)) * ((phi-1)/(m+1)) +
              exp(-delta) * ((2/(m+1))^phi) * log(2/(m+1)) * ((phi-1)/(m+1)) -
              exp(-2*delta) * ((2/(m+1))^phi) * log(2/(m+1)) * ((phi-1)/(m+1)) +
              exp(-delta) * ((2/(m+1))^phi) * (1/(m+1)) -
              exp(-2*delta) * ((2/(m+1))^phi) * (1/(m+1))
    )/( (1-exp(-delta))^2) 
    #-----------------------------------------------------------------#
    # Expected value of the second-order derivatives
    #-----------------------------------------------------------------#
    # m
    ES <- - (phi/((m+1)^2)) + ( (phi-1)/((m+1)^2)) * EY + ((phi-1)/(m+1))*(phi/(m+1)) * EY - h1
    # Expected value of the cross partial derivative
    ED <- (1/(m+1)) - (1/(m+1)) * EY - ((phi-1)/(m+1)) * EYlog - h2
    # phi
    EV <- (1/((phi-1)^2)) + ((phi-1)/phi) * EYlog2 + (2/(phi^2))* EYlog - (2/(phi^3)) * EY -  h3
    #-----------------------------------------------------------------#
    H = diag(as.vector(ES))%*%(mT^2)
    C = mT%*%diag(as.vector(ED))%*%mC
    Q = diag(as.vector(EV))%*%(mC^2)
    
    Kbb <- t(X)%*%(H)%*%X
    Kbn <- t(X)%*%(C)%*%Z
    Knn <- t(Z)%*%(Q)%*%Z
    
    #-----------------------------------------------------------------#
    Ktheta<-rbind(
      cbind(Kbb,(Kbn)),
      cbind(t(Kbn),Knn)
    )
    
    return(Ktheta)
  }
  ##------------------------ Initial Value ----------------------------------#
  ynew <- link1$linkfun( (y+1)/2)
  ajuste1 <- lm(ynew ~ X+0)
  mqo=ajuste1$coef # Ordinary Least Squares
  s2 <- 2
   
  ajuste2 = as.numeric(rep(s2,q))
  start2 = c(as.numeric(ajuste1$coefficients), ajuste2)
  #-----------------------------------------------------------------#
  # Model fitting
  fit2 <- try(nlminb(start2,fn2,control=list(rel.tol=1e-6)),silent = T)
  fit1 <- try(optim(start2, fn2, method = "Nelder-Mead", hessian = T,
                control = list(fnscale = 1, maxit = 3000, reltol = 1e-6)), silent = T)
  
  if (fit2$convergence != 0)
    warning("FUNCTION DID NOT CONVERGE!")
  #-----------------------------------------------------------------#
  # R^2
  #-----------------------------------------------------------------#
  # Likelihood function without regression structure
  lvero1 <- function(par){
    alpha = par[1]
    lambda = par[2]
    ll <- suppressWarnings(
      log(alpha) - log(lambda) + (alpha-1) * log((y+1)/lambda) - (((y+1)/lambda)^alpha) - log(  1 - exp(-((2/lambda)^alpha)) )
    )
    sum(-ll)
  }
  # nlminb()
  fitLvero1 <- try(nlminb(c(2,2),lvero1,control=list(rel.tol=1e-6)),silent = T)
  # optim()
  fitLvero2 <- try(optim(c(2,2), lvero1, method = "Nelder-Mead", hessian = T,
                         control = list(fnscale = 1, maxit = 3000, reltol = 1e-6)), silent = T)
  # COX & SNELL
  R2_nlminb<- 1-exp(-(2/n)*(-fit2$objective-fitLvero1$objective))
  R2_optim<- 1-exp(-(2/n)*(-fit1$value-fitLvero2$value))
  #-----------------------------------------------------------------#
  z <- list()
  z$X=X
  z$Z=Z
  z$conv <- fit2$convergence
  z$conv_optm <- fit1$conv
  coef <- (fit2$par)[1:(ncol(Z)+ncol(X))]
  coef_optm <- (fit1$par)[1:(ncol(Z)+ncol(X))]
  names(coef) <- c(c(paste("beta",1:ncol( as.matrix(X) ),sep="")),
                   c(paste("nu",1:ncol( as.matrix(Z) ),sep="")))
  names(coef_optm) <- c(c(paste("beta",1:ncol( as.matrix(X) ),sep="")),
                        c(paste("nu",1:ncol( as.matrix(Z) ),sep="")))
  z$coef <- coef
  z$coef_optm <- coef_optm
  #----------------------------------------------------------------------------#
  # Covariance
  #----------------------------------------------------------------------------#
  # Inverse of Fisher's matrix
  #-----------------------------------------------------------------#
  if((link.mu=="identity") & (link.phi=="identity")){
    z$fisher2 <- Fisher2(fit2$par, y) # nlminb - Fisher
    z$fisher1 <- fit1$hessian # optim
    
    vcov2 <- (1/n)*solve(z$fisher2) # nlminb
    vcov1 <- solve(z$fisher1) # optim
    #
    z$vcov2 <- vcov2
    z$vcov1 <- vcov1
    #
    stderror2 <- sqrt(diag(vcov2))
    stderror1 <- sqrt(diag(vcov1))
    # Standard Error
    z$stderror2 <- stderror2
    z$stderror1 <- stderror1
    # Wald
    z$zstat2 <- abs(z$coef/stderror2)
    z$zstat1 <- abs(z$coef/stderror1)
    z$pvalues2 <- 2*(1 - pnorm(z$zstat2) )
    z$pvalues1 <- 2*(1 - pnorm(z$zstat1) )
  }else{
    # With regression structure - optim() and nlminb()
    z$fisher2 <- Fisher(fit2$par, data = y) # nlminb - Expected Fisher Information
    z$fisher1 <- fit1$hessian # optim
    #
    vcov2 <- solve(z$fisher2) # nlminb
    vcov1 <- solve(z$fisher1) # optim
    #
    z$vcov2 <- vcov2
    z$vcov1 <- vcov1
    #
    stderror2 <- sqrt(diag(vcov2))
    stderror1 <- sqrt(diag(vcov1))
    #
    z$stderror2 <- stderror2
    z$stderror1 <- stderror1
    #
    z$zstat2 <- abs(z$coef/stderror2)
    z$zstat1 <- abs(z$coef/stderror1)
    z$pvalues2 <- 2*(1 - pnorm(z$zstat2) )
    z$pvalues1 <- 2*(1 - pnorm(z$zstat1) )
  }
  #-----------------------------------------------------------------#
  z$Par_est<-fit2$par # nlminb()
  beta <- coef[1:ncol(X)] 
  nu<- coef[(ncol(X)+1):length(coef)] 
  #
  eta1hat <- X%*%as.matrix(beta)
  if(link.mu == "identity"){
    mhat <-link1$linkinv(eta1hat)
  }else{
    mhat <- 2 * link1$linkinv(eta1hat) - 1
  }
  #mhat
  eta2hat <- Z%*%as.matrix(nu)
  phihat <- link2$linkinv(eta2hat)
  #
  alphahat=phihat
  #
  z$alpha=alphahat
  z$fitted <- mhat
  z$eta1hat <- eta1hat
  z$eta2hat<- eta2hat
  z$lambda <- (z$fitted+1)*((z$alpha/(z$alpha-1))^(1/z$alpha))
  #
  z$serie <- y
  z$WT <- names(coef)
  z$l.phi<-link.phi
  z$l.mu<-link.mu
  #-----------------------------------------------------------------#
  # optm
  z$Par_optim<-fit1$par # optim()
  beta_optm <- coef_optm[1:ncol(X)] 
  nu_optm<- coef_optm[(ncol(X)+1):length(coef_optm)] 
  #
  eta1hat_optm <- X%*%as.matrix(beta_optm)
  if(link.mu == "identity"){
    mhat_optm <-link1$linkinv(eta1hat_optm)
  }else{
    mhat_optm <- 2 * link1$linkinv(eta1hat_optm) - 1
  }
  #mhat
  eta2hat_optm <- Z%*%as.matrix(nu_optm)
  phihat_optm <- link2$linkinv(eta2hat_optm)
  #
  alphahat_optm=phihat_optm
  #
  z$alpha_optm=alphahat_optm
  z$fitted_optm <- mhat_optm
  z$eta1hat_optm <- eta1hat_optm
  z$eta2hat_optm<- eta2hat_optm
  z$lambda_optm <- (z$fitted_optm+1)*((z$alpha_optm/(z$alpha_optm-1))^(1/z$alpha_optm))
  #
  z$serie <- y
  z$WT <- names(coef_optm)
  z$l.phi<-link.phi
  z$l.mu<-link.mu
  #-----------------------------------------------------------------#
  # Residuals
  #-----------------------------------------------------------------#
  # Randomized quantile residual
  z$resid1 = qnorm( pWeibull(as.numeric(y), alfa = z$alpha, lamb = z$lambda))
  z$resid_optm = qnorm( pWeibull(as.numeric(y), alfa = z$alpha_optm, lamb = z$lambda_optm))
  #-----------------------------------------------------------------#
  # Result: nlminb
  z$loglik2 <- fit2$objective
  z$counts2 <- as.numeric(fit2$evaluations[1] )
  z$aic2 <- -2*z$loglik2+2*(length(nu)+length(beta))
  z$bic2 <- -2*z$loglik2+log(n)*(length(nu)+length(beta))
  z$r2_nlminb <- R2_nlminb
  #-----------------------------------------------------------------#
  # Result: optim
  z$loglik1 <- fit1$value
  z$counts1 <- as.numeric(fit1$counts[1] )
  z$aic1 <- -2*z$loglik1+2*(length(nu)+length(beta))
  z$bic1 <- -2*z$loglik1+log(n)*(length(nu)+length(beta))
  z$r2_optim <- R2_optim
  z$formula <-formula
  #-----------------------------------------------------------------#
  model_nlminb <- cbind(round(z$coef,4),
                        round(z$stderror2,4),
                        round(z$zstat2,4),
                        round(z$pvalues2,4)
  )
  model_optim <- cbind(round(z$coef_optm,4),
                       round(z$stderror1,4),
                       round(z$zstat1,4),
                       round(z$pvalues1,4)
  )
  colnames(model_nlminb)<-c("Estimate nlminb",
                            "Std. Error",
                            "z value",
                            "Pr(>|z|)"
  )
  colnames(model_optim)<-c("Estimate Optm",
                           "Std. Error",
                           "z value",
                           "Pr(>|z|)"
  )
  z$model_nlminb <- model_nlminb
  z$model_optm <- model_optim
  
  if(output==1){
    print(model_nlminb)
    print(" ",quote=F)
    print(c("Log-likelihood:",round(z$loglik2,4)),quote=F)
    print(c("Number of iterations in Nlminb:",z$counts2),quote=F)
    print(c("AIC:",round(z$aic2,4)," BIC:",round(z$bic2,4)),quote=F)
    print("Residuals:",quote=F)
    print(summary(z$resid1))
    print(c("R-squared:",round(z$r2_nlminb,4)),quote=F)
    #
    print(" ",quote=F)
    #
    print(model_optim)
    print(" ",quote=F)
    print(c("Log-likelihood:",round(z$loglik1,4)),quote=F)
    print(c("Number of iterations in optim:",z$counts1),quote=F)
    print(c("AIC:",round(z$aic1,4)," BIC:",round(z$bic1,4)),quote=F)
    print("Residuals:",quote=F)
    print(summary(z$resid_optm))
    print(c("R-squared:",round(z$r2_optim,4)),quote=F)
  }
  return(z)
  
}

#------------------------------------------------------------------------------#
# GRAPHICS
diag.WT.fit <- function(model,sim=100,conf=.95, data, r) {
  alfa1 <-(1-conf)/2 
  #-----------------------------------------------------------------#
  formula.model = model$formula
  form = Formula(formula.model)
  f = model.frame(form, data)
  X<-model.part(form,data=f, rhs = 1) # X
  Z<-model.part(form,data=f, rhs = 2) # Z
  #-----------------------------------------------------------------#
  alpha= model$alpha
  lambda=model$lambda
  yfitted=model$fitted
  
  alpha_optm=model$alpha_optm
  lambda_optm=model$lambda_optm
  res_optm = model$resid_optm
  
  y <-model$serie
  n <-length(y)
  res <- model$resid1 # residuals
  e <- matrix(0,n,sim)
  e1 <- numeric(n)
  e2 <- numeric(n)
  #-----------------------------------------------------------------#
  # quantis teoricos
  w = 1:n
  q = (w+n-1/8)/(2*n+0.5)
  qteo = qnorm(q)
  #-----------------------------------------------------------------#
  i<-1
  while(i<=sim) {
    ynew <- rWtruncated(Nsim =  n, alfa = alpha, lamb = lambda)
    formula_new = update(form, ynew ~ .)
    data_new <- data[,-1] |>
      mutate(ynew = ynew, .before = 1)
    New_formula = as.formula(formula_new)
    fit = try(WeibullT_reg.fit(New_formula,link.mu=model$l.mu,link.phi=model$l.phi, output = 0,data_new), silent = T)
    
    if(class(fit) != "try-error"){
      if(((fit$conv==0)  & (sum(is.na(fit$resid1))==0)   )) {
        ti <- fit$resid1
        eo <- sort(abs( fit$resid1 )) # Sorts the absolute values of the quantile residuals
        e[,i] <- eo
        i<- i+1
        print(i)
      } 
    }
  }
  # Confidence bands
  l1 = apply(e,1,min)
  l2 = apply(e,1,max)
  #
  for(i in 1:n) {
    eo <- sort(e[i,])
    e1[i] <- quantile(eo,alfa1)
    e2[i] <- quantile(eo,1-alfa1)
  }
  
  med <- apply(e,1,median)
  qq <- qnorm((n+1:n-1/8)/(2*n+0.5))
  #-----------------------------------------------------------------#
  # ggplot construction:
  # x - quantile of the normal distribution, given by the variable qteo
  # lower - minimum of the simulated values
  # upper - maximum of the simulated values
  # median - median of the simulated values
  # residuals - randomized quantile residuals (RQ)
  hnp_cores <- function(X) {
    out <- X[["residuals"]] < X[["lower"]] | X[["residuals"]] > X[["upper"]]
    c("Dentro", "Fora")[out + 1L]
  }
  hnp_texto <- function(X) {
    out <- X[["residuals"]] < X[["lower"]] | X[["residuals"]] > X[["upper"]]
    n <- sum(out)
    txt1 <- sprintf("Total de pontos: %d", nrow(X))
    txt2 <- sprintf("Pontos fora do envelope: %d (%2.4g%%)", n, 100*n/nrow(X))
    list(Total = txt1, Fora = txt2)
  }
  
  if(r==1){
    data.graf = data.frame(x = qteo,lower = e1,upper = e2,median = med, residuals = sort(abs(res)))
  }
  if(r==2){
    data.graf = data.frame(x = qteo,lower = l1,upper = l2,median = med, residuals = sort(abs(res)))
  }
  data.graf$cores <- hnp_cores(data.graf)
  data.graf_texto <- hnp_texto(data.graf)
  
  env <- ggplot(data = data.graf, aes(x)) +
    geom_point(aes(y = residuals, color = cores), show.legend = FALSE) +
    geom_text(x = 0, y = 4, label = data.graf_texto$Total, hjust = 0) +
    geom_text(x = 0, y = 3.7, label = data.graf_texto$Fora, hjust = 0) +
    #
    geom_ribbon(aes(ymin = lower, ymax = upper),alpha = 0.5
    )+
    geom_line(aes(y = median), linetype = 2, col="blue") + 
    scale_color_manual(values = c(Dentro = "black", Fora = "red")) +
    coord_cartesian(ylim = c(0, 6  )) +
    ylab("Resíduos quantílicos aleatorizados") +
    xlab("Quantis Teóricos N(0,1)")  +
    theme_bw()
  print(env)
  #-----------------------------------------------------------------#
  t<-seq(-5,n+6,by=1)
  j<-seq(n)
  
  t<-seq(-5,n+6,by=1)
  
  mar_b<-2.5
  mar_e<-2.5
  mar_c<-0.5
  mar_d<-0.5
  dist_text<-1.5
  dist_tick<-0.5
  #
  resid_ind = data.frame(indice = j, residuos = res)
  graf_ResidInd = ggplot(resid_ind, aes(indice, residuos)) + geom_point() +
    theme_bw() + geom_hline(yintercept = -3, col = "darkred", lty = 2)+
    geom_hline(yintercept = 3, col = "darkred", lty = 2) + ylab("Resíduos quantílicos aleatorizados") +
    xlab("Índice")
  print(graf_ResidInd)
  #-----------------------------------------------------------------#
  par(mar=c(mar_b, mar_e, mar_c, mar_d)) 
  par(mgp=c(dist_text, dist_tick, 0))
  densidade<-density(res)
  plot(densidade,ylab="Densidade", xlab="x",main=" ", ylim=c(0, 0.5),cex.lab=1.3)
  lines(densidade$x,dnorm(densidade$x),lty=2)
  legend("topleft",c("Densidade estimada","Normal padrão"),#pch=vpch,
         pt.bg="white", lty=c(1,2), bty="n")
  
  graf_fitted <- data.frame(res = res)
  #
  graf_fitted2 <- melt(graf_fitted)
  dens <- ggplot(graf_fitted2, aes(x = value, col=variable)) + geom_density(aes(value, colour = "Densidade estimada")) +
    stat_function(fun = dnorm, args = list(mean = mean(graf_fitted$res), sd = sd(graf_fitted$res)), aes(colour = "Normal Padrão"), 
                  lwd=0.5, geom="line", position="identity" ) + 
    labs(x = "x", y = "Densidade", colour = "") +
    scale_x_continuous(limits = c(-3,3)) +
    theme_bw() +
    theme(
      legend.key.size = unit(0.5, 'cm'), # change legend key size
      legend.key.height = unit(0.4, 'cm'), #change legend key height
      legend.key.width = unit(0.4, 'cm'), #change legend key width
      legend.title = element_text(size=8), #change legend title font size
      legend.text = element_text(size=7),
      plot.title = element_text(size = 11, family = "Tahoma", face = "bold"),
      text = element_text(size = 8, family = "Tahoma"),
      legend.position = c(0.08,0.99),
      legend.justification = c(0.08,0.99),
      legend.direction="vertical",
      legend.box="horizontal",
      legend.box.just = c("top"),
      legend.background = element_rect(size = 0.2, fill = "white", colour = "black"),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
    ) +
    scale_color_manual(values = c("black", "darkblue"), 
                       guide = guide_legend(override.aes = list(
                         linetype = c("solid", "solid")),
                         title = NULL)
    )
  print(dens)
  #-----------------------------------------------------------------#
  graf_fitted <- data.frame(y = y, yfitted = yfitted)
  obs_fitted <-  ggplot(graf_fitted, aes(x = y, y = yfitted)) + geom_point() +
    geom_abline(intercept=0, slope=1, color = "darkred") +
    labs(x='Valores observados', y='Valores ajustados', title='Ajustado vs. Observados') +
    theme_bw()
  print(obs_fitted)
  #-----------------------------------------------------------------#
}

