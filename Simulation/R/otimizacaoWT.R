getwd()
setwd("C:/Users/giord/Dropbox/Orientação Giordanna/Simulações/nlminb")
#-------------------------------------------------------------------------------
library(tidyverse)
source("random_weibullT.R")
set.seed(2)
#beta<-c(-1, 0.5)
#nu<-c(1,0.5)
#n=50
#R=1000
#link.mu = "logit"
#link.phi="log"
#-------------------------------------------------------------------------------
teste_otimizacao <- function(n,beta,nu,link.mu,link.phi,R){
  ##----------------------- Funcao de ligacao para moda --------------------------
  if(any(link.mu == c("logit","probit","cloglog"))) {
    stats1 <- make.link(link.mu)
  }else stop(paste(link.mu, "link not available, available links are",
                   "\"log\" and \"identity\""))
  
  link1 <- structure(list(link = link.mu, 
                          linkfun = stats1$linkfun,
                          linkinv = stats1$linkinv, 
                          mu.eta = stats1$mu.eta, 
                          diflink = function(t) 1/(stats1$mu.eta(stats1$linkfun(t)))
  ))
  ##----------------------- Funcao de ligacao para dispersao ---------------------
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
  #-----------------------------------------------------------------------------
  derivada_gama <- function(beta, nu, X, Z, a, c){
    eta.m <- X%*%as.matrix(beta)
    eta.phi <- Z%*%as.matrix(nu)
    
    m <- 2*link1$linkinv(eta.m)-1
    phi <- link2$linkinv(eta.phi)
    
    funcao <- function(mode, prec, a1, c1){
      function(w){
        delta <- ((2/(mode+1))^prec) * ((prec-1)/prec)
        result <- (w * delta)^(a1-1) * ((log(w*delta))^c1) * exp(-w*delta) * delta
        return(result)
      }
    }
    # limites de integracao - LI (limite inferior) e LS (limite superior)
    LI= rep(0,length(m))
    LS = rep(1,length(m))
    L <- list(m = m, phi = phi, LI = LI, LS = LS)
    L2 <- L |> transpose()
    #
    unlist(lapply(L2, function(vars){
      integrate(funcao(mode = vars$m,prec =  vars$phi, a1 = a, c1 = c), lower =  vars$LI, upper =  vars$LS)$value
    }))
  }
  # ------------------------------- Funcao de Verossimilhanca --------------------
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
  # ------------------------ Matriz de Informacao de Fisher ----------------------
  Fisher <- function(theta, data){
    beta <- theta[1:p]
    nu<- theta[(p+1):length(theta)]
    #
    eta.m <- X%*%as.matrix(beta) # eta_1i = g1(m_i*) = X^T * beta
    m <- 2*link1$linkinv(eta.m)-1 # m_i = 2 * g1^(-1) (eta.m) - 1
    eta.phi <- Z%*%as.matrix(nu) # eta_2i = g2(phi_i) = Z^{T} * nu
    phi <- link2$linkinv(eta.phi) # phi_i = g2^(-1) (eta.phi)
    
    #-----------------------------------------------------------------------
    mT <- diag(as.vector( 2 * (1/link1$diflink((m+1)/2) ) )) # (dm_i/d eta_1i = 2/ g_1^{'} (m*)) 
    mC <- diag(as.vector( link2$mu.eta(eta.phi)   )  ) # (dphi_i/eta_2i) = 1/g_{2}^{'} (phi)
    delta <- ((2/(m+1))^phi) * ((phi-1)/phi)
    #-----------------------------------------------------------------------
    # Resultados importantes
    #-----------------------------------------------------------------------
    EY =  ((phi/(phi-1))/(1- exp(-delta))) * derivada_gama(beta, nu, X, Z, 2, 0)
    EYlog = (   (1/(phi-1))/ (1 - exp(-delta)))*( derivada_gama(beta, nu, X, Z, 2, 1) + log(phi/(phi-1)) * derivada_gama(beta, nu, X, Z, 2, 0) )
    EYlog2 = (  (  (1/phi)*(1/(phi-1)) )/(1- exp(-delta))  ) * ( 2 * log(phi/(phi-1)) * derivada_gama(beta, nu, X, Z, 2, 1) + (  log(phi/(phi-1))  )^(2) * derivada_gama(beta, nu, X, Z, 2, 0) + derivada_gama(beta, nu, X, Z, 2, 2) )
    #-----------------------------------------------------------------------
    # derivada de segunda ordem em funcao de m
    h1 <- (   exp(-delta) * (2/(m+1))^(2*phi) * ( (phi-1)/(m+1) )^(2)-
                exp(-delta) * (2/(m+1))^(phi) * (phi/(m+1)) * ((phi-1)/(m+1)) +
                exp(-2*delta) * (2/(m+1))^(phi) * (phi/(m+1)) * ((phi-1)/(m+1)) -
                exp(-delta) * (2/(m+1))^(phi) * ((phi-1)/((m+1)^2)) +
                exp(-2*delta) * (2/(m+1))^(phi) * ((phi-1)/((m+1)^2))
    )/ ( (1-exp(-delta))^2 )
    # derivada de segunda ordem em funcao de phi
    h3 <- -((
      -exp(-delta) * (2/(m+1))^(2*phi) * ((phi-1)/phi)^(2) * (log( 2/(m+1)))^(2) -
        2 * exp(-delta) * (2/(m+1))^(2*phi) * ((phi-1)/phi) * (1/(phi^2)) * log(2/(m+1)) +
        exp(-delta) * (2/(m+1))^(phi) * ((phi-1)/phi) * (log(2/(m+1)) )^(2) -
        exp(-2*delta) * (2/(m+1))^(phi) * ((phi-1)/phi) * (log(2/(m+1)) )^(2) +
        2 * exp(-delta) * (2/(m+1))^(phi) * (1/(phi^2)) * log(2/(m+1)) -
        2 * exp(-2*delta) * (2/(m+1))^(phi) * (1/(phi^2)) * log(2/(m+1)) -
        exp(-delta) * (2/(m+1))^(2*phi) * (1/(phi^2))^(2) -
        exp(-delta) * (2/(m+1))^(phi) * (2/(phi^3)) +
        exp(-2*delta) * (2/(m+1))^(phi) * (2/(phi^3))
    )/( ( 1-exp(-delta))^2 ) )
    # derivada cruzada
    h2 <- ( -exp(-delta) * (2/(m+1))^(2*phi) * log(2/(m+1)) * ((phi-1)/phi) * ((phi-1)/(m+1)) -
              exp(-delta) * (2/(m+1))^(2*phi) * (1/(phi^2)) * ((phi-1)/(m+1)) +
              exp(-delta) * (2/(m+1))^(phi) * log(2/(m+1)) * ((phi-1)/(m+1)) -
              exp(-2*delta) * (2/(m+1))^(phi) * log(2/(m+1)) * ((phi-1)/(m+1)) +
              exp(-delta) * (2/(m+1))^(phi) * (1/(m+1)) -
              exp(-2*delta) * (2/(m+1))^(phi) * (1/(m+1))
    )/( (1-exp(-delta))^2) 
    # ( - Esperanca das derivadas de segunda ordem)
    # eperanca em funcao de m
    ES <- as.vector(
      - (phi/((m+1)^2)) + ( (phi-1)/((m+1)^2)) * EY + ((phi-1)/(m+1))*(phi/(m+1)) * EY - h1
    )
    # esperanca da derivada cruzada
    ED <- as.vector(
      (1/(m+1)) - (1/(m+1)) * EY - ((phi-1)/(m+1)) * EYlog - h2
    )
    # esperanca de phi
    EV <- as.vector(
      (1/((phi-1)^2)) + ((phi-1)/phi) * EYlog2 + (2/(phi^2))* EYlog - (2/(phi^3)) * EY -  h3
    )
    #-----------------------------------------------------------------------
    H = diag(as.vector(ES))%*%(mT^2)
    C = mT%*%diag(as.vector(ED))%*%mC
    Q = diag(as.vector(EV))%*%(mC^2)
    
    Kbb <- t(X)%*%(H)%*%X
    Kbn <- t(X)%*%(C)%*%Z
    Knn <- t(Z)%*%(Q)%*%Z
    
    #---------------------------------------------------------------------------
    Ktheta<-rbind(
      cbind(Kbb,(Kbn)),
      cbind(t(Kbn),Knn)
    )
    
    return(Ktheta)
  }
  #-----------------------------------------------------------------------------
  vtheta=c(beta, nu)
  p <- length(beta)
  q <- length(nu)
  k=nTheta <- length(vtheta)
  mreplicas2 = matrix(0, nrow=R,ncol=k)
  contIC90 = contIC95 = contIC99 <-  matrix(0, nrow=R,ncol=k)
  z_c90 <-qnorm(1-0.1/2)
  z_c95 <-qnorm(1-0.05/2)
  z_c99 <-qnorm(1-0.01/2)
  #-----------------------------------------------------------------------------
  l=1 # inicia
  while (l <= R) {
    if(p==1){
      X<-cbind(rep(1,n))
    }else{
      X <- cbind(1, replicate(p-1, runif(n)))
    }
    
    if(q==1){
      Z = cbind(rep(1,n))
    }else{
      Z = cbind(1, replicate(q-1, runif(n)))
    }
    y=random_weibullT(n, beta, nu, X, Z, link.mu, link.phi)
    #hist(y)
    #---------------------------------------------------------------------------
    if(max(abs(y))!=1){
      ##------------------------ Chute Inicial----------------------------------
      ynew <- link1$linkfun( (y+1)/2)
      
      ajuste1 <- lm(ynew ~ X+0)
      mqo=ajuste1$coef # minimos quadrados ordinarios
      s2 <- var(ynew)
      ajuste2 = as.numeric(rep(s2,q)) 
      start = c(as.numeric(ajuste1$coefficients), ajuste2)
      #-------------------------------------------------------------------------
      fit2 <- try(nlminb(start, fn2, hessian = T, control = list(rel.tol = 1e-6) ), silent = T)
      
      TesteInvF <- try(solve(Fisher(fit2$par, data = y)), silent = T)
      if( sum(class(TesteInvF) != "try-error") == 2 ){
        if((class(fit2) != "try-error") ){
          if(   (fit2$convergence == 0  &  (sum(is.na(sqrt(diag(solve(Fisher(fit2$par,data=y))))) ) == 0)) ){
            mreplicas2[l,]<- fit2$par
            H1 <- sqrt(diag(solve(Fisher(fit2$par, data = y))) )
            
            IC90<-cbind(fit2$par-z_c90*H1,fit2$par+z_c90*H1)
            IC95<-cbind(fit2$par-z_c95*H1,fit2$par+z_c95*H1)
            IC99<-cbind(fit2$par-z_c99*H1,fit2$par+z_c99*H1)
            
            contIC90[l,] <- as.numeric(between(vtheta,IC90[,1],IC90[,2]))
            contIC95[l,] <- as.numeric(between(vtheta,IC95[,1],IC95[,2]))
            contIC99[l,] <- as.numeric(between(vtheta,IC99[,1],IC99[,2]))
            
            l <- l+1
            #print(l)
            if((100*l/R)%%10 == 0)
              print(c(100*l/R, "%"), quote = F)
          }
        }
      } # fim 2 if
    } # fim 1 if
  }  # fim do while
  
  res <- list()
  res$nlminb<-mreplicas2
  vies_nlminb <- 100*((apply(mreplicas2,2,mean)-vtheta)/vtheta)
  media_nlminb <- apply(mreplicas2,2,mean)
  ViesR <- cbind(vies_nlminb)
  media <- cbind(media_nlminb)
  res$ViesR<-ViesR
  res$mediaMC <- media
  
  TC90<-apply(contIC90,2,mean) # taxa de aceitacao
  TC95<-apply(contIC95,2,mean) # taxa de aceitacao
  TC99<-apply(contIC99,2,mean) # taxa de aceitacao
  TC <- cbind(TC90, TC95, TC99)
  res$TC <- TC
  return(res)
  
}
#set.seed(2)
#a1 <- teste_otimizacao(20,c(-1, 0.5),c(1,0.5),link.mu = "logit",link.phi = "log", 10000)
#a2 <- teste_otimizacao(30,c(-1, 0.5),c(1,0.5),link.mu = "logit",link.phi = "log", 10000)
#a3 <- teste_otimizacao(40,c(-1, 0.5),c(1,0.5),link.mu = "logit",link.phi = "log", 10000)
#a4 <- teste_otimizacao(50,c(-1, 0.5),c(1,0.5),link.mu = "logit",link.phi = "log", 10000)
#a5 <- teste_otimizacao(100,c(-1, 0.5),c(1,0.5),link.mu = "logit",link.phi = "log", 5000)
#a6 <- teste_otimizacao(200,c(-1, 0.5),c(1,0.5),link.mu = "logit",link.phi = "log", 5000)
#a7 <- teste_otimizacao(500,c(-1, 0.5),c(1,0.5),link.mu = "logit",link.phi = "log", 1000)

#a1
#a2
#a3
#a4
#a5
#a6
#a7

#save.image("a_logit_log_Simulacao_R10000")

