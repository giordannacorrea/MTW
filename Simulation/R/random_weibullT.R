
random_weibullT <- function(n, beta, nu, X, Z, link.mu,link.phi){
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
  ###---------------------------------------------------------------------------
  r_W <- function(alpha, lambda, N){
    u <- runif(N)
    y <- -1 + lambda*(  -log(1 - u*(1 - exp(  -(2/lambda)^alpha ))) )^(1/alpha)
    y
  }
  # funcao para gerar os dados da distrbuicao weibull
  simu.WeibullT <- function(n, beta, nu, X, Z){
    eta = X%*%as.matrix(beta)
    if(link.mu == "identity"){
      m <-link1$linkinv(eta)
    }else{
      m <-2*link1$linkinv(eta)-1
    }
    eta2=Z%*%as.matrix(nu)
    phi= link2$linkinv(eta2)
    lamb = (m+1)*(phi/(phi-1))^(1/phi)
    y <- r_W(phi,lamb, n)
    
    return(y)
  }
  ###---------------------------------------------------------------------------
  p <- length(beta)
  q <- length(nu)
  ###---------------------------------------------------------------------------
  y=simu.WeibullT(n, beta, nu, X, Z)
  return(y)
}

#-------------------------------------------------------------------------------
#set.seed(2)
#link.mu = "logit"
#link.phi = "log"
#y1 <- random_weibullT(30,c(-1,0.5),c(1.4,0.4),"logit","log")
#hist(y1)
#-------------------------------------------------------------------------------
#set.seed(2)
#y2 <- random_weibullT(30,c(-1,0.5),c(1.2,1),"logit","log")
#hist(y2)
#-------------------------------------------------------------------------------
#set.seed(2)
#y3 <- random_weibullT(30,c(-1,0.5),c(1.7,0.3),"logit","log")
#hist(y3)
#-------------------------------------------------------------------------------
#set.seed(2)
#y4 <- random_weibullT(30,c(-0.8, 1),c(1,0.6),"logit","log")
#hist(y4)
#-------------------------------------------------------------------------------
#set.seed(2)
#y5 <- random_weibullT(30,c(0.7, -1.5),c(1.2,0.6),"logit","log")
#hist(y5)
#-------------------------------------------------------------------------------
#set.seed(2)
#y6 <- random_weibullT(30,c(1, -1.2),c(1.5,0.4),"logit","log")
#hist(y6)
#-------------------------------------------------------------------------------
#set.seed(2)
#y7 <- random_weibullT(30,c(0.6, -1),c(1.5,0.4),"logit","log")
#hist(y7)
#-------------------------------------------------------------------------------
#set.seed(2)
#y8 <- random_weibullT(30,c(0.6, -1),c(1.5,0.5),"logit","log")
#hist(y8)
#-------------------------------------------------------------------------------
#set.seed(2)
#y9<- random_weibullT(30,c(0.5, -1.5),c(1.5,0.5),"logit","log")
#hist(y9)
#-------------------------------------------------------------------------------
#set.seed(2)
#y10<- random_weibullT(30,c(2, -1),c(1.5,2),"logit","log")
#hist(y10)
#-------------------------------------------------------------------------------
#set.seed(2)
#y11<- random_weibullT(30,c(-1.5, 0.5),c(0.5,1),"logit","log")
#hist(y11)
#-------------------------------------------------------------------------------
#set.seed(2)
#y12<- random_weibullT(30,c(0.5, -1),c(1.5,0.5),"logit","log")
#hist(y12)


