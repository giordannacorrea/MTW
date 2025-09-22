comp.models <-  function(formula, link.mu, link.phi, data){
  #----------------------------------------------------------------------------#
  # FORMULA
  #----------------------------------------------------------------------------#
  form = Formula(formula)
  model = model.frame(form,data)
  y<-model.response(model) # variavel resposta do modelo
  X<-model.matrix(form,data=model, rhs = 1) # X
  Z<-model.matrix(form,data=model, rhs = 2) # Z 
  n <- length(y)
  p = ncol(X)
  q = ncol(Z)
  #----------------------------------------------------------------------------#
  # FUNCAO DE LIGACAO PARA A MODA
  #----------------------------------------------------------------------------#
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
  #----------------------------------------------------------------------------#
  # FUNCAO DE LIGACAO PARA PHI
  #----------------------------------------------------------------------------#
  if(any(link.phi == c("log","sqrt","identity", "1/mu^2"))) {
    stats2 <- make.link(link.phi)
  }else stop(paste(link.mu, "link not available, available links are",
                   "\"logit\", \"probit\". \"cauchit\" and \"cloglog\"" ))
  
  link2 <- structure(list(link = link.phi, 
                          linkfun = stats2$linkfun,
                          linkinv = stats2$linkinv, 
                          mu.eta = stats2$mu.eta, 
                          diflink = function(t) 1/(stats2$mu.eta(stats1$linkfun(t)))
  ))
  
  #Weibull  
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
    return(log.likelihood)
  }
  #-------------------------------------------------------------------------------
  # BETA
  fn_beta <- function(theta){
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
    m <- (m+1)/2
    
    log.likelihood <- suppressWarnings(sum(
      m * phi * log((y+1)/2) + (1-m)* phi * log(1- ((y+1)/2) ) - lbeta( m*phi+1, (1-m)*phi +1 )-log(2)
    ))
    return(log.likelihood)
  }
  #-------------------------------------------------------------------------------
  # KUMARASWAMY
  fn_kum <- function(theta){
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
    m <- (m+1)/2
    
    log.likelihood <- suppressWarnings(sum(
      log(  1 + (phi-1)*(m^(-phi))  ) + (phi-1)*log((y+1)/2) + (  (phi^(-1)) * (1+  (m^(-phi)) * (phi-1) )  -1) * log(1 - ((y+1)/2)^phi )-log(2)
    ))
    return(log.likelihood)
  }
  #-------------------------------------------------------------------------------
  # UNIT-GAMMA
  fn_ugamma <- function(theta){
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
    m <- (m+1)/2
    
    log.likelihood <- suppressWarnings(sum(
      phi *  log( (1-phi + log(m)) * (log(m)^(-1))) + (   (1+log(m)-phi) * (log(m)^(-1))  -1) * log( (y+1)/2) + (phi-1) * log(-log( (y+1)/2 )) - lgamma(phi) - log(2)
    ))
    return(log.likelihood)
  }
  #-------------------------------------------------------------------------------
  # unit gompertz
  fn_ugomp <- function(theta){
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
    m <- (m+1)/2
    log.likelihood <- suppressWarnings(sum(
      phi * log(m) + log(phi+1) - (phi+1) * log((y+1)/2) - (phi^(-1)) * (m^phi) * (phi+1) * ( ( ((y+1)/2)^(-phi) )  -1) - log(2)
    ))
    return(log.likelihood)
  }
  #### Estimação
  #-------------------------------------------------------------------------------
  # CHUTE INICIAL
  ynew <- link1$linkfun( (y+1)/2)
  ajuste1 <- lm(ynew ~ X+0) 
  mqo=ajuste1$coef # minimos quadrados ordinarios
  s2 <- 2
  ajuste2 = as.numeric(rep(s2,q))
  chute_inicial1 <- c(1.4110,7.9977,-2.5494,15.1587,-14.2304)
  chute_inicial2 = c(as.numeric(ajuste1$coefficients), ajuste2)
  #-------------------------------------------------------------------------------
  #Weibull
  fitWTM <- maxLik(logLik = fn2, start = chute_inicial2)
  ajusteWTM=summary(fitWTM)
  coef <- (ajusteWTM$estimate)[1:(ncol(Z)+ncol(X))]
  beta <- coef[1:ncol(X)]
  nu<- coef[(ncol(X)+1):length(coef)]
  AICw=-2*(-fitWTM$maximum)+2*(length(nu)+length(beta));AICw
  BICw = -2*(-fitWTM$maximum)+log(n)*(length(nu)+length(beta));BICw
  #-------------------------------------------------------------------------------
  #beta
  fitBETA <- maxLik(fn_beta,start = chute_inicial2)
  ajusteBETA=summary(fitBETA)
  coef <- (fitBETA$estimate)[1:(ncol(Z)+ncol(X))]
  beta <- coef[1:ncol(X)] 
  nu<- coef[(ncol(X)+1):length(coef)]
  AICb=-2*(-fitBETA$maximum)+2*(length(nu)+length(beta));AICb
  BICb=-2*(-fitBETA$maximum)+log(n)*(length(nu)+length(beta));BICb
  #-------------------------------------------------------------------------------
  #Kumaraswamy
  fitKUMA <- maxLik(fn_kum,start=chute_inicial1)
  ajusteKUMA=summary(fitKUMA)
  coef <- (fitKUMA$estimate)[1:(ncol(Z)+ncol(X))]
  beta <- coef[1:ncol(X)] 
  nu<- coef[(ncol(X)+1):length(coef)]
  AICk=-2*(-fitKUMA$maximum)+2*(length(nu)+length(beta));AICk
  BICk=-2*(-fitKUMA$maximum)+log(n)*(length(nu)+length(beta));BICk
  
  #-------------------------------------------------------------------------------
  #Unit Gompertz
  
  fitGZ <- maxLik(fn_ugomp, start = chute_inicial1)
  ajusteGZ = summary(fitGZ)
  coef <- (fitGZ$estimate)[1:(ncol(Z)+ncol(X))]
  beta <- coef[1:ncol(X)] 
  nu<- coef[(ncol(X)+1):length(coef)]
  AICgz=-2*(-fitGZ$maximum)+2*(length(nu)+length(beta));AICgz
  BICgz=-2*(-fitGZ$maximum)+log(n)*(length(nu)+length(beta));BICgz
  #-------------------------------------------------------------------------------
  #Unit gama
  
  fitGA <- maxLik(fn_ugamma,start = chute_inicial2)
  ajusteGA = summary(fitGA)
  coef <- (fitGA$estimate)[1:(ncol(Z)+ncol(X))]
  beta <- coef[1:ncol(X)] 
  nu<- coef[(ncol(X)+1):length(coef)]
  AICga=-2*(-fitGA$maximum)+2*(length(nu)+length(beta));AICga
  BICga=-2*(-fitGA$maximum)+log(n)*(length(nu)+length(beta));BICga
  #-------------------------------------------------------------------------------
  #Densidades
  dWeibullT <- function(y,m,phi){
    num <- (phi/(m+1))*((phi-1)/phi)*(( (y+1)/(m+1))^(phi-1))*exp(-((phi-1)/phi)*(( (y+1)/(m+1))^phi))
    den <- 1-exp( -((2/(m+1))^phi)*((phi-1)/phi))
    result <- num/den
    return(result)
  }
  #-------------------------------------------------------------------------------
  # BETA
  dbeta_V <- function(y,m,phi){
    m <- (m+1)/2
    c1 <- ( ((y+1)/2)^(m*phi) )* ( (1-((y+1)/2) )^((1-m)*phi) )
    c2 <- beta(m*phi+1,(1-m)*phi+1)
    result <- (1/2)*(c1/c2)
    return(result)
  }
  #-------------------------------------------------------------------------------
  # KUMARASWAMY
  dkuma_V <- function(y,m,phi){
    m <- (m+1)/2
    c1 <- (  1+(phi-1)*(m^(-phi))  )*(  ((y+1)/2)^(phi-1) )
    c2 <- (1-((y+1)/2)^phi)
    c3<-   ((phi^(-1))*(  1+(m^(-phi))*(phi-1)   )-1)
    result <- (1/2)*c1*(c2^c3)
    return(result)
  }
  #-------------------------------------------------------------------------------
  # GOMPERTZ
  dgomp_V <- function(y,m,phi){
    m <- (m+1)/2
    c1 <- (m^phi)*(phi+1)*( ((y+1)/2)^(-phi-1))
    c2 <- exp( (-phi^(-1))*(m^phi)*(phi+1)* (  ((y+1)/2)^(-phi)-1)  )
    result <- (1/2)*c1*c2
    return(result)
  }
  # UNIT GAMMA
  dUgamma_V <- function(y,m,phi){
    m <- (m+1)/2
    c1 <- (( (1 - phi  + log(m)  )*(log(m)^(-1))    )^phi) * ( ((y+1)/2)^(  (1-phi+log(m))*(log(m)^(-1))  -1) )
    c2 <- ((-log((y+1)/2))^(phi-1))
    c3 <- gamma(phi)
    result <- (1/2)*((c1*c2)/c3)
    return(result)
  }
  # LOG DA UNIT GAMA
  dUgamma_V2 <- function(y,m,phi){
    m <- (m+1)/2
    c1<- -lgamma(phi)+phi*log( (1+log(m)-phi)/log(m)) - log(2)
    c2<- ((1-phi)/log(m))*log((y+1)/2)+(phi-1)*log(-log((y+1)/2))
    result <- (c1+c2)
    return(result)
  }
  #-----------------------------------------------------------------------------
  # TESTE VUONG
  #-------------------------------------------------------------------------------
  f_weibul <- dWeibullT(y,2*link1$linkinv(X%*%fitWTM$estimate[1:p])-1,link2$linkinv(Z%*%fitWTM$estimate[-c(1:p)])   )
  f_beta   <- dbeta_V(y,2*link1$linkinv(X%*%fitBETA$estimate[1:p])-1,link2$linkinv(Z%*%fitBETA$estimate[-c(1:p)]))
  f_kuma   <- dkuma_V(y,2*link1$linkinv(X%*%fitKUMA$estimate[1:p])-1,link2$linkinv(Z%*%fitKUMA$estimate[-c(1:p)])   )
  f_gomp   <- dgomp_V(y,2*link1$linkinv(X%*%fitGZ$estimate[1:p])-1,link2$linkinv(Z%*%fitGZ$estimate[-c(1:p)])   )
  
  #Vuong Weibull Gama unit
  lgamunit<- dUgamma_V2(y,2*link1$linkinv(X%*%fitGA$estimate[1:p])-1,link2$linkinv(Z%*%fitGA$estimate[-c(1:p)]))
  w2 <- (1/n)*sum((log(f_weibul)-lgamunit)^2)- ((1/n)*sum(log(f_weibul)-lgamunit))^2
  T_w_g <- (1/(w2*sqrt(n)))*sum(log(f_weibul)-lgamunit)
  
  T_w_g
  voung_test <- function(f,g){
    n <- length(f)
    w2 <- mean( log(f/g)^2) - (mean(log(f/g)))^2
    temp = (1/(w2*sqrt(n))  )*sum(log(f/g))
    return(temp)
  }
  #Apenas para a gama unitária
  voung_test2 <- function(f,g){
    n <- length(f)
    w2 <- mean( (log(f)-g)^2) - (mean(log(f)-g))^2
    temp = (1/(w2*sqrt(n))  )*sum(log(f)-g)
    return(temp)
  }
  # COMAPARACAO DE TODOS OS MODELOS
  # COMPARACAO DA WEIBULL
  WB <- voung_test(f_weibul, f_beta) # Weibull ganhou
  WK <- voung_test(f_weibul, f_kuma) # Weibull ganhou
  WGz<-voung_test(f_weibul, f_gomp) # Weibull ganhou
  WGa<-voung_test2(f_weibul, lgamunit) # Weibull ganhou
  # COMPARACOES COM A BETA
  BK <- voung_test(f_beta, f_kuma)
  BGz <- voung_test(f_beta, f_gomp)
  BGa <- voung_test2(f_beta, lgamunit)
  # Comparacoes com a Kumaraswam
  KGz <- voung_test(f_kuma, f_gomp)
  KGa <- voung_test2(f_kuma, lgamunit)
  # Comparacoes com a Gompertz unitaria
  GzGa <- voung_test2(f_gomp, lgamunit)
  #-------------------------------------------------------------------------------
  # lof(f(y)) sem estrutura de regrassao para calculo do R^{2} Cox-Snell
  #-------------------------------------------------------------------------------
  #CHUTE INICAL
  #-------------------------------------------------------------------------------
  # WEIBULL
  lveroWTM <- function(par){
    alpha = par[1]
    lambda = par[2]
    ll <- suppressWarnings(
      log(alpha) - log(lambda) + (alpha-1) * log((y+1)/lambda) - (((y+1)/lambda)^alpha) - log(  1 - exp(-((2/lambda)^alpha)) )
    )
    sum(ll)
  }
  fitLveroWTM <- maxLik(lveroWTM,start = c(2,2),method = "BFGS")
  summary(fitLveroWTM)
  # COX-SNELL
  R2_WTM<- 1-exp(-(2/n)*(-(-fitWTM$maximum)-(-fitLveroWTM$maximum)))
  #-------------------------------------------------------------------------------
  # BETA
  #-------------------------------------------------------------------------------
  # a, b > 0, Y = 2X - 1; X = (Y+1)/2; |J|=1/2
  lveroBETA <- function(par){
    a = par[1]
    b = par[2]
    ll <- suppressWarnings(
      (a-1) * log((y+1)/2) + (b-1)*log(1-((y+1)/2)) - lbeta(a,b) - log(2)
    )
    sum(ll)
  }
  fitLveroBETA <- maxLik(lveroBETA,start = c(2,5), method = "BFGS")
  summary(fitLveroBETA)
  # COX-SNELL
  R2_BETA<- 1-exp(-(2/n)*(-(-fitBETA$maximum)-(-fitLveroBETA$maximum)))
  #-------------------------------------------------------------------------------
  # KUMARASWAMY
  # a, b > 0, Y = 2X - 1; X = (Y+1)/2; |J|=1/2
  lveroKUMA <- function(par){
    a = par[1]
    b = par[2]
    ll <- suppressWarnings(
      log(a*b) + (a-1)*log((y+1)/2) + (b-1)*log(1-( ((y+1)/2)^a)) - log(2)
    )
    sum(ll)
  }
  fitLveroKUMA <- maxLik(lveroKUMA,start = c(2,7), method = "BFGS")
  summary(fitLveroKUMA)
  # COX-SNELL
  R2_KUMA<- 1-exp(-(2/n)*(-(-fitKUMA$maximum)-(-fitLveroKUMA$maximum)))
  #-------------------------------------------------------------------------------
  # GAMA UNITARIA
  # a, b > 0, Y = 2X - 1; X = (Y+1)/2; |J|=1/2
  lveroGA <- function(par){
    a = par[1]
    b = par[2]
    ll <- suppressWarnings(
      a * log(b) - lgamma(a) + (b-1)*log((y+1)/2) + (a-1)*log(-log((y+1)/2)) - log(2)
    )
    sum(ll)
  }
  fitLveroGA <- maxLik(lveroGA,start = c(20,30), method = "BFGS")
  summary(fitLveroGA)
  # COX-SNELL
  R2_GA<- 1-exp(-(2/n)*(-(-fitGA$maximum)-(-fitLveroGA$maximum)))
  #-------------------------------------------------------------------------------
  # GOMPERTZ UNITARIA
  # a, b > 0, Y = 2X - 1; X = (Y+1)/2; |J|=1/2
  lveroGZ <- function(par){
    a = par[1]
    b = par[2]
    ll <- suppressWarnings(
      log(a*b) - (b+1)*log((y+1)/2) - a * ( ((y+1)/2)^(-b)  -1 ) - log(2)
    )
    sum(ll)
  }
  fitLveroGZ <- maxLik(lveroGZ,start = c(2,2),method = "BFGS")
  summary(fitLveroGZ)
  # COX-SNELL
  R2_GZ<- 1-exp(-(2/n)*(-(-fitGZ$maximum)-(-fitLveroGZ$maximum)))
  #-------------------------------------------------------------------------------
  # Residuos
  #-------------------------------------------------------------------------------
  # Resíduo Weibull Truncada
  m_WTM <- 2*link1$linkinv(X%*%fitWTM$estimate[1:p])-1
  phi_WTM <- link2$linkinv(Z%*%fitWTM$estimate[-c(1:p)])
  a <- alpha <- phi_WTM
  b <- lambda <- (m_WTM + 1)*(phi_WTM/(phi_WTM-1))^(1/phi_WTM)
  
  k1 <- ((y+1)/b)^a
  k2 <- (2/b)^a
  num <- 1 - exp(-k1)
  den <- 1 - exp(-k2)
  FWT <- num/den
  re_WT <- qnorm(FWT)
  hist(re_WT,prob=T)
  curve(dnorm(x),add=T, col = "darkred")
  WTM_g <- hnp(re_WT,pch=16,halfnormal=T)
  x_weibull <- WTM_g$x
  lower_weibull <- WTM_g$lower
  upper_weibull <- WTM_g$upper
  med_weibull <- WTM_g$median
  res_weibull <- WTM_g$residuals
  data_resWTM <- data.frame(x = x_weibull, lower = lower_weibull, upper = upper_weibull, median = med_weibull, residuals = res_weibull)
  weibull_hnp <- ggplot(data = data_resWTM, aes(x)) +
    geom_point(aes(y = residuals)) +
    geom_line(aes(y = lower)) +
    geom_line(aes(y = upper)) +
    geom_line(aes(y = median), linetype = "dashed", col = "darkred") +
    xlab("Quantil Teórico") + 
    ylab("Resíduos") +
    theme_bw()
  print(weibull_hnp)
  #-------------------------------------------------------------------------------
  #Resíduo kumaraswamy
  m_kuma <-link1$linkinv(X%*%fitKUMA$estimate[1:p])
  phi_kuma <- link2$linkinv(Z%*%fitKUMA$estimate[-c(1:p)])
  b <- (phi_kuma^(-1))*(1+(m_kuma^(-phi_kuma))*(phi_kuma-1)  )
  a <- phi_kuma
  re_kuma <- qnorm(1-(1-((y+1)/2)^a)^b)
  hist(re_kuma,prob=T)
  curve(dnorm(x),add=T,  col = "darkred", lwd = 2)
  kuma_g <- hnp(re_kuma,pch=16,halfnormal=T)
  x_kuma <- kuma_g$x
  lower_kuma <- kuma_g$lower
  upper_kuma <- kuma_g$upper
  med_kuma <- kuma_g$median
  res_kuma <- kuma_g$residuals
  data_reskuma <- data.frame(x = x_kuma, lower = lower_kuma, upper = upper_kuma, median = med_kuma, residuals = res_kuma)
  kuma_hnp <- ggplot(data = data_reskuma, aes(x)) +
    geom_point(aes(y = residuals)) +
    geom_line(aes(y = lower)) +
    geom_line(aes(y = upper)) +
    geom_line(aes(y = median), linetype = "dashed", col = "darkred") +
    xlab("Quantil Teórico") + 
    ylab("Resíduos") +
    theme_bw()
  print(kuma_hnp)
  #-------------------------------------------------------------------------------
  # gompertz
  m_ugomp <- link1$linkinv(X%*%fitGZ$estimate[1:p])
  phi_ugomp <- link2$linkinv(Z%*%fitGZ$estimate[-c(1:p)])
  b <- phi_ugomp
  a <- (phi_ugomp^(-1)) * (m_ugomp^(phi_ugomp)) * (phi_ugomp+1)
  re_ugomp <- qnorm(  exp(- a * ( (((y+1)/2)^(-b)) -1 )  )      )
  hist(re_ugomp,prob=T)
  curve(dnorm(x),add=T, col = "darkred", lwd = 2)
  gomp_g <- hnp(re_ugomp,halfnormal=T,pch=16, how.many.out = T, print.on = T)
  x_gomp <- gomp_g$x
  lower_gomp <- gomp_g$lower
  upper_gomp <- gomp_g$upper
  med_gomp <- gomp_g$median
  res_gomp <- gomp_g$residuals
  data_resgomp <- data.frame(x = x_gomp, lower = lower_gomp, upper = upper_gomp, median = med_gomp, residuals = res_gomp)
  gomp_hnp <- ggplot(data = data_resgomp, aes(x)) +
    geom_point(aes(y = residuals)) +
    geom_line(aes(y = lower)) +
    geom_line(aes(y = upper)) +
    geom_line(aes(y = median), linetype = "dashed", col = "darkred") +
    xlab("Quantil Teórico") + 
    ylab("Resíduos") +
    theme_bw()
  print(gomp_hnp)
  #-------------------------------------------------------------------------------
  # Beta
  m_beta <- link1$linkinv(X%*%fitBETA$estimate[1:p])
  phi_beta <- link2$linkinv(Z%*%fitBETA$estimate[-c(1:p)])
  a <- 1+m_beta*phi_beta
  b <- 1+phi_beta*(1-m_beta)
  res_beta <- qnorm(pbeta((y+1)/2,a,b))
  hist(res_beta, probability = T)
  curve(dnorm(x),add=T, col = "darkred", lwd = 2)
  #hnp(res_beta,pch=16,halfnormal = F)
  beta_g <- hnp(res_beta,pch=16,halfnormal = T)
  x_beta <- beta_g$x
  lower_beta <- beta_g$lower
  upper_beta <- beta_g$upper
  med_beta <- beta_g$median
  res_beta <- beta_g$residuals
  data_res <- data.frame(x = x_beta, lower = lower_beta, upper = upper_beta, median = med_beta, residuals = res_beta)
  beta_hnp <- ggplot(data = data_res, aes(x)) +
    geom_point(aes(y = residuals)) +
    geom_line(aes(y = lower)) +
    geom_line(aes(y = upper)) +
    geom_line(aes(y = median), linetype = "dashed", col = "darkred") +
    xlab("Quantil Teórico") + 
    ylab("Resíduos") +
    theme_bw()
  print(beta_hnp)
  #-------------------------------------------------------------------------------
  # Unit gama
  m_gama <- link1$linkinv(X%*%fitGA$estimate[1:p])
  phi_gama <- link2$linkinv(Z%*%fitGA$estimate[-c(1:p)])
  a <-phi_gama
  b <- (1+log(m_gama)-phi_gama)*(log(m_gama)^(-1))
  FUga <- pgamma(b * (-log((y+1)/2)),a)
  res_Uga <- qnorm(FUga)
  hist(res_Uga, probability = T)
  curve(dnorm(x),add=T, col = "darkred", lwd = 2)
  gama_g <- hnp(res_Uga,pch=16,halfnormal = T, how.many.out = T, print.on = T)
  x_gama <- gama_g$x
  lower_gama <- gama_g$lower
  upper_gama <- gama_g$upper
  med_gama <- gama_g$median
  res_gama <- gama_g$residuals
  data_resGA <- data.frame(x = x_gama, lower = lower_gama, upper = upper_gama, median = med_gama, residuals = res_gama)
  gama_hnp <- ggplot(data = data_resGA, aes(x)) +
    geom_point(aes(y = residuals)) +
    geom_line(aes(y = lower)) +
    geom_line(aes(y = upper)) +
    geom_line(aes(y = median), linetype = "dashed", col = "darkred") +
    xlab("Quantil Teórico") + 
    ylab("Resíduos") +
    theme_bw()
  print(gama_hnp)
  #-----------------------------------------------------------------------------
  z <- list()
  #-----------------------------------------------------------------------------
  Tab<-rbind(c(R2_GZ, R2_KUMA, R2_WTM, R2_GA, R2_BETA),
             c(AICgz, AICk, AICw, AICga, AICb),
             c(BICgz, BICk, BICw, BICga, BICb))
  Tab
  colnames(Tab)<-c("Gompertz","Kuma","WeibulTM","U-Gama","Beta")
  rownames(Tab) <- c("R2","AIC", "BIC")
  Tab
  mresults<-cbind(c(fitGZ$estimate,  R2_GZ,AICgz,BICgz),
                  c(fitKUMA$estimate,R2_KUMA,AICk, BICk),
                  c(fitWTM$estimate, R2_WTM,AICw, BICw),
                  c(fitGA$estimate,  R2_GA,AICga,BICga),
                  c(fitBETA$estimate,R2_BETA,AICb, BICb))
  #nomes das coluna
  rownames(mresults)<-c(c(paste("beta",1:ncol( as.matrix(X) ),sep="")), c(paste("nu",1:ncol( as.matrix(Z) ),sep="")),
                        "R2","AIC", "BIC")
  #nomes das linhas
  colnames(mresults)<-c("U-Gompertz",  "Kumaraswamy","WeibulTM","U-Gama", "Beta")
  
  #---------------------------------
  # P-VALOR
  pvalueWB <-(1-pnorm(abs(WB)))
  pvalueWK <-(1-pnorm(abs(WK)))
  pvalueWGz<-(1-pnorm(abs(WGz)))
  pvalueWGa<-(1-pnorm(abs(WGa)))
  #---------------------------------
  pvalueBK <-(1-pnorm(abs(BK)))
  pvalueBGz<-(1-pnorm(abs(BGz)))
  pvalueBGa<-(1-pnorm(abs(BGa)))
  #---------------------------------
  pvalueKGz<-(1-pnorm(abs(KGz)))
  pvalueKGa<-(1-pnorm(abs(KGa)))
  #---------------------------------
  pvalueGzGa<-(1-pnorm(abs(GzGa)))
  #---------------------------------
  comp = cbind(c(WB, WK, WGz, WGa,
                 BK, BGz, BGa,
                 KGz, KGa,
                 GzGa),
               round(c(pvalueWB, pvalueWK, pvalueWGz, pvalueWGa,
                       pvalueBK, pvalueBGz, pvalueBGa,
                       pvalueKGz, pvalueKGa,
                       pvalueGzGa),4)
  )
  colnames(comp)<-c("Vuong", "p-valor")
  rownames(comp) <- c("Weibull versus Beta","Weibull versus Kuma", "Weibull versus UGompertz", "Weibull versus UGama",
                      "Beta versus Kuma", "Beta versus UGompertz", "Beta versus UGama",
                      "Kuma versus Ugompertz", "Kuma versus UGama",
                      "Ugompertz versus UGama")
  
  z$resumo <- mresults
  z$WTM <- summary(fitWTM)
  z$Beta <- summary(fitBETA)
  z$kuma <- summary(fitKUMA)
  z$UGa <- summary(fitGA)
  z$Ugz <- summary(fitGZ)
  z$comparacao <- comp
  
  return(z)
  #-------------------------------------------------------------------------------
}