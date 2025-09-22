#-------------------------------------------------------------------------------
# Truncated Weibull distribution functions
#-------------------------------------------------------------------------------

# Probability Density Function (PDF)
dWeibullT <- function(y,m,phi){
  num <- (phi/(m+1))*((phi-1)/phi)*(( (y+1)/(m+1))^(phi-1))*exp(-((phi-1)/phi)*(( (y+1)/(m+1))^phi))
  den <- 1-exp( -((2/(m+1))^phi)*((phi-1)/phi))
  result <- num/den
  return(result)
}
#------------------------------------------------------------------------------#
# Simulation function for the Truncated Weibull distribution
rWtruncated <- function(Nsim,alfa,lamb){
  u <- runif(Nsim)
  y <- -1 + lamb*(  -log(1 - u*(1 - exp(  -(2/lamb)^alfa ))) )^(1/alfa)
  y
}
#------------------------------------------------------------------------------#
# Simulator for parameters M and PHI
rWeibullT <- function(par, Nsim){
  m <- par[1]
  phi <- par[2]
  
  u <- runif(Nsim)
  f <- -1 + ( (m+1)*(phi/(phi-1))^(1/phi)) * (-log(1 - u*(1 - exp(  -(2/(m+1))^(phi) * ((phi-1)/phi)))))^(1/phi)
  return(f)
}
#------------------------------------------------------------------------------#
# Cumulative Distribution
pWeibullT <- function(x, m, phi){
  k1 <- ((phi-1)/phi) * (((x+1)/(m+1))^phi)
  k2 <- ((phi-1)/phi) * ((2/(m+1))^phi)
  num <- 1 - exp(-k1)
  den <- 1 - exp(-k2)
  result <- num/den
  return(result)
}
#-------------------------------------------------------------------------------
pWeibull <- function(x, alfa, lamb){
  num = 1 - exp(-(((x+1)/lamb)^alfa))
  den = 1 - exp(-((2/lamb)^alfa))
  resultado = num/den
  return(resultado)
}
#-------------------------------------------------------------------------------
Hessian <- function(par, x){
  m <- par[1]
  phi <- par[2]
  #
  n <- length(x)
  #
  k <- ((2/(m+1))^(phi) * ((phi-1)/phi))
  # derivada m
  h1 <- (   exp(-k) * (2/(m+1))^(2*phi) * ( ((phi-1)/(m+1)) )^(2)-
               exp(-k) * (2/(m+1))^(phi) * (phi/(m+1)) * ((phi-1)/(m+1)) +
               exp(-2*k) * (2/(m+1))^(phi) * (phi/(m+1)) * ((phi-1)/(m+1)) -
               exp(-k) * (2/(m+1))^(phi) * ((phi-1)/((m+1)^2)) +
               exp(-2*k) * (2/(m+1))^(phi) * ((phi-1)/((m+1)^2))
  )/ ((1-exp(-k))^2)
  h11 <- n * (phi/((m+1)^2)) - ((phi-1)/((m+1)^2)) * sum( ((x+1)/(m+1))^(phi)) - ((phi-1)/(m+1)) * (phi/(m+1)) * sum( ((x+1)/(m+1))^(phi) ) + n * h1
  # derivada phi
  h2 <- -((
    -exp(-k) * (2/(m+1))^(2*phi) * ((phi-1)/phi)^(2) * (log( 2/(m+1)))^(2) -
      2 * exp(-k) * (2/(m+1))^(2*phi) * ((phi-1)/phi) * (1/(phi^2)) * log(2/(m+1)) +
      exp(-k) * (2/(m+1))^(phi) * ((phi-1)/phi) * (log(2/(m+1)) )^(2) -
      exp(-2*k) * (2/(m+1))^(phi) * ((phi-1)/phi) * (log(2/(m+1)) )^(2) +
      2 * exp(-k) * (2/(m+1))^(phi) * (1/(phi^2)) * log(2/(m+1)) -
      2 * exp(-2*k) * (2/(m+1))^(phi) * (1/(phi^2)) * log(2/(m+1)) -
      exp(-k) * (2/(m+1))^(2*phi) * (1/(phi^2))^(2) -
      exp(-k) * (2/(m+1))^(phi) * (2/(phi^3)) +
      exp(-2*k) * (2/(m+1))^(phi) * (2/(phi^3))
  )/((1-exp(-k))^2))
  h22 <- -(n/((phi-1)^2)) + (2/(phi^3)) * sum( ((x+1)/(m+1))^(phi)  ) - 2 * (1/(phi^2)) * sum( ((x+1)/(m+1))^(phi) * log(  ((x+1)/(m+1))  )) -
    ((phi-1)/phi) * sum(  ((x+1)/(m+1))^(phi) * (log(   ((x+1)/(m+1))  ) )^(2)    ) + n* h2
  # derivada cruzada
  hc <- ( -exp(-k) * (2/(m+1))^(2*phi) * log(2/(m+1)) * ((phi-1)/phi) * ((phi-1)/(m+1)) -
            exp(-k) * (2/(m+1))^(2*phi) * (1/(phi^2)) * ((phi-1)/(m+1)) +
            exp(-k) * (2/(m+1))^(phi) * log(2/(m+1)) * ((phi-1)/(m+1)) -
            exp(-2*k) * (2/(m+1))^(phi) * log(2/(m+1)) * ((phi-1)/(m+1)) +
            exp(-k) * (2/(m+1))^(phi) * (1/(m+1)) -
            exp(-2*k) * (2/(m+1))^(phi) * (1/(m+1))
  )/((1-exp(-k))^2)
  h12 = h21 <- - (n/(m+1)) + (1/(m+1)) * sum( ((x+1)/(m+1))^(phi) ) + ((phi-1)/(m+1)) * sum(  ((x+1)/(m+1))^(phi) * log(  ((x+1)/(m+1)))  ) + n * hc
  
  
  H<-matrix(0,2,2)
  
  H[1,1] <- h11
  H[2,2] <- h22
  H[1,2] <- h12
  H[2,1] <- h21
  
  return(H)
}
#-------------------------------------------------------------------------------
gama_derivada <- function(m, phi, a, c){
  delta = ((2/(m+1))^phi)*((phi-1)/phi)
  h<- list()
  g<- numeric()
  for(i in 1:length(m)){
    h[[i]] = function(t){
      (t^(a-1) )* (log(t)^c) * exp(-t)
    }
    g[i] = integrate(h[[i]], lower=0, upper=delta[i])$value
  }
  return(g)
}
# derivative of the gamma function
derivada_gama <- function(m, phi, a, c){
  funcao <- function(mode, prec, a1, c1){
    delta <- ((2/(mode+1))^prec) * ((prec-1)/prec)
    function(w){
      result <- ((w * delta)^(a1-1)) * ((log(w*delta))^c1) * exp(-w*delta) * delta
      return(result)
    }
  }
  # limits of integration - LI (lower limit) and LS (upper limit)
  LI=rep(0,length(m)) 
  LS = rep(1,length(m))
  L <- list(m = m, phi = phi, LI = LI, LS = LS)
  L2 <- L |> transpose()
  #
  unlist(lapply(L2, function(vars){
    integrate(funcao(vars$m, vars$phi, a1 = a, c1 = c), vars$LI, vars$LS)$value
  } ))
}
#-------------------------------------------------------------------------------
# Fisher Information
Fisher2 <- function(par, x){
  m <- par[1]
  phi <- par[2]
  k <-  ( (2/(m+1))^phi) * ((phi-1)/phi)
  
  # Expected value calculation
  EY =  ((phi/(phi-1))/(1- exp(-k))) * derivada_gama(m,phi, 2, 0)
  EYlog = (   (1/(phi-1))/ (1 - exp(-k)))*( derivada_gama(m, phi, 2, 1) + log(phi/(phi-1)) * derivada_gama(m,phi, 2, 0) )
  EYlog2 = (  (  (1/phi)*(1/(phi-1)) )/(1- exp(-k))  ) * ( 2 * log(phi/(phi-1)) * derivada_gama(m,phi, 2, 1) + (  log(phi/(phi-1))  )^(2) * derivada_gama(m,phi, 2, 0) + derivada_gama(m,phi, 2, 2) )
  
  #-----------------------------------------------------------------------
  # Important results
  #-----------------------------------------------------------------------
  # derivative with respect to m
  h1 <- (   exp(-k) * (2/(m+1))^(2*phi) * ( (phi-1)/(m+1) )^(2)-
              exp(-k) * (2/(m+1))^(phi) * (phi/(m+1)) * ((phi-1)/(m+1)) +
              exp(-2*k) * (2/(m+1))^(phi) * (phi/(m+1)) * ((phi-1)/(m+1)) -
              exp(-k) * (2/(m+1))^(phi) * ((phi-1)/((m+1)^2)) +
              exp(-2*k) * (2/(m+1))^(phi) * ((phi-1)/((m+1)^2))
  )/ ( (1-exp(-k))^2 )
  #
  I11 <-  - (phi/((m+1)^2)) + ( (phi-1)/((m+1)^2)) * EY + ((phi-1)/(m+1))*(phi/(m+1)) * EY - h1
  # derivative with respect to phi
  h2 <- -((
    -exp(-k) * (2/(m+1))^(2*phi) * ((phi-1)/phi)^(2) * (log( 2/(m+1)))^(2) -
      2 * exp(-k) * (2/(m+1))^(2*phi) * ((phi-1)/phi) * (1/(phi^2)) * log(2/(m+1)) +
      exp(-k) * (2/(m+1))^(phi) * ((phi-1)/phi) * (log(2/(m+1)) )^(2) -
      exp(-2*k) * (2/(m+1))^(phi) * ((phi-1)/phi) * (log(2/(m+1)) )^(2) +
      2 * exp(-k) * (2/(m+1))^(phi) * (1/(phi^2)) * log(2/(m+1)) -
      2 * exp(-2*k) * (2/(m+1))^(phi) * (1/(phi^2)) * log(2/(m+1)) -
      exp(-k) * (2/(m+1))^(2*phi) * (1/(phi^2))^(2) -
      exp(-k) * (2/(m+1))^(phi) * (2/(phi^3)) +
      exp(-2*k) * (2/(m+1))^(phi) * (2/(phi^3))
  )/((1-exp(-k))^(2)))
  #
  I22 <- (1/((phi-1)^2)) + ((phi-1)/phi) * EYlog2 + (2/(phi^2))* EYlog - (2/(phi^3)) * EY -  h2
  # mixed partial derivative
  hc <- ( -exp(-k) * (2/(m+1))^(2*phi) * log(2/(m+1)) * ((phi-1)/phi) * ((phi-1)/(m+1)) -
            exp(-k) * (2/(m+1))^(2*phi) * (1/(phi^2)) * ((phi-1)/(m+1)) +
            exp(-k) * (2/(m+1))^(phi) * log(2/(m+1)) * ((phi-1)/(m+1)) -
            exp(-2*k) * (2/(m+1))^(phi) * log(2/(m+1)) * ((phi-1)/(m+1)) +
            exp(-k) * (2/(m+1))^(phi) * (1/(m+1)) -
            exp(-2*k) * (2/(m+1))^(phi) * (1/(m+1))
  )/((1-exp(-k))^2)   
  #
  I12 = I21 <- (1/(m+1)) - (1/(m+1)) * EY - ((phi-1)/(m+1)) * EYlog - hc 
  
  I<-matrix(0,2,2)
  
  I[1,1] <- I11
  I[2,2] <- I22
  I[1,2] <- I12
  I[2,1] <- I21
  
  return(I)
}