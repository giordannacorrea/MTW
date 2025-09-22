library(Stat2Data)
library(hnp)
library(Formula)
library(tidyverse)

graficoWT <- function(model, link.mu, link.phi, data){
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
  #-------------------------------------------------------------------------------
  # Resíduo Weibull Truncada
  #m_WTM <- 2*link1$linkinv(X%*%model$Par_est[1:p])-1
  #phi_WTM <- link2$linkinv(Z%*%model$Par_est[-c(1:p)])
  #a <- alpha <- phi_WTM
  #b <- lambda <- (m_WTM + 1)*(phi_WTM/(phi_WTM-1))^(1/phi_WTM)
  lambda=model$lambda
  alpha = model$alpha
  y <- model$serie
  k1 <- ((y+1)/lambda)^alpha
  k2 <- (2/lambda)^alpha
  num <- 1 - exp(-k1)
  den <- 1 - exp(-k2)
  FWT <- num/den
  re_WT <- qnorm(FWT)
  hist(re_WT,prob=T)
  curve(dnorm(x),add=T, col = "darkred")
  WTM_g <- hnp(re_WT,pch=16,halfnormal=T, how.many.out = T, print.on = T)
  x_weibull <- WTM_g$x
  lower_weibull <- WTM_g$lower
  upper_weibull <- WTM_g$upper
  med_weibull <- WTM_g$median
  res_weibull <- WTM_g$residuals
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
  data_resWTM <- data.frame(x = x_weibull, lower = lower_weibull, upper = upper_weibull, median = med_weibull, residuals = res_weibull)
  data_resWTM$cores <- hnp_cores(data_resWTM)
  data_resWTM_texto <- hnp_texto(data_resWTM)
  
  weibull_hnp <- ggplot(data = data_resWTM, aes(x)) +
    geom_point(aes(y = residuals, color = cores), show.legend = F) +
    geom_text(x = 0, y = 3, label = data_resWTM_texto$Total, hjust = 0) +
    geom_text(x = 0, y = 2.8, label = data_resWTM_texto$Fora, hjust = 0) +
    scale_color_manual(values = c(Dentro = "black", Fora = "red"))+
    geom_line(aes(y = lower)) +
    geom_line(aes(y = upper)) +
    geom_line(aes(y = median), linetype = "dashed", col = "darkred") +
    xlab("Quantil Teórico") + 
    ylab("Resíduos") +
    theme_bw()
  
  print(weibull_hnp)
}