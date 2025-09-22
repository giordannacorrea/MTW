source("Simulation/R/otimizacaoWT.R")
source("Simulation/R/random_weibullT.R")

# Cenario B
#beta<-c(0.5, -1, -0.5)
#nu <- c(1.5,0.5)


set.seed(2)
e1 <- teste_otimizacao(20,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "logit",link.phi = "log", 10000)
e2 <- teste_otimizacao(30,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "logit",link.phi = "log", 10000)
e3 <- teste_otimizacao(40,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "logit",link.phi = "log", 10000)
e4 <- teste_otimizacao(50,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "logit",link.phi = "log", 10000)
e5 <- teste_otimizacao(100,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "logit",link.phi = "log", 10000)
e6 <- teste_otimizacao(200,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "logit",link.phi = "log", 10000)
e7 <- teste_otimizacao(500,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "logit",link.phi = "log", 10000)

e1
e2
e3
e4
e5
e6
e7
 
save.image("eNEW_logit_log_Simulacao_R10000.RData")
#-----------------------------------------------------------------------------

set.seed(2)
ee1 <- teste_otimizacao(20,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "probit",link.phi = "log", 10000)
ee2 <- teste_otimizacao(30,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "probit",link.phi = "log", 10000)
ee3 <- teste_otimizacao(40,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "probit",link.phi = "log", 10000)
ee4 <- teste_otimizacao(50,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "probit",link.phi = "log", 10000)
ee5 <- teste_otimizacao(100,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "probit",link.phi = "log", 10000)
ee6 <- teste_otimizacao(200,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "probit",link.phi = "log", 10000)
ee7 <- teste_otimizacao(500,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "probit",link.phi = "log", 10000)

ee1
ee2
ee3
ee4
ee5
ee6
ee7

save.image("e_probit_log_Simulacao_R10000.RData")
set.seed(2)
eee1 <- teste_otimizacao(20,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "cloglog",link.phi = "log", 10000)
eee2 <- teste_otimizacao(30,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "cloglog",link.phi = "log", 10000)
eee3 <- teste_otimizacao(40,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "cloglog",link.phi = "log", 10000)
eee4 <- teste_otimizacao(50,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "cloglog",link.phi = "log", 10000)
eee5 <- teste_otimizacao(100,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "cloglog",link.phi = "log", 10000)
eee6 <- teste_otimizacao(200,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "cloglog",link.phi = "log", 10000)
eee7 <- teste_otimizacao(500,c(0.5, -1, -0.5),c(1.5,0.5),link.mu = "cloglog",link.phi = "log", 10000)

eee1
eee2
eee3
eee4
eee5
eee6
eee7

save.image("e_cloglog_log_Simulacao_R10000.RData")