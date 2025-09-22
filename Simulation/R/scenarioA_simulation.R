source("Simulation/R/otimizacaoWT.R")
source("Simulation/R/random_weibullT.R")

set.seed(2)
a1 <- teste_otimizacao(20,c(-1, 0.5),c(1,0.5),link.mu = "logit",link.phi = "log", 10000)
a2 <- teste_otimizacao(30,c(-1, 0.5),c(1,0.5),link.mu = "logit",link.phi = "log", 10000)
a3 <- teste_otimizacao(40,c(-1, 0.5),c(1,0.5),link.mu = "logit",link.phi = "log", 10000)
a4 <- teste_otimizacao(50,c(-1, 0.5),c(1,0.5),link.mu = "logit",link.phi = "log", 10000)
a5 <- teste_otimizacao(100,c(-1, 0.5),c(1,0.5),link.mu = "logit",link.phi = "log", 10000)
a6 <- teste_otimizacao(200,c(-1, 0.5),c(1,0.5),link.mu = "logit",link.phi = "log", 10000)
a7 <- teste_otimizacao(500,c(-1, 0.5),c(1,0.5),link.mu = "logit",link.phi = "log", 10000)

a1
a2
a3
a4
a5
a6
a7

save.image("aTeste_logit_log_Simulacao_R10000.RData")

set.seed(2)
aa1 <- teste_otimizacao(20,c(-1, 0.5),c(1,0.5),link.mu = "probit",link.phi = "log", 10000)
aa2 <- teste_otimizacao(30,c(-1, 0.5),c(1,0.5),link.mu = "probit",link.phi = "log", 10000)
aa3 <- teste_otimizacao(40,c(-1, 0.5),c(1,0.5),link.mu = "probit",link.phi = "log", 10000)
aa4 <- teste_otimizacao(50,c(-1, 0.5),c(1,0.5),link.mu = "probit",link.phi = "log", 10000)
aa5 <- teste_otimizacao(100,c(-1, 0.5),c(1,0.5),link.mu = "probit",link.phi = "log", 10000)
aa6 <- teste_otimizacao(200,c(-1, 0.5),c(1,0.5),link.mu = "probit",link.phi = "log", 10000)
aa7 <- teste_otimizacao(500,c(-1, 0.5),c(1,0.5),link.mu = "probit",link.phi = "log", 10000)

aa1
aa2
aa3
aa4
aa5
aa6
aa7

save.image("aTeste_probit_log_Simulacao_R10000.RData")
set.seed(2)
aaa1 <- teste_otimizacao(20,c(-1, 0.5),c(1,0.5),link.mu = "cloglog",link.phi = "log", 10000)
aaa2 <- teste_otimizacao(30,c(-1, 0.5),c(1,0.5),link.mu = "cloglog",link.phi = "log", 10000)
aaa3 <- teste_otimizacao(40,c(-1, 0.5),c(1,0.5),link.mu = "cloglog",link.phi = "log", 10000)
aaa4 <- teste_otimizacao(50,c(-1, 0.5),c(1,0.5),link.mu = "cloglog",link.phi = "log", 10000)
aaa5 <- teste_otimizacao(100,c(-1, 0.5),c(1,0.5),link.mu = "cloglog",link.phi = "log", 10000)
aaa6 <- teste_otimizacao(200,c(-1, 0.5),c(1,0.5),link.mu = "cloglog",link.phi = "log", 10000)
aaa7 <- teste_otimizacao(500,c(-1, 0.5),c(1,0.5),link.mu = "cloglog",link.phi = "log", 10000)

aaa1
aaa2
aaa3
aaa4
aaa5
aaa6
aaa7

save.image("aTeste_cloglog_log_Simulacao_R10000.RData")