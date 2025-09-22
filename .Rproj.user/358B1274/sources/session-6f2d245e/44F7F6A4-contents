# Load required packages
library(VGAM)
library(tidyverse)
library(hnp)
library(betareg)
library(Stat2Data)
library(usmap)
library(PNWColors)
library(car)
library(ggthemes)
library(MASS)
library(reshape2)
library(reshape)
library(maxLik)
library(here)
library(usmapdata)

# Load custom functions (order matters for dependencies)
source(here("USA_election", "R", "funcoesWeibullT.R"))
source(here("USA_election", "R", "WeibullT_Regression.R"))
source(here("USA_election", "R", "comp.models.R"))
source(here("USA_election", "R", "graficoWT.R"))

# Load example dataset
data(Election16)
data(statepop) # Used for usmap visualizations

# The dataset Election16, used in the analysis, has been updated.
Election16$Dem.Rep <- c(
  -17.3,-17.2,-1.4,-6.9, 15.5,-1.4,10.5, 5.8, 0.7,-3.9,
  20.8,-30, 10.8,-6, -3.3,-13.1,-2.8,-1.6,-3.7, 15,
  18.6, 4, 0.5,-5.2,-7.9, -14.4,-8.7,-4.1,-8.8,10.2,
  11.6, 18.4, 0.2,-20.8, 0.2,-13.4,5.8, 3, 19.4,-12.5,
  -15.8,-10.3,-5.5,-29.9, 21.7, -1.6,9.1,-6.1,-0.7,-31.8)

# Check updates
head(Election16$Dem.Rep)
head(Election16)


# Attach data for simplified variable access 
attach(Election16)
#------------------------------------------------------------------------------
# Plot construction - data set containing graph variables ("fips", "abbr", "full")

statepop2 <- statepop |> 
  filter(abbr != "DC")# remove the "District of Columbia"
#------------------------------------------------------------------------------
# Merging the datasets by state code

eleicoes <- dplyr::right_join(statepop2, Election16, by = c("abbr" = "Abr")) |>
  mutate(DDR = Dem.Rep / 100) |>                            # scaling the response variable
  mutate(TrumpWin = as.factor(TrumpWin)) |>                 # converting binary variable to a factor
  mutate(Status = ifelse(TrumpWin == "1", "Won", "Lost"),
         .after = "TrumpWin") |>                            # creates a status for Plot 2
  mutate(statusDDR = ifelse(DDR > "0", "Won", "Lost")) # classification based on DDR
names(eleicoes)

#------------------------------------------------------------------------------#
# Plot 1 - United States Election (2016), states where Trump Won/lost.

trump_win <- eleicoes |>
  mutate(state = str_trim(State, side = "right")) |>
  dplyr::select(fips, abbr, full, Status)

TRUMP16 = usmap::plot_usmap(
  data = trump_win,
  values = "Status",
  labels = T,
  label_color = "white",
  exclude = c("DC")
) +
  scale_fill_manual(values = c("Won" = "darkred",
                               "Lost" = "darkblue"),
                    name = "Trump") +
  theme(
    legend.position = "right",
    panel.background = element_rect(colour = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(title = "") +
  xlab("") +
  ylab("")
TRUMP16

#ggsave("USA_election/figures/Trump2016.png", width = 10, height = 7, dpi = 300)
#------------------------------------------------------------------------------#
# Plot 2 - Response variable

grafDDR <- plot_usmap(data = eleicoes,
                      values = "DDR",
                      labels = TRUE,
                      exclude = c("DC")
                      ) +
  scale_fill_gradientn(colours = pnw_palette("Bay", n = 50, type = "continuous"),
                       name = "DDR",
                       limits = c(-0.35, 0.35)) +
  theme(legend.position = "right",
        panel.background = element_rect(colour = "black")) + ylab("") +
  xlab("")

grafDDR
#ggsave("USA_election/figures/USmap_2016.png", width = 10, height = 7, dpi = 300)
#------------------------------------------------------------------------------#
# Summary of the response variable - Democratic and Republican Difference (DDR)
y = DDR <- eleicoes$DDR
summary(y)

m0 <- 0.0201
phi0 <-8.5809

hist <- data.frame(DDR = Election16$Dem.Rep/100)

ggplot(hist,aes(x=DDR))+xlim(-0.5,0.4)+
  theme_classic()+
  theme(text = element_text(size = 15),
        panel.border = element_rect(linetype = "solid", fill = NA))+
  geom_histogram(aes(y =..density..),bins=7,
                 colour = 1, fill = "white")+
  geom_function(fun=dWeibullT,args=list(phi=phi0, m=m0))+
  geom_function(fun=dWeibullT,args=list(phi=phi0, m=m0),lwd=1.5)+
  ylab("Density")
#ggsave("USA_election/figures/histDDR.pdf", width = 8, height = 6)
#-------------------------------------------------------------------------------
# Covariates
RPC<-Income/10000 # Per capita income (scaled by 10,000)

PGE<-HS # Percentage of high school graduates

PGR<-BA # Percentage of college graduates (bachelor's degree)

PDA<-Adv # Percentage with advanced degrees
#------------------------------------------------------------------------------#
#Response variable
y<-eleicoes$DDR
n=length(y)

# Data frame for model fitting
dados_ajuste <- data.frame(
  Y3 = y,
  X1 = PGR / 100, # Proportion of college graduates
  X2 = PDA / 100, # Proportion with advanced degrees
  X3 = PGE / 100, # Proportion of high school graduates
  X4 = RPC        # Per capita income 
)
#-------------------------------------------------------------------------------
# Association between the response variable and the explanatory variables RPC, PGE, PDA, and PGR.

pairs(dados_ajuste,pch=16)

#-------------------------------------------------------------------------------
# Descriptive analysis

fit = WeibullT_reg.fit(
  formula = Y3 ~ 1 |
    1,
  link.mu = "identity",
  link.phi = "identity",
  output = 1,
  dados_ajuste
)

# Histogram of the response variable with overlaid density curves

ggplot(data = dados_ajuste, aes(x =  Y3)) +
  geom_histogram(
    aes(y = ..density..),
    binwidth = 0.1,
    fill = "darkgray",
    col = "white"
  ) +
  stat_function(
    fun = dWeibullT,
    args = list(fit$Par_est[1], fit$Par_est[2]),
    color = "darkred"
  ) +
  xlab("y") +
  ylab("Densidade") +
  scale_x_continuous(breaks = seq(from = -0.4, to = 0.35, by = 0.1),
                     limits = c(-0.4, 0.38)) +
  ggtitle("") +
  theme_bw()
#-------------------------------------------------------------------------------
# Fitted Models
#-------------------------------------------------------------------------------
# Model 1: logit (first model)
# Y3 ~ PGR + PDA + PGE + log(RPC) | PDA + PGE

fit1 <- WeibullT_reg.fit(formula = Y3 ~ X1+X2+X3+log(X4)|X2+X3,
                         output = 0,
                         link.mu = "logit",
                         link.phi = "log",
                         dados_ajuste)
fit1$model_nlminb
round(fit1$r2_nlminb,4)
round(fit1$aic2,4)
round(fit1$bic2,4)
# -------------------------------------------------------------------------
# Model 2: probit 

fit2 <- WeibullT_reg.fit(formula = Y3 ~ X1+X2+X3+log(X4)|X2+X3,
                         output = 0,
                         link.mu = "probit",
                         link.phi = "log",
                         dados_ajuste)
fit2$model_nlminb
round(fit2$r2_nlminb,4)
round(fit2$aic2,4)
round(fit2$bic2,4)
# ------------------------------------------------------------------------------
# Model 3: cloglog

fit3 <- WeibullT_reg.fit(formula = Y3 ~ X1+X2+X3+log(X4)|X2+X3,
                         output = 0,
                         link.mu = "cloglog",
                         link.phi = "log",
                         dados_ajuste)
fit3$model_nlminb
round(fit3$r2_nlminb,4)
round(fit3$aic2,4)
round(fit3$bic2,4)
#------------------------------------------------------------------------------#
# The best model - Y3 ~ PDA+PGE|PGE 
#------------------------------------------------------------------------------#
# Logit

fit_logit = WeibullT_reg.fit(formula = Y3 ~X2+X3|X3,
                             output = 0,
                             "logit",
                             "log",
                             dados_ajuste)

fit_logit$coef
beta <- fit_logit$coef[1:ncol(fit_logit$X)];beta
nu<- fit_logit$coef[(ncol(fit_logit$X)+1):length(fit_logit$coef)];nu
AIC=-2*fit_logit$loglik2+2*(2+length(beta));AIC
BIC=-2*fit_logit$loglik2+log(n)*(2+length(beta));BIC
round(fit_logit$aic2,5)
round(fit_logit$bic2,5)
fit_logit$loglik2
fit_logit$model_nlminb
fit_logit$Par_est[1]
round(fit_logit$r2_nlminb,5)
# Plot
set.seed(2)
diag.WT.fit(model = fit_logit, sim = 100, conf = 0.95, data = dados_ajuste, r = 1)
#------------------------------------------------------------------------------#
# Probit

fit_probit = WeibullT_reg.fit(formula = Y3 ~ X2+X3|X3,
                              output = 0,
                              "probit",
                              "log",
                              dados_ajuste)
fit_probit$model_nlminb
fit_probit$coef
beta <- fit_probit$coef[1:ncol(fit_probit$X)] 
nu<- fit_probit$coef[(ncol(fit_probit$X)+1):length(fit_probit$coef)]
AIC=-2*fit_probit$loglik2+2*(2+length(beta));AIC
BIC=-2*fit_probit$loglik2+log(n)*(2+length(beta));BIC
round(fit_probit$aic2,5)
round(fit_probit$bic2,5)
round(fit_probit$r2_nlminb,5)
# Plot
set.seed(2)
diag.WT.fit(model = fit_probit, sim = 100, conf = 0.95, data = dados_ajuste, r = 1)
#-------------------------------------------------------------------------------
# Cloglog

fit_cloglog = WeibullT_reg.fit(formula = Y3 ~ X2+X3|X3,
                               output = 0,
                               "cloglog",
                               "log",
                               dados_ajuste)
fit_cloglog$model_nlminb
fit_cloglog$coef
beta <- fit_cloglog$coef[1:ncol(fit_cloglog$X)] 
nu<- fit_cloglog$coef[(ncol(fit_cloglog$X)+1):length(fit_cloglog$coef)]
AIC=-2*fit_cloglog$loglik2+2*(2+length(beta));AIC
BIC=-2*fit_cloglog$loglik2+log(n)*(2+length(beta));BIC
round(fit_cloglog$aic2,5)
round(fit_cloglog$bic2,5)
round(fit_cloglog$r2_nlminb,5)
# Plot
set.seed(2)
diag.WT.fit(model = fit_cloglog, sim = 100, conf = 0.95, data = dados_ajuste, r = 1)

#-------------------------------------------------------------------------------
# Vuong test for comparing the WTM, Kumaraswamy, Unit-Gompertz, 
# Unit-Gamma, and Beta models for the electoral data
#-------------------------------------------------------------------------------
# The best model selected was:

formula = Y3~X2+X3|X3
link.mu = "logit"
link.phi = "log"
data = dados_ajuste
set.seed(2)
# Comparison of models using the logit link function
fit_logit = WeibullT_reg.fit(formula = Y3 ~X2+X3|X3,
                             output = 0,
                             "logit",
                             "log",
                             dados_ajuste)
fit_logit$model_nlminb
set.seed(2)
logit_graf <- graficoWT(model = fit_logit, "logit", "log", dados_ajuste)
fit_logit$r2_nlminb

# Comparison

logito<-comp.models(formula = Y3~X2+X3|X3,
                    "logit",
                    "log",
                    dados_ajuste)

# Summary of the model adjustments for the electoral data
logito$resumo
# Vuong test
logito$comparacao

# Estimate and Stardard Error
# Beta
a <- logito$Beta
round(a$estimate,4)[,1:2]
# Kumaraswamy
b <- logito$kuma
round(b$estimate,4)[,1:2]
# Ugompertz
c <- logito$UGa
round(c$estimate,4)[,1:2]
# Ugama
d <- logito$Ugz
round(d$estimate,4)[,1:2]
