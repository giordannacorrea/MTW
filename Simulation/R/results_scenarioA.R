# Load required packages
library(here)
library(moments)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

source(here("Simulation", "R", "functions_graph_summary.R"))
#-------------------------------------------------------------------------------
# Scenario A - logit
# beta = c(-1,0.5)
# nu = c(1,0.5)
# The objects a1 through a7 are loaded into the environment when executing
load(here("Simulation", "simulated_data", "aTeste_logit_log_Simulacao_R10000.RData"))

vtheta = c(-1,0.5, 1,0.5)

summary_n20_logit <- gerar_resumo(a1, vtheta)
summary_n30_logit <- gerar_resumo(a2, vtheta)
summary_n40_logit <- gerar_resumo(a3, vtheta)
summary_n50_logit <- gerar_resumo(a4, vtheta)
summary_n100_logit <- gerar_resumo(a5, vtheta)
summary_n200_logit <- gerar_resumo(a6, vtheta)
summary_n500_logit <- gerar_resumo(a7, vtheta)
# 
saveRDS(summary_n20_logit, here("Simulation", "summary_tables", "scenarioA", "scenarioA_logit_n20.Rds"))
saveRDS(summary_n30_logit, here("Simulation", "summary_tables", "scenarioA", "scenarioA_logit_n30.Rds"))
saveRDS(summary_n40_logit, here("Simulation", "summary_tables", "scenarioA", "scenarioA_logit_n40.Rds"))
saveRDS(summary_n50_logit, here("Simulation", "summary_tables", "scenarioA", "scenarioA_logit_n50.Rds"))
saveRDS(summary_n100_logit, here("Simulation", "summary_tables", "scenarioA", "scenarioA_logit_n100.Rds"))
saveRDS(summary_n200_logit, here("Simulation", "summary_tables","scenarioA", "scenarioA_logit_n200.Rds"))
saveRDS(summary_n500_logit, here("Simulation", "summary_tables","scenarioA", "scenarioA_logit_n500.Rds"))
#
completed_scenarioA_logit <- bind_rows(
  mutate(summary_n20_logit, n = 20),
  mutate(summary_n30_logit, n = 30),
  mutate(summary_n40_logit, n = 40),
  mutate(summary_n50_logit, n = 50),
  mutate(summary_n100_logit, n = 100),
  mutate(summary_n200_logit, n = 200),
  mutate(summary_n500_logit, n = 500)
)
completed_scenarioA_logit

# Save objects
saveRDS(completed_scenarioA_logit, here("Simulation", "summary_tables", "scenarioA", "summaryA_logit_completed.Rds"))
# Read the file resume_A_completed.Rds
summary_A_logit <- readRDS("Simulation/summary_tables/scenarioA/summaryA_logit_completed.Rds")
#-------------------------------------------------------------------------------
# Figures
# Define the list of simulation objects

sim_list <- list(
  "20" = a1, "30" = a2, "40" = a3,
  "50" = a4, "100" = a5, "200" = a6, "500" = a7
)

# True values of the parameters
true_vals = c(-1,0.5, 1,0.5) 

# Labels for plot axis (use expression for Greek letters)
param_labels <- c(expression(beta[0]), expression(beta[1]),
                  expression(nu[0]), expression(nu[1]))

# Generate boxplots
boxplots <- parameter_boxplots(sim_list, true_vals, param_labels)

# Example:
beta0_logit_scenarioA <- boxplots[[1]] # beta0
beta1_logit_scenarioA <- boxplots[[2]] # beta1
nu0_logit_scenarioA <- boxplots[[3]] # nu0
nu1_logit_scenarioA <- boxplots[[4]] # nu1

ggsave(here("Simulation", "figures", "scenarioA", "beta0_logit_scenarioA.png"), beta0_logit_scenarioA, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioA", "beta1_logit_scenarioA.png"), beta1_logit_scenarioA, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioA", "nu0_logit_scenarioA.png"), nu0_logit_scenarioA, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioA", "nu1_logit_scenarioA.png"), nu1_logit_scenarioA, width = 6, height = 4)
#-------------------------------------------------------------------------------
# Scenario A - probit
# beta = c(-1,0.5)
# nu = c(1,0.5)
# The objects aa1 through aa7 are loaded into the environment when executing
load(here("Simulation", "simulated_data", "aTeste_probit_log_Simulacao_R10000.RData"))

vtheta = c(-1,0.5, 1,0.5)

summary_n20_probit <- gerar_resumo(aa1, vtheta)
summary_n30_probit <- gerar_resumo(aa2, vtheta)
summary_n40_probit <- gerar_resumo(aa3, vtheta)
summary_n50_probit <- gerar_resumo(aa4, vtheta)
summary_n100_probit <- gerar_resumo(aa5, vtheta)
summary_n200_probit <- gerar_resumo(aa6, vtheta)
summary_n500_probit <- gerar_resumo(aa7, vtheta)
# 
saveRDS(summary_n20_probit, here("Simulation", "summary_tables", "scenarioA", "scenarioA_probit_n20.Rds"))
saveRDS(summary_n30_probit, here("Simulation", "summary_tables", "scenarioA", "scenarioA_probit_n30.Rds"))
saveRDS(summary_n40_probit, here("Simulation", "summary_tables","scenarioA",  "scenarioA_probit_n40.Rds"))
saveRDS(summary_n50_probit, here("Simulation", "summary_tables", "scenarioA","scenarioA_probit_n50.Rds"))
saveRDS(summary_n100_probit, here("Simulation", "summary_tables", "scenarioA","scenarioA_probit_n100.Rds"))
saveRDS(summary_n200_probit, here("Simulation", "summary_tables", "scenarioA","scenarioA_probit_n200.Rds"))
saveRDS(summary_n500_probit, here("Simulation", "summary_tables", "scenarioA","scenarioA_probit_n500.Rds"))
#
completed_scenarioA_probit <- bind_rows(
  mutate(summary_n20_probit, n = 20),
  mutate(summary_n30_probit, n = 30),
  mutate(summary_n40_probit, n = 40),
  mutate(summary_n50_probit, n = 50),
  mutate(summary_n100_probit, n = 100),
  mutate(summary_n200_probit, n = 200),
  mutate(summary_n500_probit, n = 500)
)
completed_scenarioA_probit

# Save objects
saveRDS(completed_scenarioA_probit, here("Simulation", "summary_tables", "scenarioA", "summaryA_probit_completed.Rds"))
# Read the file 
summary_A_probit <- readRDS("Simulation/summary_tables/scenarioA/summaryA_probit_completed.Rds")
#-------------------------------------------------------------------------------
# Figures
# Define the list of simulation objects

sim_list <- list(
  "20" = aa1, "30" = aa2, "40" = aa3,
  "50" = aa4, "100" = aa5, "200" = aa6, "500" = aa7
)

# True values of the parameters
true_vals = c(-1,0.5, 1,0.5) 

# Labels for plot axis (use expression for Greek letters)
param_labels <- c(expression(beta[0]), expression(beta[1]),
                  expression(nu[0]), expression(nu[1]))

# Generate boxplots
boxplots <- parameter_boxplots(sim_list, true_vals, param_labels)

# Example:
beta0_probit_scenarioA <- boxplots[[1]] # beta0
beta1_probit_scenarioA <- boxplots[[2]] # beta1
nu0_probit_scenarioA <- boxplots[[3]] # nu0
nu1_probit_scenarioA <- boxplots[[4]] # nu1

ggsave(here("Simulation", "figures", "scenarioA", "beta0_probit_scenarioA.png"), beta0_probit_scenarioA, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioA", "beta1_probit_scenarioA.png"), beta1_probit_scenarioA, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioA", "nu0_probit_scenarioA.png"), nu0_probit_scenarioA, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioA", "nu1_probit_scenarioA.png"), nu1_probit_scenarioA, width = 6, height = 4)
#-------------------------------------------------------------------------------
# Scenario A - cloglog
# beta = c(-1,0.5)
# nu = c(1,0.5)
# The objects aaa1 through aaa7 are loaded into the environment when executing

load(here("Simulation", "simulated_data", "aTeste_cloglog_log_Simulacao_R10000.RData"))

vtheta = c(-1,0.5, 1,0.5)

summary_n20_cloglog <- gerar_resumo(aaa1, vtheta)
summary_n30_cloglog <- gerar_resumo(aaa2, vtheta)
summary_n40_cloglog <- gerar_resumo(aaa3, vtheta)
summary_n50_cloglog <- gerar_resumo(aaa4, vtheta)
summary_n100_cloglog <- gerar_resumo(aaa5, vtheta)
summary_n200_cloglog <- gerar_resumo(aaa6, vtheta)
summary_n500_cloglog <- gerar_resumo(aaa7, vtheta)
# 
saveRDS(summary_n20_cloglog, here("Simulation", "summary_tables", "scenarioA","scenarioA_cloglog_n20.Rds"))
saveRDS(summary_n30_cloglog, here("Simulation", "summary_tables", "scenarioA","scenarioA_cloglog_n30.Rds"))
saveRDS(summary_n40_cloglog, here("Simulation", "summary_tables", "scenarioA","scenarioA_cloglog_n40.Rds"))
saveRDS(summary_n50_cloglog, here("Simulation", "summary_tables", "scenarioA","scenarioA_cloglog_n50.Rds"))
saveRDS(summary_n100_cloglog, here("Simulation", "summary_tables", "scenarioA","scenarioA_cloglog_n100.Rds"))
saveRDS(summary_n200_cloglog, here("Simulation", "summary_tables", "scenarioA","scenarioA_cloglog_n200.Rds"))
saveRDS(summary_n500_cloglog, here("Simulation", "summary_tables", "scenarioA","scenarioA_cloglog_n500.Rds"))
#
completed_scenarioA_cloglog <- bind_rows(
  mutate(summary_n20_cloglog, n = 20),
  mutate(summary_n30_cloglog, n = 30),
  mutate(summary_n40_cloglog, n = 40),
  mutate(summary_n50_cloglog, n = 50),
  mutate(summary_n100_cloglog, n = 100),
  mutate(summary_n200_cloglog, n = 200),
  mutate(summary_n500_cloglog, n = 500)
)
completed_scenarioA_cloglog

# Save objects
saveRDS(completed_scenarioA_cloglog, here("Simulation", "summary_tables", "scenarioA", "summaryA_cloglog_completed.Rds"))
# Read the file 
summary_A_cloglog <- readRDS("Simulation/summary_tables/scenarioA/summaryA_cloglog_completed.Rds")
#-------------------------------------------------------------------------------
# Figures
# Define the list of simulation objects
sim_list <- list(
  "20" = aaa1, "30" = aaa2, "40" = aaa3,
  "50" = aaa4, "100" = aaa5, "200" = aaa6, "500" = aaa7
)

# True values of the parameters
true_vals = c(-1,0.5, 1,0.5) 

# Labels for plot axis (use expression for Greek letters)
param_labels <- c(expression(beta[0]), expression(beta[1]),
                  expression(nu[0]), expression(nu[1]))

# Generate boxplots
boxplots <- parameter_boxplots(sim_list, true_vals, param_labels)

# Example:
beta0_cloglog_scenarioA <- boxplots[[1]] # beta0
beta1_cloglog_scenarioA <- boxplots[[2]] # beta1
nu0_cloglog_scenarioA <- boxplots[[3]] # nu0
nu1_cloglog_scenarioA <- boxplots[[4]] # nu1

ggsave(here("Simulation", "figures", "scenarioA", "beta0_cloglog_scenarioA.png"), beta0_cloglog_scenarioA, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioA", "beta1_cloglog_scenarioA.png"), beta1_cloglog_scenarioA, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioA", "nu0_cloglog_scenarioA.png"), nu0_cloglog_scenarioA, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioA", "nu1_cloglog_scenarioA.png"), nu1_cloglog_scenarioA, width = 6, height = 4)
