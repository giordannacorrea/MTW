# Load required packages
library(here)
library(moments)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

source(here("Simulation", "R", "functions_graph_summary.R"))
#-------------------------------------------------------------------------------
# Scenario B - logit
#beta<-c(0.5, -1, -0.5)
#nu <- c(1.5,0.5)

# The objects e1 through e7 are loaded into the environment when executing
load(here("Simulation", "simulated_data", "eNEW_logit_log_Simulacao_R10000.RData"))

vtheta = c(0.5, -1, -0.5, 1.5,0.5)

summary_n20_logit <- gerar_resumo(e1, vtheta)
summary_n30_logit <- gerar_resumo(e2, vtheta)
summary_n40_logit <- gerar_resumo(e3, vtheta)
summary_n50_logit <- gerar_resumo(e4, vtheta)
summary_n100_logit <- gerar_resumo(e5, vtheta)
summary_n200_logit <- gerar_resumo(e6, vtheta)
summary_n500_logit <- gerar_resumo(e7, vtheta)
# 
saveRDS(summary_n20_logit, here("Simulation", "summary_tables", "scenarioB", "scenarioB_logit_n20.Rds"))
saveRDS(summary_n30_logit, here("Simulation", "summary_tables", "scenarioB", "scenarioB_logit_n30.Rds"))
saveRDS(summary_n40_logit, here("Simulation", "summary_tables", "scenarioB", "scenarioB_logit_n40.Rds"))
saveRDS(summary_n50_logit, here("Simulation", "summary_tables", "scenarioB", "scenarioB_logit_n50.Rds"))
saveRDS(summary_n100_logit, here("Simulation", "summary_tables", "scenarioB", "scenarioB_logit_n100.Rds"))
saveRDS(summary_n200_logit, here("Simulation", "summary_tables","scenarioB", "scenarioB_logit_n200.Rds"))
saveRDS(summary_n500_logit, here("Simulation", "summary_tables","scenarioB", "scenarioB_logit_n500.Rds"))
#
completed_scenarioB_logit <- bind_rows(
  mutate(summary_n20_logit, n = 20),
  mutate(summary_n30_logit, n = 30),
  mutate(summary_n40_logit, n = 40),
  mutate(summary_n50_logit, n = 50),
  mutate(summary_n100_logit, n = 100),
  mutate(summary_n200_logit, n = 200),
  mutate(summary_n500_logit, n = 500)
)
completed_scenarioB_logit

# Save objects
saveRDS(completed_scenarioB_logit, here("Simulation", "summary_tables", "scenarioB", "summaryB_logit_completed.Rds"))
# Read the file resume_A_completed.Rds
summary_B_logit <- readRDS("Simulation/summary_tables/scenarioB/summaryB_logit_completed.Rds")
#------------------------------------------------------------------------------
# Figures
# Define the list of simulation objects

sim_list <- list(
  "20" = e1, "30" = e2, "40" = e3,
  "50" = e4, "100" = e5, "200" = e6, "500" = e7
)

# True values of the parameters
true_vals <- c(0.5, -1.0, -0.5, 1.5, 0.5)

# Labels for plot axis (use expression for Greek letters)
param_labels <- c(expression(beta[0]), expression(beta[1]),expression(beta[2]),
                  expression(nu[0]), expression(nu[1]))

# Generate boxplots
boxplots <- parameter_boxplots(sim_list, true_vals, param_labels)

# Example:
beta0_logit_scenarioB <- boxplots[[1]] # beta0
beta1_logit_scenarioB <- boxplots[[2]] # beta1
beta2_logit_scenarioB <- boxplots[[3]] # beta2
nu0_logit_scenarioB <- boxplots[[4]] # nu0
nu1_logit_scenarioB <- boxplots[[5]] # nu1

ggsave(here("Simulation", "figures", "scenarioB", "beta0_logit_scenarioB.png"), beta0_logit_scenarioB, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioB", "beta1_logit_scenarioB.png"), beta1_logit_scenarioB, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioB", "beta2_logit_scenarioB.png"), beta2_logit_scenarioB, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioB", "nu0_logit_scenarioB.png"), nu0_logit_scenarioB, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioB", "nu1_logit_scenarioB.png"), nu1_logit_scenarioB, width = 6, height = 4)
#-------------------------------------------------------------------------------
# Scenario B - probit
# beta<-c(0.5, -1, -0.5)
# nu <- c(1.5,0.5)

# The objects ee1 through ee7 are loaded into the environment when executing
load(here("Simulation", "simulated_data", "e_probit_log_Simulacao_R10000.RData"))

vtheta = c(0.5, -1, -0.5, 1.5,0.5)

summary_n20_probit <- gerar_resumo(ee1, vtheta)
summary_n30_probit <- gerar_resumo(ee2, vtheta)
summary_n40_probit <- gerar_resumo(ee3, vtheta)
summary_n50_probit <- gerar_resumo(ee4, vtheta)
summary_n100_probit <- gerar_resumo(ee5, vtheta)
summary_n200_probit <- gerar_resumo(ee6, vtheta)
summary_n500_probit <- gerar_resumo(ee7, vtheta)
# 
saveRDS(summary_n20_probit, here("Simulation", "summary_tables", "scenarioB", "scenarioB_probit_n20.Rds"))
saveRDS(summary_n30_probit, here("Simulation", "summary_tables", "scenarioB", "scenarioB_probit_n30.Rds"))
saveRDS(summary_n40_probit, here("Simulation", "summary_tables","scenarioB",  "scenarioB_probit_n40.Rds"))
saveRDS(summary_n50_probit, here("Simulation", "summary_tables", "scenarioB","scenarioB_probit_n50.Rds"))
saveRDS(summary_n100_probit, here("Simulation", "summary_tables", "scenarioB","scenarioB_probit_n100.Rds"))
saveRDS(summary_n200_probit, here("Simulation", "summary_tables", "scenarioB","scenarioB_probit_n200.Rds"))
saveRDS(summary_n500_probit, here("Simulation", "summary_tables", "scenarioB","scenarioB_probit_n500.Rds"))
#
completed_scenarioB_probit <- bind_rows(
  mutate(summary_n20_probit, n = 20),
  mutate(summary_n30_probit, n = 30),
  mutate(summary_n40_probit, n = 40),
  mutate(summary_n50_probit, n = 50),
  mutate(summary_n100_probit, n = 100),
  mutate(summary_n200_probit, n = 200),
  mutate(summary_n500_probit, n = 500)
)
completed_scenarioB_probit

# Save objects
saveRDS(completed_scenarioB_probit, here("Simulation", "summary_tables", "scenarioB", "summaryB_probit_completed.Rds"))
# Read the file 
summary_B_probit <- readRDS("Simulation/summary_tables/scenarioB/summaryB_probit_completed.Rds")
#------------------------------------------------------------------------------
# Figures
# Define the list of simulation objects
sim_list <- list(
  "20" = ee1, "30" = ee2, "40" = ee3,
  "50" = ee4, "100" = ee5, "200" = ee6, "500" = ee7
)

# True values of the parameters
true_vals <- c(0.5, -1.0, -0.5, 1.5, 0.5)

# Labels for plot axis (use expression for Greek letters)
param_labels <- c(expression(beta[0]), expression(beta[1]),expression(beta[2]),
                  expression(nu[0]), expression(nu[1]))

# Generate boxplots
boxplots <- parameter_boxplots(sim_list, true_vals, param_labels)

# Example:
beta0_probit_scenarioB <- boxplots[[1]] # beta0
beta1_probit_scenarioB <- boxplots[[2]] # beta1
beta2_probit_scenarioB <- boxplots[[3]] # beta2
nu0_probit_scenarioB <- boxplots[[4]] # nu0
nu1_probit_scenarioB <- boxplots[[5]] # nu1

ggsave(here("Simulation", "figures", "scenarioB", "beta0_probit_scenarioB.png"), beta0_probit_scenarioB, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioB", "beta1_probit_scenarioB.png"), beta1_probit_scenarioB, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioB", "beta2_probit_scenarioB.png"), beta2_probit_scenarioB, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioB", "nu0_probit_scenarioB.png"), nu0_probit_scenarioB, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioB", "nu1_probit_scenarioB.png"), nu1_probit_scenarioB, width = 6, height = 4)
#------------------------------------------------------------------------------
# Scenario B - cloglog
# beta<-c(0.5, -1, -0.5)
# nu <- c(1.5,0.5)

# The objects eee1 through eee7 are loaded into the environment when executing

load(here("Simulation", "simulated_data", "e_cloglog_log_Simulacao_R10000.RData"))

vtheta = c(0.5, -1, -0.5, 1.5,0.5)

summary_n20_cloglog <- gerar_resumo(eee1, vtheta)
summary_n30_cloglog <- gerar_resumo(eee2, vtheta)
summary_n40_cloglog <- gerar_resumo(eee3, vtheta)
summary_n50_cloglog <- gerar_resumo(eee4, vtheta)
summary_n100_cloglog <- gerar_resumo(eee5, vtheta)
summary_n200_cloglog <- gerar_resumo(eee6, vtheta)
summary_n500_cloglog <- gerar_resumo(eee7, vtheta)
# 
saveRDS(summary_n20_cloglog, here("Simulation", "summary_tables", "scenarioB","scenarioB_cloglog_n20.Rds"))
saveRDS(summary_n30_cloglog, here("Simulation", "summary_tables", "scenarioB","scenarioB_cloglog_n30.Rds"))
saveRDS(summary_n40_cloglog, here("Simulation", "summary_tables", "scenarioB","scenarioB_cloglog_n40.Rds"))
saveRDS(summary_n50_cloglog, here("Simulation", "summary_tables", "scenarioB","scenarioB_cloglog_n50.Rds"))
saveRDS(summary_n100_cloglog, here("Simulation", "summary_tables", "scenarioB","scenarioB_cloglog_n100.Rds"))
saveRDS(summary_n200_cloglog, here("Simulation", "summary_tables", "scenarioB","scenarioB_cloglog_n200.Rds"))
saveRDS(summary_n500_cloglog, here("Simulation", "summary_tables", "scenarioB","scenarioB_cloglog_n500.Rds"))
#
completed_scenarioB_cloglog <- bind_rows(
  mutate(summary_n20_cloglog, n = 20),
  mutate(summary_n30_cloglog, n = 30),
  mutate(summary_n40_cloglog, n = 40),
  mutate(summary_n50_cloglog, n = 50),
  mutate(summary_n100_cloglog, n = 100),
  mutate(summary_n200_cloglog, n = 200),
  mutate(summary_n500_cloglog, n = 500)
)
completed_scenarioB_cloglog

# Save objects
saveRDS(completed_scenarioB_cloglog, here("Simulation", "summary_tables", "scenarioB", "summaryB_cloglog_completed.Rds"))
# Read the file 
summary_B_cloglog <- readRDS("Simulation/summary_tables/scenarioB/summaryB_cloglog_completed.Rds")
#-------------------------------------------------------------------------------
# Figures

# Define the list of simulation objects
sim_list <- list(
  "20" = eee1, "30" = eee2, "40" = eee3,
  "50" = eee4, "100" = eee5, "200" = eee6, "500" = eee7
)

# True values of the parameters
true_vals <- c(0.5, -1.0, -0.5, 1.5, 0.5)

# Labels for plot axis (use expression for Greek letters)
param_labels <- c(expression(beta[0]), expression(beta[1]),expression(beta[2]),
                  expression(nu[0]), expression(nu[1]))

# Generate boxplots
boxplots <- parameter_boxplots(sim_list, true_vals, param_labels)

# Example:
beta0_cloglog_scenarioB <- boxplots[[1]] # beta0
beta1_cloglog_scenarioB <- boxplots[[2]] # beta1
beta2_cloglog_scenarioB <- boxplots[[3]] # beta2
nu0_cloglog_scenarioB <- boxplots[[4]] # nu0
nu1_cloglog_scenarioB <- boxplots[[5]] # nu1

ggsave(here("Simulation", "figures", "scenarioB", "beta0_cloglog_scenarioB.png"), beta0_cloglog_scenarioB, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioB", "beta1_cloglog_scenarioB.png"), beta1_cloglog_scenarioB, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioB", "beta2_cloglog_scenarioB.png"), beta2_cloglog_scenarioB, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioB", "nu0_cloglog_scenarioB.png"), nu0_cloglog_scenarioB, width = 6, height = 4)
ggsave(here("Simulation", "figures", "scenarioB", "nu1_cloglog_scenarioB.png"), nu1_cloglog_scenarioB, width = 6, height = 4)
