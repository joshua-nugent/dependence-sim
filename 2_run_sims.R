library(tidyverse)
library(tmle)
library(ltmle)
library(dbarts)
library(doParallel)
library(foreach)
library(sandwich)
library(lme4)
library(gee)
source("0_utils.R")
source("1_sim_and_fit.R")
source("Stage2_cluster.R")

# IRGTT medium / big sample size
# Balanced medium / big sample size
# Unbalanced medium / big (but not singletons)
#                                   #irgtt  #unbal  #bal
params <- cbind.data.frame(tac =  c(25, 25, 25, 25, 25, 25),
                           tacs = c( 5,  5, 10, 10,  5,  5),
                           cac =  c(50, 50, 25, 25, 25, 25),
                           cacs = c( 1,  1,  5,  5,  5,  5),
                           te =   c( 0,.25,  0,.25,  0,.25))
set.seed(11)
for(i in 1:nrow(params)){
  registerDoParallel(cores = 8)
  nsims <- 1000
  SLL <- c("SL.mean", "SL.glm", "SL.earth", "SL.glm.interaction", "SL.glmnet")
  results <- data.frame(foreach(j = 1:nsims, .combine = rbind, .errorhandling = "remove", .verbose = F) %dopar% {
    run_sim(treatment_arm_clusters = params[i,"tac"], treatment_arm_cluster_size = params[i,"tacs"],
            control_arm_clusters = params[i,"cac"], control_arm_cluster_size = params[i,"cacs"],
            txt_eff = params[i,"te"], SLL = SLL,
            informative_cluster_size = F)
  })
  stopImplicitCluster()
  fname <- paste0("results/results",
                 "_tac", params[i,"tac"], "_tacs", params[i,"tacs"],
                 "_cac", params[i,"cac"], "_cacs", params[i,"cacs"],
                 "_te", params[i,"te"],
                 ".rds")
  saveRDS(results, file = fname)
}





