library(tidyverse)
library(tmle)
library(ltmle)
library(dbarts)
library(doParallel)
library(foreach)
library(sandwich)
library(lme4)
library(gee)
library(lmerTest)

run_sim <- function(treatment_arm_cluster_size = 5, treatment_arm_clusters = 7,
                    control_arm_cluster_size = 1, control_arm_clusters = 35,
                    txt_eff = .25, ranef_sd = .35, covar1_coef = 1, seed = NA,
                    informative_cluster_size = F, SLL = c("SL.mean", "SL.glm", "SL.earth", "SL.glm.interaction")){
  
  dat <- gen_all(treatment_arm_clusters = treatment_arm_clusters,
                 treatment_arm_cluster_size = treatment_arm_cluster_size,
                 control_arm_clusters = control_arm_clusters,
                 control_arm_cluster_size = control_arm_cluster_size,
                 covar1_coef = covar1_coef,
                 ranef_sd = ranef_sd,
                 seed = seed,
                 txt_eff = txt_eff,
                 informative_cluster_size = informative_cluster_size)
  
  n <- nrow(dat)
  M <- J <- length(unique(dat$cluster_id))
  
  # Individual level effect, aggregating the IC accounting for cluster size - this is \Psi^I
  ind_agg_ind_weight <- suppressWarnings(suppressMessages(
    ltmle(data = dat %>% select(all_of(c("covar0", "covar1", "A", "Y"))),
          Ynodes = "Y", Anodes = "A",
          id = dat$cluster_id,
          variance.method = "ic",
          SL.library = SLL, abar = list(1,0))))
  summary_ind_agg_ind_weight <- summary(ind_agg_ind_weight)
  
  # GEE 
  gee <- gee(Y ~ A + covar0 + covar1, id = dat$id, corstr = "exchangeable",
             data = dat)
  gee_pe <- summary(gee)$coefficients["A","Estimate"]
  gee_se <- summary(gee)$coefficients["A","Robust S.E."]
  gee_lo <- gee_pe - qnorm(0.975) * gee_se
  gee_hi <- gee_pe + qnorm(0.975) * gee_se
  
    
  # GLMM
  # Note using Nugent paper to find optimal DF for LMM in this situation
  if(control_arm_cluster_size == 1){
    glmm <- lmer(Y ~ A + covar0 + covar1 +
                   (0 + A | cluster_id), data = dat, REML = T)
  } else {
    glmm <- lmer(Y ~ A + covar0 + covar1 +
                   (1 | cluster_id), data = dat, REML = T)
  }
  glmm_pe <- summary(glmm)$coefficients["A","Estimate"]
  glmm_se <- summary(glmm)$coefficients["A","Std. Error"]
  satterthwaite_df <- summary(glmm)$coefficients["A","df"]
  glmm_lo <- glmm_pe - qt(0.975, df = satterthwaite_df) * glmm_se
  glmm_hi <- glmm_pe + qt(0.975, df = satterthwaite_df) * glmm_se
  
  
  
  return(cbind.data.frame(treatment_arm_clusters = treatment_arm_clusters,
                          treatment_arm_cluster_size = treatment_arm_cluster_size,
                          control_arm_clusters = control_arm_clusters,
                          control_arm_cluster_size = control_arm_cluster_size,
                          covar1_coef = covar1_coef,
                          seed = seed,
                          ranef_sd = ranef_sd,
                          txt_eff = txt_eff,
                          informative_cluster_size = informative_cluster_size,
                          ind_agg_ind_weight_ltmle_pe = summary_ind_agg_ind_weight$effect.measures$ATE$estimate,
                          ind_agg_ind_weight_ltmle_lo = summary_ind_agg_ind_weight$effect.measures$ATE$CI[1],
                          ind_agg_ind_weight_ltmle_hi = summary_ind_agg_ind_weight$effect.measures$ATE$CI[2],
                          gee_pe = gee_pe,
                          gee_lo = gee_lo,
                          gee_hi = gee_hi,
                          glmm_pe = glmm_pe,
                          glmm_lo = glmm_lo,
                          glmm_hi = glmm_hi,
                          unit_mean_diff = dat$unit_mean_diff[1],
                          cluster_mean_diff = dat$cluster_mean_diff[1],
                          ind_weighted_cluster_mean_diff = dat$ind_weighted_cluster_mean_diff[1]
  ))
}








