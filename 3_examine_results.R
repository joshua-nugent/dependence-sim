library(tidyverse)

params <- cbind.data.frame(tac =  c(25, 25, 25, 25, 25, 25),
                           tacs = c( 5,  5, 10, 10,  5,  5),
                           cac =  c(50, 50, 25, 25, 25, 25),
                           cacs = c( 1,  1,  5,  5,  5,  5),
                           te =   c( 0,.25,  0,.25,  0,.25))

for(i in 1:nrow(params)){
  fname <- paste0("results/results",
                  "_tac", params[i,"tac"], "_tacs", params[i,"tacs"],
                  "_cac", params[i,"cac"], "_cacs", params[i,"cacs"],
                  "_te", params[i,"te"], 
                  ".rds")
  results <- readRDS(file = fname)
  print("===========================================")
  print(fname)
  
  print(paste0("Bias ltmle: ", round(mean(results$ind_agg_ind_weight_ltmle_pe - results$unit_mean_diff), digits = 3)))
  print(paste0("Bias gee  : ", round(mean(results$gee_pe - results$unit_mean_diff), digits = 3)))
  print(paste0("Bias glmm : ", round(mean(results$glmm_pe - results$unit_mean_diff), digits = 3)))

  print(paste0("SD of point estimates ltmle: ", round(sd(results$ind_agg_ind_weight_ltmle_pe), digits = 3)))
  print(paste0("SD of point estimates gee  : ", round(sd(results$gee_pe), digits = 3)))
  print(paste0("SD of point estimates glmm : ", round(sd(results$glmm_pe), digits = 3)))

  print(paste0("ltmle coverage: ", round(mean(results$ind_agg_ind_weight_ltmle_lo < results$unit_mean_diff &
               results$ind_agg_ind_weight_ltmle_hi > results$unit_mean_diff), digits = 3)))
  print(paste0("  gee coverage: ", round(mean(results$gee_lo < results$unit_mean_diff &
               results$gee_hi > results$unit_mean_diff), digits = 3)))
  print(paste0(" glmm coverage: ", round(mean(results$glmm_lo < results$unit_mean_diff &
               results$glmm_hi > results$unit_mean_diff), digits = 3)))
  if(params[i,"te"] == 0){
    print(paste0("ltmle TIE rate = ", 1 - round(mean(results$ind_agg_ind_weight_ltmle_lo < results$unit_mean_diff &
                                                 results$ind_agg_ind_weight_ltmle_hi > results$unit_mean_diff), digits = 3)))
    print(paste0("  gee TIE rate = ", 1 - round(mean(results$gee_lo < results$unit_mean_diff &
                                                 results$gee_hi > results$unit_mean_diff), digits = 3)))
    print(paste0(" glmm TIE rate = ", 1 - round(mean(results$glmm_lo < results$unit_mean_diff &
                                                 results$glmm_hi > results$unit_mean_diff), digits = 3)))
  } else {
    print(paste0("ltmle power = ", round(mean(results$ind_agg_ind_weight_ltmle_lo > 0), digits = 3)))
    print(paste0("  gee power = ", round(mean(results$gee_lo > 0), digits = 3)))
    print(paste0(" glmm power = ", round(mean(results$glmm_lo > 0), digits = 3)))
  }
  print(paste0("Mean CI width ltmle: ", round(mean(results$ind_agg_ind_weight_ltmle_hi - results$ind_agg_ind_weight_ltmle_lo), digits = 3)))
  print(paste0("Mean CI width gee  : ", round(mean(results$gee_hi - results$gee_lo), digits = 3)))
  print(paste0("Mean CI width glmm : ", round(mean(results$glmm_hi - results$glmm_lo), digits = 3)))
}



