#waic comparison and save the ELPD resutls
library(loo)
library(rstan)
options(mc.cores = 6)# 
nsubj = 33
lbd_list = c(1, 15, 40, 70, 110)
lbd_b_list = c(.01, 1, 10, 20)
lbd_list_rv = c(1, 15, 40, 70, 110)
lbd_b_list_rv = c(.01, 1, 10, 20)

ref.ind <- 20


for (subjs in 1:nsubj) {
        loo_results = list()
        k = 1
        for (i in 1: length(lbd_list)){
            for (j in 1: length(lbd_b_list)){
                    for (m in 1:length(lbd_b_rv)){
                            for (n in 1:length(lbd_b_list_rv)){
                                    fits <- readRDS(file =paste0("sym_registration_Subject_",
                                                                 subjs,
                                                                 "_Ref_subj_",
                                                                 ref.ind, 
                                                                 "_lbdT_", lbd_list[i],
                                                                 "_lbdB_",lbd_b_list[j],
                                                                 "_lbdT_rv_", lbd_list_rv[m],
                                                                 "_lbdB_rv_",lbd_b_list_rv[n],
                                                                 ".rds"))
                                    fit.summary = rstan::summary(fits)
                                    rhat = max(fit.summary$summary[, "Rhat"])
                                    if (rhat>1.01){
                                            loo_results[[k]] = NA
                                    }else{
                                            loo_results[[k]] <-loo::loo(fits, K = 10, parameter_name = "log_lik", cores = 10)
                                    }
                                    k = k+1     
                            }
                    }
            }
            print(c(subjs,i,rhat))
        }
        
        elpds <- rep(NA, length(lbd_list)*length(lbd_b_list))
        k = 1
        for (i in 1:length(lbd_list)) {
                for (j in 1:length(lbd_b_list)) {
                        for (m in 1:length(lbd_b_rv)) {
                                for (n in 1:length(lbd_b_list_rv)) {
                                        if (!is.na(loo_results[[k]]))
                                                elpds[k] = loo_results[[k]]$estimates[1, 1]
                                        k = k + 1
                                }
                        }
                }
        }
        saveRDS(loo_results, file = paste0("../results/waic/sym_loo_subject_",subjs,".rds"))
        saveRDS(elpds, file = paste0("../results/waic/sym_elpds_subject_",subjs,".rds"))
}
