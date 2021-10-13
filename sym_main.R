# fit stan model ----------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(rstan)
library(geoR)
library(R.matlab)
# load data and priors  ---------------------------------------------------
painDat = readRDS(file = paste0("data/Paindata_Trial_30by30_round.rds"))
ini.prior.procrustes <- R.matlab::readMat(paste0("data/symmetric_priors_.mat"))$trsall
ini.imgs <- R.matlab::readMat(paste0("data/symmetric_IMG_priors_.mat"))

ref.mat <- readMat("data/reference_3pks.mat")$ref 
ref.dim <- dim(ref.mat)
ref.nrow <- ref.dim[1]
ref.ncol <- ref.dim[2]
mean.coord <- 15
ref.coord <- 
  expand.grid(x = 10:24,
              y = 13:25)-mean.coord

# prepare for HPC computing -----------------------------------------------

# temp = commandArgs(trailingOnly = TRUE)
temp = c(1, 4, 1, 4, 20, 7)

rv_lbdid = as.integer(temp[1])
rv_lbd_bid = as.integer(temp[2])
lbd_bid = as.integer(temp[3])
lbdid = as.integer(temp[4])
refid = as.integer(temp[5])
subjs = as.integer(temp[6])

ref.ind = subjlist[refid]
subjlist = 1:33

lbd_list = c(1, 15, 40, 70, 110)
lbd_b_list = c(.01, 1, 10, 20)
lbd_list_rv = c(1, 15, 40, 70, 110)
lbd_b_list_rv = c(.01, 1, 10, 20)

subjs = subjlist[as.integer(temp[1])]
top = toplist[1]
lbds = lbd_list[lbdid]
lbd_b = lbd_b_list[lbd_bid]
lbds_rv = lbd_list_rv[rv_lbdid]
lbd_b_rv= lbd_b_list_rv[rv_lbd_bid]

X <- painDat[,subjs]
X_rv = painDat[, refid]
nrow = length(unique(painDat$x))
ncol = length(unique(painDat$y))
coord <- painDat[,c("x","y")]
ref <- as.vector(ref.mat)

# estimate the spatial decay parameter, rho
fit.mle <-  geoR::likfit(data=as.vector(X),coords=coord/nrow,
                         fix.nugget=T,
                         nugget = 0,
                         cov.model="exponential",
                         ini = c(var(as.vector(X)), 1),
                         message = F)
rhos <- 1/fit.mle$phi
# prepare the parameters for prior distribution for all subjects
ini.prior <- data.frame(log_sigma_x = log(ini.prior.procrustes[,,1]$scaling[,1]),
                        log_sigma_y = log(ini.prior.procrustes[,,1]$scaling[,2]),
                        rota = ini.prior.procrustes[,,1]$rotation[,2],
                        theta_x = ini.prior.procrustes[,,1]$b[,1],
                        theta_y = ini.prior.procrustes[,,1]$b[,2],
                        phi = rep(5e-2, 33),
                        b = ini.prior.procrustes[,,1]$scale,
                        log_sigma_x_rv = -log(ini.prior.procrustes[,,1]$scaling[,1]),
                        log_sigma_y_rv = -log(ini.prior.procrustes[,,1]$scaling[,2]),
                        rota_rv = -ini.prior.procrustes[,,1]$rotation[,2],
                        theta_x_rv = -ini.prior.procrustes[,,1]$b[,1],
                        theta_y_rv = -ini.prior.procrustes[,,1]$b[,2],
                        b_rv = ini.prior.procrustes[,,1]$scale.inv)
# data input for STAN fitting
data2d.input = list(N = length(X), 
                    NROW = nrow,
                    N_ref = length(ref),
                    N_rv = nrow(ini.imgs[[1]][[subj]][[1]]),
                    ref = as.vector(ref),
                    X = X, 
                    ref_rv = as.vector(ini.imgs[[1]][[subj]][[1]][,3]),
                    X_rv = painDat[, subj], 
                    coord = t(coord),
                    coord_ref = t(as.matrix(ref.coord)),
                    coord_rv = t(as.matrix(ini.imgs[[1]][[subj]][[1]][,1:2])),
                    lim_var = 5,
                    lim_theta = c(10,10),
                    lim_log_sigma = c(2,2),
                    limRota = 2,
                    rho = rhos,
                    b0 = ini.prior$b[subjs],
                    tx0 = ini.prior$theta_x[subjs],
                    ty0 = ini.prior$theta_y[subjs],
                    lsx0 = ini.prior$log_sigma_x[subjs],
                    lsy0 = ini.prior$log_sigma_y[subjs],
                    r0 = ini.prior$rota[subjs],
                    phi0 = ini.prior$phi[subjs],
                    b0_rv= ini.prior$b_rv[subjs],
                    tx0_rv = ini.prior$theta_x_rv[subjs],
                    ty0_rv = ini.prior$theta_y_rv[subjs],
                    lsx0_rv = ini.prior$log_sigma_x_rv[subjs],
                    lsy0_rv = ini.prior$log_sigma_y_rv[subjs],
                    r0_rv = ini.prior$rota_rv[subjs],
                    lbd_T = lbds,
                    lbd_b = lbd_b,
                    lbd_T_rv = lbds_rv,
                    lbd_b_rv = lbd_b_rv
                    )
# initialize subject-specifc parameters
init = list(b = ini.prior$b[subjs],
            theta_x = ini.prior$theta_x[subjs],
            theta_y = ini.prior$theta_y[subjs],
            log_sigma_x = ini.prior$log_sigma_x[subjs],
            log_sigma_y = ini.prior$log_sigma_y[subjs],
            rota = ini.prior$rota[subjs],
            b_rv = ini.prior$b_rv[subjs],
            theta_x_rv = ini.prior$theta_x_rv[subjs],
            theta_y_rv = ini.prior$theta_y_rv[subjs],
            log_sigma_x_rv = ini.prior$log_sigma_x_rv[subjs],
            log_sigma_y_rv = ini.prior$log_sigma_y_rv[subjs],
            rota_rv = ini.prior$rota_rv[subjs])
myinit = list(list(theta_x = init$theta_x + 1,
                   theta_y = init$theta_y - 1, 
                   log_sigma_x = init$log_sigma_x+.01, 
                   log_sigma_y = init$log_sigma_y+.01, 
                   rota = init$rota, 
                   theta_x_rv = init$theta_x_rv,
                   theta_y_rv = init$theta_y_rv, 
                   log_sigma_x_rv = init$log_sigma_x_rv, 
                   log_sigma_y_rv = init$log_sigma_y_rv, 
                   rota_rv = init$rota_rv, 
                   b = init$b,
                   b_rv = init$b_rv,
                   phi = 1e-3
                   ),
              list(theta_x = init$theta_x - 1,
                   theta_y = init$theta_y - 1, 
                   log_sigma_x = init$log_sigma_x+.01, 
                   log_sigma_y = init$log_sigma_y+.01, 
                   rota = init$rota,
                   theta_x_rv = init$theta_x_rv,
                   theta_y_rv = init$theta_y_rv, 
                   log_sigma_x_rv = init$log_sigma_x_rv, 
                   log_sigma_y_rv = init$log_sigma_y_rv, 
                   rota_rv = init$rota_rv, 
                   b = init$b,
                   b_rv = init$b_rv,
                   phi = 1e-3
                   ),
              list(theta_x = init$theta_x + 1,
                   theta_y = init$theta_y + 1, 
                   log_sigma_x = init$log_sigma_x+.01, 
                   log_sigma_y = init$log_sigma_y+.01, 
                   rota = init$rota,
                   theta_x_rv = init$theta_x_rv,
                   theta_y_rv = init$theta_y_rv , 
                   log_sigma_x_rv = init$log_sigma_x_rv, 
                   log_sigma_y_rv = init$log_sigma_y_rv, 
                   rota_rv = init$rota_rv, 
                   b = init$b,
                   b_rv = init$b_rv,
                   phi = 1e-3),
              list(theta_x = init$theta_x - 1,
                   theta_y = init$theta_y + 1, 
                   log_sigma_x = init$log_sigma_x+.01, 
                   log_sigma_y = init$log_sigma_y+.01, 
                   rota = init$rota,
                   theta_x_rv = init$theta_x_rv,
                   theta_y_rv = init$theta_y_rv, 
                   log_sigma_x_rv = init$log_sigma_x_rv, 
                   log_sigma_y_rv = init$log_sigma_y_rv, 
                   rota_rv = init$rota_rv, 
                   b = init$b,
                   b_rv = init$b_rv,
                   phi = 1e-3))

# sampling by STAN ---------------------
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

stanmodel = stanc(file = "code/stan/sym_sqloss_2D.stan", verbose = T)
stanmodel = stan_model(model_code = stanmodel$model_code)

fit.painData = sampling(stanmodel,
                      data = data2d.input, 
                      chains = 4, 
                      pars = c("theta_x", "theta_y", "sigma_x", "sigma_y", "rota", "b",
                               "theta_x_rv", "theta_y_rv", "sigma_x_rv", "sigma_y_rv", "rota_rv", "b_rv", "phi", "lp__"),
                      iter = 1e4, 
                      seed = 1,
                      init = myinit,
		                  thin = 1)
  saveRDS(fit.painData, 
      file = paste0("sym_registration_Subject_",
                    subjs,
                    "_Ref_subj_",
                    ref.ind, 
	                  "_lbdT_", data2d.input$lbd_T,
                    "_lbdB_",data2d.input$lbd_b,
                    "_lbdT_rv_", data2d.input$lbd_T_rv,
                    "_lbdB_rv_",data2d.input$lbd_b_rv,
                    ".rds"))
gc()
