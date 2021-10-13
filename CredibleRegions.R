# credible intervals ------------------------------------------------------
source("supp_2D.R")
# library(gplots)
library(dplyr)
library(alphahull)
library(rstan)

# library(igraph)
vertices <- as.matrix(sweep(data.frame(x = c(rep(10, 12), rep(25, 12), 10:25,       10:25),
                                       y = c(13:24,       13:24,       rep(13, 16), rep(24, 16))), 
                            1, c(15,15), "-"))
c0 = c(0,0)
nsubj = 33
A = array(dim = c(2,2,nsubj))
b = array(dim = c(nsubj,2))
subjs = 1
ref.ind = 20
mean_coord = 15

lbd_list = c(1, 15, 40, 70, 110)
lbd_b_list = c(.01, 1, 10, 20)
lbd_list_rv = c(1, 15, 40, 70, 110)
lbd_b_list_rv = c(.01, 1, 10, 20)

M = length(lbd_list)
N = length(lbd_b_list)
P = length(lbd_list_rv)
Q = length(lbd_b_list_rv)

indmat = array(1:(M*N*P*Q), dim = c(M, N, P, Q))

for (subjs in 1:33){
  print(subjs)
  elpds <- readRDS(paste0("../results/waic/sym_elpds_subject_",subjs,".rds"))
  elpds[elpds > min(elpds, na.rm = T) + 20] = NaN
  indbT = which(indmat==which.max(elpds), arr.ind = T)
  if( abs(diff(range(elpds,na.rm = T))) < 2) lbds <- lbd_list[!is.na(elpds)][1]
  fits <- readRDS(file =paste0("sym_registration_Subject_",
                               subjs,
                               "_Ref_subj_",
                               ref.ind, 
                               "_lbdT_", lbd_list[indbT[1]],
                               "_lbdB_",lbd_b_list[indbT[2]],
                               "_lbdT_rv_", lbd_list_rv[indbT[3]],
                               "_lbdB_rv_",lbd_b_list_rv[indbT[4]],
                               ".rds"))
  
  samples.par <- extract(model.selected, permuted = F)
  samples.par.chain <- samples.par[,1,]
  samples.par.median <- apply(samples.par.chain,2, median) %>% as.list()
  
  A[,,subjs] = diag(c(samples.par.median$sigma_x, samples.par.median$sigma_y)) %*% t(rotate(samples.par.median$rota))
  b[subjs, ] = c(samples.par.median$theta_x, samples.par.median$theta_y)
  vert.array <- array(dim = c(dim(samples.par)[1], dim(vertices)[1],2))
  vert.array.inv <- vert.array
  vert.median <- (A[,,subjs] %*% t(vertices) + b[subjs,] + as.vector(A[,,subjs] %*% as.matrix(c0))) + mean_coord
  
  for (ind.samples in 1: dim(samples.par)[1]){
    samples <- samples.par.chain[ind.samples, ]
    As = diag(c(samples["sigma_x"], samples["sigma_y"])) %*% t(rotate(samples["rota"]))
    vert.array[ind.samples,,] <- 
      t(As %*% t(vertices) +
          c(samples["theta_x"], samples["theta_y"])+
          as.vector(As %*% as.matrix(c0)) + mean_coord)
  }
  
  samp.pars <- samples.par.chain[,c("theta_x", "theta_y", "sigma_x", "sigma_y", "rota")] %>% as.data.frame()
  samp.pars.T <- apply(samp.pars, 1, function(r){
    A = diag(c(r["sigma_x"], r["sigma_y"])) %*% t(rotate(r["rota"]))
    return(c(as.vector(A), r[1:2]))
  }) %>% t()
  # 
  elist <- seq(.1, 2, by = .01)
  prop <- elist
  
  for (e in 1:length(elist)){
    T.db <- dbscan::dbscan(samp.pars.T, eps = elist[e])
    prop[e] <- max(tabulate(T.db$cluster)/dim(samples.par)[1])
  }
  e.select <- which.min(abs(prop-.95))
  T.db <- dbscan::dbscan(samp.pars.T, eps = elist[e.select])
  indicators <- T.db$cluster==which.max(tabulate(T.db$cluster))
  
  vert.array.inside <- matrix(vert.array[indicators,,], sum(indicators) * dim(vertices)[1], 2)
  vert.array.inside.unique <- unique(vert.array.inside)
  
  vert.ashape <- ashape(vert.array.inside.unique, alpha = .5)
  
  R.matlab::writeMat(paste0("../data/CIs/sym_edges_", subjs, ".mat"), edges = vert.ashape$edges[,3:4])
}
R.matlab::writeMat(paste0("../data/sym_posters_30by30_round_3pks_realdata.mat"), A = A, b = b)

