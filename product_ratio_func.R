#product-ratio method

setwd("~/Dropbox/Finished Work (sharefolder)/cs_joint_densities/code/product-ratio_method")
source('ER_GR.R')
source('density_norm.R')

library(pracma)
library(ftsa)
library(flexmix)
mape = ftsa:::mape

# product-ratio method
## INPUTS:
#  dat_1, dat_2: n_year * n_age
#  fh: forecasting horizion
#  fmethod: forecasting method for the component scores
## OUTPUTS:
#  fore_count_1, fore_count_2: p*fh
PR_forecast = function(dat_1,dat_2,fh,fmethod){
  
  n_year = nrow(dat_1)
  n_age = ncol(dat_1)
  year_index = rownames(dat_1)
  age_index = colnames(dat_1)
  dat_1_center = sweep(dat_1, 1, apply(dat_1, 1, sum), "/")
  dat_2_center = sweep(dat_2, 1, apply(dat_2, 1, sum), "/")
  alpha_x_1 = alpha_x_2 = vector("numeric", n_age)
  for (ik in 1:n_age) {
    alpha_x_1[ik] = geometric.mean(dat_1_center[, ik])
    alpha_x_2[ik] = geometric.mean(dat_2_center[, ik])
  }
  f_x_t_1 = f_x_t_2 = matrix(NA, n_year, n_age)
  for (ik in 1:n_year) {
    f_x_t_1[ik, ] = (dat_1[ik, ]/alpha_x_1)/sum(dat_1[ik, ]/alpha_x_1)
    f_x_t_2[ik, ] = (dat_2[ik, ]/alpha_x_2)/sum(dat_2[ik, ]/alpha_x_2)
  }
  g_t_1 = g_t_2 = vector("numeric", n_year)
  h_x_t_1 = h_x_t_2 = matrix(NA, n_year, n_age)
  for (ik in 1:n_year) {
    g_t_1[ik] = geometric.mean(f_x_t_1[ik, ])
    h_x_t_1[ik, ] = log(f_x_t_1[ik, ]/g_t_1[ik])
    g_t_2[ik] = geometric.mean(f_x_t_2[ik, ])
    h_x_t_2[ik, ] = log(f_x_t_2[ik, ]/g_t_2[ik])
  }
  colnames(h_x_t_1) = colnames(h_x_t_2) = age_index
  rownames(h_x_t_1) = rownames(h_x_t_2) = year_index
  
  # product-ratio method applied to the unconstrained beta (i.e, beta_j \equiv log(f_j) in Hyndman(2013))
  log_prod = (h_x_t_1 + h_x_t_2)/2
  log_ratio = (h_x_t_1 - h_x_t_2)/2
  
  dum_prod = ER_GR(data = log_prod)             
  ncomp_prod = max(dum_prod$k_ER, dum_prod$k_GR)
  
  dum_ratio = ER_GR(data = log_ratio)
  ncomp_ratio = max(dum_ratio$k_ER, dum_ratio$k_GR)
  rm(dum_prod); rm(dum_ratio)
  
  est_long_run_cov_prod = long_run_covariance_estimation(dat = t(log_prod))            # n_age * n_age
  eigen_prod = eigen(est_long_run_cov_prod, symmetric = TRUE)
  # eigen_prod = eigen(t(log_prod)%*%log_prod/(n_year*n_age), symmetric = TRUE)
  basis_prod = matrix(eigen_prod$vectors[,1:ncomp_prod],nrow=n_age,ncol=ncomp_prod)       # n_age * ncomp_prod
  score_prod = t(basis_prod) %*% t(log_prod)                                              # ncomp_prod * n_year
  recon_prod = basis_prod %*% score_prod                                                  # n_age * n_year
  
  
  est_long_run_cov_ratio = long_run_covariance_estimation(dat = t(log_ratio))
  eigen_ratio = eigen(est_long_run_cov_ratio, symmetric = TRUE)
  # eigen_ratio = eigen(t(log_ratio)%*%log_ratio/(n_year*n_age), symmetric = TRUE)
  basis_ratio = matrix(eigen_ratio$vectors[,1:ncomp_ratio],nrow=n_age,ncol=ncomp_ratio)
  score_ratio = t(basis_ratio) %*% t(log_ratio)
  recon_ratio = basis_ratio %*% score_ratio

  # reconstructed beta_F and beta_M
  recon_beta1 = recon_prod + recon_ratio                                                   # n_age * n_year
  recon_beta2 = recon_prod - recon_ratio                                                   # n_age * n_year
  
  # forecasts of principal component scores
  score_fore_prod = matrix(NA, ncomp_prod, fh)
  for (ik in 1:ncomp_prod){
    if (fmethod == "RWF_no_drift") {
      score_fore_prod[ik,] = rwf(as.numeric(as.numeric(score_prod[ik,])), h = fh, drift = FALSE)$mean
    }
    else if (fmethod == "RWF_drift") {
      score_fore_prod[ik,] = rwf(as.numeric(as.numeric(score_prod[ik,])), h = fh, drift = TRUE)$mean
    }
    else if (fmethod == "ETS") {
      score_fore_prod[ik,] = forecast(ets(as.numeric(score_prod[ik,])), h = fh)$mean
    }
    else if (fmethod == "ARIMA") {
      score_fore_prod[ik,] = forecast(auto.arima(as.numeric(score_prod[ik,])), h = fh)$mean
    }
    else {
      warning("Univariate time-series forecasting method is not listed.")
    }
    
  }
  
  score_fore_ratio = matrix(NA, ncomp_ratio, fh)
  for (ik in 1:ncomp_ratio){
    if (fmethod == "RWF_no_drift") {
      score_fore_ratio[ik,] = rwf(as.numeric(as.numeric(score_ratio[ik,])), h = fh, drift = FALSE)$mean
    }
    else if (fmethod == "RWF_drift") {
      score_fore_ratio[ik,] = rwf(as.numeric(as.numeric(score_ratio[ik,])), h = fh, drift = TRUE)$mean
    }
    else if (fmethod == "ETS") {
      score_fore_ratio[ik,] = forecast(ets(as.numeric(score_ratio[ik,])), h = fh)$mean
    }
    else if (fmethod == "ARIMA") {
      score_fore_ratio[ik,] = forecast(auto.arima(as.numeric(score_ratio[ik,])), h = fh)$mean
    }
    else {
      warning("Univariate time-series forecasting method is not listed.")
    }
  }
  
  colnames(score_fore_prod) = colnames(score_fore_ratio) = 1:fh
  
  # obtain forecasts in log-scale
  fore_log_prod = basis_prod %*% score_fore_prod                                           # n_age * fh
  fore_log_ratio = basis_ratio %*% score_fore_ratio
  
  fore_beta1 = fore_log_prod + fore_log_ratio
  fore_beta2 = fore_log_prod - fore_log_ratio
  
  # obtain forecasts in original scale
  f_x_t_star_recon_1 = f_x_t_star_recon_2 = d_x_t_star_recon_1 = d_x_t_star_recon_2 = matrix(NA, 
                                                                                             n_age, n_year)
  for (ik in 1:n_year) {
    f_x_t_star_recon_1[, ik] = exp(recon_beta1[, ik])/sum(exp(recon_beta1[, ik]))
    f_x_t_star_recon_2[, ik] = exp(recon_beta2[, ik])/sum(exp(recon_beta2[, ik]))
    d_x_t_star_recon_1[, ik] = (f_x_t_star_recon_1[, ik] * 
                                  alpha_x_1)/sum(f_x_t_star_recon_1[, ik] * alpha_x_1)
    d_x_t_star_recon_2[, ik] = (f_x_t_star_recon_2[, ik] * 
                                  alpha_x_2)/sum(f_x_t_star_recon_2[, ik] * alpha_x_2)
  }
  R2_1 = 1 - sum((t(d_x_t_star_recon_1) * 10^5 - dat_1)^2)/sum((dat_1 - 
                                                                  colMeans(dat_1))^2)
  R2_2 = 1 - sum((t(d_x_t_star_recon_2) * 10^5 - dat_2)^2)/sum((dat_2 - 
                                                                  colMeans(dat_2))^2)
  f_x_t_star_fore_1 = f_x_t_star_fore_2 = d_x_t_star_fore_1 = d_x_t_star_fore_2 = matrix(NA, 
                                                                                         n_age, fh)
  for (ik in 1:fh) {
    f_x_t_star_fore_1[, ik] = exp(fore_beta1[, ik])/sum(exp(fore_beta1[, ik]))
    f_x_t_star_fore_2[, ik] = exp(fore_beta2[, ik])/sum(exp(fore_beta2[, ik]))
    d_x_t_star_fore_1[, ik] = (f_x_t_star_fore_1[, ik] * 
                                 alpha_x_1)/sum((f_x_t_star_fore_1[, ik] * alpha_x_1))
    d_x_t_star_fore_2[, ik] = (f_x_t_star_fore_2[, ik] * 
                                 alpha_x_2)/sum((f_x_t_star_fore_2[, ik] * alpha_x_2))
  }
  colnames(d_x_t_star_fore_1) = colnames(d_x_t_star_fore_2) = 1:fh
  rownames(d_x_t_star_fore_1) = rownames(d_x_t_star_fore_2) = age_index
  
  return(list(R2_1 = R2_1, R2_2 = R2_2, 
              recon_1 = d_x_t_star_recon_1 * 10^5, recon_2 = d_x_t_star_recon_2 * 10^5, 
              fore_count_1 = d_x_t_star_fore_1 * 10^5, fore_count_2 = d_x_t_star_fore_2 * 10^5, 
              alpha_x_1 = alpha_x_1, alpha_x_2 = alpha_x_2, 
              h_x_t_1 = h_x_t_1, h_x_t_2 = h_x_t_2))
  
}

#expanding window forecast
PR_eval <- function(dat_1,dat_2, fh,fmethod){
  # sub1, sub2: T *n
  # fh: forecast horizon
  # fmethod: forecasting method
  
  Tt = dim(dat_1)[1]
  p = ncol(dat_1)
  train = Tt-31
  
  fore_h30_1 = fore_h30_2 = matrix(NA, p, (31 - fh))                
  
  fore_h30_1_KLdiv = fore_h30_2_KLdiv = matrix(NA, (31 - fh), 2)
  
  fore_h30_1_density_norm = fore_h30_1_mape = fore_h30_1_JSdiv = fore_h30_1_JSdiv_geo = 
    fore_h30_2_density_norm = fore_h30_2_mape = fore_h30_2_JSdiv = fore_h30_2_JSdiv_geo =  vector("numeric", (31 - fh))
  
  for(ik in 1:(31 - fh))
  {
    
    PR = PR_forecast(dat_1 = dat_1[1:(train+ik),], dat_2 = dat_2[1:(train+ik),], fh = fh,fmethod=fmethod)
    fore_h30_1[,ik] = PR$fore_count_1[,fh]
    fore_h30_2[,ik] = PR$fore_count_2[,fh]
    
    # MAPE 
    
    fore_h30_1_mape[ik] = mape(forecast = t(fore_h30_1[,ik]), true = dat_1[(train+fh+ik),])        
    fore_h30_2_mape[ik] = mape(forecast = t(fore_h30_2[,ik]), true = dat_2[(train+fh+ik),])        
    
    # density_norm
    
    fore_h30_1_density_norm[ik] = density_norm(d1 = t(fore_h30_1[,ik]), d2 = dat_1[(train+fh+ik),], 
                                                     time_vector = 1:111)
    fore_h30_2_density_norm[ik] = density_norm(d1 = t(fore_h30_2[,ik]), d2 = dat_2[(train+fh+ik),], 
                                               time_vector = 1:111)
    
    # KLdiv
    
    dat = cbind(true = dat_1[(train+fh+ik),], forecast = as.numeric(t(fore_h30_1[,ik])))
    colnames(dat) = c("True", "Estimate")
    fore_h30_1_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
    rm(dat)
    
    dat = cbind(true = dat_2[(train+fh+ik),], forecast = as.numeric(t(fore_h30_2[,ik])))
    colnames(dat) = c("True", "Estimate")
    fore_h30_2_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
    rm(dat)
    
    # JSdiv(simple)
    
    dat = cbind(true = dat_1[(train+fh+ik),], forecast = as.numeric(t(fore_h30_1[,ik])))
    M = rowMeans(dat)
    P_M = cbind(dat_1[(train+fh+ik),], M)
    E_M = cbind(as.numeric(t(fore_h30_1[,ik])), M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    fore_h30_1_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M)
    
    dat = cbind(true = dat_2[(train+fh+ik),], forecast = as.numeric(t(fore_h30_2[,ik])))
    M = rowMeans(dat)
    P_M = cbind(dat_2[(train+fh+ik),], M)
    E_M = cbind(as.numeric(t(fore_h30_2[,ik])), M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    fore_h30_2_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M)
    
    # JSdiv(geo)
    
    dat = cbind(true = dat_1[(train+fh+ik),], forecast = as.numeric(t(fore_h30_1[,ik])))
    M = apply(dat, 1, geometric.mean)
    P_M = cbind(dat_1[(train+fh+ik),], M)
    E_M = cbind(as.numeric(t(fore_h30_1[,ik])), M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    fore_h30_1_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M)
    
    dat = cbind(true = dat_2[(train+fh+ik),], forecast = as.numeric(t(fore_h30_2[,ik])))
    M = apply(dat, 1, geometric.mean)
    P_M = cbind(dat_2[(train+fh+ik),], M)
    E_M = cbind(as.numeric(t(fore_h30_2[,ik])), M)
    colnames(E_M) = colnames(P_M) = c("True", "M")
    fore_h30_2_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
    rm(M); rm(P_M); rm(E_M)
    
  }
  
 
  
  ###############################################################
  # take mean over the number of years in the forecasting period
  ###############################################################
  
  return(list(PR_1_mape = mean(fore_h30_1_mape), 
              PR_2_mape = mean(fore_h30_2_mape),
              
              PR_1_density_norm = mean(fore_h30_1_density_norm),
              PR_2_density_norm = mean(fore_h30_2_density_norm),
              
              
              PR_1_KLdiv = mean(fore_h30_1_KLdiv),
              PR_2_KLdiv = mean(fore_h30_2_KLdiv),
              
              
              PR_1_JSdiv = mean(fore_h30_1_JSdiv),
              PR_2_JSdiv = mean(fore_h30_2_JSdiv),
              
              PR_1_JSdiv_geo = mean(fore_h30_1_JSdiv_geo),
              PR_2_JSdiv_geo = mean(fore_h30_2_JSdiv_geo)))
              
}











# cohort result
setwd("~/Dropbox/Finished Work (sharefolder)/cs_joint_densities/data/cohort_lt")
female.int = t(as.matrix(read.table('uk_female_pop_complete.txt')))
male.int = t(as.matrix(read.table('uk_male_pop_complete.txt')))

#replace the 0 death count with 0.01
death0.f = which(is.na(female.int) | female.int==0)
death0.m = which(is.na(male.int) | male.int==0)

EW_female_pop = replace(female.int,death0.f,0.01) #total is 130 years
EW_male_pop = replace(male.int,death0.m,0.01)


PR_1_mape = PR_1_density_norm = PR_1_KLdiv = PR_1_JSdiv = PR_1_JSdiv_geo =
  PR_2_err = PR_2_mape = PR_2_density_norm = PR_2_KLdiv = PR_2_JSdiv = PR_2_JSdiv_geo= vector("numeric", 30)

for(ik in 1:30)
{
  print(ik)
  dum = PR_eval(dat_1 = EW_female_pop, dat_2 = EW_male_pop, fh = ik,fmethod='ARIMA')
  
  ## 1st population
  
  # MAPE
  PR_1_mape[ik] = dum$PR_1_mape
  PR_1_density_norm[ik] = dum$PR_1_density_norm
  PR_1_KLdiv[ik] = dum$PR_1_KLdiv
  PR_1_JSdiv[ik] = dum$PR_1_JSdiv
  PR_1_JSdiv_geo[ik] =dum$PR_1_JSdiv_geo
    
  PR_2_mape[ik] = dum$PR_2_mape
  PR_2_density_norm[ik] = dum$PR_2_density_norm
  PR_2_KLdiv[ik] = dum$PR_2_KLdiv
  PR_2_JSdiv[ik] = dum$PR_2_JSdiv
  PR_2_JSdiv_geo[ik] =dum$PR_2_JSdiv_geo 
}

round(mean(PR_1_mape),3)
round(mean(PR_1_KLdiv),3)
round(mean(PR_1_JSdiv),3)
round(mean(PR_1_JSdiv),3)

round(mean(PR_2_mape),3)
round(mean(PR_2_KLdiv),3)
round(mean(PR_2_JSdiv),3)
round(mean(PR_2_JSdiv_geo),3)


# period result
setwd("~/Dropbox/Finished Work (sharefolder)/cs_joint_densities/data/period_lt")


female_qx = t(matrix(read.table("EW_lt_female_death.txt", header = TRUE)[,4], 111, 178))
male_qx   = t(matrix(read.table("EW_lt_male_death.txt", header = TRUE)[,4], 111, 178))

n_col = ncol(female_qx)
n_row = nrow(female_qx)
female_pop = male_pop = matrix(NA, n_row, n_col)
for(ij in 1:n_row)
{
  start_pop_female = start_pop_male = 100000
  for(ik in 1:n_col)
  {
    female_pop[ij,ik] = female_qx[ij,ik] * start_pop_female
    start_pop_female = start_pop_female - female_pop[ij,ik]
    
    male_pop[ij,ik] = male_qx[ij,ik] * start_pop_male
    start_pop_male = start_pop_male - male_pop[ij,ik]
    rm(ik)
  }
  rm(ij)
}
colnames(female_pop) = colnames(male_pop) = 0:110
rownames(female_pop) = rownames(male_pop) = 1841:2018

EW_female_pop_period = t(female_pop)
EW_male_pop_period = t(male_pop)
colnames(EW_female_pop_period) = colnames(EW_male_pop_period) = 1841:2018
rownames(EW_female_pop_period) = rownames(EW_male_pop_period) = 0:110 


PR_1_mape_period = PR_1_density_norm_period = PR_1_KLdiv_period = PR_1_JSdiv_period = PR_1_JSdiv_geo_period =
  PR_2_err_period = PR_2_mape_period = PR_2_density_norm_period = PR_2_KLdiv_period = PR_2_JSdiv_period = PR_2_JSdiv_geo_period = vector("numeric", 30)

for(ik in 1:30)
{
  print(ik)
  dum = PR_eval(dat_1 = t(EW_female_pop_period), dat_2 = t(EW_male_pop_period), fh = ik,fmethod='ARIMA')
  
  ## 1st population
  
  # MAPE
  PR_1_mape_period[ik] = dum$PR_1_mape
  PR_1_density_norm_period[ik] = dum$PR_1_density_norm
  PR_1_KLdiv_period[ik] = dum$PR_1_KLdiv
  PR_1_JSdiv_period[ik] = dum$PR_1_JSdiv
  PR_1_JSdiv_geo_period[ik] =dum$PR_1_JSdiv_geo
  
  PR_2_mape_period[ik] = dum$PR_2_mape
  PR_2_density_norm_period[ik] = dum$PR_2_density_norm
  PR_2_KLdiv_period[ik] = dum$PR_2_KLdiv
  PR_2_JSdiv_period[ik] = dum$PR_2_JSdiv
  PR_2_JSdiv_geo_period[ik] =dum$PR_2_JSdiv_geo 
}

round(mean(PR_1_mape_period),3)
# mean(PR_1_density_norm_period)
round(mean(PR_1_KLdiv_period),3)
round(mean(PR_1_JSdiv_period),3)
round(mean(PR_1_JSdiv_geo_period),3)

round(mean(PR_2_mape_period),3)
# mean(PR_2_density_norm_period)
round(mean(PR_2_KLdiv_period),3)
round(mean(PR_2_JSdiv_period),3)
round(mean(PR_2_JSdiv_geo_period),3)



