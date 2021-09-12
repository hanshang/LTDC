################
# CoDa function
################

require(psych)
require(ftsa)


R_square_fit=function (dat_1, dat_2, fh, modeling_method = c("FTS", 
                                                "MFTS", "MLFTS"), forecasting_method = c("RWF_no_drift", 
                                                                                         "RWF_drift", "ETS", "ARIMA")) 
{
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
    f_x_t_1[ik, ] = (dat_1[ik, ]/alpha_x_1)/sum(dat_1[ik, 
                                                      ]/alpha_x_1)
    f_x_t_2[ik, ] = (dat_2[ik, ]/alpha_x_2)/sum(dat_2[ik, 
                                                      ]/alpha_x_2)
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
  if (modeling_method == "FTS") {
    est_long_run_cov = long_run_covariance_estimation(dat = t(h_x_t_1))
    est_eigen_decomp = eigen(est_long_run_cov, symmetric = TRUE)
    lambda_val = est_eigen_decomp$values[which(est_eigen_decomp$values > 
                                                 0)]
    k_max = length(which(lambda_val >= mean(lambda_val)))
    tau = 1/log(max(lambda_val[1], length(lambda_val)))
    eigen_val_ratio = vector("numeric", k_max)
    for (ik in 1:k_max) {
      eigen_val_ratio[ik] = lambda_val[ik + 1]/lambda_val[ik] * 
        ifelse(lambda_val[ik]/lambda_val[1] >= tau, 1, 
               0) + ifelse(lambda_val[ik]/lambda_val[1] < 
                             tau, 1, 0)
    }
    ncomp_1 = which.min(eigen_val_ratio)
    rm(tau)
    rm(eigen_val_ratio)
    rm(k_max)
    rm(lambda_val)
    rm(est_long_run_cov)
    rm(est_eigen_decomp)
    est_long_run_cov = long_run_covariance_estimation(dat = t(h_x_t_2))
    est_eigen_decomp = eigen(est_long_run_cov, symmetric = TRUE)
    lambda_val = est_eigen_decomp$values[which(est_eigen_decomp$values > 
                                                 0)]
    k_max = length(which(lambda_val >= mean(lambda_val)))
    tau = 1/log(max(lambda_val[1], length(lambda_val)))
    eigen_val_ratio = vector("numeric", k_max)
    for (ik in 1:k_max) {
      eigen_val_ratio[ik] = lambda_val[ik + 1]/lambda_val[ik] * 
        ifelse(lambda_val[ik]/lambda_val[1] >= tau, 1, 
               0) + ifelse(lambda_val[ik]/lambda_val[1] < 
                             tau, 1, 0)
    }
    ncomp_2 = which.min(eigen_val_ratio)
    rm(tau)
    rm(eigen_val_ratio)
    rm(k_max)
    rm(lambda_val)
    rm(est_long_run_cov)
    rm(est_eigen_decomp)
    SVD_decomp_1 = svd(h_x_t_1)
    SVD_decomp_2 = svd(h_x_t_2)
    basis_1 = SVD_decomp_1$v[, 1:ncomp_1]
    basis_2 = SVD_decomp_2$v[, 1:ncomp_2]
    score_1 = t(basis_1) %*% t(h_x_t_1)
    score_2 = t(basis_2) %*% t(h_x_t_2)
    recon_1 = basis_1 %*% score_1
    recon_2 = basis_2 %*% score_2
    colnames(score_1) = colnames(score_2) = year_index
    colnames(recon_1) = colnames(recon_2) = year_index
    score_fore_1 = matrix(NA, ncomp_1, fh)
    score_fore_2 = matrix(NA, ncomp_2, fh)
    for (ik in 1:ncomp_1) {
      if (forecasting_method == "RWF_no_drift") {
        score_fore_1[ik, ] = rwf(as.numeric(score_1[ik, 
                                                    ]), h = fh, drift = FALSE)$mean
      }
      else if (forecasting_method == "RWF_drift") {
        score_fore_1[ik, ] = rwf(as.numeric(score_1[ik, 
                                                    ]), h = fh, drift = TRUE)$mean
      }
      else if (forecasting_method == "ETS") {
        score_fore_1[ik, ] = forecast(ets(as.numeric(score_1[ik, 
                                                             ])), h = fh)$mean
      }
      else if (forecasting_method == "ARIMA") {
        score_fore_1[ik, ] = forecast(auto.arima(as.numeric(score_1[ik, 
                                                                    ])), h = fh)$mean
      }
      else {
        warning("Univariate time-series forecasting method is not listed.")
      }
    }
    for (ik in 1:ncomp_2) {
      if (forecasting_method == "RWF_no_drift") {
        score_fore_2[ik, ] = rwf(as.numeric(score_2[ik, 
                                                    ]), h = fh, drift = FALSE)$mean
      }
      else if (forecasting_method == "RWF_drift") {
        score_fore_2[ik, ] = rwf(as.numeric(score_2[ik, 
                                                    ]), h = fh, drift = TRUE)$mean
      }
      else if (forecasting_method == "ETS") {
        score_fore_2[ik, ] = forecast(ets(as.numeric(score_2[ik, 
                                                             ])), h = fh)$mean
      }
      else if (forecasting_method == "ARIMA") {
        score_fore_2[ik, ] = forecast(auto.arima(as.numeric(score_2[ik, 
                                                                    ])), h = fh)$mean
      }
      else {
        warning("Univariate time-series forecasting method is not listed.")
      }
    }
    colnames(score_fore_1) = colnames(score_fore_2) = 1:fh
    fore_val_1 = basis_1 %*% score_fore_1
    fore_val_2 = basis_2 %*% score_fore_2
  }
  else if (modeling_method == "MFTS") {
    data_comb = list()
    data_comb[[1]] = t(h_x_t_1)
    data_comb[[2]] = t(h_x_t_2)
    rowmeans_object = sd_object = decenter_object = list()
    for (ik in 1:2) {
      rowmeans_object[[ik]] = rowMeans(data_comb[[ik]], 
                                       na.rm = TRUE)
      sd_object[[ik]] = apply(data_comb[[ik]], 1, sd, na.rm = TRUE)
      decenter_object[[ik]] = t(scale(t(data_comb[[ik]]), 
                                      center = TRUE, scale = TRUE))
    }
    comb_object = do.call(rbind, decenter_object)
    est_long_run_cov = long_run_covariance_estimation(dat = comb_object)
    est_eigen_decomp = eigen(est_long_run_cov, symmetric = TRUE)
    lambda_val = est_eigen_decomp$values[which(est_eigen_decomp$values > 
                                                 0)]
    k_max = length(which(lambda_val >= mean(lambda_val)))
    tau = 1/log(max(lambda_val[1], length(lambda_val)))
    eigen_val_ratio = vector("numeric", k_max)
    for (ik in 1:k_max) {
      eigen_val_ratio[ik] = lambda_val[ik + 1]/lambda_val[ik] * 
        ifelse(lambda_val[ik]/lambda_val[1] >= tau, 1, 
               0) + ifelse(lambda_val[ik]/lambda_val[1] < 
                             tau, 1, 0)
    }
    ncomp = which.min(eigen_val_ratio)
    rm(tau)
    rm(eigen_val_ratio)
    rm(k_max)
    rm(lambda_val)
    rm(est_long_run_cov)
    rm(est_eigen_decomp)
    SVD_decomp = svd(t(comb_object))
    basis = as.matrix(SVD_decomp$v[, 1:ncomp])
    score = t(basis) %*% comb_object
    recon = (basis %*% score) * do.call(c, sd_object) + do.call(c, 
                                                                rowmeans_object)
    recon_1 = recon[1:n_age, ]
    recon_2 = recon[(n_age + 1):(n_age * 2), ]
    score_fore = matrix(NA, ncomp, fh)
    for (ik in 1:ncomp) {
      if (forecasting_method == "RWF_no_drift") {
        score_fore[ik, ] = rwf(as.numeric(score[ik, ]), 
                               h = fh, drift = FALSE)$mean
      }
      else if (forecasting_method == "RWF_drift") {
        score_fore[ik, ] = rwf(as.numeric(score[ik, ]), 
                               h = fh, drift = TRUE)$mean
      }
      else if (forecasting_method == "ETS") {
        score_fore[ik, ] = forecast(ets(as.numeric(score[ik, 
                                                         ])), h = fh)$mean
      }
      else if (forecasting_method == "ARIMA") {
        score_fore[ik, ] = forecast(auto.arima(as.numeric(score[ik, 
                                                                ])), h = fh)$mean
      }
      else {
        warning("Univariate time-series forecasting method is not listed.")
      }
    }
    colnames(score_fore) = 1:fh
    fore_val_1 = as.matrix(((basis %*% score_fore) * do.call(c, 
                                                             sd_object) + do.call(c, rowmeans_object))[1:n_age, 
                                                                                                       ])
    fore_val_2 = as.matrix(((basis %*% score_fore) * do.call(c, 
                                                             sd_object) + do.call(c, rowmeans_object))[(n_age + 
                                                                                                          1):(n_age * 2), ])
  }
  else if (modeling_method == "MLFTS") {
    h_x_t_ave = (h_x_t_1 + h_x_t_2)/2
    est_long_run_cov = long_run_covariance_estimation(dat = t(h_x_t_ave))
    est_eigen_decomp = eigen(est_long_run_cov, symmetric = TRUE)
    lambda_val = est_eigen_decomp$values[which(est_eigen_decomp$values > 
                                                 0)]
    k_max = length(which(lambda_val >= mean(lambda_val)))
    tau = 1/log(max(lambda_val[1], length(lambda_val)))
    eigen_val_ratio = vector("numeric", k_max)
    for (ik in 1:k_max) {
      eigen_val_ratio[ik] = lambda_val[ik + 1]/lambda_val[ik] * 
        ifelse(lambda_val[ik]/lambda_val[1] >= tau, 1, 
               0) + ifelse(lambda_val[ik]/lambda_val[1] < 
                             tau, 1, 0)
    }
    ncomp_ave = which.min(eigen_val_ratio)
    rm(tau)
    rm(eigen_val_ratio)
    rm(k_max)
    rm(lambda_val)
    rm(est_long_run_cov)
    rm(est_eigen_decomp)
    SVD_decomp_ave = svd(h_x_t_ave)
    basis_ave = SVD_decomp_ave$v[, 1:ncomp_ave]
    basis_ave_svd_value = SVD_decomp_ave$d[1:ncomp_ave]
    score_ave = matrix(t(basis_ave) %*% t(h_x_t_ave), ncomp_ave, 
    )
    recon_ave = basis_ave %*% score_ave
    h_x_t_resi_1 = h_x_t_1 - t(recon_ave)
    h_x_t_resi_2 = h_x_t_2 - t(recon_ave)
    est_long_run_cov = long_run_covariance_estimation(dat = t(h_x_t_resi_1))
    est_eigen_decomp = eigen(est_long_run_cov, symmetric = TRUE)
    lambda_val = est_eigen_decomp$values[which(est_eigen_decomp$values > 
                                                 0)]
    k_max = length(which(lambda_val >= mean(lambda_val)))
    tau = 1/log(max(lambda_val[1], length(lambda_val)))
    eigen_val_ratio = vector("numeric", k_max)
    for (ik in 1:k_max) {
      eigen_val_ratio[ik] = lambda_val[ik + 1]/lambda_val[ik] * 
        ifelse(lambda_val[ik]/lambda_val[1] >= tau, 1, 
               0) + ifelse(lambda_val[ik]/lambda_val[1] < 
                             tau, 1, 0)
    }
    ncomp_resi_1 = which.min(eigen_val_ratio)
    rm(tau)
    rm(eigen_val_ratio)
    rm(k_max)
    rm(lambda_val)
    rm(est_long_run_cov)
    rm(est_eigen_decomp)
    est_long_run_cov = long_run_covariance_estimation(dat = t(h_x_t_resi_2))
    est_eigen_decomp = eigen(est_long_run_cov, symmetric = TRUE)
    lambda_val = est_eigen_decomp$values[which(est_eigen_decomp$values > 
                                                 0)]
    k_max = length(which(lambda_val >= mean(lambda_val)))
    tau = 1/log(max(lambda_val[1], length(lambda_val)))
    eigen_val_ratio = vector("numeric", k_max)
    for (ik in 1:k_max) {
      eigen_val_ratio[ik] = lambda_val[ik + 1]/lambda_val[ik] * 
        ifelse(lambda_val[ik]/lambda_val[1] >= tau, 1, 
               0) + ifelse(lambda_val[ik]/lambda_val[1] < 
                             tau, 1, 0)
    }
    ncomp_resi_2 = which.min(eigen_val_ratio)
    rm(tau)
    rm(eigen_val_ratio)
    rm(k_max)
    rm(lambda_val)
    rm(est_long_run_cov)
    rm(est_eigen_decomp)
    SVD_decomp_resi_1 = svd(h_x_t_resi_1)
    basis_resi_1 = SVD_decomp_resi_1$v[, 1:ncomp_resi_1]
    basis_resi_1_svd_value = SVD_decomp_resi_1$d[1:ncomp_resi_1]
    score_resi_1 = matrix(t(basis_resi_1) %*% t(h_x_t_resi_1), 
                          ncomp_resi_1, )
    recon_resi_1 = basis_resi_1 %*% score_resi_1
    SVD_decomp_resi_2 = svd(h_x_t_resi_2)
    basis_resi_2 = SVD_decomp_resi_2$v[, 1:ncomp_resi_2]
    basis_resi_2_svd_value = SVD_decomp_resi_2$d[1:ncomp_resi_2]
    score_resi_2 = matrix(t(basis_resi_2) %*% t(h_x_t_resi_2), 
                          ncomp_resi_2, )
    recon_resi_2 = basis_resi_2 %*% score_resi_2
    recon_1 = recon_ave + recon_resi_1
    recon_2 = recon_ave + recon_resi_2
    score_ave_fore = matrix(NA, ncomp_ave, fh)
    score_resi_fore_1 = matrix(NA, ncomp_resi_1, fh)
    score_resi_fore_2 = matrix(NA, ncomp_resi_2, fh)
    for (ik in 1:ncomp_ave) {
      if (forecasting_method == "RWF_no_drift") {
        score_ave_fore[ik, ] = rwf(as.numeric(score_ave[ik, 
                                                        ]), h = fh, drift = FALSE)$mean
      }
      else if (forecasting_method == "RWF_drift") {
        score_ave_fore[ik, ] = rwf(as.numeric(score_ave[ik, 
                                                        ]), h = fh, drift = TRUE)$mean
      }
      else if (forecasting_method == "ETS") {
        score_ave_fore[ik, ] = forecast(ets(as.numeric(score_ave[ik, 
                                                                 ])), h = fh)$mean
      }
      else if (forecasting_method == "ARIMA") {
        score_ave_fore[ik, ] = forecast(auto.arima(as.numeric(score_ave[ik, 
                                                                        ])), h = fh)$mean
      }
      else {
        warning("Univariate time-series forecasting method is not listed.")
      }
    }
    for (ik in 1:ncomp_resi_1) {
      if (forecasting_method == "RWF_no_drift") {
        score_resi_fore_1[ik, ] = rwf(as.numeric(score_resi_1[ik, 
                                                              ]), h = fh, drift = FALSE)$mean
      }
      else if (forecasting_method == "RWF_drift") {
        score_resi_fore_1[ik, ] = rwf(as.numeric(score_resi_1[ik, 
                                                              ]), h = fh, drift = TRUE)$mean
      }
      else if (forecasting_method == "ETS") {
        score_resi_fore_1[ik, ] = forecast(ets(as.numeric(score_resi_1[ik, 
                                                                       ])), h = fh)$mean
      }
      else if (forecasting_method == "ARIMA") {
        score_resi_fore_1[ik, ] = forecast(auto.arima(as.numeric(score_resi_1[ik, 
                                                                              ])), h = fh)$mean
      }
      else {
        warning("Univariate time-series forecasting method is not listed.")
      }
    }
    for (ik in 1:ncomp_resi_2) {
      if (forecasting_method == "RWF_no_drift") {
        score_resi_fore_2[ik, ] = rwf(as.numeric(score_resi_2[ik, 
                                                              ]), h = fh, drift = FALSE)$mean
      }
      else if (forecasting_method == "RWF_drift") {
        score_resi_fore_2[ik, ] = rwf(as.numeric(score_resi_2[ik, 
                                                              ]), h = fh, drift = TRUE)$mean
      }
      else if (forecasting_method == "ETS") {
        score_resi_fore_2[ik, ] = forecast(ets(as.numeric(score_resi_2[ik, 
                                                                       ])), h = fh)$mean
      }
      else if (forecasting_method == "ARIMA") {
        score_resi_fore_2[ik, ] = forecast(auto.arima(as.numeric(score_resi_2[ik, 
                                                                              ])), h = fh)$mean
      }
      else {
        warning("Univariate time-series forecasting method is not listed.")
      }
    }
    colnames(score_ave_fore) = colnames(score_resi_fore_1) = colnames(score_resi_fore_2) = 1:fh
    fore_val_1 = basis_ave %*% score_ave_fore + basis_resi_1 %*% 
      score_resi_fore_1
    fore_val_2 = basis_ave %*% score_ave_fore + basis_resi_2 %*% 
      score_resi_fore_2
  }
  else {
    warning("Modeling_method needs to be selected from the provided list.")
  }
  f_x_t_star_recon_1 = f_x_t_star_recon_2 = d_x_t_star_recon_1 = d_x_t_star_recon_2 = matrix(NA, 
                                                                                             n_age, n_year)
  for (ik in 1:n_year) {
    f_x_t_star_recon_1[, ik] = exp(recon_1[, ik])/sum(exp(recon_1[, 
                                                                  ik]))
    f_x_t_star_recon_2[, ik] = exp(recon_2[, ik])/sum(exp(recon_2[, 
                                                                  ik]))
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
    f_x_t_star_fore_1[, ik] = exp(fore_val_1[, ik])/sum(exp(fore_val_1[, 
                                                                       ik]))
    f_x_t_star_fore_2[, ik] = exp(fore_val_2[, ik])/sum(exp(fore_val_2[, 
                                                                       ik]))
    d_x_t_star_fore_1[, ik] = (f_x_t_star_fore_1[, ik] * 
                                 alpha_x_1)/sum((f_x_t_star_fore_1[, ik] * alpha_x_1))
    d_x_t_star_fore_2[, ik] = (f_x_t_star_fore_2[, ik] * 
                                 alpha_x_2)/sum((f_x_t_star_fore_2[, ik] * alpha_x_2))
  }
  colnames(d_x_t_star_fore_1) = colnames(d_x_t_star_fore_2) = 1:fh
  rownames(d_x_t_star_fore_1) = rownames(d_x_t_star_fore_2) = age_index
  if (modeling_method == "MLFTS") {
    return(list(R2_1 = R2_1, R2_2 = R2_2, recon_1 = d_x_t_star_recon_1 * 
                  10^5, recon_2 = d_x_t_star_recon_2 * 10^5, fore_count_1 = d_x_t_star_fore_1 * 
                  10^5, fore_count_2 = d_x_t_star_fore_2 * 10^5, alpha_x_1 = alpha_x_1, 
                alpha_x_2 = alpha_x_2, h_x_t_1 = h_x_t_1, h_x_t_2 = h_x_t_2, 
                basis_ave_svd_value = basis_ave_svd_value, basis_resi_1_svd_value = basis_resi_1_svd_value, 
                basis_resi_2_svd_value = basis_resi_2_svd_value, 
                ncomp_ave = ncomp_ave, ncomp_resi_1 = ncomp_resi_1, 
                ncomp_resi_2 = ncomp_resi_2))
  }
  else {
    return(list(R2_1 = R2_1, R2_2 = R2_2, recon_1 = d_x_t_star_recon_1 * 
                  10^5, recon_2 = d_x_t_star_recon_2 * 10^5, fore_count_1 = d_x_t_star_fore_1 * 
                  10^5, fore_count_2 = d_x_t_star_fore_2 * 10^5, alpha_x_1 = alpha_x_1, 
                alpha_x_2 = alpha_x_2, h_x_t_1 = h_x_t_1, h_x_t_2 = h_x_t_2))
  }
}

