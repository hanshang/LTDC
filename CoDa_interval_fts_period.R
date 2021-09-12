##################
# load R packages
##################

load.packages(c("psych", "ftsa", "tseries", "sandwich"))
require(psych)
require(ftsa)
require(tseries)
require(sandwich)

#################
# interval score
#################

interval_score <- function(holdout, lb, ub, alpha)
{
    lb_ind = ifelse(holdout < lb, 1, 0)
    ub_ind = ifelse(holdout > ub, 1, 0)
    score = (ub - lb) + 2/alpha * ((lb - holdout) * lb_ind + (holdout - ub) * ub_ind)
    cover = 1 - (length(which(lb_ind == 1)) + length(which(ub_ind == 1)))/length(holdout)
    cpd = abs(cover - (1 - alpha))
    return(c(mean(score), cpd))
}

#####################
# Interval forecasts
#####################

CoDa_recon_int <- function(dat_1, dat_2, fore_method = c("ets", "arima", "rwf", "rw"), fh,  
                           no_boot, alpha)
{
    n_year = nrow(dat_1)
    n_age = ncol(dat_1)
    year_index = rownames(dat_1)
    age_index = colnames(dat_1)
    
    # standardize life table death to sum to 1
    
    dat_1_center = sweep(dat_1, 1, apply(dat_1, 1, sum), "/")
    dat_2_center = sweep(dat_2, 1, apply(dat_2, 1, sum), "/")
    
    alpha_x_1 = alpha_x_2 = vector("numeric", n_age)
    for(ik in 1:n_age)
    {
        alpha_x_1[ik] = geometric.mean(dat_1_center[,ik])
        alpha_x_2[ik] = geometric.mean(dat_2_center[,ik])
    }
    
    f_x_t_1 = f_x_t_2 = matrix(NA, n_year, n_age)
    for(ik in 1:n_year)
    {
        f_x_t_1[ik,] = (dat_1[ik,]/alpha_x_1)/sum(dat_1[ik,]/alpha_x_1)
        f_x_t_2[ik,] = (dat_2[ik,]/alpha_x_2)/sum(dat_2[ik,]/alpha_x_2)
    }
    
    g_t_1 = g_t_2 = vector("numeric", n_year)
    h_x_t_1 = h_x_t_2 = matrix(NA, n_year, n_age)
    for(ik in 1:n_year)
    {
        g_t_1[ik] = geometric.mean(f_x_t_1[ik,])
        h_x_t_1[ik,] = log(f_x_t_1[ik,]/g_t_1[ik])
        
        g_t_2[ik] = geometric.mean(f_x_t_2[ik,])
        h_x_t_2[ik,] = log(f_x_t_2[ik,]/g_t_2[ik])
    }
    colnames(h_x_t_1) = colnames(h_x_t_2) = age_index
    rownames(h_x_t_1) = rownames(h_x_t_2) = year_index
 
    dum_1 = ER_GR(data = h_x_t_1)
    dum_2 = ER_GR(data = h_x_t_2)
    ncomp_1 = max(dum_1$k_ER, dum_1$k_GR)
    ncomp_2 = max(dum_2$k_ER, dum_2$k_GR)
    rm(dum_1); rm(dum_2)
    
    SVD_decomp_1 = svd(h_x_t_1) # 1st component only
    SVD_decomp_2 = svd(h_x_t_2)
    
    basis_1 = SVD_decomp_1$v[,1:ncomp_1]
    basis_2 = SVD_decomp_2$v[,1:ncomp_2]
    
    score_1 = t(basis_1) %*% t(h_x_t_1)
    score_2 = t(basis_2) %*% t(h_x_t_2)
    
    recon_1 = as.matrix(basis_1) %*% score_1
    recon_2 = as.matrix(basis_2) %*% score_2
    colnames(score_1) = colnames(score_2) = year_index
    colnames(recon_1) = colnames(recon_2) = year_index
    
    resi_1 = t(h_x_t_1) - recon_1
    resi_2 = t(h_x_t_2) - recon_2
    
    # determine in-sample forecast error
    
    olivia_1 = matrix(NA, ncomp_1, fh)
    olivia_2 = matrix(NA, ncomp_2, fh)
    if(fore_method == "ets")
    {
        for(ij in 1:ncomp_1)
        {
            olivia_1[ij,] = forecast(ets(score_1[ij,]), h = fh)$mean
        }
        for(ij in 1:ncomp_2)
        {
            olivia_2[ij,] = forecast(ets(score_2[ij,]), h = fh)$mean
        }
    }
    else if(fore_method == "arima")
    {
        for(ij in 1:ncomp_1)
        {
            olivia_1[ij,] = forecast(auto.arima(score_1[ij,]), h = fh)$mean
        }
        for(ij in 1:ncomp_2)
        {
            olivia_2[ij,] = forecast(auto.arima(score_2[ij,]), h = fh)$mean
        }
    }
    else if(fore_method == "rwf")
    {
        for(ij in 1:ncomp_1)
        {
            olivia_1[ij,] = rwf(score_1[ij,], h = fh, drift = TRUE)$mean
        }
        for(ij in 1:ncomp_2)
        {
            olivia_2[ij,] = rwf(score_2[ij,], h = fh, drift = TRUE)$mean
        }
    }
    else if(fore_method == "rw")
    {
        for(ij in 1:ncomp_1)
        {
            olivia_1[ij,] = rwf(score_1[ij,], h = fh, drift = FALSE)$mean
        }
        for(ij in 1:ncomp_2)
        {
            olivia_2[ij,] = rwf(score_2[ij,], h = fh, drift = FALSE)$mean
        }
    }
    else
    {
        warning("Forecasting method is not on the list.")
    }
    forerr_1 = matrix(NA, (n_year - ncomp_1 - fh + 1), ncomp_1)
    for(i in fh:(n_year - ncomp_1))
    {
        k = i + (ncomp_1 - fh)
        fore = matrix(NA, 1, ncomp_1)
        if(fore_method == "ets")
        {
            for(j in 1:ncomp_1)
            {
                fore[,j] = forecast(ets(score_1[j,1:k]), h = fh)$mean[fh]
            }
        }
        else if(fore_method == "arima")
        {
            for(j in 1:ncomp_1)
            {
                fore[,j] = forecast(auto.arima(score_1[j,1:k]), h = fh)$mean[fh]
            }
        }
        else if(fore_method == "rwf")
        {
            if(k <= 2)
            {
                for(j in 1:ncomp_1)
                {
                    fore[,j] = score_1[j,k]
                }
            }
            if(k > 2)
            {
                for(j in 1:ncomp_1)
                {
                    fore[,j] = rwf(score_1[j,1:k], h = fh, drift = TRUE)$mean[fh]
                }
            }
        }
        else if(fore_method == "rw")
        {
            if(k == 1)
            {
                for(j in 1:ncomp_1)
                {
                    fore[,j] = score_1[j,1]
                }
            }
            if(k > 1)
            {
                for(j in 1:ncomp_1)
                {
                    fore[,j] = rwf(score_1[j,1:k], h = fh, drift = FALSE)$mean[fh]
                }
            }
        }
        else
        {
            warning("Forecasting method is not on the list.")
        }
        forerr_1[i - fh + 1,] = score_1[, k + fh] - fore
    }
    
    forerr_2 = matrix(NA, (n_year - ncomp_2 - fh + 1), ncomp_2)
    for(i in fh:(n_year - ncomp_2))
    {
        k = i + (ncomp_2 - fh)
        fore = matrix(NA, 1, ncomp_2)
        if(fore_method == "ets")
        {
            for(j in 1:ncomp_2)
            {
                fore[,j] = forecast(ets(score_2[j,1:k]), h = fh)$mean[fh]
            }
        }
        else if(fore_method == "arima")
        {
            for(j in 1:ncomp_2)
            {
                fore[,j] = forecast(auto.arima(score_2[j,1:k]), h = fh)$mean[fh]
            }
        }
        else if(fore_method == "rwf")
        {
            if(k <= 2)
            {
                for(j in 1:ncomp_2)
                {
                    fore[,j] = score_2[j,k]
                }
            }
            if(k > 2)
            {
                for(j in 1:ncomp_2)
                {
                    fore[,j] = rwf(score_2[j,1:k], h = fh, drift = TRUE)$mean[fh]
                }
            }
        }
        else if(fore_method == "rw")
        {
            if(k == 1)
            {
                for(j in 1:ncomp_2)
                {
                    fore[,j] = score_2[j,1]
                }
            }
            if(k > 1)
            {
                for(j in 1:ncomp_2)
                {
                    fore[,j] = rwf(score_2[j,1:k], h = fh, drift = FALSE)$mean[fh]
                }
            }
        }
        else
        {
            warning("Forecasting method is not on the list.")
        }
        forerr_2[i - fh + 1,] = score_2[, k + fh] - fore
    }
    
    # bootstrapping residuals
    
    q_1 = q_2 = array(NA, dim = c(n_age, no_boot, no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:n_age)
        {
            for(k in 1:no_boot)
            {
                q_1[i,,k,j] = sample(resi_1[i,], size = no_boot, replace = TRUE)
                q_2[i,,k,j] = sample(resi_2[i,], size = no_boot, replace = TRUE)
            }
        }
    }
    rm(i); rm(j); rm(k)
    
    # bootstrapping PC score errors
    
    ny_1 = array(NA, dim = c(ncomp_1, no_boot, fh))
    ny_2 = array(NA, dim = c(ncomp_2, no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:ncomp_1)
        {
            ny_1[i,,j] = sample(forerr_1[,i], size = no_boot, replace = TRUE)
        }
        rm(i)
        for(i in 1:ncomp_2)
        {
            ny_2[i,,j] = sample(forerr_2[,i], size = no_boot, replace = TRUE)
        }
        rm(i)
    }
    rm(j)
    
    # adding the PC score error to the predicted score
    
    oli_1 = array(rep(olivia_1, no_boot * fh), dim = c(ncomp_1, no_boot, fh))
    fo_1 = array(NA, dim = c(ncomp_1, no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:no_boot)
        {
            fo_1[,i,j] = oli_1[,i,j] + ny_1[,i,j]
        }
    }
    rm(i); rm(j)
    
    oli_2 = array(rep(olivia_2, no_boot * fh), dim = c(ncomp_2, no_boot, fh))
    fo_2 = array(NA, dim = c(ncomp_2, no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:no_boot)
        {
            fo_2[,i,j] = oli_2[,i,j] + ny_2[,i,j]
        }
    }
    rm(i); rm(j)
    
    # construct bootstrapped samples
    
    pred_1 = pred_2 = array(NA, dim = c(n_age, no_boot, no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:no_boot)
        {
            for(k in 1:no_boot)
            {
                pred_1[,i,k,j] = as.matrix(basis_1) %*% fo_1[,i,j] + q_1[,i,k,j]
                pred_2[,i,k,j] = as.matrix(basis_2) %*% fo_2[,i,j] + q_2[,i,k,j]
            }
        }
    }
    rm(i); rm(j); rm(k)
    
    pred_1_resize = pred_2_resize = array(NA, dim = c(n_age, no_boot * no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:no_boot)
        {
            pred_1_resize[, (((i-1)*no_boot+1):(i*no_boot)), ] = pred_1[,i,,j]
            pred_2_resize[, (((i-1)*no_boot+1):(i*no_boot)), ] = pred_2[,i,,j]
        }
    }
    rm(i); rm(j)
    
    # transform back
    
    f_x_t_star_fore_1 = d_x_t_star_fore_1 = 
    f_x_t_star_fore_2 = d_x_t_star_fore_2 = array(NA, dim = c(n_age, no_boot * no_boot, fh))
    for(iw in 1:fh)
    {
        for(ij in 1:(no_boot * no_boot))
        {
            f_x_t_star_fore_1[,ij,iw] = exp(pred_1_resize[,ij,iw])/sum(exp(pred_1_resize[,ij,iw]))
            d_x_t_star_fore_1[,ij,iw] = (f_x_t_star_fore_1[,ij,iw] * alpha_x_1)/sum((f_x_t_star_fore_1[,ij,iw] * alpha_x_1))

            f_x_t_star_fore_2[,ij,iw] = exp(pred_2_resize[,ij,iw])/sum(exp(pred_2_resize[,ij,iw]))
            d_x_t_star_fore_2[,ij,iw] = (f_x_t_star_fore_2[,ij,iw] * alpha_x_2)/sum((f_x_t_star_fore_2[,ij,iw] * alpha_x_2))
        }
    }
    rm(iw); rm(ij)
    d_boot_1 = d_x_t_star_fore_1 * (10^5)
    d_boot_2 = d_x_t_star_fore_2 * (10^5)
    return(list(PI_1 = apply(d_boot_1, c(1, 3), quantile, c(alpha/2, 1-alpha/2)),
                PI_2 = apply(d_boot_2, c(1, 3), quantile, c(alpha/2, 1-alpha/2))))
}

##############
# h = 1 to 30
##############

CoDa_int <- function(fh, method_fore, number_boot, sig)
{
    fore_h30_1 = fore_h30_2 = array(NA, dim = c(2, nrow(EW_female_pop), (31 - fh)))
    for(ik in 1:(31 - fh))
    {
        dum = CoDa_recon_int(dat_1 = t(EW_female_pop)[1:(147+ik),], 
                             dat_2 = t(EW_male_pop)[1:(147+ik),], 
                             fore_method = method_fore,
                             fh = fh, no_boot = number_boot, alpha = sig)
        fore_h30_1[,,ik] = dum$PI_1[,,fh]  
        fore_h30_2[,,ik] = dum$PI_2[,,fh]  
        rm(ik); rm(dum)
    }
    return(list(fore_h30_1 = fore_h30_1, fore_h30_2 = fore_h30_2))
}

####################
# CoDa (traditional)
#####################

## ARIMA

CoDa_int_array_ARIMA_80 = list()
for(iwk in 1:30)
{
    CoDa_int_array_ARIMA_80[[iwk]] = CoDa_int(fh = iwk, method_fore = "arima", number_boot = 100, sig = 0.2) 
    print(iwk); rm(iwk)
}


CoDa_int_array_ARIMA_95 = list()
for(iwk in 1:30)
{
    CoDa_int_array_ARIMA_95[[iwk]] = CoDa_int(fh = iwk, method_fore = "arima", number_boot = 100, sig = 0.05)
    print(iwk); rm(iwk)
}





# compute mean interval score (80%) for females and males

CoDa_ECP_female_ARIMA_80 = CoDa_ECP_male_ARIMA_80 = vector("numeric", 30)
CoDa_mean_score_female_ARIMA_80 = CoDa_mean_score_male_ARIMA_80 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    CoDa_ECP_female_ARIMA_80[ik] = 1 - (length(which(EW_female_pop[,(148+ik):178] < (CoDa_int_array_ARIMA_80[[ik]])$fore_h30_1[1,,])) + 
                                        length(which(EW_female_pop[,(148+ik):178] > (CoDa_int_array_ARIMA_80[[ik]])$fore_h30_1[2,,])))/length(EW_female_pop[,(148+ik):178])

    CoDa_ECP_male_ARIMA_80[ik] = 1 - (length(which(EW_male_pop[,(148+ik):178] < (CoDa_int_array_ARIMA_80[[ik]])$fore_h30_1[1,,])) + 
                                      length(which(EW_male_pop[,(148+ik):178] > (CoDa_int_array_ARIMA_80[[ik]])$fore_h30_1[2,,])))/length(EW_male_pop[,(148+ik):178])
    
    CoDa_mean_score_female_ARIMA_80[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                 lb = (CoDa_int_array_ARIMA_80[[ik]])$fore_h30_1[1,,],
                                                 ub = (CoDa_int_array_ARIMA_80[[ik]])$fore_h30_1[2,,], alpha = 0.2)    
    
    CoDa_mean_score_male_ARIMA_80[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                     lb = (CoDa_int_array_ARIMA_80[[ik]])$fore_h30_2[1,,],
                                                     ub = (CoDa_int_array_ARIMA_80[[ik]])$fore_h30_2[2,,], alpha = 0.2)    
}

# compute mean interval score (95%)

CoDa_ECP_female_ARIMA_95 = CoDa_ECP_male_ARIMA_95 = vector("numeric", 30)
CoDa_mean_score_female_ARIMA_95 = CoDa_mean_score_male_ARIMA_95 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    CoDa_ECP_female_ARIMA_95[ik] = 1 - (length(which(EW_female_pop[,(148+ik):178] < (CoDa_int_array_ARIMA_95[[ik]])$fore_h30_1[1,,])) + 
                                            length(which(EW_female_pop[,(148+ik):178] > (CoDa_int_array_ARIMA_95[[ik]])$fore_h30_1[2,,])))/length(EW_female_pop[,(148+ik):178])
    
    CoDa_ECP_male_ARIMA_95[ik] = 1 - (length(which(EW_male_pop[,(148+ik):178] < (CoDa_int_array_ARIMA_95[[ik]])$fore_h30_1[1,,])) + 
                                          length(which(EW_male_pop[,(148+ik):178] > (CoDa_int_array_ARIMA_95[[ik]])$fore_h30_1[2,,])))/length(EW_male_pop[,(148+ik):178])
    
    CoDa_mean_score_female_ARIMA_95[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                 lb = (CoDa_int_array_ARIMA_95[[ik]])$fore_h30_1[1,,],
                                                 ub = (CoDa_int_array_ARIMA_95[[ik]])$fore_h30_1[2,,], alpha = 0.05)    
    
    CoDa_mean_score_male_ARIMA_95[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                lb = (CoDa_int_array_ARIMA_95[[ik]])$fore_h30_2[1,,],
                                                ub = (CoDa_int_array_ARIMA_95[[ik]])$fore_h30_2[2,,], alpha = 0.05)
}
colnames(CoDa_mean_score_female_ARIMA_80) = colnames(CoDa_mean_score_male_ARIMA_80) = 
colnames(CoDa_mean_score_female_ARIMA_95) = colnames(CoDa_mean_score_male_ARIMA_95) = c("Mean interval score", "CPD")

rownames(CoDa_mean_score_female_ARIMA_80) = rownames(CoDa_mean_score_male_ARIMA_80) = 
rownames(CoDa_mean_score_female_ARIMA_95) = rownames(CoDa_mean_score_male_ARIMA_95) = 1:30
    
round(colMeans(CoDa_mean_score_female_ARIMA_80), 4) #  829.7658 0.2877 
round(colMeans(CoDa_mean_score_female_ARIMA_95), 4) # 1267.7368 0.1637
round(colMeans(CoDa_mean_score_male_ARIMA_80), 4)   # 1830.2730 0.3388
round(colMeans(CoDa_mean_score_male_ARIMA_95), 4)   # 3954.0509 0.2413

## ETS

CoDa_int_array_ETS_80 = list()
for(iwk in 1:30)
{
    CoDa_int_array_ETS_80[[iwk]] = CoDa_int(fh = iwk, method_fore = "ets", number_boot = 100, sig = 0.2)
    print(iwk); rm(iwk)
}

CoDa_int_array_ETS_95 = list()
for(iwk in 1:30)
{
    CoDa_int_array_ETS_95[[iwk]] = CoDa_int(fh = iwk, method_fore = "ets", number_boot = 100, sig = 0.05)
    print(iwk); rm(iwk)   
}

# compute mean interval score

CoDa_mean_score_female_ETS_80 = CoDa_mean_score_male_ETS_80 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    CoDa_mean_score_female_ETS_80[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                        lb = (CoDa_int_array_ETS_80[[ik]])$fore_h30_1[1,,],
                                                        ub = (CoDa_int_array_ETS_80[[ik]])$fore_h30_1[2,,], alpha = 0.2)    
    
    CoDa_mean_score_male_ETS_80[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                      lb = (CoDa_int_array_ETS_80[[ik]])$fore_h30_2[1,,],
                                                      ub = (CoDa_int_array_ETS_80[[ik]])$fore_h30_2[2,,], alpha = 0.2)    
}

CoDa_mean_score_female_ETS_95 = CoDa_mean_score_male_ETS_95 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    CoDa_mean_score_female_ETS_95[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                        lb = (CoDa_int_array_ETS_95[[ik]])$fore_h30_1[1,,],
                                                        ub = (CoDa_int_array_ETS_95[[ik]])$fore_h30_1[2,,], alpha = 0.05)    
    
    CoDa_mean_score_male_ETS_95[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                      lb = (CoDa_int_array_ETS_95[[ik]])$fore_h30_2[1,,],
                                                      ub = (CoDa_int_array_ETS_95[[ik]])$fore_h30_2[2,,], alpha = 0.05)
}

colnames(CoDa_mean_score_female_ETS_80) = colnames(CoDa_mean_score_male_ETS_80) = 
colnames(CoDa_mean_score_female_ETS_95) = colnames(CoDa_mean_score_male_ETS_95) = c("Mean interval score", "CPD")
rownames(CoDa_mean_score_female_ETS_80) = rownames(CoDa_mean_score_male_ETS_80) = 
rownames(CoDa_mean_score_female_ETS_95) = rownames(CoDa_mean_score_male_ETS_95) = 1:30

round(colMeans(CoDa_mean_score_female_ETS_80), 4) # 718.4675 0.2466
round(colMeans(CoDa_mean_score_female_ETS_95), 4) # 931.7810 0.1067
round(colMeans(CoDa_mean_score_male_ETS_80), 4)   # 1958.6489 0.2917
round(colMeans(CoDa_mean_score_male_ETS_95), 4)   # 4437.7861 0.2336


## RWD

CoDa_int_array_RWD_80 = list()
for(iwk in 1:20)
{
    CoDa_int_array_RWD_80[[iwk]] = CoDa_int(fh = iwk, method_fore = "rwf", number_boot = 100, sig = 0.2)
    print(iwk); rm(iwk)
}

CoDa_int_array_RWD_95 = list()
for(iwk in 1:30)
{
    CoDa_int_array_RWD_95[[iwk]] = CoDa_int(fh = iwk, method_fore = "rwf", number_boot = 100, sig = 0.05)
    print(iwk); rm(iwk)   
}

# compute mean interval score

CoDa_mean_score_female_RWD_80 = CoDa_mean_score_male_RWD_80 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    CoDa_mean_score_female_RWD_80[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                        lb = (CoDa_int_array_RWD_80[[ik]])$fore_h30_1[1,,],
                                                        ub = (CoDa_int_array_RWD_80[[ik]])$fore_h30_1[2,,], alpha = 0.2)    
    
    CoDa_mean_score_male_RWD_80[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                      lb = (CoDa_int_array_RWD_80[[ik]])$fore_h30_2[1,,],
                                                      ub = (CoDa_int_array_RWD_80[[ik]])$fore_h30_2[2,,], alpha = 0.2)    
}

CoDa_mean_score_female_RWD_95 = CoDa_mean_score_male_RWD_95 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    CoDa_mean_score_female_RWD_95[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                        lb = (CoDa_int_array_RWD_95[[ik]])$fore_h30_1[1,,],
                                                        ub = (CoDa_int_array_RWD_95[[ik]])$fore_h30_1[2,,], alpha = 0.05)    
    
    CoDa_mean_score_male_RWD_95[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                      lb = (CoDa_int_array_RWD_95[[ik]])$fore_h30_2[1,,],
                                                      ub = (CoDa_int_array_RWD_95[[ik]])$fore_h30_2[2,,], alpha = 0.05)
}

colnames(CoDa_mean_score_female_RWD_80) = colnames(CoDa_mean_score_male_RWD_80) = 
colnames(CoDa_mean_score_female_RWD_95) = colnames(CoDa_mean_score_male_RWD_95) = c("Mean interval score", "CPD")
rownames(CoDa_mean_score_female_RWD_80) = rownames(CoDa_mean_score_male_RWD_80) = 
rownames(CoDa_mean_score_female_RWD_95) = rownames(CoDa_mean_score_male_RWD_95) = 1:30

round(colMeans(CoDa_mean_score_female_RWD_80), 4) # 797.0984 0.1633
round(colMeans(CoDa_mean_score_female_RWD_95), 4) # 994.4774 0.0341
round(colMeans(CoDa_mean_score_male_RWD_80), 4)   # 1781.9027 0.3077
round(colMeans(CoDa_mean_score_male_RWD_95), 4)   # 3393.8619 0.2013


## RW

CoDa_int_array_RW_80 = list()
for(iwk in 1:30)
{
    CoDa_int_array_RW_80[[iwk]] = CoDa_int(fh = iwk, method_fore = "rw", number_boot = 100, sig = 0.2)
    print(iwk); rm(iwk)
}

CoDa_int_array_RW_95 = list()
for(iwk in 1:30)
{
    CoDa_int_array_RW_95[[iwk]] = CoDa_int(fh = iwk, method_fore = "rw", number_boot = 100, sig = 0.05)
    print(iwk); rm(iwk)   
}

# compute mean interval score

CoDa_mean_score_female_RW_80 = CoDa_mean_score_male_RW_80 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    CoDa_mean_score_female_RW_80[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                       lb = (CoDa_int_array_RW_80[[ik]])$fore_h30_1[1,,],
                                                       ub = (CoDa_int_array_RW_80[[ik]])$fore_h30_1[2,,], alpha = 0.2)    
    
    CoDa_mean_score_male_RW_80[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                     lb = (CoDa_int_array_RW_80[[ik]])$fore_h30_2[1,,],
                                                     ub = (CoDa_int_array_RW_80[[ik]])$fore_h30_2[2,,], alpha = 0.2)    
}

CoDa_mean_score_female_RW_95 = CoDa_mean_score_male_RW_95 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    CoDa_mean_score_female_RW_95[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                       lb = (CoDa_int_array_RW_95[[ik]])$fore_h30_1[1,,],
                                                       ub = (CoDa_int_array_RW_95[[ik]])$fore_h30_1[2,,], alpha = 0.05)    
    
    CoDa_mean_score_male_RW_95[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                     lb = (CoDa_int_array_RW_95[[ik]])$fore_h30_2[1,,],
                                                     ub = (CoDa_int_array_RW_95[[ik]])$fore_h30_2[2,,], alpha = 0.05)
}

colnames(CoDa_mean_score_female_RW_80) = colnames(CoDa_mean_score_male_RW_80) = 
colnames(CoDa_mean_score_female_RW_95) = colnames(CoDa_mean_score_male_RW_95) = c("Mean interval score", "CPD")
rownames(CoDa_mean_score_female_RW_80) = rownames(CoDa_mean_score_male_RW_80) = 
rownames(CoDa_mean_score_female_RW_95) = rownames(CoDa_mean_score_male_RW_95) = 1:30

round(colMeans(CoDa_mean_score_female_RW_80), 4) # 926.1453  0.1183 
round(colMeans(CoDa_mean_score_female_RW_95), 4) # 1437.3608 0.0703 
round(colMeans(CoDa_mean_score_male_RW_80), 4)   # 1942.4707 0.2836
round(colMeans(CoDa_mean_score_male_RW_95), 4)   # 4306.5804 0.2294

