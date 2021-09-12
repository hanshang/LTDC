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

#############################
# Interval forecasts (MLFTS)
#############################

CoDa_recon_int_mlfts <- function(dat_1, dat_2, fore_method = c("ets", "arima", "rwf", "rw"), 
                                 fh, no_boot = 100, alpha)
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
    
    # construct average of two populations
    
    h_x_t_ave = (h_x_t_1 + h_x_t_2)/2
    dum = ER_GR(h_x_t_ave)
    ncomp_ave = max(dum$k_ER, dum$k_GR)
    rm(dum)
    
    SVD_decomp_ave = svd(h_x_t_ave)
    basis_ave = SVD_decomp_ave$v[,1:ncomp_ave]
    score_ave = matrix(t(basis_ave) %*% t(h_x_t_ave), ncomp_ave, )
    recon_ave = basis_ave %*% score_ave
    resi_ave = t(h_x_t_ave) - recon_ave
    
    # obtain sex-specific residuals
    
    h_x_t_resi_1 = h_x_t_1 - t(recon_ave)
    h_x_t_resi_2 = h_x_t_2 - t(recon_ave)
    
    dum_resi_1 = ER_GR(h_x_t_resi_1)
    ncomp_resi_1 = max(dum_resi_1$k_ER, dum_resi_1$k_GR)
    
    dum_resi_2 = ER_GR(h_x_t_resi_2)
    ncomp_resi_2 = max(dum_resi_2$k_ER, dum_resi_2$k_GR)
    
    # 1st population
    
    SVD_decomp_resi_1 = svd(h_x_t_resi_1)
    basis_resi_1 = SVD_decomp_resi_1$v[,1:ncomp_resi_1]
    score_resi_1 = matrix(t(basis_resi_1) %*% t(h_x_t_resi_1), ncomp_resi_1, )
    recon_resi_1 = basis_resi_1 %*% score_resi_1
    resi_resi_1  = t(h_x_t_resi_1) - recon_resi_1
    
    # 2nd population
    
    SVD_decomp_resi_2 = svd(h_x_t_resi_2)
    basis_resi_2 = SVD_decomp_resi_2$v[,1:ncomp_resi_2]
    score_resi_2 = matrix(t(basis_resi_2) %*% t(h_x_t_resi_2), ncomp_resi_2, )
    recon_resi_2 = basis_resi_2 %*% score_resi_2
    resi_resi_2  = t(h_x_t_resi_2) - recon_resi_2
    
    # reconstruction of original curves for both populations
    
    recon_1 = recon_ave + recon_resi_1
    recon_2 = recon_ave + recon_resi_2
    
    ##################
    # forecast scores
    ##################
    
    # ave
    
    olivia_ave = matrix(NA, ncomp_ave, fh)
    if(fore_method == "ets")
    {
        for(ij in 1:ncomp_ave)
        {
            olivia_ave[ij,] = forecast(ets(score_ave[ij,]), h = fh)$mean
        }
    }
    else if(fore_method == "arima")
    {
        for(ij in 1:ncomp_ave)
        {
            olivia_ave[ij,] = forecast(auto.arima(score_ave[ij,]), h = fh)$mean
        }
    }
    else if(fore_method == "rwf")
    {
        for(ij in 1:ncomp_ave)
        {
            olivia_ave[ij,] = rwf(score_ave[ij,], h = fh, drift = TRUE)$mean
        }
    }
    else if(fore_method == "rw")
    {
        for(ij in 1:ncomp_ave)
        {
            olivia_ave[ij,] = rwf(score_ave[ij,], h = fh, drift = FALSE)$mean
        }
    }
    else
    {
        warning("Forecasting method is not on the list.")
    }
    rm(ij)
    
    # resi_1
    
    olivia_resi_1 = matrix(NA, ncomp_resi_1, fh)
    if(fore_method == "ets")
    {
        for(ij in 1:ncomp_resi_1)
        {
            olivia_resi_1[ij,] = forecast(ets(score_resi_1[ij,]), h = fh)$mean
        }
    }
    else if(fore_method == "arima")
    {
        for(ij in 1:ncomp_resi_1)
        {
            olivia_resi_1[ij,] = forecast(auto.arima(score_resi_1[ij,]), h = fh)$mean
        }
    }
    else if(fore_method == "rwf")
    {
        for(ij in 1:ncomp_resi_1)
        {
            olivia_resi_1[ij,] = rwf(score_resi_1[ij,], h = fh, drift = TRUE)$mean
        }
    }
    else if(fore_method == "rw")
    {
        for(ij in 1:ncomp_resi_1)
        {
            olivia_resi_1[ij,] = rwf(score_resi_1[ij,], h = fh, drift = FALSE)$mean
        }
    }
    else
    {
        warning("Forecasting method is not on the list.")
    }
    rm(ij)
    
    # resi_2
    
    olivia_resi_2 = matrix(NA, ncomp_resi_2, fh)
    if(fore_method == "ets")
    {
        for(ij in 1:ncomp_resi_2)
        {
            olivia_resi_2[ij,] = forecast(ets(score_resi_2[ij,]), h = fh)$mean
        }
    }
    else if(fore_method == "arima")
    {
        for(ij in 1:ncomp_resi_2)
        {
            olivia_resi_2[ij,] = forecast(auto.arima(score_resi_2[ij,]), h = fh)$mean
        }
    }
    else if(fore_method == "rwf")
    {
        for(ij in 1:ncomp_resi_2)
        {
            olivia_resi_2[ij,] = rwf(score_resi_2[ij,], h = fh, drift = TRUE)$mean
        }
    }
    else if(fore_method == "rw")
    {
        for(ij in 1:ncomp_resi_2)
        {
            olivia_resi_2[ij,] = rwf(score_resi_2[ij,], h = fh, drift = FALSE)$mean
        }
    }
    else
    {
        warning("Forecasting method is not on the list.")
    }
    rm(ij)
    
    ############################
    # determine in-sample error
    ############################
    
    forerr_ave = matrix(NA, (n_year - ncomp_ave - fh + 1), ncomp_ave)
    for(i in fh:(n_year - ncomp_ave))
    {
        k = i + (ncomp_ave - fh)
        fore = matrix(NA, 1, ncomp_ave)
        if(fore_method == "ets")
        {
            for(j in 1:ncomp_ave)
            {
                fore[,j] = forecast(ets(score_ave[j,1:k]), h = fh)$mean[fh]
            }
        }
        else if(fore_method == "arima")
        {
            for(j in 1:ncomp_ave)
            {
                fore[,j] = forecast(auto.arima(score_ave[j,1:k]), h = fh)$mean[fh]
            }
        }
        else if(fore_method == "rwf")
        {
            if(k <= 2)
            {
                for(j in 1:ncomp_ave)
                {
                    fore[,j] = score_ave[j,k]
                }
            }
            if(k > 2)
            {
                for(j in 1:ncomp_ave)
                {
                    fore[,j] = rwf(score_ave[j,1:k], h = fh, drift = TRUE)$mean[fh]
                }
            }
        }
        else if(fore_method == "rw")
        {
            if(k == 1)
            {
                for(j in 1:ncomp_ave)
                {
                    fore[,j] = score_ave[j,1]
                }
            }
            if(k > 1)
            {
                for(j in 1:ncomp_ave)
                {
                    fore[,j] = rwf(score_ave[j,1:k], h = fh, drift = FALSE)$mean[fh]
                }
            }
        }
        else
        {
            warning("Forecasting method is not on the list.")
        }
        forerr_ave[i - fh + 1,] = score_ave[, k + fh] - fore
    }
    rm(fore)
    
    # resi_1
    
    forerr_resi_1 = matrix(NA, (n_year - ncomp_resi_1 - fh + 1), ncomp_resi_1)
    for(i in fh:(n_year - ncomp_resi_1))
    {
        k = i + (ncomp_resi_1 - fh)
        fore = matrix(NA, 1, ncomp_resi_1)
        if(fore_method == "ets")
        {
            for(j in 1:ncomp_resi_1)
            {
                fore[,j] = forecast(ets(score_resi_1[j,1:k]), h = fh)$mean[fh]
            }
        }
        else if(fore_method == "arima")
        {
            for(j in 1:ncomp_resi_1)
            {
                fore[,j] = forecast(auto.arima(score_resi_1[j,1:k]), h = fh)$mean[fh]
            }
        }
        else if(fore_method == "rwf")
        {
            if(k <= 2)
            {
                for(j in 1:ncomp_resi_1)
                {
                    fore[,j] = score_resi_1[j,k]
                }
            }
            if(k > 2)
            {
                for(j in 1:ncomp_resi_1)
                {
                    fore[,j] = rwf(score_resi_1[j,1:k], h = fh, drift = TRUE)$mean[fh]
                }
            }
        }
        else if(fore_method == "rw")
        {
            if(k == 1)
            {
                for(j in 1:ncomp_resi_1)
                {
                    fore[,j] = score_resi_1[j,1]
                }
            }
            if(k > 1)
            {
                for(j in 1:ncomp_resi_1)
                {
                    fore[,j] = rwf(score_resi_1[j,1:k], h = fh, drift = FALSE)$mean[fh]
                }
            }
        }
        else
        {
            warning("Forecasting method is not on the list.")
        }
        forerr_resi_1[i - fh + 1,] = score_resi_1[, k + fh] - fore
    }
    rm(fore)
    
    # resi_2
    
    forerr_resi_2 = matrix(NA, (n_year - ncomp_resi_2 - fh + 1), ncomp_resi_2)
    for(i in fh:(n_year - ncomp_resi_2))
    {
        k = i + (ncomp_resi_2 - fh)
        fore = matrix(NA, 1, ncomp_resi_2)
        if(fore_method == "ets")
        {
            for(j in 1:ncomp_resi_2)
            {
                fore[,j] = forecast(ets(score_resi_2[j,1:k]), h = fh)$mean[fh]
            }
        }
        else if(fore_method == "arima")
        {
            for(j in 1:ncomp_resi_2)
            {
                fore[,j] = forecast(auto.arima(score_resi_2[j,1:k]), h = fh)$mean[fh]
            }
        }
        else if(fore_method == "rwf")
        {
            if(k <= 2)
            {
                for(j in 1:ncomp_resi_2)
                {
                    fore[,j] = score_resi_2[j,k]
                }
            }
            if(k > 2)
            {
                for(j in 1:ncomp_resi_2)
                {
                    fore[,j] = rwf(score_resi_2[j,1:k], h = fh, drift = TRUE)$mean[fh]
                }
            }
        }
        else if(fore_method == "rw")
        {
            if(k == 1)
            {
                for(j in 1:ncomp_resi_2)
                {
                    fore[,j] = score_resi_2[j,1]
                }
            }
            if(k > 1)
            {
                for(j in 1:ncomp_resi_2)
                {
                    fore[,j] = rwf(score_resi_2[j,1:k], h = fh, drift = FALSE)$mean[fh]
                }
            }
        }
        else
        {
            warning("Forecasting method is not on the list.")
        }
        forerr_resi_2[i - fh + 1,] = score_resi_2[, k + fh] - fore
    }
    rm(fore)

    ##########################
    # bootstrapping residuals
    ##########################
    
    # ave
    
    q_ave = array(NA, dim = c(n_age, no_boot, no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:n_age)
        {
            for(k in 1:no_boot)
            {
                q_ave[i,,k,j] = sample(resi_ave[i,], size = no_boot, replace = TRUE)
            }
        }
    }
    rm(i); rm(j); rm(k)
    
    # resi_1
    
    q_resi_1 = array(NA, dim = c(n_age, no_boot, no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:n_age)
        {
            for(k in 1:no_boot)
            {
                q_resi_1[i,,k,j] = sample(resi_resi_1[i,], size = no_boot, replace = TRUE)
            }
        }
    }
    rm(i); rm(j); rm(k)
    
    # resi_2
    
    q_resi_2 = array(NA, dim = c(n_age, no_boot, no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:n_age)
        {
            for(k in 1:no_boot)
            {
                q_resi_2[i,,k,j] = sample(resi_resi_2[i,], size = no_boot, replace = TRUE)
            }
        }
    }
    rm(i); rm(j); rm(k)
    
    ################################
    # bootstrapping PC score errors
    ################################
    
    ny_ave = array(NA, dim = c(ncomp_ave, no_boot, fh))
    ny_resi_1 = array(NA, dim = c(ncomp_resi_1, no_boot, fh))
    ny_resi_2 = array(NA, dim = c(ncomp_resi_2, no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:ncomp_ave)
        {
            ny_ave[i,,j] = sample(forerr_ave[,i], size = no_boot, replace = TRUE)
        }
        rm(i)
        for(i in 1:ncomp_resi_1)
        {
            ny_resi_1[i,,j] = sample(forerr_resi_1[,i], size = no_boot, replace = TRUE)
        }
        rm(i)
        for(i in 1:ncomp_resi_2)
        {
            ny_resi_2[i,,j] = sample(forerr_resi_2[,i], size = no_boot, replace = TRUE)
        }
        rm(i)
    }
    rm(j)
    
    ###################################################
    # adding the PC score error to the predicted score
    ###################################################
    
    # ave
    
    oli_ave = array(rep(olivia_ave, no_boot * fh), dim = c(ncomp_ave, no_boot, fh))
    fo_ave = array(NA, dim = c(ncomp_ave, no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:no_boot)
        {
            fo_ave[,i,j] = oli_ave[,i,j] + ny_ave[,i,j]
        }
    }
    rm(i); rm(j)
    
    # resi_1
    
    oli_resi_1 = array(rep(olivia_resi_1, no_boot * fh), dim = c(ncomp_resi_1, no_boot, fh))
    fo_resi_1 = array(NA, dim = c(ncomp_resi_1, no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:no_boot)
        {
            fo_resi_1[,i,j] = oli_resi_1[,i,j] + ny_resi_1[,i,j]
        }
    }
    rm(i); rm(j)
    
    # resi_2
    
    oli_resi_2 = array(rep(olivia_resi_2, no_boot * fh), dim = c(ncomp_resi_2, no_boot, fh))
    fo_resi_2 = array(NA, dim = c(ncomp_resi_2, no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:no_boot)
        {
            fo_resi_2[,i,j] = oli_resi_2[,i,j] + ny_resi_2[,i,j]
        }
    }
    rm(i); rm(j)
    
    #################################
    # construct bootstrapped samples
    #################################
    
    pred_ave = pred_resi_1 = pred_resi_2 = array(NA, dim = c(n_age, no_boot, no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:no_boot)
        {
            for(k in 1:no_boot)
            {
                pred_ave[,i,k,j] = as.matrix(basis_ave) %*% fo_ave[,i,j] + q_ave[,i,k,j]
                pred_resi_1[,i,k,j] = as.matrix(basis_resi_1) %*% fo_resi_1[,i,j] + q_resi_1[,i,k,j]
                pred_resi_2[,i,k,j] = as.matrix(basis_resi_2) %*% fo_resi_2[,i,j] + q_resi_2[,i,k,j]
            }
        }
    }
    rm(i); rm(j); rm(k)
    
    pred_ave_resize = pred_resi_1_resize = pred_resi_2_resize = array(NA, dim = c(n_age, no_boot * no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:no_boot)
        {
            pred_ave_resize[, (((i-1)*no_boot+1):(i*no_boot)),] = pred_ave[,i,,j]
            pred_resi_1_resize[,(((i-1)*no_boot+1):(i*no_boot)),] = pred_resi_1[,i,,j]
            pred_resi_2_resize[,(((i-1)*no_boot+1):(i*no_boot)),] = pred_resi_2[,i,,j]
        }
    }
    rm(i); rm(j)
    
    fore_val_1 = pred_ave_resize + pred_resi_1_resize
    fore_val_2 = pred_ave_resize + pred_resi_2_resize
    
    #################
    # transform back
    #################
    
    f_x_t_star_fore_1 = d_x_t_star_fore_1 =
    f_x_t_star_fore_2 = d_x_t_star_fore_2 = array(NA, dim = c(n_age, no_boot * no_boot, fh))
    for(iw in 1:fh)
    {
        for(ij in 1:(no_boot * no_boot))
        {
            f_x_t_star_fore_1[,ij,iw] = exp(fore_val_1[,ij,iw])/sum(exp(fore_val_1[,ij,iw]))
            d_x_t_star_fore_1[,ij,iw] = (f_x_t_star_fore_1[,ij,iw] * alpha_x_1)/sum((f_x_t_star_fore_1[,ij,iw] * alpha_x_1))
            
            f_x_t_star_fore_2[,ij,iw] = exp(fore_val_2[,ij,iw])/sum(exp(fore_val_2[,ij,iw]))
            d_x_t_star_fore_2[,ij,iw] = (f_x_t_star_fore_2[,ij,iw] * alpha_x_2)/sum((f_x_t_star_fore_2[,ij,iw] * alpha_x_2))
        }
    }
    d_boot_1 = d_x_t_star_fore_1 * (10^5)
    d_boot_2 = d_x_t_star_fore_2 * (10^5)
    rm(iw); rm(ij); rm(fore_val_1); rm(fore_val_2); rm(f_x_t_star_fore_1); rm(f_x_t_star_fore_2)
    rm(d_x_t_star_fore_1); rm(d_x_t_star_fore_2)
    return(list(d_boot_1 = d_boot_1, d_boot_2 = d_boot_2,
                PI_1 = apply(d_boot_1, c(1,3), quantile, c(alpha/2, 1-alpha/2)),
                PI_2 = apply(d_boot_2, c(1,3), quantile, c(alpha/2, 1-alpha/2))))
}
    
##############
# h = 1 to 20
##############

CoDa_int_mlfts <- function(fh, method_fore, number_boot, sig)
{
    fore_h30_1 = fore_h30_2 = array(NA, dim = c(2, nrow(EW_female_pop), (31 - fh)))
    for(ik in 1:(31 - fh))
    {
        dum = CoDa_recon_int_mlfts(dat_1 = t(EW_female_pop)[1:(147+ik),],
                                   dat_2 = t(EW_male_pop)[1:(147+ik),],
                                   fore_method = method_fore,
                                   fh = fh, no_boot = number_boot, alpha = sig)
        fore_h30_1[,,ik] = dum$PI_1[,,fh]
        fore_h30_2[,,ik] = dum$PI_2[,,fh]
        rm(ik); rm(dum)
    }
    return(list(fore_h30_1 = fore_h30_1, fore_h30_2 = fore_h30_2))
}
 
#######
# CoDa
#######

## ARIMA

CoDa_int_mlfts_ARIMA_80 = list()
for(iwk in 1:30)
{
    CoDa_int_mlfts_ARIMA_80[[iwk]] = CoDa_int_mlfts(fh = iwk, method_fore = "arima", number_boot = 100, sig = 0.2)
    print(iwk); rm(iwk)
}

CoDa_int_mlfts_ARIMA_95 = list()
for(iwk in 1:30)
{
    CoDa_int_mlfts_ARIMA_95[[iwk]] = CoDa_int_mlfts(fh = iwk, method_fore = "arima", number_boot = 100, sig = 0.05)
    print(iwk); rm(iwk)
}

# compute mean interval score (80%) for females and males

CoDa_mean_score_mlfts_female_ARIMA_80 = CoDa_mean_score_mlfts_male_ARIMA_80 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    CoDa_mean_score_mlfts_female_ARIMA_80[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                          lb = (CoDa_int_mlfts_ARIMA_80[[ik]])$fore_h30_1[1,,],
                                                          ub = (CoDa_int_mlfts_ARIMA_80[[ik]])$fore_h30_1[2,,], alpha = 0.2)    
    
    CoDa_mean_score_mlfts_male_ARIMA_80[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                        lb = (CoDa_int_mlfts_ARIMA_80[[ik]])$fore_h30_2[1,,],
                                                        ub = (CoDa_int_mlfts_ARIMA_80[[ik]])$fore_h30_2[2,,], alpha = 0.2)    
    rm(ik)
}

# compute mean interval score (95%)

CoDa_mean_score_mlfts_female_ARIMA_95 = CoDa_mean_score_mlfts_male_ARIMA_95 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    CoDa_mean_score_mlfts_female_ARIMA_95[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                          lb = (CoDa_int_mlfts_ARIMA_95[[ik]])$fore_h30_1[1,,],
                                                          ub = (CoDa_int_mlfts_ARIMA_95[[ik]])$fore_h30_1[2,,], alpha = 0.05)    
    
    CoDa_mean_score_mlfts_male_ARIMA_95[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                        lb = (CoDa_int_mlfts_ARIMA_95[[ik]])$fore_h30_2[1,,],
                                                        ub = (CoDa_int_mlfts_ARIMA_95[[ik]])$fore_h30_2[2,,], alpha = 0.05)
    rm(ik)
}
colnames(CoDa_mean_score_mlfts_female_ARIMA_80) = colnames(CoDa_mean_score_mlfts_male_ARIMA_80) = 
colnames(CoDa_mean_score_mlfts_female_ARIMA_95) = colnames(CoDa_mean_score_mlfts_male_ARIMA_95) = c("Mean interval score", "CPD")

rownames(CoDa_mean_score_mlfts_female_ARIMA_80) = rownames(CoDa_mean_score_mlfts_male_ARIMA_80) = 
rownames(CoDa_mean_score_mlfts_female_ARIMA_95) = rownames(CoDa_mean_score_mlfts_male_ARIMA_95) = 1:30

round(colMeans(CoDa_mean_score_mlfts_female_ARIMA_80), 4) #  650.2248 0.0939
round(colMeans(CoDa_mean_score_mlfts_female_ARIMA_95), 4) #  866.9047 0.0451
round(colMeans(CoDa_mean_score_mlfts_male_ARIMA_80), 4)   #  941.8934 0.1391
round(colMeans(CoDa_mean_score_mlfts_male_ARIMA_95), 4)   # 1221.6110 0.0450

############################
# Figures in the Appendices
############################

plot(1:30, CoDa_mean_score_male_RWD_80[,1], type = "l", col = 1, lty = 1, xlab = "", ylab = "",
     ylim = c(350, 2600), main = "Male")
lines(1:30, CoDa_mfts_mean_score_male_RWD_80[,1], col = 2, lty = 2)
lines(1:30, CoDa_mean_score_mlfts_male_RWD_80[,1], col = 4, lty = 4)

plot(1:30, CoDa_mean_score_male_RWD_95[,1], type = "l", col = 1, lty = 1, xlab = "Forecast horizon",
     ylab = "", ylim = c(550, 5600), main = "Male")
lines(1:30, CoDa_mfts_mean_score_male_RWD_95[,1], col = 2, lty = 2)
lines(1:30, CoDa_mean_score_mlfts_male_RWD_95[,1], col = 4, lty = 4)

## ETS

CoDa_int_mlfts_ETS_80 = list()
for(iwk in 1:30)
{
    CoDa_int_mlfts_ETS_80[[iwk]] = CoDa_int_mlfts(fh = iwk, method_fore = "ets", number_boot = 100, sig = 0.2)
    print(iwk); rm(iwk)
}

CoDa_int_mlfts_ETS_95 = list()
for(iwk in 1:30)
{
    CoDa_int_mlfts_ETS_95[[iwk]] = CoDa_int_mlfts(fh = iwk, method_fore = "ets", number_boot = 100, sig = 0.05)
    print(iwk); rm(iwk)
}

# compute mean interval score (80%) for females and males

CoDa_mean_score_mlfts_female_ETS_80 = CoDa_mean_score_mlfts_male_ETS_80 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    CoDa_mean_score_mlfts_female_ETS_80[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                              lb = (CoDa_int_mlfts_ETS_80[[ik]])$fore_h30_1[1,,],
                                                              ub = (CoDa_int_mlfts_ETS_80[[ik]])$fore_h30_1[2,,], alpha = 0.2)    
    
    CoDa_mean_score_mlfts_male_ETS_80[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                            lb = (CoDa_int_mlfts_ETS_80[[ik]])$fore_h30_2[1,,],
                                                            ub = (CoDa_int_mlfts_ETS_80[[ik]])$fore_h30_2[2,,], alpha = 0.2)    
}

# compute mean interval score (95%)

CoDa_mean_score_mlfts_female_ETS_95 = CoDa_mean_score_mlfts_male_ETS_95 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    CoDa_mean_score_mlfts_female_ETS_95[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                              lb = (CoDa_int_mlfts_ETS_95[[ik]])$fore_h30_1[1,,],
                                                              ub = (CoDa_int_mlfts_ETS_95[[ik]])$fore_h30_1[2,,], alpha = 0.05)    
    
    CoDa_mean_score_mlfts_male_ETS_95[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                            lb = (CoDa_int_mlfts_ETS_95[[ik]])$fore_h30_2[1,,],
                                                            ub = (CoDa_int_mlfts_ETS_95[[ik]])$fore_h30_2[2,,], alpha = 0.05)
}
colnames(CoDa_mean_score_mlfts_female_ETS_80) = colnames(CoDa_mean_score_mlfts_male_ETS_80) = 
colnames(CoDa_mean_score_mlfts_female_ETS_95) = colnames(CoDa_mean_score_mlfts_male_ETS_95) = c("Mean interval score", "CPD")

rownames(CoDa_mean_score_mlfts_female_ETS_80) = rownames(CoDa_mean_score_mlfts_male_ETS_80) = 
rownames(CoDa_mean_score_mlfts_female_ETS_95) = rownames(CoDa_mean_score_mlfts_male_ETS_95) = 1:30

round(colMeans(CoDa_mean_score_mlfts_female_ETS_80), 4) #  701.0976 0.1392 
round(colMeans(CoDa_mean_score_mlfts_female_ETS_95), 4) #  949.3477 0.0460
round(colMeans(CoDa_mean_score_mlfts_male_ETS_80), 4)   # 1009.0984 0.1192
round(colMeans(CoDa_mean_score_mlfts_male_ETS_95), 4)   # 1330.4392 0.0468

## RWD

CoDa_int_mlfts_RWD_80 = list()
for(iwk in 1:30)
{
    CoDa_int_mlfts_RWD_80[[iwk]] = CoDa_int_mlfts(fh = iwk, method_fore = "rwf", number_boot = 100, sig = 0.2)
    print(iwk); rm(iwk)
}

CoDa_int_mlfts_RWD_95 = list()
for(iwk in 1:30)
{
    CoDa_int_mlfts_RWD_95[[iwk]] = CoDa_int_mlfts(fh = iwk, method_fore = "rwf", number_boot = 100, sig = 0.05)
    print(iwk); rm(iwk)
}

# compute mean interval score (80%) for females and males

CoDa_mean_score_mlfts_female_RWD_80 = CoDa_mean_score_mlfts_male_RWD_80 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    CoDa_mean_score_mlfts_female_RWD_80[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                              lb = (CoDa_int_mlfts_RWD_80[[ik]])$fore_h30_1[1,,],
                                                              ub = (CoDa_int_mlfts_RWD_80[[ik]])$fore_h30_1[2,,], alpha = 0.2)    
    
    CoDa_mean_score_mlfts_male_RWD_80[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                            lb = (CoDa_int_mlfts_RWD_80[[ik]])$fore_h30_2[1,,],
                                                            ub = (CoDa_int_mlfts_RWD_80[[ik]])$fore_h30_2[2,,], alpha = 0.2)    
}

# compute mean interval score (95%)

CoDa_mean_score_mlfts_female_RWD_95 = CoDa_mean_score_mlfts_male_RWD_95 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    CoDa_mean_score_mlfts_female_RWD_95[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                              lb = (CoDa_int_mlfts_RWD_95[[ik]])$fore_h30_1[1,,],
                                                              ub = (CoDa_int_mlfts_RWD_95[[ik]])$fore_h30_1[2,,], alpha = 0.05)    
    
    CoDa_mean_score_mlfts_male_RWD_95[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                            lb = (CoDa_int_mlfts_RWD_95[[ik]])$fore_h30_2[1,,],
                                                            ub = (CoDa_int_mlfts_RWD_95[[ik]])$fore_h30_2[2,,], alpha = 0.05)
}
colnames(CoDa_mean_score_mlfts_female_RWD_80) = colnames(CoDa_mean_score_mlfts_male_RWD_80) = 
colnames(CoDa_mean_score_mlfts_female_RWD_95) = colnames(CoDa_mean_score_mlfts_male_RWD_95) = c("Mean interval score", "CPD")

rownames(CoDa_mean_score_mlfts_female_RWD_80) = rownames(CoDa_mean_score_mlfts_male_RWD_80) = 
rownames(CoDa_mean_score_mlfts_female_RWD_95) = rownames(CoDa_mean_score_mlfts_male_RWD_95) = 1:30

round(colMeans(CoDa_mean_score_mlfts_female_RWD_80), 4) #  655.5298 0.1947 
round(colMeans(CoDa_mean_score_mlfts_female_RWD_95), 4) # 1056.2860 0.050
round(colMeans(CoDa_mean_score_mlfts_male_RWD_80), 4)   #  738.5429 0.1348
round(colMeans(CoDa_mean_score_mlfts_male_RWD_95), 4)   # 1077.238  0.0490

## RW

CoDa_int_mlfts_RW_80 = list()
for(iwk in 1:30)
{
    CoDa_int_mlfts_RW_80[[iwk]] = CoDa_int_mlfts(fh = iwk, method_fore = "rw", number_boot = 100, sig = 0.2)
    print(iwk); rm(iwk)
}

CoDa_int_mlfts_RW_95 = list()
for(iwk in 1:30)
{
    CoDa_int_mlfts_RW_95[[iwk]] = CoDa_int_mlfts(fh = iwk, method_fore = "rw", number_boot = 100, sig = 0.05)
    print(iwk); rm(iwk)
}

# compute mean interval score (80%) for females and males

CoDa_mean_score_mlfts_female_RW_80 = CoDa_mean_score_mlfts_male_RW_80 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    CoDa_mean_score_mlfts_female_RW_80[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                             lb = (CoDa_int_mlfts_RW_80[[ik]])$fore_h30_1[1,,],
                                                             ub = (CoDa_int_mlfts_RW_80[[ik]])$fore_h30_1[2,,], alpha = 0.2)    
    
    CoDa_mean_score_mlfts_male_RW_80[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                           lb = (CoDa_int_mlfts_RW_80[[ik]])$fore_h30_2[1,,],
                                                           ub = (CoDa_int_mlfts_RW_80[[ik]])$fore_h30_2[2,,], alpha = 0.2)    
}

# compute mean interval score (95%)

CoDa_mean_score_mlfts_female_RW_95 = CoDa_mean_score_mlfts_male_RW_95 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    CoDa_mean_score_mlfts_female_RW_95[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                             lb = (CoDa_int_mlfts_RW_95[[ik]])$fore_h30_1[1,,],
                                                             ub = (CoDa_int_mlfts_RW_95[[ik]])$fore_h30_1[2,,], alpha = 0.05)    
    
    CoDa_mean_score_mlfts_male_RW_95[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                           lb = (CoDa_int_mlfts_RW_95[[ik]])$fore_h30_2[1,,],
                                                           ub = (CoDa_int_mlfts_RW_95[[ik]])$fore_h30_2[2,,], alpha = 0.05)
}
colnames(CoDa_mean_score_mlfts_female_RW_80) = colnames(CoDa_mean_score_mlfts_male_RW_80) = 
colnames(CoDa_mean_score_mlfts_female_RW_95) = colnames(CoDa_mean_score_mlfts_male_RW_95) = c("Mean interval score", "CPD")

rownames(CoDa_mean_score_mlfts_female_RW_80) = rownames(CoDa_mean_score_mlfts_male_RW_80) = 
rownames(CoDa_mean_score_mlfts_female_RW_95) = rownames(CoDa_mean_score_mlfts_male_RW_95) = 1:30

round(colMeans(CoDa_mean_score_mlfts_female_RW_80), 4) # 711.6352  0.1410
round(colMeans(CoDa_mean_score_mlfts_female_RW_95), 4) # 918.5135  0.0450
round(colMeans(CoDa_mean_score_mlfts_male_RW_80), 4)   # 1000.4340 0.1220
round(colMeans(CoDa_mean_score_mlfts_male_RW_95), 4)   # 1327.6690 0.0462

