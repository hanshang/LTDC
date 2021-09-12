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

###########################
# Interval forecasts (MFTS)
###########################

# dat_1: first data set (n by p)
# dat_2: second data set (n by p)
# fore_method: univariate time series forecasting method
# fh: forecast horizon
# no_boot: number of bootstrap replications
# alpha: level of significance

CoDa_recon_int_mfts <- function(dat_1, dat_2, fore_method = c("ets", "arima", "rwf", "rw"), fh,
                                no_boot = 100, alpha)
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
    
    data_comb = list()
    data_comb[[1]] = t(h_x_t_1)
    data_comb[[2]] = t(h_x_t_2)
    
    rowmeans_object = sd_object = decenter_object = list()
    for(ik in 1:2)
    {
        # compute mean and standard deviation functions
        rowmeans_object[[ik]] = rowMeans(data_comb[[ik]], na.rm = TRUE)
        sd_object[[ik]] = apply(data_comb[[ik]], 1, sd, na.rm = TRUE)
        
        # de-center functional data
        decenter_object[[ik]] = t(scale(t(data_comb[[ik]]), center = TRUE, scale = TRUE))
    }
    comb_object = do.call(rbind, decenter_object)
    
    # select ncomp by ER and GR
    
    dum = ER_GR(t(comb_object))
    ncomp = max(dum$k_ER, dum$k_GR)
    rm(dum)
    
    SVD_decomp = svd(t(comb_object))
    basis = as.matrix(SVD_decomp$v[,1:ncomp])
    score = t(basis) %*% comb_object
    
    recon = (basis %*% score) * do.call(c, sd_object) + do.call(c, rowmeans_object)
    resi = rbind(data_comb[[1]], data_comb[[2]]) - recon
        
    # determine in-sample forecast error
    
    olivia = matrix(NA, ncomp, fh)
    if(fore_method == "ets")
    {
        for(ij in 1:ncomp)
        {
            olivia[ij,] = forecast(ets(score[ij,]), h = fh)$mean
        }
    }
    else if(fore_method == "arima")
    {
        for(ij in 1:ncomp)
        {
            olivia[ij,] = forecast(auto.arima(score[ij,]), h = fh)$mean
        }
    }
    else if(fore_method == "rwf")
    {
        for(ij in 1:ncomp)
        {
            olivia[ij,] = rwf(score[ij,], h = fh, drift = TRUE)$mean
        }
    }
    else if(fore_method == "rw")
    {
        for(ij in 1:ncomp)
        {
            olivia[ij,] = rwf(score[ij,], h = fh, drift = FALSE)$mean
        }
    }
    else
    {
        warning("Forecasting method is not on the list.")
    }
    rm(ij)
    forerr = matrix(NA, (n_year - ncomp - fh + 1), ncomp)
    for(i in fh:(n_year - ncomp))
    {
        k = i + (ncomp - fh)
        fore = matrix(NA, 1, ncomp)
        if(fore_method == "ets")
        {
            for(j in 1:ncomp)
            {
                fore[,j] = forecast(ets(score[j,1:k]), h = fh)$mean[fh]
            }
        }
        else if(fore_method == "arima")
        {
            for(j in 1:ncomp)
            {
                fore[,j] = forecast(auto.arima(score[j,1:k]), h = fh)$mean[fh]
            }
        }
        else if(fore_method == "rwf")
        {
            if(k <= 2)
            {
                for(j in 1:ncomp)
                {
                    fore[,j] = score[j,k]
                }
            }
            if(k > 2)
            {
                for(j in 1:ncomp)
                {
                    fore[,j] = rwf(score[j,1:k], h = fh, drift = TRUE)$mean[fh]
                }
            }
        }
        else if(fore_method == "rw")
        {
            if(k == 1)
            {
                for(j in 1:ncomp)
                {
                    fore[,j] = score[j,1]
                }
            }
            if(k > 1)
            {
                for(j in 1:ncomp)
                {
                    fore[,j] = rwf(score[j,1:k], h = fh, drift = FALSE)$mean[fh]
                }
            }
        }
        else
        {
            warning("Forecasting method is not on the list.")
        }
        forerr[i - fh + 1,] = score[, k + fh] - fore
    }
    
    # bootstrapping residuals
    
    q = array(NA, dim = c(n_age * 2, no_boot, no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:(n_age * 2))
        {
            for(k in 1:no_boot)
            {
                q[i,,k,j] = sample(resi[i,], size = no_boot, replace = TRUE)
            }
        }
    }
    rm(i); rm(j); rm(k)
    
    # bootstrapping PC score vectors
    
    ny = array(NA, dim = c(ncomp, no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:ncomp)
        {
            ny[i,,j] = sample(forerr[,i], size = no_boot, replace = TRUE)
        }
        rm(i)
    }
    rm(j)
    
    # adding the PC score error to the predicted score
    
    oli = array(rep(olivia, no_boot * fh), dim = c(ncomp, no_boot, fh))
    fo = array(NA, dim = c(ncomp, no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:no_boot)
        {
            fo[,i,j] = oli[,i,j] + ny[,i,j]
        }
    }
    rm(i); rm(j)
    
    # construct bootstrapped samples
    
    pred = array(NA, dim = c(n_age*2, no_boot, no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:no_boot)
        {
            for(k in 1:no_boot)
            {
                pred[,i,k,j] = (as.matrix(basis) %*% fo[,i,j]) + as.matrix(q[,i,k,j])
            }
        }
    }
    rm(i); rm(j); rm(k)
    
    pred_resize = array(NA, dim = c(n_age*2, no_boot * no_boot, fh))
    for(j in 1:fh)
    {
        for(i in 1:no_boot)
        {
            pred_resize[, (((i-1)*no_boot+1):(i*no_boot)),] = pred[,i,,j]
        }
    }
    rm(i); rm(j)
    
    pred_resize_transform = array(NA, dim = c(n_age * 2, no_boot * no_boot, fh))
    for(iwk in 1:(no_boot * no_boot))
    {
        pred_resize_transform[,iwk,] = pred_resize[,iwk,] * do.call(c, sd_object) + do.call(c, rowmeans_object)
    }
    
    pred_1_resize = pred_2_resize = array(NA, dim = c(n_age, no_boot * no_boot, fh))
    if(fh == 1)
    {
        pred_1_resize[,,1] = pred_resize_transform[1:n_age,,]
        pred_2_resize[,,1] = pred_resize_transform[(n_age + 1):(n_age * 2),,]
    }
    else
    {
        pred_1_resize = pred_resize_transform[1:n_age,,]
        pred_2_resize = pred_resize_transform[(n_age + 1):(n_age * 2),,]
    }
    rm(pred_resize_transform); rm(pred_resize)
    
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
    return(list(d_boot_1 = d_boot_1, d_boot_2 = d_boot_2,
                PI_1 = apply(d_boot_1, c(1, 3), quantile, c(alpha/2, 1-alpha/2)),
                PI_2 = apply(d_boot_2, c(1, 3), quantile, c(alpha/2, 1-alpha/2))))
}

##############
# h = 1 to 20
##############

CoDa_int_mfts <- function(fh, method_fore, sig)
{
    fore_h30_1 = fore_h30_2 = array(NA, dim = c(2, nrow(EW_female_pop), (31 - fh)))
    for(ik in 1:(31 - fh))
    {
        dum = CoDa_recon_int_mfts(dat_1 = t(EW_female_pop)[1:(147+ik),],
                                  dat_2 = t(EW_male_pop)[1:(147+ik),],
                                  fore_method = method_fore,
                                  fh = fh, no_boot = 100, alpha = sig)
        fore_h30_1[,,ik] = dum$PI_1[,,fh]
        fore_h30_2[,,ik] = dum$PI_2[,,fh]
        rm(ik); rm(dum)
    }
    return(list(fore_h30_1 = fore_h30_1, fore_h30_2 = fore_h30_2))
}

##############
# CoDa female
##############

## ARIMA

CoDa_int_mfts_array_ARIMA_80 = list()
for(iwk in 1:30)
{
    CoDa_int_mfts_array_ARIMA_80[[iwk]] = CoDa_int_mfts(fh = iwk, method_fore = "arima", sig = 0.2)
    print(iwk); rm(iwk)
}

CoDa_int_mfts_array_ARIMA_95 = list()
for(iwk in 1:30)
{
    CoDa_int_mfts_array_ARIMA_95[[iwk]] = CoDa_int_mfts(fh = iwk, method_fore = "arima", sig = 0.05)
    print(iwk); rm(iwk)
}

# compute mean interval score and empirical coverage probability deviance

CoDa_mfts_mean_score_female_ARIMA_80 = CoDa_mfts_mean_score_male_ARIMA_80 = 
CoDa_mfts_mean_score_female_ARIMA_95 = CoDa_mfts_mean_score_male_ARIMA_95 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    # female
    
    CoDa_mfts_mean_score_female_ARIMA_80[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                          lb = (CoDa_int_mfts_array_ARIMA_80[[ik]])$fore_h30_1[1,,],
                                                          ub = (CoDa_int_mfts_array_ARIMA_80[[ik]])$fore_h30_1[2,,], alpha = 0.2)    
    
    # male
    
    CoDa_mfts_mean_score_male_ARIMA_80[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                        lb = (CoDa_int_mfts_array_ARIMA_80[[ik]])$fore_h30_2[1,,],
                                                        ub = (CoDa_int_mfts_array_ARIMA_80[[ik]])$fore_h30_2[2,,], alpha = 0.2)    

    # female
    
    CoDa_mfts_mean_score_female_ARIMA_95[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                        lb = (CoDa_int_mfts_array_ARIMA_95[[ik]])$fore_h30_1[1,,],
                                                        ub = (CoDa_int_mfts_array_ARIMA_95[[ik]])$fore_h30_1[2,,], alpha = 0.05)    

    # male
    
    CoDa_mfts_mean_score_male_ARIMA_95[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                        lb = (CoDa_int_mfts_array_ARIMA_95[[ik]])$fore_h30_2[1,,],
                                                        ub = (CoDa_int_mfts_array_ARIMA_95[[ik]])$fore_h30_2[2,,], alpha = 0.05)    
}
colnames(CoDa_mfts_mean_score_female_ARIMA_80) = colnames(CoDa_mfts_mean_score_male_ARIMA_80) = 
colnames(CoDa_mfts_mean_score_female_ARIMA_95) = colnames(CoDa_mfts_mean_score_male_ARIMA_95) = c("Mean interval score", "CPD")

rownames(CoDa_mfts_mean_score_female_ARIMA_80) = rownames(CoDa_mfts_mean_score_male_ARIMA_80) = 
rownames(CoDa_mfts_mean_score_female_ARIMA_95) = rownames(CoDa_mfts_mean_score_male_ARIMA_95) = 1:30

round(colMeans(CoDa_mfts_mean_score_female_ARIMA_80), 4) #  421.5240 0.0889 
round(colMeans(CoDa_mfts_mean_score_female_ARIMA_95), 4) #  657.4738 0.0372 
round(colMeans(CoDa_mfts_mean_score_male_ARIMA_80), 4)   # 2022.3470 0.3308
round(colMeans(CoDa_mfts_mean_score_male_ARIMA_95), 4)   # 4731.3844 0.2927

## ETS

CoDa_int_mfts_array_ETS_80 = list()
for(iwk in 1:30)
{
    CoDa_int_mfts_array_ETS_80[[iwk]] = CoDa_int_mfts(fh = iwk, method_fore = "ets", sig = 0.2)
    print(iwk); rm(iwk)
}

CoDa_int_mfts_array_ETS_95 = list()
for(iwk in 1:30)
{
    CoDa_int_mfts_array_ETS_95[[iwk]] = CoDa_int_mfts(fh = iwk, method_fore = "ets", sig = 0.05)
    print(iwk); rm(iwk)
}

# compute mean interval score and empirical coverage probability deviance

CoDa_mfts_mean_score_female_ETS_80 = CoDa_mfts_mean_score_male_ETS_80 = 
CoDa_mfts_mean_score_female_ETS_95 = CoDa_mfts_mean_score_male_ETS_95 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    # female
    
    CoDa_mfts_mean_score_female_ETS_80[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                               lb = (CoDa_int_mfts_array_ETS_80[[ik]])$fore_h30_1[1,,],
                                                               ub = (CoDa_int_mfts_array_ETS_80[[ik]])$fore_h30_1[2,,], alpha = 0.2)    
    
    # male
    
    CoDa_mfts_mean_score_male_ETS_80[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                             lb = (CoDa_int_mfts_array_ETS_80[[ik]])$fore_h30_2[1,,],
                                                             ub = (CoDa_int_mfts_array_ETS_80[[ik]])$fore_h30_2[2,,], alpha = 0.2)    

    # female
    
    CoDa_mfts_mean_score_female_ETS_95[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                             lb = (CoDa_int_mfts_array_ETS_95[[ik]])$fore_h30_1[1,,],
                                                             ub = (CoDa_int_mfts_array_ETS_95[[ik]])$fore_h30_1[2,,], alpha = 0.05)
    
    # male
    
    CoDa_mfts_mean_score_male_ETS_95[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                           lb = (CoDa_int_mfts_array_ETS_95[[ik]])$fore_h30_2[1,,],
                                                           ub = (CoDa_int_mfts_array_ETS_95[[ik]])$fore_h30_2[2,,], alpha = 0.05)
}

colnames(CoDa_mfts_mean_score_female_ETS_80) = colnames(CoDa_mfts_mean_score_male_ETS_80) = 
colnames(CoDa_mfts_mean_score_female_ETS_95) = colnames(CoDa_mfts_mean_score_male_ETS_95) = c("Mean interval score", "CPD")

rownames(CoDa_mfts_mean_score_female_ETS_80) = rownames(CoDa_mfts_mean_score_male_ETS_80) = 
rownames(CoDa_mfts_mean_score_female_ETS_95) = rownames(CoDa_mfts_mean_score_male_ETS_95) = 1:30

round(colMeans(CoDa_mfts_mean_score_female_ETS_80), 4) #  480.0395 0.0799
round(colMeans(CoDa_mfts_mean_score_female_ETS_95), 4) #  676.0784 0.0409
round(colMeans(CoDa_mfts_mean_score_male_ETS_80), 4)   # 2120.3771 0.3094
round(colMeans(CoDa_mfts_mean_score_male_ETS_95), 4)   # 5435.2535 0.2819


## RWD

CoDa_int_mfts_array_RWD_80 = list()
for(iwk in 1:30)
{  
    CoDa_int_mfts_array_RWD_80[[iwk]] = CoDa_int_mfts(fh = iwk, method_fore = "rwf", sig = 0.2)
    print(iwk); rm(iwk)
}

CoDa_int_mfts_array_RWD_95 = list()
for(iwk in 1:30)
{
    CoDa_int_mfts_array_RWD_95[[iwk]] = CoDa_int_mfts(fh = iwk, method_fore = "rwf", sig = 0.05)
    print(iwk); rm(iwk)
}

# compute mean interval score and empirical coverage probability deviance

CoDa_mfts_mean_score_female_RWD_80 = CoDa_mfts_mean_score_male_RWD_80 = 
CoDa_mfts_mean_score_female_RWD_95 = CoDa_mfts_mean_score_male_RWD_95 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    # female
    
    CoDa_mfts_mean_score_female_RWD_80[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                             lb = (CoDa_int_mfts_array_RWD_80[[ik]])$fore_h30_1[1,,],
                                                             ub = (CoDa_int_mfts_array_RWD_80[[ik]])$fore_h30_1[2,,], alpha = 0.2)    
    
    # male
    
    CoDa_mfts_mean_score_male_RWD_80[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                           lb = (CoDa_int_mfts_array_RWD_80[[ik]])$fore_h30_2[1,,],
                                                           ub = (CoDa_int_mfts_array_RWD_80[[ik]])$fore_h30_2[2,,], alpha = 0.2)    
    
    # female
    
    CoDa_mfts_mean_score_female_RWD_95[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                             lb = (CoDa_int_mfts_array_RWD_95[[ik]])$fore_h30_1[1,,],
                                                             ub = (CoDa_int_mfts_array_RWD_95[[ik]])$fore_h30_1[2,,], alpha = 0.05)
    
    # male
    
    CoDa_mfts_mean_score_male_RWD_95[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                           lb = (CoDa_int_mfts_array_RWD_95[[ik]])$fore_h30_2[1,,],
                                                           ub = (CoDa_int_mfts_array_RWD_95[[ik]])$fore_h30_2[2,,], alpha = 0.05)
}

colnames(CoDa_mfts_mean_score_female_RWD_80) = colnames(CoDa_mfts_mean_score_male_RWD_80) = 
colnames(CoDa_mfts_mean_score_female_RWD_95) = colnames(CoDa_mfts_mean_score_male_RWD_95) = c("Mean interval score", "CPD")

rownames(CoDa_mfts_mean_score_female_RWD_80) = rownames(CoDa_mfts_mean_score_male_RWD_80) = 
rownames(CoDa_mfts_mean_score_female_RWD_95) = rownames(CoDa_mfts_mean_score_male_RWD_95) = 1:30

round(colMeans(CoDa_mfts_mean_score_female_RWD_80), 4) #  443.6864 0.1581
round(colMeans(CoDa_mfts_mean_score_female_RWD_95), 4) #  796.8239 0.0485
round(colMeans(CoDa_mfts_mean_score_male_RWD_80), 4)   # 1503.8608 0.2151
round(colMeans(CoDa_mfts_mean_score_male_RWD_95), 4)   # 2344.6733 0.1479


## RW

CoDa_int_mfts_array_RW_80 = list()
for(iwk in 1:30)
{
    CoDa_int_mfts_array_RW_80[[iwk]] = CoDa_int_mfts(fh = iwk, method_fore = "rw", sig = 0.2)
    print(iwk); rm(iwk)
}

CoDa_int_mfts_array_RW_95 = list()
for(iwk in 1:30)
{
    CoDa_int_mfts_array_RW_95[[iwk]] = CoDa_int_mfts(fh = iwk, method_fore = "rw", sig = 0.05)
    print(iwk); rm(iwk)
}

# compute mean interval score and empirical coverage probability deviance

CoDa_mfts_mean_score_female_RW_80 = CoDa_mfts_mean_score_male_RW_80 = 
CoDa_mfts_mean_score_female_RW_95 = CoDa_mfts_mean_score_male_RW_95 = matrix(NA, 30, 2)
for(ik in 1:30)
{
    # female
    
    CoDa_mfts_mean_score_female_RW_80[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                             lb = (CoDa_int_mfts_array_RW_80[[ik]])$fore_h30_1[1,,],
                                                             ub = (CoDa_int_mfts_array_RW_80[[ik]])$fore_h30_1[2,,], alpha = 0.2)    
    
    # male
    
    CoDa_mfts_mean_score_male_RW_80[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                           lb = (CoDa_int_mfts_array_RW_80[[ik]])$fore_h30_2[1,,],
                                                           ub = (CoDa_int_mfts_array_RW_80[[ik]])$fore_h30_2[2,,], alpha = 0.2)    
    
    # female
    
    CoDa_mfts_mean_score_female_RW_95[ik,] = interval_score(holdout = EW_female_pop[,(148+ik):178],
                                                             lb = (CoDa_int_mfts_array_RW_95[[ik]])$fore_h30_1[1,,],
                                                             ub = (CoDa_int_mfts_array_RW_95[[ik]])$fore_h30_1[2,,], alpha = 0.05)
    
    # male
    
    CoDa_mfts_mean_score_male_RW_95[ik,] = interval_score(holdout = EW_male_pop[,(148+ik):178],
                                                           lb = (CoDa_int_mfts_array_RW_95[[ik]])$fore_h30_2[1,,],
                                                           ub = (CoDa_int_mfts_array_RW_95[[ik]])$fore_h30_2[2,,], alpha = 0.05)
}

colnames(CoDa_mfts_mean_score_female_RW_80) = colnames(CoDa_mfts_mean_score_male_RW_80) = 
colnames(CoDa_mfts_mean_score_female_RW_95) = colnames(CoDa_mfts_mean_score_male_RW_95) = c("Mean interval score", "CPD")

rownames(CoDa_mfts_mean_score_female_RW_80) = rownames(CoDa_mfts_mean_score_male_RW_80) = 
rownames(CoDa_mfts_mean_score_female_RW_95) = rownames(CoDa_mfts_mean_score_male_RW_95) = 1:30

round(colMeans(CoDa_mfts_mean_score_female_RW_80), 4) #  459.7545 0.1121 
round(colMeans(CoDa_mfts_mean_score_female_RW_95), 4) #  673.5936 0.0452 
round(colMeans(CoDa_mfts_mean_score_male_RW_80), 4)   # 2037.3201 0.2731
round(colMeans(CoDa_mfts_mean_score_male_RW_95), 4)   # 5095.6282 0.2608

############################
# Figures in the Appendices
############################

plot(1:30, CoDa_mean_score_female_ARIMA_80[,1], col = 1, lty = 1, type = "l", xlab = "",
     ylab = expression(bar(S)[0.2](h)), ylim = c(200, 1050), main = "Female")
lines(1:30, CoDa_mfts_mean_score_female_ARIMA_80[,1], col = 2, lty = 2)
lines(1:30, CoDa_mean_score_mlfts_female_ARIMA_80[,1], col = 4, lty = 4)
legend("topleft", c("Univariate functional time series method", "Multivariate functional time series method", "Multilevel functional time series method"), 
       col = c(1,2,4), lty = c(1,2,4), cex=0.7)


plot(1:30, CoDa_mean_score_female_ARIMA_95[,1], col = 1, lty = 1, type = "l", xlab = "Forecast horizon",
     ylab = expression(bar(S)[0.05](h)), ylim = c(400, 1600), main = "Female")
lines(1:30, CoDa_mfts_mean_score_female_ARIMA_95[,1], col = 2, lty = 2)
lines(1:30, CoDa_mean_score_mlfts_female_ARIMA_95[,1], col = 4, lty = 4)

