#################
# load R package
#################

require(psych)
require(flexmix)
require(ftsa)
require(meboot)
require(pracma)
mape = ftsa:::mape

source("ER_GR.R")
source("CoDa_model_fitting.R")

# dat_1: female data set (n x p)
# dat_2: male data set (n x p)
# fh: forecast horizon from one to 30
# fmethod: forecasting method = c("ARIMA", "ETS", "RW", "RWD")

CoDa_point <- function(dat_1, dat_2, fh, fmethod = c("ARIMA", "ETS", "RWF_no_drift", "RWF_drift"))
{
    # dat: original data matrix (n by p)
    # fh: forecast horizon
    # fmethod: forecasting method
    
    fore_h30_FTS_1 = fore_h30_MFTS_1 = fore_h30_MLFTS_1 = 
    fore_h30_FTS_2 = fore_h30_MFTS_2 = fore_h30_MLFTS_2 = matrix(NA, ncol(dat_1), (31 - fh))
    
    fore_h30_FTS_1_KLdiv = fore_h30_FTS_2_KLdiv = 
    fore_h30_MFTS_1_KLdiv = fore_h30_MFTS_2_KLdiv = 
    fore_h30_MLFTS_1_KLdiv = fore_h30_MLFTS_2_KLdiv = 
    RW_1_KLdiv = RW_2_KLdiv = 
    RWD_1_KLdiv = RWD_2_KLdiv = matrix(NA, (31 - fh), 2)
        
    fore_h30_FTS_1_density_norm = fore_h30_FTS_1_mape = fore_h30_FTS_1_JSdiv = fore_h30_FTS_1_JSdiv_geo = 
    fore_h30_FTS_2_density_norm = fore_h30_FTS_2_mape = fore_h30_FTS_2_JSdiv = fore_h30_FTS_2_JSdiv_geo =    
        
    fore_h30_MFTS_1_density_norm = fore_h30_MFTS_1_mape = fore_h30_MFTS_1_JSdiv = fore_h30_MFTS_1_JSdiv_geo = 
    fore_h30_MFTS_2_density_norm = fore_h30_MFTS_2_mape = fore_h30_MFTS_2_JSdiv = fore_h30_MFTS_2_JSdiv_geo =     
        
    fore_h30_MLFTS_1_density_norm = fore_h30_MLFTS_1_mape = fore_h30_MLFTS_1_JSdiv = fore_h30_MLFTS_1_JSdiv_geo = 
    fore_h30_MLFTS_2_density_norm = fore_h30_MLFTS_2_mape = fore_h30_MLFTS_2_JSdiv = fore_h30_MLFTS_2_JSdiv_geo = 
    
    RW_1_density_norm = RW_1_mape = RW_1_JSdiv = RW_1_JSdiv_geo = 
    RW_2_density_norm = RW_2_mape = RW_2_JSdiv = RW_2_JSdiv_geo = 
        
    RWD_1_density_norm = RWD_1_mape = RWD_1_JSdiv = RWD_1_JSdiv_geo = 
    RWD_2_density_norm = RWD_2_mape = RWD_2_JSdiv = RWD_2_JSdiv_geo = vector("numeric", (31 - fh))
    for(ik in 1:(31 - fh))
    {
        #######
        ## FTS 
        #######
        
        R2_female = R_square_fit(dat_1 = dat_1[1:(147+ik),], dat_2 = dat_2[1:(147+ik),],
                                 fh = fh, modeling_method = "FTS",
                                 forecasting_method = fmethod)
        fore_h30_FTS_1[,ik] = R2_female$fore_count_1[,fh]
        
        # MAPE (1st population)
        
        fore_h30_FTS_1_mape[ik] = mape(forecast = t(fore_h30_FTS_1[,ik]), true = dat_1[(147+fh+ik),])
        
        # density_norm
        
        fore_h30_FTS_1_density_norm[ik] = density_norm(d1 = t(fore_h30_FTS_1[,ik]), d2 = dat_1[(147+fh+ik),], 
                                                       time_vector = 1:111)
        
        # KLdiv
        
        dat = cbind(true = dat_1[(147+fh+ik),], forecast = as.numeric(t(fore_h30_FTS_1[,ik])))
        colnames(dat) = c("True", "Estimate")
        fore_h30_FTS_1_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)        
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_1[(147+fh+ik),], forecast = as.numeric(t(fore_h30_FTS_1[,ik])))
        M = rowMeans(dat)
        P_M = cbind(dat_1[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_FTS_1[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_FTS_1_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_1[(147+fh+ik),], forecast = as.numeric(t(fore_h30_FTS_1[,ik])))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_1[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_FTS_1[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_FTS_1_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        

        # MAPE (2nd population)
        
        fore_h30_FTS_2[,ik] = R2_female$fore_count_2[,fh]
        fore_h30_FTS_2_mape[ik] = mape(forecast = t(fore_h30_FTS_2[,ik]), true = dat_2[(147+fh+ik),])
        
        # density_norm
        
        fore_h30_FTS_2_density_norm[ik] = density_norm(d1 = t(fore_h30_FTS_2[,ik]), d2 = dat_2[(147+fh+ik),], 
                                                       time_vector = 1:111)
        
        # KLdiv
        
        dat = cbind(true = dat_2[(147+fh+ik),], forecast = as.numeric(t(fore_h30_FTS_2[,ik])))
        colnames(dat) = c("True", "Estimate")
        fore_h30_FTS_2_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_2[(147+fh+ik),], forecast = as.numeric(t(fore_h30_FTS_2[,ik])))
        M = rowMeans(dat)
        P_M = cbind(dat_2[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_FTS_2[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_FTS_2_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_2[(147+fh+ik),], forecast = as.numeric(t(fore_h30_FTS_2[,ik])))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_2[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_FTS_2[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_FTS_2_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        rm(R2_female)
        
        ########
        ## MFTS
        ########
        
        R2_female = R_square_fit(dat_1 = dat_1[1:(147+ik),], dat_2 = dat_2[1:(147+ik),], 
                                 fh = fh, modeling_method = "MFTS",
                                 forecasting_method = fmethod)
        fore_h30_MFTS_1[,ik] = R2_female$fore_count_1[,fh]
        
        # MAPE (1st population)
        
        fore_h30_MFTS_1_mape[ik] = mape(forecast = t(fore_h30_MFTS_1[,ik]), true = dat_1[(147+fh+ik),])
        
        # density_norm
        
        fore_h30_MFTS_1_density_norm[ik] = density_norm(d1 = t(fore_h30_MFTS_1[,ik]), d2 = dat_1[(147+fh+ik),], 
                                                        time_vector = 1:111)
        
        # KLdiv
        
        dat = cbind(true = dat_1[(147+fh+ik),], forecast = as.numeric(t(fore_h30_MFTS_1[,ik])))
        colnames(dat) = c("True", "Estimate")      
        fore_h30_MFTS_1_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_1[(147+fh+ik),], forecast = as.numeric(t(fore_h30_MFTS_1[,ik])))
        M = rowMeans(dat)
        P_M = cbind(dat_1[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_MFTS_1[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")        
        fore_h30_MFTS_1_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_1[(147+fh+ik),], forecast = as.numeric(t(fore_h30_MFTS_1[,ik])))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_1[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_MFTS_1[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_MFTS_1_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        fore_h30_MFTS_2[,ik] = R2_female$fore_count_2[,fh]
        
        # MAPE (2nd population)
        
        fore_h30_MFTS_2_mape[ik] = mape(forecast = t(fore_h30_MFTS_2[,ik]), true = dat_2[(147+fh+ik),])
        
        # density_norm
        
        fore_h30_MFTS_2_density_norm[ik] = density_norm(d1 = t(fore_h30_MFTS_2[,ik]), d2 = dat_2[(147+fh+ik),], 
                                                        time_vector = 1:111)
        
        # KLdiv
        
        dat = cbind(true = dat_2[(147+fh+ik),], forecast = as.numeric(t(fore_h30_MFTS_2[,ik])))
        colnames(dat) = c("True", "Estimate")
        fore_h30_MFTS_2_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_2[(147+fh+ik),], forecast = as.numeric(t(fore_h30_MFTS_2[,ik])))
        M = rowMeans(dat)
        P_M = cbind(dat_2[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_MFTS_2[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_MFTS_2_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_2[(147+fh+ik),], forecast = as.numeric(t(fore_h30_MFTS_2[,ik])))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_2[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_MFTS_2[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_MFTS_2_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        rm(R2_female)
        
        #########
        ## MLFTS
        #########
        
        R2_female = R_square_fit(dat_1 = dat_1[1:(147+ik),], dat_2 = dat_2[1:(147+ik),], 
                                 fh = fh, modeling_method = "MLFTS",
                                 forecasting_method = fmethod)
        fore_h30_MLFTS_1[,ik] = R2_female$fore_count_1[,fh]
        
        # MAPE (1st population)
        
        fore_h30_MLFTS_1_mape[ik] = mape(forecast = t(fore_h30_MLFTS_1[,ik]), true = dat_1[(147+fh+ik),])
        
        # density_norm
        
        fore_h30_MLFTS_1_density_norm[ik] = density_norm(d1 = t(fore_h30_MLFTS_1[,ik]), d2 = dat_1[(147+fh+ik),], 
                                                         time_vector = 1:111)
        
        # KLdiv
        
        dat = cbind(true = dat_1[(147+fh+ik),], forecast = as.numeric(t(fore_h30_MLFTS_1[,ik])))
        colnames(dat) = c("True", "Estimate")
        fore_h30_MLFTS_1_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_1[(147+fh+ik),], forecast = as.numeric(t(fore_h30_MLFTS_1[,ik])))
        M = rowMeans(dat)
        P_M = cbind(dat_1[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_MLFTS_1[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_MLFTS_1_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_1[(147+fh+ik),], forecast = as.numeric(t(fore_h30_MLFTS_1[,ik])))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_1[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_MLFTS_1[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_MLFTS_1_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        
        fore_h30_MLFTS_2[,ik] = R2_female$fore_count_2[,fh]
        
        # MAPE (2nd population)
        
        fore_h30_MLFTS_2_mape[ik] = mape(forecast = t(fore_h30_MLFTS_2[,ik]), true = dat_2[(147+fh+ik),])
        
        # density_norm
        
        fore_h30_MLFTS_2_density_norm[ik] = density_norm(d1 = t(fore_h30_MLFTS_2[,ik]), d2 = dat_2[(147+fh+ik),], 
                                                         time_vector = 1:111)
        
        # KLdiv
        
        dat = cbind(true = dat_2[(147+fh+ik),], forecast = as.numeric(t(fore_h30_MLFTS_2[,ik])))
        colnames(dat) = c("True", "Estimate")
        fore_h30_MLFTS_2_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_2[(147+fh+ik),], forecast = as.numeric(t(fore_h30_MLFTS_2[,ik])))
        M = rowMeans(dat)
        P_M = cbind(dat_2[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_MLFTS_2[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_MLFTS_2_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_2[(147+fh+ik),], forecast = as.numeric(t(fore_h30_MLFTS_2[,ik])))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_2[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_MLFTS_2[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_MLFTS_2_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        rm(R2_female)
        
        ####################
        ## CoDa (naive RWD)
        ####################
        
        R2_female = naive_fun(dat = dat_1[1:(147+ik),], fh = fh, drift_term = TRUE, level_ci = 80)$fore_count[,fh]
        
        # MAPE (1st population)
        
        RWD_1_mape[ik] = mape(forecast = t(R2_female), true = dat_1[(147+fh+ik),])        
        
        # density_norm
        
        RWD_1_density_norm[ik] = density_norm(d1 = t(R2_female), d2 = dat_1[(147+fh+ik),], time_vector = 1:111)
        
        # KLdiv
        
        dat = cbind(true = dat_1[(147+fh+ik),], forecast = as.numeric(t(R2_female)))
        colnames(dat) = c("True", "Estimate")
        RWD_1_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_1[(147+fh+ik),], forecast = as.numeric(t(R2_female)))
        M = rowMeans(dat)
        P_M = cbind(dat_1[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(R2_female)), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        RWD_1_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M); rm(dat)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_1[(147+fh+ik),], forecast = as.numeric(t(R2_female)))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_1[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(R2_female)), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        RWD_1_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M); rm(dat); rm(R2_female)
        
        
        R2_female = naive_fun(dat = dat_2[1:(147+ik),], fh = fh, drift_term = TRUE, level_ci = 80)$fore_count[,fh]
        
        # MAPE (2nd population)
        
        RWD_2_mape[ik] = mape(forecast = t(R2_female), true = dat_2[(147+fh+ik),])        
        
        # density_norm
        
        RWD_2_density_norm[ik] = density_norm(d1 = t(R2_female), d2 = dat_2[(147+fh+ik),], time_vector = 1:111)
        
        # KLdiv
        
        dat = cbind(true = dat_2[(147+fh+ik),], forecast = as.numeric(t(R2_female)))
        colnames(dat) = c("True", "Estimate")
        RWD_2_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_2[(147+fh+ik),], forecast = as.numeric(t(R2_female)))
        M = rowMeans(dat)
        P_M = cbind(dat_2[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(R2_female)), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        RWD_2_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M); rm(dat)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_2[(147+fh+ik),], forecast = as.numeric(t(R2_female)))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_2[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(R2_female)), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        RWD_2_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M); rm(dat); rm(R2_female)
        
        ##################
        ## CoDa (naive RW)
        ##################
        
        R2_female = naive_fun(dat = dat_1[1:(147+ik),], fh = fh, drift_term = FALSE, level_ci = 80)$fore_count[,fh]
        
        # MAPE (1st population)
        
        RW_1_mape[ik] = mape(forecast = t(R2_female), true = dat_1[(147+fh+ik),])        
        
        # density_norm
        
        RW_1_density_norm[ik] = density_norm(d1 = t(R2_female), d2 = dat_1[(147+fh+ik),], time_vector = 1:111)
        
        # KLdiv
        
        dat = cbind(true = dat_1[(147+fh+ik),], forecast = as.numeric(t(R2_female)))
        colnames(dat) = c("True", "Estimate")
        RW_1_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_1[(147+fh+ik),], forecast = as.numeric(t(R2_female)))
        M = rowMeans(dat)
        P_M = cbind(dat_1[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(R2_female)), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        RW_1_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M); rm(dat)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_1[(147+fh+ik),], forecast = as.numeric(t(R2_female)))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_1[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(R2_female)), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        RW_1_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M); rm(dat); rm(R2_female)
        
        
        R2_female = naive_fun(dat = dat_2[1:(147+ik),], fh = fh, drift_term = FALSE, level_ci = 80)$fore_count[,fh]
        
        # MAPE (2nd population)
        
        RW_2_mape[ik] = mape(forecast = t(R2_female), true = dat_2[(147+fh+ik),])        
        
        # density_norm
        
        RW_2_density_norm[ik] = density_norm(d1 = t(R2_female), d2 = dat_2[(147+fh+ik),], time_vector = 1:111)
        
        # KLdiv
        
        dat = cbind(true = dat_2[(147+fh+ik),], forecast = as.numeric(t(R2_female)))
        colnames(dat) = c("True", "Estimate")
        RW_2_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_2[(147+fh+ik),], forecast = as.numeric(t(R2_female)))
        M = rowMeans(dat)
        P_M = cbind(dat_2[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(R2_female)), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        RW_2_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M); rm(dat)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_2[(147+fh+ik),], forecast = as.numeric(t(R2_female)))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_2[(147+fh+ik),], M)
        E_M = cbind(as.numeric(t(R2_female)), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        RW_2_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M); rm(dat); rm(R2_female)
        rm(ik)
    }
    
    ###############################################################
    # take mean over the number of years in the forecasting period
    ###############################################################
    
    return(list(RW_1_err = mean(RW_1_mape), 
                RWD_1_err = mean(RWD_1_mape),
                CoDa_FTS_1_err = mean(fore_h30_FTS_1_mape),
                CoDa_MFTS_1_err = mean(fore_h30_MFTS_1_mape),
                CoDa_MLFTS_1_err = mean(fore_h30_MLFTS_1_mape),
                
                RW_1_density_norm = mean(RW_1_density_norm),
                RWD_1_density_norm = mean(RWD_1_density_norm),
                CoDa_FTS_1_density_norm = mean(fore_h30_FTS_1_density_norm),
                CoDa_MFTS_1_density_norm = mean(fore_h30_MFTS_1_density_norm),
                CoDa_MLFTS_1_density_norm = mean(fore_h30_MLFTS_1_density_norm),
                
                RW_1_KLdiv = mean(RW_1_KLdiv),
                RWD_1_KLdiv = mean(RWD_1_KLdiv),
                CoDa_FTS_1_KLdiv = mean(fore_h30_FTS_1_KLdiv),
                CoDa_MFTS_1_KLdiv = mean(fore_h30_MFTS_1_KLdiv),                
                CoDa_MLFTS_1_KLdiv = mean(fore_h30_MLFTS_1_KLdiv),
                
                RW_1_JSdiv = mean(RW_1_JSdiv),
                RWD_1_JSdiv = mean(RWD_1_JSdiv),
                CoDa_FTS_1_JSdiv = mean(fore_h30_FTS_1_JSdiv),
                CoDa_MFTS_1_JSdiv = mean(fore_h30_MFTS_1_JSdiv),
                CoDa_MLFTS_1_JSdiv = mean(fore_h30_MLFTS_1_JSdiv),
                
                RW_1_JSdiv_geo = mean(RW_1_JSdiv_geo),
                RWD_1_JSdiv_geo = mean(RWD_1_JSdiv_geo),
                CoDa_FTS_1_JSdiv_geo = mean(fore_h30_FTS_1_JSdiv_geo),
                CoDa_MFTS_1_JSdiv_geo = mean(fore_h30_MFTS_1_JSdiv_geo),
                CoDa_MLFTS_1_JSdiv_geo = mean(fore_h30_MLFTS_1_JSdiv_geo),
                
                RW_2_err = mean(RW_2_mape), 
                RWD_2_err = mean(RWD_2_mape),
                CoDa_FTS_2_err = mean(fore_h30_FTS_2_mape),
                CoDa_MFTS_2_err = mean(fore_h30_MFTS_2_mape),
                CoDa_MLFTS_2_err = mean(fore_h30_MLFTS_2_mape),
                
                RW_2_density_norm = mean(RW_2_density_norm),
                RWD_2_density_norm = mean(RWD_2_density_norm),
                CoDa_FTS_2_density_norm = mean(fore_h30_FTS_2_density_norm),
                CoDa_MFTS_2_density_norm = mean(fore_h30_MFTS_2_density_norm),
                CoDa_MLFTS_2_density_norm = mean(fore_h30_MLFTS_2_density_norm),
                
                RW_2_KLdiv = mean(RW_2_KLdiv),
                RWD_2_KLdiv = mean(RWD_2_KLdiv),
                CoDa_FTS_2_KLdiv = mean(fore_h30_FTS_2_KLdiv),
                CoDa_MFTS_2_KLdiv = mean(fore_h30_MFTS_2_KLdiv),                
                CoDa_MLFTS_2_KLdiv = mean(fore_h30_MLFTS_2_KLdiv),
                
                RW_2_JSdiv = mean(RW_2_JSdiv),
                RWD_2_JSdiv = mean(RWD_2_JSdiv),
                CoDa_FTS_2_JSdiv = mean(fore_h30_FTS_2_JSdiv),
                CoDa_MFTS_2_JSdiv = mean(fore_h30_MFTS_2_JSdiv),
                CoDa_MLFTS_2_JSdiv = mean(fore_h30_MLFTS_2_JSdiv),
                
                RW_2_JSdiv_geo = mean(RW_2_JSdiv_geo),
                RWD_2_JSdiv_geo = mean(RWD_2_JSdiv_geo),
                CoDa_FTS_2_JSdiv_geo = mean(fore_h30_FTS_2_JSdiv_geo),
                CoDa_MFTS_2_JSdiv_geo = mean(fore_h30_MFTS_2_JSdiv_geo),
                CoDa_MLFTS_2_JSdiv_geo = mean(fore_h30_MLFTS_2_JSdiv_geo)))
}

#########
## ARIMA 
#########

RW_1_point_ARIMA_err = RWD_1_point_ARIMA_err = CoDa_FTS_1_point_ARIMA_err = 
CoDa_MFTS_1_point_ARIMA_err = CoDa_MLFTS_1_point_ARIMA_err = 

RW_1_point_ARIMA_density_norm = RWD_1_point_ARIMA_density_norm = CoDa_FTS_1_point_ARIMA_density_norm = 
CoDa_MFTS_1_point_ARIMA_density_norm = CoDa_MLFTS_1_point_ARIMA_density_norm =     
    
RW_1_point_ARIMA_KLdiv = RWD_1_point_ARIMA_KLdiv = CoDa_FTS_1_point_ARIMA_KLdiv = 
CoDa_MFTS_1_point_ARIMA_KLdiv = CoDa_MLFTS_1_point_ARIMA_KLdiv = 
    
RW_1_point_ARIMA_JSdiv = RWD_1_point_ARIMA_JSdiv = CoDa_FTS_1_point_ARIMA_JSdiv = 
CoDa_MFTS_1_point_ARIMA_JSdiv = CoDa_MLFTS_1_point_ARIMA_JSdiv = 
    
RW_1_point_ARIMA_JSdiv_geo = RWD_1_point_ARIMA_JSdiv_geo = CoDa_FTS_1_point_ARIMA_JSdiv_geo = 
CoDa_MFTS_1_point_ARIMA_JSdiv_geo = CoDa_MLFTS_1_point_ARIMA_JSdiv_geo = 
    
RW_2_point_ARIMA_err = RWD_2_point_ARIMA_err = CoDa_FTS_2_point_ARIMA_err = 
CoDa_MFTS_2_point_ARIMA_err = CoDa_MLFTS_2_point_ARIMA_err = 
    
RW_2_point_ARIMA_density_norm = RWD_2_point_ARIMA_density_norm = CoDa_FTS_2_point_ARIMA_density_norm = 
CoDa_MFTS_2_point_ARIMA_density_norm = CoDa_MLFTS_2_point_ARIMA_density_norm = 
    
RW_2_point_ARIMA_KLdiv = RWD_2_point_ARIMA_KLdiv = CoDa_FTS_2_point_ARIMA_KLdiv = 
CoDa_MFTS_2_point_ARIMA_KLdiv = CoDa_MLFTS_2_point_ARIMA_KLdiv = 
    
RW_2_point_ARIMA_JSdiv = RWD_2_point_ARIMA_JSdiv = CoDa_FTS_2_point_ARIMA_JSdiv = 
CoDa_MFTS_2_point_ARIMA_JSdiv = CoDa_MLFTS_2_point_ARIMA_JSdiv = 
    
RW_2_point_ARIMA_JSdiv_geo = RWD_2_point_ARIMA_JSdiv_geo = CoDa_FTS_2_point_ARIMA_JSdiv_geo = 
CoDa_MFTS_2_point_ARIMA_JSdiv_geo = CoDa_MLFTS_2_point_ARIMA_JSdiv_geo = vector("numeric", 30)
for(ik in 1:30)
{
    dum = CoDa_point(dat_1 = t(EW_female_pop), dat_2 = t(EW_male_pop), fh = ik, fmethod = "ARIMA")
    
    ## 1st population
    
    # MAPE
    
    RW_1_point_ARIMA_err[ik] = dum$RW_1_err
    RWD_1_point_ARIMA_err[ik] = dum$RWD_1_err
    CoDa_FTS_1_point_ARIMA_err[ik] = dum$CoDa_FTS_1_err
    CoDa_MFTS_1_point_ARIMA_err[ik] = dum$CoDa_MFTS_1_err
    CoDa_MLFTS_1_point_ARIMA_err[ik] = dum$CoDa_MLFTS_1_err
    
    # density_norm
    
    RW_1_point_ARIMA_density_norm[ik] = dum$RW_1_density_norm
    RWD_1_point_ARIMA_density_norm[ik] = dum$RWD_1_density_norm
    CoDa_FTS_1_point_ARIMA_density_norm[ik] = dum$CoDa_FTS_1_density_norm
    CoDa_MFTS_1_point_ARIMA_density_norm[ik] = dum$CoDa_MFTS_1_density_norm
    CoDa_MLFTS_1_point_ARIMA_density_norm[ik] = dum$CoDa_MLFTS_1_density_norm
    
    # KL-div
    
    RW_1_point_ARIMA_KLdiv[ik] = dum$RW_1_KLdiv
    RWD_1_point_ARIMA_KLdiv[ik] = dum$RWD_1_KLdiv
    CoDa_FTS_1_point_ARIMA_KLdiv[ik] = dum$CoDa_FTS_1_KLdiv
    CoDa_MFTS_1_point_ARIMA_KLdiv[ik] = dum$CoDa_MFTS_1_KLdiv
    CoDa_MLFTS_1_point_ARIMA_KLdiv[ik] = dum$CoDa_MLFTS_1_KLdiv
    
    # JSdiv_simple
    
    RW_1_point_ARIMA_JSdiv[ik] = dum$RW_1_JSdiv
    RWD_1_point_ARIMA_JSdiv[ik] = dum$RWD_1_JSdiv
    CoDa_FTS_1_point_ARIMA_JSdiv[ik] = dum$CoDa_FTS_1_JSdiv
    CoDa_MFTS_1_point_ARIMA_JSdiv[ik] = dum$CoDa_MFTS_1_JSdiv
    CoDa_MLFTS_1_point_ARIMA_JSdiv[ik] = dum$CoDa_MLFTS_1_JSdiv
    
    # JSdiv_geo
    
    RW_1_point_ARIMA_JSdiv_geo[ik] = dum$RW_1_JSdiv_geo
    RWD_1_point_ARIMA_JSdiv_geo[ik] = dum$RWD_1_JSdiv_geo
    CoDa_FTS_1_point_ARIMA_JSdiv_geo[ik] = dum$CoDa_FTS_1_JSdiv_geo
    CoDa_MFTS_1_point_ARIMA_JSdiv_geo[ik] = dum$CoDa_MFTS_1_JSdiv_geo
    CoDa_MLFTS_1_point_ARIMA_JSdiv_geo[ik] = dum$CoDa_MLFTS_1_JSdiv_geo
    
    ## 2nd population
    
    # MAPE
    
    RW_2_point_ARIMA_err[ik] = dum$RW_2_err
    RWD_2_point_ARIMA_err[ik] = dum$RWD_2_err
    CoDa_FTS_2_point_ARIMA_err[ik] = dum$CoDa_FTS_2_err
    CoDa_MFTS_2_point_ARIMA_err[ik] = dum$CoDa_MFTS_2_err
    CoDa_MLFTS_2_point_ARIMA_err[ik] = dum$CoDa_MLFTS_2_err
    
    # density_norm
    
    RW_2_point_ARIMA_density_norm[ik] = dum$RW_2_density_norm
    RWD_2_point_ARIMA_density_norm[ik] = dum$RWD_2_density_norm
    CoDa_FTS_2_point_ARIMA_density_norm[ik] = dum$CoDa_FTS_2_density_norm
    CoDa_MFTS_2_point_ARIMA_density_norm[ik] = dum$CoDa_MFTS_2_density_norm
    CoDa_MLFTS_2_point_ARIMA_density_norm[ik] = dum$CoDa_MLFTS_2_density_norm
    
    # KL-div
    
    RW_2_point_ARIMA_KLdiv[ik] = dum$RW_2_KLdiv
    RWD_2_point_ARIMA_KLdiv[ik] = dum$RWD_2_KLdiv
    CoDa_FTS_2_point_ARIMA_KLdiv[ik] = dum$CoDa_FTS_2_KLdiv
    CoDa_MFTS_2_point_ARIMA_KLdiv[ik] = dum$CoDa_MFTS_2_KLdiv
    CoDa_MLFTS_2_point_ARIMA_KLdiv[ik] = dum$CoDa_MLFTS_2_KLdiv
    
    # JSdiv_simple
    
    RW_2_point_ARIMA_JSdiv[ik] = dum$RW_2_JSdiv
    RWD_2_point_ARIMA_JSdiv[ik] = dum$RWD_2_JSdiv
    CoDa_FTS_2_point_ARIMA_JSdiv[ik] = dum$CoDa_FTS_2_JSdiv
    CoDa_MFTS_2_point_ARIMA_JSdiv[ik] = dum$CoDa_MFTS_2_JSdiv
    CoDa_MLFTS_2_point_ARIMA_JSdiv[ik] = dum$CoDa_MLFTS_2_JSdiv
    
    # JSdiv_geo
    
    RW_2_point_ARIMA_JSdiv_geo[ik] = dum$RW_2_JSdiv_geo
    RWD_2_point_ARIMA_JSdiv_geo[ik] = dum$RWD_2_JSdiv_geo
    CoDa_FTS_2_point_ARIMA_JSdiv_geo[ik] = dum$CoDa_FTS_2_JSdiv_geo
    CoDa_MFTS_2_point_ARIMA_JSdiv_geo[ik] = dum$CoDa_MFTS_2_JSdiv_geo
    CoDa_MLFTS_2_point_ARIMA_JSdiv_geo[ik] = dum$CoDa_MLFTS_2_JSdiv_geo
    print(ik); rm(dum); rm(ik)
}

#################
# 1st population
#################

CoDa_FTS_1_ARIMA_MAPE = round(mean(CoDa_FTS_1_point_ARIMA_err), 4)   # 29.8999
CoDa_MFTS_1_ARIMA_MAPE = round(mean(CoDa_MFTS_1_point_ARIMA_err), 4)  # 18.9661 (28.3678)
CoDa_MLFTS_1_ARIMA_MAPE = round(mean(CoDa_MLFTS_1_point_ARIMA_err), 4) # 22.9876 (23.8481)
RW_1_ARIMA_MAPE = round(mean(RW_1_point_ARIMA_err), 4)         # 34.2087
RWD_1_ARIMA_MAPE = round(mean(RWD_1_point_ARIMA_err), 4)        # 17.7396

round(mean(CoDa_FTS_1_point_ARIMA_density_norm), 4)   # 1805.588
round(mean(CoDa_MFTS_1_point_ARIMA_density_norm), 4)  # 2213.518
round(mean(CoDa_MLFTS_1_point_ARIMA_density_norm), 4) # 1100.288
round(mean(RW_1_point_ARIMA_density_norm), 4)         # 1707.906
round(mean(RWD_1_point_ARIMA_density_norm), 4)        # 688.1755

CoDa_FTS_1_ARIMA_KLdiv = round(mean(CoDa_FTS_1_point_ARIMA_KLdiv), 4)   # 0.0413
CoDa_MFTS_1_ARIMA_KLdiv = round(mean(CoDa_MFTS_1_point_ARIMA_KLdiv), 4)  # 0.0209 (0.0533)
CoDa_MLFTS_1_ARIMA_KLdiv = round(mean(CoDa_MLFTS_1_point_ARIMA_KLdiv), 4) # 0.0335 (0.0371)
RW_1_ARIMA_KLdiv = round(mean(RW_1_point_ARIMA_KLdiv), 4)         # 0.0415
RWD_1_ARIMA_KLdiv = round(mean(RWD_1_point_ARIMA_KLdiv), 4)        # 0.0170

CoDa_FTS_1_ARIMA_JSdiv = round(mean(CoDa_FTS_1_point_ARIMA_JSdiv), 4)   # 0.0102
CoDa_MFTS_1_ARIMA_JSdiv = round(mean(CoDa_MFTS_1_point_ARIMA_JSdiv), 4)  # 0.0052 (0.0131)
CoDa_MLFTS_1_ARIMA_JSdiv = round(mean(CoDa_MLFTS_1_point_ARIMA_JSdiv), 4) # 0.0083 (0.0092)
RW_1_ARIMA_JSdiv = round(mean(RW_1_point_ARIMA_JSdiv), 4)         # 0.0102
RWD_1_ARIMA_JSdiv = round(mean(RWD_1_point_ARIMA_JSdiv), 4)        # 0.0042

CoDa_FTS_1_ARIMA_JSdiv_geo = round(mean(CoDa_FTS_1_point_ARIMA_JSdiv_geo), 4)   # 0.0103
CoDa_MFTS_1_ARIMA_JSdiv_geo = round(mean(CoDa_MFTS_1_point_ARIMA_JSdiv_geo), 4)  # 0.0052 (0.0133)
CoDa_MLFTS_1_ARIMA_JSdiv_geo = round(mean(CoDa_MLFTS_1_point_ARIMA_JSdiv_geo), 4) # 0.0084 (0.0093)
RW_1_ARIMA_JSdiv_geo = round(mean(RW_1_point_ARIMA_JSdiv_geo), 4)         # 0.0103
RWD_1_ARIMA_JSdiv_geo = round(mean(RWD_1_point_ARIMA_JSdiv_geo), 4)        # 0.0042

# 2nd population

CoDa_FTS_2_ARIMA_MAPE = round(mean(CoDa_FTS_2_point_ARIMA_err), 4)   # 42.7374
CoDa_MFTS_2_ARIMA_MAPE = round(mean(CoDa_MFTS_2_point_ARIMA_err), 4)  # 38.3179 (43.8327)
CoDa_MLFTS_2_ARIMA_MAPE = round(mean(CoDa_MLFTS_2_point_ARIMA_err), 4) # 27.4993 (32.5176)
RW_2_ARIMA_MAPE = round(mean(RW_2_point_ARIMA_err), 4)         # 43.1510 
RWD_2_ARIMA_MAPE = round(mean(RWD_2_point_ARIMA_err), 4)        # 26.2370

round(mean(CoDa_FTS_2_point_ARIMA_density_norm), 4)   # 5480.242
round(mean(CoDa_MFTS_2_point_ARIMA_density_norm), 4)  # 7802.926
round(mean(CoDa_MLFTS_2_point_ARIMA_density_norm), 4) # 2587.587
round(mean(RW_2_point_ARIMA_density_norm), 4)         # 4069.460
round(mean(RWD_2_point_ARIMA_density_norm), 4)        # 2296.811

CoDa_FTS_2_ARIMA_KLdiv = round(mean(CoDa_FTS_2_point_ARIMA_KLdiv), 4)   # 0.1304
CoDa_MFTS_2_ARIMA_KLdiv = round(mean(CoDa_MFTS_2_point_ARIMA_KLdiv), 4)  # 0.0991 (0.1549)
CoDa_MLFTS_2_ARIMA_KLdiv = round(mean(CoDa_MLFTS_2_point_ARIMA_KLdiv), 4) # 0.0549 (0.0630)
RW_2_ARIMA_KLdiv = round(mean(RW_2_point_ARIMA_KLdiv), 4)         # 0.0856
RWD_2_ARIMA_KLdiv = round(mean(RWD_2_point_ARIMA_KLdiv), 4)        # 0.0629

CoDa_FTS_2_ARIMA_JSdiv = round(mean(CoDa_FTS_2_point_ARIMA_JSdiv), 4)   # 0.0314
CoDa_MFTS_2_ARIMA_JSdiv = round(mean(CoDa_MFTS_2_point_ARIMA_JSdiv), 4)  # 0.0241 (0.0369)
CoDa_MLFTS_2_ARIMA_JSdiv = round(mean(CoDa_MLFTS_2_point_ARIMA_JSdiv), 4) # 0.0135 (0.0154)
RW_2_ARIMA_JSdiv = round(mean(RW_2_point_ARIMA_JSdiv), 4)         # 0.0208
RWD_2_ARIMA_JSdiv = round(mean(RWD_2_point_ARIMA_JSdiv), 4)        # 0.0154

CoDa_FTS_2_ARIMA_JSdiv_geo = round(mean(CoDa_FTS_2_point_ARIMA_JSdiv_geo), 4)   # 0.0326
CoDa_MFTS_2_ARIMA_JSdiv_geo = round(mean(CoDa_MFTS_2_point_ARIMA_JSdiv_geo), 4)  # 0.0247 (0.0388)
CoDa_MLFTS_2_ARIMA_JSdiv_geo = round(mean(CoDa_MLFTS_2_point_ARIMA_JSdiv_geo), 4) # 0.0137 (0.0157)
RW_2_ARIMA_JSdiv_geo = round(mean(RW_2_point_ARIMA_JSdiv_geo), 4)         # 0.0214
RWD_2_ARIMA_JSdiv_geo = round(mean(RWD_2_point_ARIMA_JSdiv_geo), 4)        # 0.0157

#######
## ETS
#######

CoDa_FTS_1_point_ETS_err = CoDa_MFTS_1_point_ETS_err = CoDa_MLFTS_1_point_ETS_err =
CoDa_FTS_1_point_ETS_density_norm = CoDa_MFTS_1_point_ETS_density_norm = CoDa_MLFTS_1_point_ETS_density_norm =
CoDa_FTS_1_point_ETS_KLdiv = CoDa_MFTS_1_point_ETS_KLdiv = CoDa_MLFTS_1_point_ETS_KLdiv =     
CoDa_FTS_1_point_ETS_JSdiv = CoDa_MFTS_1_point_ETS_JSdiv = CoDa_MLFTS_1_point_ETS_JSdiv =     
CoDa_FTS_1_point_ETS_JSdiv_geo = CoDa_MFTS_1_point_ETS_JSdiv_geo = CoDa_MLFTS_1_point_ETS_JSdiv_geo =     
    
CoDa_FTS_2_point_ETS_err = CoDa_MFTS_2_point_ETS_err = CoDa_MLFTS_2_point_ETS_err = 
CoDa_FTS_2_point_ETS_density_norm = CoDa_MFTS_2_point_ETS_density_norm = CoDa_MLFTS_2_point_ETS_density_norm =
CoDa_FTS_2_point_ETS_KLdiv = CoDa_MFTS_2_point_ETS_KLdiv = CoDa_MLFTS_2_point_ETS_KLdiv = 
CoDa_FTS_2_point_ETS_JSdiv = CoDa_MFTS_2_point_ETS_JSdiv = CoDa_MLFTS_2_point_ETS_JSdiv = 
CoDa_FTS_2_point_ETS_JSdiv_geo = CoDa_MFTS_2_point_ETS_JSdiv_geo = CoDa_MLFTS_2_point_ETS_JSdiv_geo = vector("numeric", 30)
for(iwk in 1:30)
{
    dum = CoDa_point(dat_1 = t(EW_female_pop), dat_2 = t(EW_male_pop), fh = iwk, fmethod = "ETS")
    
    ## 1st population
    
    # MAPE
    
    CoDa_FTS_1_point_ETS_err[iwk] = dum$CoDa_FTS_1_err
    CoDa_MFTS_1_point_ETS_err[iwk] = dum$CoDa_MFTS_1_err
    CoDa_MLFTS_1_point_ETS_err[iwk] = dum$CoDa_MLFTS_1_err

    # density_norm
    
    CoDa_FTS_1_point_ETS_density_norm[iwk] = dum$CoDa_FTS_1_density_norm
    CoDa_MFTS_1_point_ETS_density_norm[iwk] = dum$CoDa_MFTS_1_density_norm
    CoDa_MLFTS_1_point_ETS_density_norm[iwk] = dum$CoDa_MLFTS_1_density_norm
    
    # KLdiv
    
    CoDa_FTS_1_point_ETS_KLdiv[iwk] = dum$CoDa_FTS_1_KLdiv
    CoDa_MFTS_1_point_ETS_KLdiv[iwk] = dum$CoDa_MFTS_1_KLdiv
    CoDa_MLFTS_1_point_ETS_KLdiv[iwk] = dum$CoDa_MLFTS_1_KLdiv
    
    # JSdiv
    
    CoDa_FTS_1_point_ETS_JSdiv[iwk] = dum$CoDa_FTS_1_JSdiv
    CoDa_MFTS_1_point_ETS_JSdiv[iwk] = dum$CoDa_MFTS_1_JSdiv
    CoDa_MLFTS_1_point_ETS_JSdiv[iwk] = dum$CoDa_MLFTS_1_JSdiv
    
    # JSdiv_geo
    
    CoDa_FTS_1_point_ETS_JSdiv_geo[iwk] = dum$CoDa_FTS_1_JSdiv_geo
    CoDa_MFTS_1_point_ETS_JSdiv_geo[iwk] = dum$CoDa_MFTS_1_JSdiv_geo
    CoDa_MLFTS_1_point_ETS_JSdiv_geo[iwk] = dum$CoDa_MLFTS_1_JSdiv_geo
    
    ## 2nd population
    
    # MAPE
    
    CoDa_FTS_2_point_ETS_err[iwk] = dum$CoDa_FTS_2_err
    CoDa_MFTS_2_point_ETS_err[iwk] = dum$CoDa_MFTS_2_err
    CoDa_MLFTS_2_point_ETS_err[iwk] = dum$CoDa_MLFTS_2_err
    
    # density_norm
    
    CoDa_FTS_2_point_ETS_density_norm[iwk] = dum$CoDa_FTS_2_density_norm
    CoDa_MFTS_2_point_ETS_density_norm[iwk] = dum$CoDa_MFTS_2_density_norm
    CoDa_MLFTS_2_point_ETS_density_norm[iwk] = dum$CoDa_MLFTS_2_density_norm
    
    # KLdiv
    
    CoDa_FTS_2_point_ETS_KLdiv[iwk] = dum$CoDa_FTS_2_KLdiv
    CoDa_MFTS_2_point_ETS_KLdiv[iwk] = dum$CoDa_MFTS_2_KLdiv
    CoDa_MLFTS_2_point_ETS_KLdiv[iwk] = dum$CoDa_MLFTS_2_KLdiv
    
    # JSdiv
    
    CoDa_FTS_2_point_ETS_JSdiv[iwk] = dum$CoDa_FTS_2_JSdiv
    CoDa_MFTS_2_point_ETS_JSdiv[iwk] = dum$CoDa_MFTS_2_JSdiv
    CoDa_MLFTS_2_point_ETS_JSdiv[iwk] = dum$CoDa_MLFTS_2_JSdiv
    
    # JSdiv_geo 
    
    CoDa_FTS_2_point_ETS_JSdiv_geo[iwk] = dum$CoDa_FTS_2_JSdiv_geo
    CoDa_MFTS_2_point_ETS_JSdiv_geo[iwk] = dum$CoDa_MFTS_2_JSdiv_geo
    CoDa_MLFTS_2_point_ETS_JSdiv_geo[iwk] = dum$CoDa_MLFTS_2_JSdiv_geo
    print(iwk); rm(dum); rm(iwk)
}

# 1st population

CoDa_FTS_1_ETS_MAPE = round(mean(CoDa_FTS_1_point_ETS_err), 4)   # 33.4047
CoDa_MFTS_1_ETS_MAPE = round(mean(CoDa_MFTS_1_point_ETS_err), 4)  # 32.2143 (36.077)
CoDa_MLFTS_1_ETS_MAPE = round(mean(CoDa_MLFTS_1_point_ETS_err), 4) # 32.8372 (33.7274)

round(mean(CoDa_FTS_1_point_ETS_density_norm), 4)   # 2010.598
round(mean(CoDa_MFTS_1_point_ETS_density_norm), 4)  # 3756.545 
round(mean(CoDa_MLFTS_1_point_ETS_density_norm), 4) # 1902.170

CoDa_FTS_1_ETS_KLdiv = round(mean(CoDa_FTS_1_point_ETS_KLdiv), 4)   # 0.0319
CoDa_MFTS_1_ETS_KLdiv = round(mean(CoDa_MFTS_1_point_ETS_KLdiv), 4)  # 0.0442 (0.0833)
CoDa_MLFTS_1_ETS_KLdiv = round(mean(CoDa_MLFTS_1_point_ETS_KLdiv), 4) # 0.0539 (0.0577)

CoDa_FTS_1_ETS_JSdiv = round(mean(CoDa_FTS_1_point_ETS_JSdiv), 4)   # 0.0079 
CoDa_MFTS_1_ETS_JSdiv = round(mean(CoDa_MFTS_1_point_ETS_JSdiv), 4)  # 0.0109 (0.0202)
CoDa_MLFTS_1_ETS_JSdiv = round(mean(CoDa_MLFTS_1_point_ETS_JSdiv), 4) # 0.0133 (0.0142)

CoDa_FTS_1_ETS_JSdiv_geo = round(mean(CoDa_FTS_1_point_ETS_JSdiv_geo), 4)   # 0.0080
CoDa_MFTS_1_ETS_JSdiv_geo = round(mean(CoDa_MFTS_1_point_ETS_JSdiv_geo), 4)  # 0.0110 (0.0208)
CoDa_MLFTS_1_ETS_JSdiv_geo = round(mean(CoDa_MLFTS_1_point_ETS_JSdiv_geo), 4) # 0.0134 (0.0144)

# 2nd population

CoDa_FTS_2_ETS_MAPE = round(mean(CoDa_FTS_2_point_ETS_err), 4)   # 47.2228
CoDa_MFTS_2_ETS_MAPE = round(mean(CoDa_MFTS_2_point_ETS_err), 4)  # 52.6093 (53.8414)
CoDa_MLFTS_2_ETS_MAPE = round(mean(CoDa_MLFTS_2_point_ETS_err), 4) # 39.2112 (38.8475)

round(mean(CoDa_FTS_2_point_ETS_density_norm), 4)   # 7166.643
round(mean(CoDa_MFTS_2_point_ETS_density_norm), 4)  # 10175.370
round(mean(CoDa_MLFTS_2_point_ETS_density_norm), 4) # 3894.925

CoDa_FTS_2_ETS_KLdiv = round(mean(CoDa_FTS_2_point_ETS_KLdiv), 4)   # 0.1531
CoDa_MFTS_2_ETS_KLdiv = round(mean(CoDa_MFTS_2_point_ETS_KLdiv), 4)  # 0.1229 (0.1804)
CoDa_MLFTS_2_ETS_KLdiv = round(mean(CoDa_MLFTS_2_point_ETS_KLdiv), 4) # 0.0789 (0.0865)

CoDa_FTS_2_ETS_JSdiv = round(mean(CoDa_FTS_2_point_ETS_JSdiv), 4)   # 0.0366 
CoDa_MFTS_2_ETS_JSdiv = round(mean(CoDa_MFTS_2_point_ETS_JSdiv), 4)  # 0.0296 (0.0426)
CoDa_MLFTS_2_ETS_JSdiv = round(mean(CoDa_MLFTS_2_point_ETS_JSdiv), 4) # 0.0192 (0.021)

CoDa_FTS_2_ETS_JSdiv_geo = round(mean(CoDa_FTS_2_point_ETS_JSdiv_geo), 4)   # 0.0383
CoDa_MFTS_2_ETS_JSdiv_geo = round(mean(CoDa_MFTS_2_point_ETS_JSdiv_geo), 4)  # 0.0307 (0.0452)
CoDa_MLFTS_2_ETS_JSdiv_geo = round(mean(CoDa_MLFTS_2_point_ETS_JSdiv_geo), 4) # 0.0197 (0.0216)

################
## RWF_no_drift
################

CoDa_FTS_1_point_RW_err = CoDa_MFTS_1_point_RW_err = CoDa_MLFTS_1_point_RW_err = 
CoDa_FTS_1_point_RW_density_norm = CoDa_MFTS_1_point_RW_density_norm = CoDa_MLFTS_1_point_RW_density_norm = 
CoDa_FTS_1_point_RW_KLdiv = CoDa_MFTS_1_point_RW_KLdiv = CoDa_MLFTS_1_point_RW_KLdiv =     
CoDa_FTS_1_point_RW_JSdiv = CoDa_MFTS_1_point_RW_JSdiv = CoDa_MLFTS_1_point_RW_JSdiv =     
CoDa_FTS_1_point_RW_JSdiv_geo = CoDa_MFTS_1_point_RW_JSdiv_geo = CoDa_MLFTS_1_point_RW_JSdiv_geo =     
    
CoDa_FTS_2_point_RW_err = CoDa_MFTS_2_point_RW_err = CoDa_MLFTS_2_point_RW_err = 
CoDa_FTS_2_point_RW_density_norm = CoDa_MFTS_2_point_RW_density_norm = CoDa_MLFTS_2_point_RW_density_norm = 
CoDa_FTS_2_point_RW_KLdiv = CoDa_MFTS_2_point_RW_KLdiv = CoDa_MLFTS_2_point_RW_KLdiv = 
CoDa_FTS_2_point_RW_JSdiv = CoDa_MFTS_2_point_RW_JSdiv = CoDa_MLFTS_2_point_RW_JSdiv = 
CoDa_FTS_2_point_RW_JSdiv_geo = CoDa_MFTS_2_point_RW_JSdiv_geo = CoDa_MLFTS_2_point_RW_JSdiv_geo = vector("numeric", 30)
for(ik in 1:30)
{
    dum = CoDa_point(dat_1 = t(EW_female_pop), dat_2 = t(EW_male_pop), fh = ik, fmethod = "RWF_no_drift")
    
    ## 1st population
    
    # MAPE
    
    CoDa_FTS_1_point_RW_err[ik] = dum$CoDa_FTS_1_err
    CoDa_MFTS_1_point_RW_err[ik] = dum$CoDa_MFTS_1_err
    CoDa_MLFTS_1_point_RW_err[ik] = dum$CoDa_MLFTS_1_err

    # density_norm
    
    CoDa_FTS_1_point_RW_density_norm[ik] = dum$CoDa_FTS_1_density_norm
    CoDa_MFTS_1_point_RW_density_norm[ik] = dum$CoDa_MFTS_1_density_norm
    CoDa_MLFTS_1_point_RW_density_norm[ik] = dum$CoDa_MLFTS_1_density_norm
    
    # KLdiv

    CoDa_FTS_1_point_RW_KLdiv[ik] = dum$CoDa_FTS_1_KLdiv
    CoDa_MFTS_1_point_RW_KLdiv[ik] = dum$CoDa_MFTS_1_KLdiv
    CoDa_MLFTS_1_point_RW_KLdiv[ik] = dum$CoDa_MLFTS_1_KLdiv
    
    # JSdiv
    
    CoDa_FTS_1_point_RW_JSdiv[ik] = dum$CoDa_FTS_1_JSdiv
    CoDa_MFTS_1_point_RW_JSdiv[ik] = dum$CoDa_MFTS_1_JSdiv
    CoDa_MLFTS_1_point_RW_JSdiv[ik] = dum$CoDa_MLFTS_1_JSdiv
    
    # JSdiv_geo
    
    CoDa_FTS_1_point_RW_JSdiv_geo[ik] = dum$CoDa_FTS_1_JSdiv_geo
    CoDa_MFTS_1_point_RW_JSdiv_geo[ik] = dum$CoDa_MFTS_1_JSdiv_geo
    CoDa_MLFTS_1_point_RW_JSdiv_geo[ik] = dum$CoDa_MLFTS_1_JSdiv_geo
    
    ## 2nd population
    
    # MAPE
    
    CoDa_FTS_2_point_RW_err[ik] = dum$CoDa_FTS_2_err
    CoDa_MFTS_2_point_RW_err[ik] = dum$CoDa_MFTS_2_err
    CoDa_MLFTS_2_point_RW_err[ik] = dum$CoDa_MLFTS_2_err
    
    # density_norm
    
    CoDa_FTS_2_point_RW_density_norm[ik] = dum$CoDa_FTS_2_density_norm
    CoDa_MFTS_2_point_RW_density_norm[ik] = dum$CoDa_MFTS_2_density_norm
    CoDa_MLFTS_2_point_RW_density_norm[ik] = dum$CoDa_MLFTS_2_density_norm
    
    # KLdiv
    
    CoDa_FTS_2_point_RW_KLdiv[ik] = dum$CoDa_FTS_2_KLdiv
    CoDa_MFTS_2_point_RW_KLdiv[ik] = dum$CoDa_MFTS_2_KLdiv
    CoDa_MLFTS_2_point_RW_KLdiv[ik] = dum$CoDa_MLFTS_2_KLdiv
    
    # JSdiv
    
    CoDa_FTS_2_point_RW_JSdiv[ik] = dum$CoDa_FTS_2_JSdiv
    CoDa_MFTS_2_point_RW_JSdiv[ik] = dum$CoDa_MFTS_2_JSdiv
    CoDa_MLFTS_2_point_RW_JSdiv[ik] = dum$CoDa_MLFTS_2_JSdiv
    
    # JSdiv_geo
    
    CoDa_FTS_2_point_RW_JSdiv_geo[ik] = dum$CoDa_FTS_2_JSdiv_geo
    CoDa_MFTS_2_point_RW_JSdiv_geo[ik] = dum$CoDa_MFTS_2_JSdiv_geo
    CoDa_MLFTS_2_point_RW_JSdiv_geo[ik] = dum$CoDa_MLFTS_2_JSdiv_geo
    print(ik); rm(dum); rm(ik)
}

# 1st population

CoDa_FTS_1_RW_MAPE = round(mean(CoDa_FTS_1_point_RW_err), 4)   # 31.4225
CoDa_MFTS_1_RW_MAPE = round(mean(CoDa_MFTS_1_point_RW_err), 4)  # 33.9531 (36.9431)
CoDa_MLFTS_1_RW_MAPE = round(mean(CoDa_MLFTS_1_point_RW_err), 4) # 34.7617 (36.7725)

round(mean(CoDa_FTS_1_point_RW_density_norm), 4)   # 2842.658
round(mean(CoDa_MFTS_1_point_RW_density_norm), 4)  # 3912.211
round(mean(CoDa_MLFTS_1_point_RW_density_norm), 4) # 1921.358

CoDa_FTS_1_RW_KLdiv = round(mean(CoDa_FTS_1_point_RW_KLdiv), 4)   # 0.0695
CoDa_MFTS_1_RW_KLdiv = round(mean(CoDa_MFTS_1_point_RW_KLdiv), 4)  # 0.0412 (0.0868)
CoDa_MLFTS_1_RW_KLdiv = round(mean(CoDa_MLFTS_1_point_RW_KLdiv), 4) # 0.0539 (0.0584)

CoDa_FTS_1_RW_JSdiv = round(mean(CoDa_FTS_1_point_RW_JSdiv), 4)   # 0.0170
CoDa_MFTS_1_RW_JSdiv = round(mean(CoDa_MFTS_1_point_RW_JSdiv), 4)  # 0.0102 (0.0211)
CoDa_MLFTS_1_RW_JSdiv = round(mean(CoDa_MLFTS_1_point_RW_JSdiv), 4) # 0.0133 (0.0144)

CoDa_FTS_1_RW_JSdiv_geo = round(mean(CoDa_FTS_1_point_RW_JSdiv_geo), 4)   # 0.0174
CoDa_MFTS_1_RW_JSdiv_geo = round(mean(CoDa_MFTS_1_point_RW_JSdiv_geo), 4)  # 0.0103 (0.0217)
CoDa_MLFTS_1_RW_JSdiv_geo = round(mean(CoDa_MLFTS_1_point_RW_JSdiv_geo), 4) # 0.0134 (0.0146)

# 2nd population

CoDa_FTS_2_RW_MAPE = round(mean(CoDa_FTS_2_point_RW_err), 4)   # 47.0765
CoDa_MFTS_2_RW_MAPE = round(mean(CoDa_MFTS_2_point_RW_err), 4)  # 54.9510 (55.3208)
CoDa_MLFTS_2_RW_MAPE = round(mean(CoDa_MLFTS_2_point_RW_err), 4) # 39.7968 (39.1694)

round(mean(CoDa_FTS_2_point_RW_density_norm), 4)   # 7087.29
round(mean(CoDa_MFTS_2_point_RW_density_norm), 4)  # 10461.29
round(mean(CoDa_MLFTS_2_point_RW_density_norm), 4) # 3761.872

CoDa_FTS_2_RW_KLdiv = round(mean(CoDa_FTS_2_point_RW_KLdiv), 4)   # 0.1518
CoDa_MFTS_2_RW_KLdiv = round(mean(CoDa_MFTS_2_point_RW_KLdiv), 4)  # 0.1177 (0.1838)
CoDa_MLFTS_2_RW_KLdiv = round(mean(CoDa_MLFTS_2_point_RW_KLdiv), 4) # 0.0777 (0.0860)

CoDa_FTS_2_RW_JSdiv = round(mean(CoDa_FTS_2_point_RW_JSdiv), 4)   # 0.0363
CoDa_MFTS_2_RW_JSdiv = round(mean(CoDa_MFTS_2_point_RW_JSdiv), 4)  # 0.0285 (0.0433)
CoDa_MLFTS_2_RW_JSdiv = round(mean(CoDa_MLFTS_2_point_RW_JSdiv), 4) # 0.0189 (0.0209)

CoDa_FTS_2_RW_JSdiv_geo = round(mean(CoDa_FTS_2_point_RW_JSdiv_geo), 4)   # 0.0380
CoDa_MFTS_2_RW_JSdiv_geo = round(mean(CoDa_MFTS_2_point_RW_JSdiv_geo), 4)  # 0.0293 (0.0461)
CoDa_MLFTS_2_RW_JSdiv_geo = round(mean(CoDa_MLFTS_2_point_RW_JSdiv_geo), 4) # 0.0194 (0.0215)

#############
## RWF_drift
#############

CoDa_FTS_1_point_RWD_err = CoDa_MFTS_1_point_RWD_err = CoDa_MLFTS_1_point_RWD_err = 
CoDa_FTS_1_point_RWD_density_norm = CoDa_MFTS_1_point_RWD_density_norm = CoDa_MLFTS_1_point_RWD_density_norm = 
CoDa_FTS_1_point_RWD_KLdiv = CoDa_MFTS_1_point_RWD_KLdiv = CoDa_MLFTS_1_point_RWD_KLdiv =     
CoDa_FTS_1_point_RWD_JSdiv = CoDa_MFTS_1_point_RWD_JSdiv = CoDa_MLFTS_1_point_RWD_JSdiv =     
CoDa_FTS_1_point_RWD_JSdiv_geo = CoDa_MFTS_1_point_RWD_JSdiv_geo = CoDa_MLFTS_1_point_RWD_JSdiv_geo =     
    
CoDa_FTS_2_point_RWD_err = CoDa_MFTS_2_point_RWD_err = CoDa_MLFTS_2_point_RWD_err = 
CoDa_FTS_2_point_RWD_density_norm = CoDa_MFTS_2_point_RWD_density_norm = CoDa_MLFTS_2_point_RWD_density_norm = 
CoDa_FTS_2_point_RWD_KLdiv = CoDa_MFTS_2_point_RWD_KLdiv = CoDa_MLFTS_2_point_RWD_KLdiv = 
CoDa_FTS_2_point_RWD_JSdiv = CoDa_MFTS_2_point_RWD_JSdiv = CoDa_MLFTS_2_point_RWD_JSdiv = 
CoDa_FTS_2_point_RWD_JSdiv_geo = CoDa_MFTS_2_point_RWD_JSdiv_geo = CoDa_MLFTS_2_point_RWD_JSdiv_geo = vector("numeric", 30)
for(iwk in 1:30)
{
    dum = CoDa_point(dat_1 = t(EW_female_pop), dat_2 = t(EW_male_pop), fh = iwk, fmethod = "RWF_drift")
    
    ## 1st population
    
    # MAPE
    
    CoDa_FTS_1_point_RWD_err[iwk] = dum$CoDa_FTS_1_err
    CoDa_MFTS_1_point_RWD_err[iwk] = dum$CoDa_MFTS_1_err
    CoDa_MLFTS_1_point_RWD_err[iwk] = dum$CoDa_MLFTS_1_err
    
    # density_norm
    
    CoDa_FTS_1_point_RWD_density_norm[iwk] = dum$CoDa_FTS_1_density_norm
    CoDa_MFTS_1_point_RWD_density_norm[iwk] = dum$CoDa_MFTS_1_density_norm
    CoDa_MLFTS_1_point_RWD_density_norm[iwk] = dum$CoDa_MLFTS_1_density_norm
    
    # KLdiv
    
    CoDa_FTS_1_point_RWD_KLdiv[iwk] = dum$CoDa_FTS_1_KLdiv
    CoDa_MFTS_1_point_RWD_KLdiv[iwk] = dum$CoDa_MFTS_1_KLdiv
    CoDa_MLFTS_1_point_RWD_KLdiv[iwk] = dum$CoDa_MLFTS_1_KLdiv
    
    # JSdiv
    
    CoDa_FTS_1_point_RWD_JSdiv[iwk] = dum$CoDa_FTS_1_JSdiv
    CoDa_MFTS_1_point_RWD_JSdiv[iwk] = dum$CoDa_MFTS_1_JSdiv
    CoDa_MLFTS_1_point_RWD_JSdiv[iwk] = dum$CoDa_MLFTS_1_JSdiv
    
    # JSdiv_geo
    
    CoDa_FTS_1_point_RWD_JSdiv_geo[iwk] = dum$CoDa_FTS_1_JSdiv_geo
    CoDa_MFTS_1_point_RWD_JSdiv_geo[iwk] = dum$CoDa_MFTS_1_JSdiv_geo
    CoDa_MLFTS_1_point_RWD_JSdiv_geo[iwk] = dum$CoDa_MLFTS_1_JSdiv_geo
    
    ## 2nd population
    
    # MAPE
    
    CoDa_FTS_2_point_RWD_err[iwk] = dum$CoDa_FTS_2_err
    CoDa_MFTS_2_point_RWD_err[iwk] = dum$CoDa_MFTS_2_err
    CoDa_MLFTS_2_point_RWD_err[iwk] = dum$CoDa_MLFTS_2_err
    
    # density_norm
    
    CoDa_FTS_2_point_RWD_density_norm[iwk] = dum$CoDa_FTS_2_density_norm
    CoDa_MFTS_2_point_RWD_density_norm[iwk] = dum$CoDa_MFTS_2_density_norm
    CoDa_MLFTS_2_point_RWD_density_norm[iwk] = dum$CoDa_MLFTS_2_density_norm
    
    # KLdiv
    
    CoDa_FTS_2_point_RWD_KLdiv[iwk] = dum$CoDa_FTS_2_KLdiv
    CoDa_MFTS_2_point_RWD_KLdiv[iwk] = dum$CoDa_MFTS_2_KLdiv
    CoDa_MLFTS_2_point_RWD_KLdiv[iwk] = dum$CoDa_MLFTS_2_KLdiv
    
    # JSdiv
    
    CoDa_FTS_2_point_RWD_JSdiv[iwk] = dum$CoDa_FTS_2_JSdiv
    CoDa_MFTS_2_point_RWD_JSdiv[iwk] = dum$CoDa_MFTS_2_JSdiv
    CoDa_MLFTS_2_point_RWD_JSdiv[iwk] = dum$CoDa_MLFTS_2_JSdiv
    
    # JSdiv_geo
    
    CoDa_FTS_2_point_RWD_JSdiv_geo[iwk] = dum$CoDa_FTS_2_JSdiv_geo
    CoDa_MFTS_2_point_RWD_JSdiv_geo[iwk] = dum$CoDa_MFTS_2_JSdiv_geo
    CoDa_MLFTS_2_point_RWD_JSdiv_geo[iwk] = dum$CoDa_MLFTS_2_JSdiv_geo
    print(iwk); rm(dum); rm(iwk)
}

## 1st population

CoDa_FTS_1_RWD_MAPE = round(mean(CoDa_FTS_1_point_RWD_err), 4)   # 29.9143
CoDa_MFTS_1_RWD_MAPE = round(mean(CoDa_MFTS_1_point_RWD_err), 4)  # 18.0311 (28.3098)
CoDa_MLFTS_1_RWD_MAPE = round(mean(CoDa_MLFTS_1_point_RWD_err), 4) # 20.8329 (21.8436)

round(mean(CoDa_FTS_1_point_RWD_density_norm), 4)   # 1820.475
round(mean(CoDa_MFTS_1_point_RWD_density_norm), 4)  # 2261.276
round(mean(CoDa_MLFTS_1_point_RWD_density_norm), 4) # 880.9847

CoDa_FTS_1_RWD_KLdiv = round(mean(CoDa_FTS_1_point_RWD_KLdiv), 4)   # 0.0416
CoDa_MFTS_1_RWD_KLdiv = round(mean(CoDa_MFTS_1_point_RWD_KLdiv), 4)  # 0.0180 (0.0546)
CoDa_MLFTS_1_RWD_KLdiv = round(mean(CoDa_MLFTS_1_point_RWD_KLdiv), 4) # 0.0275 (0.0312)

CoDa_FTS_1_RWD_JSdiv = round(mean(CoDa_FTS_1_point_RWD_JSdiv), 4)   # 0.0103
CoDa_MFTS_1_RWD_JSdiv = round(mean(CoDa_MFTS_1_point_RWD_JSdiv), 4)  # 0.0045 (0.0134)
CoDa_MLFTS_1_RWD_JSdiv = round(mean(CoDa_MLFTS_1_point_RWD_JSdiv), 4) # 0.0068 (0.0077)

CoDa_FTS_1_RWD_JSdiv_geo = round(mean(CoDa_FTS_1_point_RWD_JSdiv_geo), 4)   # 0.0104
CoDa_MFTS_1_RWD_JSdiv_geo = round(mean(CoDa_MFTS_1_point_RWD_JSdiv_geo), 4)  # 0.0045 (0.0137)
CoDa_MLFTS_1_RWD_JSdiv_geo = round(mean(CoDa_MLFTS_1_point_RWD_JSdiv_geo), 4) # 0.0069 (0.0078)

## 2nd population

CoDa_FTS_2_RWD_MAPE = round(mean(CoDa_FTS_2_point_RWD_err), 4)   # 43.1748
CoDa_MFTS_2_RWD_MAPE = round(mean(CoDa_MFTS_2_point_RWD_err), 4)  # 38.0957 (44.1254)
CoDa_MLFTS_2_RWD_MAPE = round(mean(CoDa_MLFTS_2_point_RWD_err), 4) # 27.0489 (32.1465)

round(mean(CoDa_FTS_2_point_RWD_density_norm), 4)   # 5330.345
round(mean(CoDa_MFTS_2_point_RWD_density_norm), 4)  # 7904.559
round(mean(CoDa_MLFTS_2_point_RWD_density_norm), 4) # 2561.725

CoDa_FTS_2_RWD_KLdiv = round(mean(CoDa_FTS_2_point_RWD_KLdiv), 4)   # 0.1269
CoDa_MFTS_2_RWD_KLdiv = round(mean(CoDa_MFTS_2_point_RWD_KLdiv), 4)  # 0.0923 (0.1560)
CoDa_MLFTS_2_RWD_KLdiv = round(mean(CoDa_MLFTS_2_point_RWD_KLdiv), 4) # 0.0529 (0.0614)

CoDa_FTS_2_RWD_JSdiv = round(mean(CoDa_FTS_2_point_RWD_JSdiv), 4)   # 0.0306
CoDa_MFTS_2_RWD_JSdiv = round(mean(CoDa_MFTS_2_point_RWD_JSdiv), 4)  # 0.0225 (0.0371)
CoDa_MLFTS_2_RWD_JSdiv = round(mean(CoDa_MLFTS_2_point_RWD_JSdiv), 4) # 0.0130 (0.0150)

CoDa_FTS_2_RWD_JSdiv_geo = round(mean(CoDa_FTS_2_point_RWD_JSdiv_geo), 4)   # 0.0317
CoDa_MFTS_2_RWD_JSdiv_geo = round(mean(CoDa_MFTS_2_point_RWD_JSdiv_geo), 4)  # 0.0230 (0.0391)
CoDa_MLFTS_2_RWD_JSdiv_geo = round(mean(CoDa_MLFTS_2_point_RWD_JSdiv_geo), 4) # 0.0132 (0.0153)

# summary

require(xtable)
results = rbind(c(CoDa_FTS_1_ARIMA_MAPE, CoDa_FTS_1_ARIMA_KLdiv, CoDa_FTS_1_ARIMA_JSdiv, CoDa_FTS_1_ARIMA_JSdiv_geo,
  CoDa_FTS_2_ARIMA_MAPE, CoDa_FTS_2_ARIMA_KLdiv, CoDa_FTS_2_ARIMA_JSdiv, CoDa_FTS_2_ARIMA_JSdiv_geo),
c(CoDa_FTS_1_ETS_MAPE, CoDa_FTS_1_ETS_KLdiv, CoDa_FTS_1_ETS_JSdiv, CoDa_FTS_1_ETS_JSdiv_geo,
  CoDa_FTS_2_ETS_MAPE, CoDa_FTS_2_ETS_KLdiv, CoDa_FTS_2_ETS_JSdiv, CoDa_FTS_2_ETS_JSdiv_geo),
c(CoDa_FTS_1_RW_MAPE, CoDa_FTS_1_RW_KLdiv, CoDa_FTS_1_RW_JSdiv, CoDa_FTS_1_RW_JSdiv_geo,
  CoDa_FTS_2_RW_MAPE, CoDa_FTS_2_RW_KLdiv, CoDa_FTS_2_RW_JSdiv, CoDa_FTS_2_RW_JSdiv_geo),
c(CoDa_FTS_1_RWD_MAPE, CoDa_FTS_1_RWD_KLdiv, CoDa_FTS_1_RWD_JSdiv, CoDa_FTS_1_RWD_JSdiv_geo,
  CoDa_FTS_2_RWD_MAPE, CoDa_FTS_2_RWD_KLdiv, CoDa_FTS_2_RWD_JSdiv, CoDa_FTS_2_RWD_JSdiv_geo),

c(CoDa_MFTS_1_ARIMA_MAPE, CoDa_MFTS_1_ARIMA_KLdiv, CoDa_MFTS_1_ARIMA_JSdiv, CoDa_MFTS_1_ARIMA_JSdiv_geo,
  CoDa_MFTS_2_ARIMA_MAPE, CoDa_MFTS_2_ARIMA_KLdiv, CoDa_MFTS_2_ARIMA_JSdiv, CoDa_MFTS_2_ARIMA_JSdiv_geo),
c(CoDa_MFTS_1_ETS_MAPE, CoDa_MFTS_1_ETS_KLdiv, CoDa_MFTS_1_ETS_JSdiv, CoDa_MFTS_1_ETS_JSdiv_geo,
  CoDa_MFTS_2_ETS_MAPE, CoDa_MFTS_2_ETS_KLdiv, CoDa_MFTS_2_ETS_JSdiv, CoDa_MFTS_2_ETS_JSdiv_geo),
c(CoDa_MFTS_1_RW_MAPE, CoDa_MFTS_1_RW_KLdiv, CoDa_MFTS_1_RW_JSdiv, CoDa_MFTS_1_RW_JSdiv_geo,
  CoDa_MFTS_2_RW_MAPE, CoDa_MFTS_2_RW_KLdiv, CoDa_MFTS_2_RW_JSdiv, CoDa_MFTS_2_RW_JSdiv_geo),
c(CoDa_MFTS_1_RWD_MAPE, CoDa_MFTS_1_RWD_KLdiv, CoDa_MFTS_1_RWD_JSdiv, CoDa_MFTS_1_RWD_JSdiv_geo,
  CoDa_MFTS_2_RWD_MAPE, CoDa_MFTS_2_RWD_KLdiv, CoDa_MFTS_2_RWD_JSdiv, CoDa_MFTS_2_RWD_JSdiv_geo),

c(CoDa_MLFTS_1_ARIMA_MAPE, CoDa_MLFTS_1_ARIMA_KLdiv, CoDa_MLFTS_1_ARIMA_JSdiv, CoDa_MLFTS_1_ARIMA_JSdiv_geo,
  CoDa_MLFTS_2_ARIMA_MAPE, CoDa_MLFTS_2_ARIMA_KLdiv, CoDa_MLFTS_2_ARIMA_JSdiv, CoDa_MLFTS_2_ARIMA_JSdiv_geo),
c(CoDa_MLFTS_1_ETS_MAPE, CoDa_MLFTS_1_ETS_KLdiv, CoDa_MLFTS_1_ETS_JSdiv, CoDa_MLFTS_1_ETS_JSdiv_geo,
  CoDa_MLFTS_2_ETS_MAPE, CoDa_MLFTS_2_ETS_KLdiv, CoDa_MLFTS_2_ETS_JSdiv, CoDa_MLFTS_2_ETS_JSdiv_geo),
c(CoDa_MLFTS_1_RW_MAPE, CoDa_MLFTS_1_RW_KLdiv, CoDa_MLFTS_1_RW_JSdiv, CoDa_MLFTS_1_RW_JSdiv_geo,
  CoDa_MLFTS_2_RW_MAPE, CoDa_MLFTS_2_RW_KLdiv, CoDa_MLFTS_2_RW_JSdiv, CoDa_MLFTS_2_RW_JSdiv_geo),
c(CoDa_MLFTS_1_RWD_MAPE, CoDa_MLFTS_1_RWD_KLdiv, CoDa_MLFTS_1_RWD_JSdiv, CoDa_MLFTS_1_RWD_JSdiv_geo,
  CoDa_MLFTS_2_RWD_MAPE, CoDa_MLFTS_2_RWD_KLdiv, CoDa_MLFTS_2_RWD_JSdiv, CoDa_MLFTS_2_RWD_JSdiv_geo),

c(RW_1_ARIMA_MAPE, RW_1_ARIMA_KLdiv, RW_1_ARIMA_JSdiv, RW_1_ARIMA_JSdiv_geo,
    RW_2_ARIMA_MAPE, RW_2_ARIMA_KLdiv, RW_2_ARIMA_JSdiv, RW_2_ARIMA_JSdiv_geo),

c(RWD_1_ARIMA_MAPE, RWD_1_ARIMA_KLdiv, RWD_1_ARIMA_JSdiv, RWD_1_ARIMA_JSdiv_geo,
  RWD_2_ARIMA_MAPE, RWD_2_ARIMA_KLdiv, RWD_2_ARIMA_JSdiv, RWD_2_ARIMA_JSdiv_geo))

rownames(results) = c("ARIMA", "ETS", "RW", "RWD", 
                      "ARIMA", "ETS", "RW", "RWD", 
                      "ARIMA", "ETS", "RW", "RWD", 
                      "RW", "RWD")



########################
# Figures in Appendices
########################

## Female (RWD)

# MAPE

plot(1:30, CoDa_FTS_1_point_RWD_err, type = "l", xlab = "", col = 1,
     ylim = c(8, 80), ylab = "MAPE", main = "Female", lty = 1)
lines(1:30, CoDa_MFTS_1_point_RWD_err, type = "l", col = 2, lty = 2)
lines(1:30, CoDa_MLFTS_1_point_RWD_err, type = "l", col = 3, lty = 3)
lines(1:30, RW_1_point_ARIMA_err, type = "l", col = 4, lty = 4)
lines(1:30, RWD_1_point_ARIMA_err, type = "l", col = 5, lty = 5)
lines(1:30, test_eval_female_err[,4], type = "l", col = 6, lty = 6)
legend("topleft", c("Univariate functional time series method", "Multivariate functional time series method", 
                    "Multilevel functional time series method", "RW", "RWD", "SES"), 
       col = 1:6, lty = 1:6, cex = 0.6)

# KLD

plot(1:30, CoDa_FTS_1_point_RWD_KLdiv, type = "l", xlab = "", col = 1,
     ylab = "Kullback-Leibler divergence", lty = 1, ylim = c(0, 0.1))
lines(1:30, CoDa_MFTS_1_point_RWD_KLdiv, type = "l", col = 2, lty = 2)
lines(1:30, CoDa_MLFTS_1_point_RWD_KLdiv, type = "l", col = 3, lty = 3)
lines(1:30, RW_1_point_ARIMA_KLdiv, type = "l", col = 4, lty = 4)
lines(1:30, RWD_1_point_ARIMA_KLdiv, type = "l", col = 5, lty = 5)
lines(1:30, test_eval_female_err[,1], type = "l", col = 6, lty = 6)

# JSD (simple mean)

plot(1:30, CoDa_FTS_1_point_RWD_JSdiv, type = "l", xlab = "", col = 1,
     ylab = "Jensen-Shannon divergence (simple mean)", lty = 1, ylim = c(0, 0.025))
lines(1:30, CoDa_MFTS_1_point_RWD_JSdiv, type = "l", col = 2, lty = 2)
lines(1:30, CoDa_MLFTS_1_point_RWD_JSdiv, type = "l", col = 3, lty = 3)
lines(1:30, RW_1_point_ARIMA_JSdiv, type = "l", col = 4, lty = 4)
lines(1:30, RWD_1_point_ARIMA_JSdiv, type = "l", col = 5, lty = 5)
lines(1:30, test_eval_female_err[,2], type = "l", col = 6, lty = 6)

# JSD (geometric mean)

plot(1:30, CoDa_FTS_1_point_RWD_JSdiv_geo, type = "l", xlab = "Forecast horizon",
     col = 1, ylab = "Jensen-Shannon divergence (geometric mean)", lty = 1, ylim = c(0, 0.025))
lines(1:30, CoDa_MFTS_1_point_RWD_JSdiv_geo, type = "l", col = 2, lty = 2)
lines(1:30, CoDa_MLFTS_1_point_RWD_JSdiv_geo, type = "l", col = 3, lty = 3)
lines(1:30, RW_1_point_ARIMA_JSdiv_geo, type = "l", col = 4, lty = 4)
lines(1:30, RWD_1_point_ARIMA_JSdiv_geo, type = "l", col = 5, lty = 5)
lines(1:30, test_eval_female_err[,3], type = "l", col = 6, lty = 6)

## Male (RWD)

# MAPE

plot(1:30, CoDa_FTS_2_point_RWD_err, type = "l", xlab = "", ylim = c(8, 80), ylab = "", main = "Male")
lines(1:30, CoDa_MFTS_2_point_RWD_err, type = "l", col = 2, lty = 2)
lines(1:30, CoDa_MLFTS_2_point_RWD_err, type = "l", col = 3, lty = 3)
lines(1:30, RW_2_point_ARIMA_err, type = "l", col = 4, lty = 4)
lines(1:30, RWD_2_point_ARIMA_err, type = "l", col = 5, lty = 5)
lines(1:30, test_eval_male_err[,4], type = "l", col = 6, lty = 6)

# KLD

plot(1:30, CoDa_FTS_2_point_RWD_KLdiv, type = "l", xlab = "", ylim = c(0, 0.28), ylab = "")
lines(1:30, CoDa_MFTS_2_point_RWD_KLdiv, type = "l", col = 2, lty = 2)
lines(1:30, CoDa_MLFTS_2_point_RWD_KLdiv, type = "l", col = 3, lty = 3)
lines(1:30, RW_2_point_ARIMA_KLdiv, type = "l", col = 4, lty = 4)
lines(1:30, RWD_2_point_ARIMA_KLdiv, type = "l", col = 5, lty = 5)
lines(1:30, test_eval_male_err[,1], type = "l", col = 6, lty = 6)

# JSD (simple mean)

plot(1:30, CoDa_FTS_2_point_RWD_JSdiv, type = "l", xlab = "", ylim = c(0, 0.065), ylab = "")
lines(1:30, CoDa_MFTS_2_point_RWD_JSdiv, type = "l", col = 2, lty = 2)
lines(1:30, CoDa_MLFTS_2_point_RWD_JSdiv, type = "l", col = 3, lty = 3)
lines(1:30, RW_2_point_ARIMA_JSdiv, type = "l", col = 4, lty = 4)
lines(1:30, RWD_2_point_ARIMA_JSdiv, type = "l", col = 5, lty = 5)
lines(1:30, test_eval_male_err[,2], type = "l", col = 6, lty = 6)

# JSD (geometric mean)

plot(1:30, CoDa_FTS_2_point_RWD_JSdiv_geo, type = "l", xlab = "Forecast horizon",
     ylim = c(0, 0.065), ylab = "")
lines(1:30, CoDa_MFTS_2_point_RWD_JSdiv_geo, type = "l", col = 2, lty = 2)
lines(1:30, CoDa_MLFTS_2_point_RWD_JSdiv_geo, type = "l", col = 3, lty = 3)
lines(1:30, RW_2_point_ARIMA_JSdiv_geo, type = "l", col = 4, lty = 4)
lines(1:30, RWD_2_point_ARIMA_JSdiv_geo, type = "l", col = 5, lty = 5)
lines(1:30, test_eval_male_err[,3], type = "l", col = 6, lty = 6)

