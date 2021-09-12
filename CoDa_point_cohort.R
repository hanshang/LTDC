#################
# load data
#################
female.int = t(as.matrix(read.table('uk_female_pop_complete.txt')))
male.int = t(as.matrix(read.table('uk_male_pop_complete.txt')))

#replace the 0 death count with 0.01
death0.f = which(is.na(female.int) | female.int==0)
death0.m = which(is.na(male.int) | male.int==0)

EW_female_pop = replace(female.int,death0.f,0.01) #total is 130 years
EW_male_pop = replace(male.int,death0.m,0.01)



#################
# load R package
#################

require(psych)
require(flexmix)
require(ftsa)
require(meboot)
require(pracma)
mape = ftsa:::mape

source("lib/CoDa_naive_fun.R")
source("lib/CoDa_model_fitting_cohort.R")
source("lib/product_ratio_func.R")



CoDa_point <- function(dat_1, dat_2, fh, fmethod = c("ARIMA", "ETS", "RWF_no_drift", "RWF_drift"))
{
    # dat: original data matrix (n by p)
    # fh: forecast horizon
    # fmethod: forecasting method
    
    Tt = dim(dat_1)[1]
    train = Tt-30-1
    
    fore_h30_FTS_1 = fore_h30_MFTS_1 = fore_h30_MLFTS_1 = 
        fore_h30_FTS_2 = fore_h30_MFTS_2 = fore_h30_MLFTS_2 = matrix(NA, ncol(dat_1), (31 - fh))
    
    fore_h30_FTS_1_KLdiv = fore_h30_FTS_2_KLdiv = 
        fore_h30_MFTS_1_KLdiv = fore_h30_MFTS_2_KLdiv = 
        fore_h30_MLFTS_1_KLdiv = fore_h30_MLFTS_2_KLdiv = 
        RW_1_KLdiv = RW_2_KLdiv = 
        RWD_1_KLdiv = RWD_2_KLdiv = matrix(NA, (31 - fh), 2)
    
    fore_h30_FTS_1_mape = fore_h30_FTS_1_JSdiv = fore_h30_FTS_1_JSdiv_geo = 
        fore_h30_FTS_2_mape = fore_h30_FTS_2_JSdiv = fore_h30_FTS_2_JSdiv_geo =    
        
       fore_h30_MFTS_1_mape = fore_h30_MFTS_1_JSdiv = fore_h30_MFTS_1_JSdiv_geo = 
        fore_h30_MFTS_2_mape = fore_h30_MFTS_2_JSdiv = fore_h30_MFTS_2_JSdiv_geo =     
        
        fore_h30_MLFTS_1_mape = fore_h30_MLFTS_1_JSdiv = fore_h30_MLFTS_1_JSdiv_geo = 
        fore_h30_MLFTS_2_mape = fore_h30_MLFTS_2_JSdiv = fore_h30_MLFTS_2_JSdiv_geo = 
        
        RW_1_mape = RW_1_JSdiv = RW_1_JSdiv_geo = 
        RW_2_mape = RW_2_JSdiv = RW_2_JSdiv_geo = 
        
        RWD_1_mape = RWD_1_JSdiv = RWD_1_JSdiv_geo = 
        RWD_2_mape = RWD_2_JSdiv = RWD_2_JSdiv_geo = vector("numeric", (31 - fh))
    for(ik in 1:(31 - fh))
    {
        #######
        ## FTS 
        #######
        
        R2_female = R_square_fit(dat_1 = dat_1[1:(train+ik),], dat_2 = dat_2[1:(train+ik),],
                                 fh = fh, modeling_method = "FTS",
                                 forecasting_method = fmethod)
        fore_h30_FTS_1[,ik] = R2_female$fore_count_1[,fh]
        
        # MAPE (1st population)
        
        fore_h30_FTS_1_mape[ik] = mape(forecast = t(fore_h30_FTS_1[,ik]), true = dat_1[(train+fh+ik),])
        
       
        # KLdiv
        
        dat = cbind(true = dat_1[(train+fh+ik),], forecast = as.numeric(t(fore_h30_FTS_1[,ik])))
        colnames(dat) = c("True", "Estimate")
        fore_h30_FTS_1_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)        
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_1[(train+fh+ik),], forecast = as.numeric(t(fore_h30_FTS_1[,ik])))
        M = rowMeans(dat)
        P_M = cbind(dat_1[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_FTS_1[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_FTS_1_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_1[(train+fh+ik),], forecast = as.numeric(t(fore_h30_FTS_1[,ik])))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_1[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_FTS_1[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_FTS_1_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        
        # MAPE (2nd population)
        
        fore_h30_FTS_2[,ik] = R2_female$fore_count_2[,fh]
        fore_h30_FTS_2_mape[ik] = mape(forecast = t(fore_h30_FTS_2[,ik]), true = dat_2[(train+fh+ik),])
        
       
        # KLdiv
        
        dat = cbind(true = dat_2[(train+fh+ik),], forecast = as.numeric(t(fore_h30_FTS_2[,ik])))
        colnames(dat) = c("True", "Estimate")
        fore_h30_FTS_2_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_2[(train+fh+ik),], forecast = as.numeric(t(fore_h30_FTS_2[,ik])))
        M = rowMeans(dat)
        P_M = cbind(dat_2[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_FTS_2[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_FTS_2_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_2[(train+fh+ik),], forecast = as.numeric(t(fore_h30_FTS_2[,ik])))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_2[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_FTS_2[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_FTS_2_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        rm(R2_female)
        
        ########
        ## MFTS
        ########
        
        R2_female = R_square_fit(dat_1 = dat_1[1:(train+ik),], dat_2 = dat_2[1:(train+ik),], 
                                 fh = fh, modeling_method = "MFTS",
                                 forecasting_method = fmethod)
        fore_h30_MFTS_1[,ik] = R2_female$fore_count_1[,fh]
        
        # MAPE (1st population)
        
        fore_h30_MFTS_1_mape[ik] = mape(forecast = t(fore_h30_MFTS_1[,ik]), true = dat_1[(train+fh+ik),])
        
        
        
        # KLdiv
        
        dat = cbind(true = dat_1[(train+fh+ik),], forecast = as.numeric(t(fore_h30_MFTS_1[,ik])))
        colnames(dat) = c("True", "Estimate")      
        fore_h30_MFTS_1_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_1[(train+fh+ik),], forecast = as.numeric(t(fore_h30_MFTS_1[,ik])))
        M = rowMeans(dat)
        P_M = cbind(dat_1[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_MFTS_1[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")        
        fore_h30_MFTS_1_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_1[(train+fh+ik),], forecast = as.numeric(t(fore_h30_MFTS_1[,ik])))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_1[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_MFTS_1[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_MFTS_1_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        fore_h30_MFTS_2[,ik] = R2_female$fore_count_2[,fh]
        
        # MAPE (2nd population)
        
        fore_h30_MFTS_2_mape[ik] = mape(forecast = t(fore_h30_MFTS_2[,ik]), true = dat_2[(train+fh+ik),])
        
        
        
        # KLdiv
        
        dat = cbind(true = dat_2[(train+fh+ik),], forecast = as.numeric(t(fore_h30_MFTS_2[,ik])))
        colnames(dat) = c("True", "Estimate")
        fore_h30_MFTS_2_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_2[(train+fh+ik),], forecast = as.numeric(t(fore_h30_MFTS_2[,ik])))
        M = rowMeans(dat)
        P_M = cbind(dat_2[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_MFTS_2[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_MFTS_2_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_2[(train+fh+ik),], forecast = as.numeric(t(fore_h30_MFTS_2[,ik])))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_2[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_MFTS_2[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_MFTS_2_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        rm(R2_female)
        
        #########
        ## MLFTS
        #########
        
        R2_female = R_square_fit(dat_1 = dat_1[1:(train+ik),], dat_2 = dat_2[1:(train+ik),], 
                                 fh = fh, modeling_method = "MLFTS",
                                 forecasting_method = fmethod)
        fore_h30_MLFTS_1[,ik] = R2_female$fore_count_1[,fh]
        
        # MAPE (1st population)
        
        fore_h30_MLFTS_1_mape[ik] = mape(forecast = t(fore_h30_MLFTS_1[,ik]), true = dat_1[(train+fh+ik),])
        
        
        
        # KLdiv
        
        dat = cbind(true = dat_1[(train+fh+ik),], forecast = as.numeric(t(fore_h30_MLFTS_1[,ik])))
        colnames(dat) = c("True", "Estimate")
        fore_h30_MLFTS_1_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_1[(train+fh+ik),], forecast = as.numeric(t(fore_h30_MLFTS_1[,ik])))
        M = rowMeans(dat)
        P_M = cbind(dat_1[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_MLFTS_1[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_MLFTS_1_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_1[(train+fh+ik),], forecast = as.numeric(t(fore_h30_MLFTS_1[,ik])))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_1[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_MLFTS_1[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_MLFTS_1_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        
        fore_h30_MLFTS_2[,ik] = R2_female$fore_count_2[,fh]
        
        # MAPE (2nd population)
        
        fore_h30_MLFTS_2_mape[ik] = mape(forecast = t(fore_h30_MLFTS_2[,ik]), true = dat_2[(train+fh+ik),])
        
       
        
        # KLdiv
        
        dat = cbind(true = dat_2[(train+fh+ik),], forecast = as.numeric(t(fore_h30_MLFTS_2[,ik])))
        colnames(dat) = c("True", "Estimate")
        fore_h30_MLFTS_2_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_2[(train+fh+ik),], forecast = as.numeric(t(fore_h30_MLFTS_2[,ik])))
        M = rowMeans(dat)
        P_M = cbind(dat_2[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_MLFTS_2[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_MLFTS_2_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_2[(train+fh+ik),], forecast = as.numeric(t(fore_h30_MLFTS_2[,ik])))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_2[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(fore_h30_MLFTS_2[,ik])), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        fore_h30_MLFTS_2_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M)
        
        rm(R2_female)
        
        ####################
        ## CoDa (naive RWD)
        ####################
        
        R2_female = naive_fun(dat = dat_1[1:(train+ik),], fh = fh, drift_term = TRUE, level_ci = 80)$fore_count[,fh]
        
        # MAPE (1st population)
        
        RWD_1_mape[ik] = mape(forecast = t(R2_female), true = dat_1[(train+fh+ik),])        
        
      
        
        # KLdiv
        
        dat = cbind(true = dat_1[(train+fh+ik),], forecast = as.numeric(t(R2_female)))
        colnames(dat) = c("True", "Estimate")
        RWD_1_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_1[(train+fh+ik),], forecast = as.numeric(t(R2_female)))
        M = rowMeans(dat)
        P_M = cbind(dat_1[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(R2_female)), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        RWD_1_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M); rm(dat)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_1[(train+fh+ik),], forecast = as.numeric(t(R2_female)))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_1[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(R2_female)), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        RWD_1_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M); rm(dat); rm(R2_female)
        
        
        R2_female = naive_fun(dat = dat_2[1:(train+ik),], fh = fh, drift_term = TRUE, level_ci = 80)$fore_count[,fh]
        
        # MAPE (2nd population)
        
        RWD_2_mape[ik] = mape(forecast = t(R2_female), true = dat_2[(train+fh+ik),])        
        
       
        
        # KLdiv
        
        dat = cbind(true = dat_2[(train+fh+ik),], forecast = as.numeric(t(R2_female)))
        colnames(dat) = c("True", "Estimate")
        RWD_2_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_2[(train+fh+ik),], forecast = as.numeric(t(R2_female)))
        M = rowMeans(dat)
        P_M = cbind(dat_2[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(R2_female)), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        RWD_2_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M); rm(dat)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_2[(train+fh+ik),], forecast = as.numeric(t(R2_female)))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_2[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(R2_female)), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        RWD_2_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M); rm(dat); rm(R2_female)
        
        ##################
        ## CoDa (naive RW)
        ##################
        
        R2_female = naive_fun(dat = dat_1[1:(train+ik),], fh = fh, drift_term = FALSE, level_ci = 80)$fore_count[,fh]
        
        # MAPE (1st population)
        
        RW_1_mape[ik] = mape(forecast = t(R2_female), true = dat_1[(train+fh+ik),])        
        
        
        
        # KLdiv
        
        dat = cbind(true = dat_1[(train+fh+ik),], forecast = as.numeric(t(R2_female)))
        colnames(dat) = c("True", "Estimate")
        RW_1_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_1[(train+fh+ik),], forecast = as.numeric(t(R2_female)))
        M = rowMeans(dat)
        P_M = cbind(dat_1[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(R2_female)), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        RW_1_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M); rm(dat)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_1[(train+fh+ik),], forecast = as.numeric(t(R2_female)))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_1[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(R2_female)), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        RW_1_JSdiv_geo[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M); rm(dat); rm(R2_female)
        
        
        R2_female = naive_fun(dat = dat_2[1:(train+ik),], fh = fh, drift_term = FALSE, level_ci = 80)$fore_count[,fh]
        
        # MAPE (2nd population)
        
        RW_2_mape[ik] = mape(forecast = t(R2_female), true = dat_2[(train+fh+ik),])        
        
        
        
        # KLdiv
        
        dat = cbind(true = dat_2[(train+fh+ik),], forecast = as.numeric(t(R2_female)))
        colnames(dat) = c("True", "Estimate")
        RW_2_KLdiv[ik,] = as.numeric(KLdiv(dat, eps = 1e-16))[2:3]
        rm(dat)
        
        # JSdiv(simple)
        
        dat = cbind(true = dat_2[(train+fh+ik),], forecast = as.numeric(t(R2_female)))
        M = rowMeans(dat)
        P_M = cbind(dat_2[(train+fh+ik),], M)
        E_M = cbind(as.numeric(t(R2_female)), M)
        colnames(E_M) = colnames(P_M) = c("True", "M")
        RW_2_JSdiv[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
        rm(M); rm(P_M); rm(E_M); rm(dat)
        
        # JSdiv(geo)
        
        dat = cbind(true = dat_2[(train+fh+ik),], forecast = as.numeric(t(R2_female)))
        M = apply(dat, 1, geometric.mean)
        P_M = cbind(dat_2[(train+fh+ik),], M)
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
    
    RW_1_point_ARIMA_KLdiv = RWD_1_point_ARIMA_KLdiv = CoDa_FTS_1_point_ARIMA_KLdiv = 
    CoDa_MFTS_1_point_ARIMA_KLdiv = CoDa_MLFTS_1_point_ARIMA_KLdiv = 
    
    RW_1_point_ARIMA_JSdiv = RWD_1_point_ARIMA_JSdiv = CoDa_FTS_1_point_ARIMA_JSdiv = 
    CoDa_MFTS_1_point_ARIMA_JSdiv = CoDa_MLFTS_1_point_ARIMA_JSdiv = 
    
    RW_1_point_ARIMA_JSdiv_geo = RWD_1_point_ARIMA_JSdiv_geo = CoDa_FTS_1_point_ARIMA_JSdiv_geo = 
    CoDa_MFTS_1_point_ARIMA_JSdiv_geo = CoDa_MLFTS_1_point_ARIMA_JSdiv_geo = 
    
    RW_2_point_ARIMA_err = RWD_2_point_ARIMA_err = CoDa_FTS_2_point_ARIMA_err = 
    CoDa_MFTS_2_point_ARIMA_err = CoDa_MLFTS_2_point_ARIMA_err = 
    
    
    RW_2_point_ARIMA_KLdiv = RWD_2_point_ARIMA_KLdiv = CoDa_FTS_2_point_ARIMA_KLdiv = 
    CoDa_MFTS_2_point_ARIMA_KLdiv = CoDa_MLFTS_2_point_ARIMA_KLdiv = 
    
    RW_2_point_ARIMA_JSdiv = RWD_2_point_ARIMA_JSdiv = CoDa_FTS_2_point_ARIMA_JSdiv = 
    CoDa_MFTS_2_point_ARIMA_JSdiv = CoDa_MLFTS_2_point_ARIMA_JSdiv = 
    
    RW_2_point_ARIMA_JSdiv_geo = RWD_2_point_ARIMA_JSdiv_geo = CoDa_FTS_2_point_ARIMA_JSdiv_geo = 
    CoDa_MFTS_2_point_ARIMA_JSdiv_geo = CoDa_MLFTS_2_point_ARIMA_JSdiv_geo = vector("numeric", 30)

for(ik in 1:30)
{
    dum = CoDa_point(dat_1 = EW_female_pop, dat_2 = EW_male_pop, fh = ik, fmethod = "ARIMA")
    
    ## 1st population
    
    # MAPE
    
    RW_1_point_ARIMA_err[ik] = dum$RW_1_err
    RWD_1_point_ARIMA_err[ik] = dum$RWD_1_err
    CoDa_FTS_1_point_ARIMA_err[ik] = dum$CoDa_FTS_1_err
    CoDa_MFTS_1_point_ARIMA_err[ik] = dum$CoDa_MFTS_1_err
    CoDa_MLFTS_1_point_ARIMA_err[ik] = dum$CoDa_MLFTS_1_err
    
    
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

round(mean(CoDa_FTS_1_point_ARIMA_err), 4)   # 46.77
round(mean(CoDa_MFTS_1_point_ARIMA_err), 4)  # 51.9685
round(mean(CoDa_MLFTS_1_point_ARIMA_err), 4) # 43.3679
round(mean(RW_1_point_ARIMA_err), 4)         # 53.9087
round(mean(RWD_1_point_ARIMA_err), 4)        # 34.5075


round(mean(CoDa_FTS_1_point_ARIMA_KLdiv), 4)   # 0.0298
round(mean(CoDa_MFTS_1_point_ARIMA_KLdiv), 4)  # 0.0288
round(mean(CoDa_MLFTS_1_point_ARIMA_KLdiv), 4) # 0.0405
round(mean(RW_1_point_ARIMA_KLdiv), 4)         # 0.0853
round(mean(RWD_1_point_ARIMA_KLdiv), 4)        # 0.0237

round(mean(CoDa_FTS_1_point_ARIMA_JSdiv), 4)   # 0.0072
round(mean(CoDa_MFTS_1_point_ARIMA_JSdiv), 4)  # 0.0069
round(mean(CoDa_MLFTS_1_point_ARIMA_JSdiv), 4) # 0.0098
round(mean(RW_1_point_ARIMA_JSdiv), 4)         # 0.0204
round(mean(RWD_1_point_ARIMA_JSdiv), 4)        # 0.0058

round(mean(CoDa_FTS_1_point_ARIMA_JSdiv_geo), 4)   # 0.0076
round(mean(CoDa_MFTS_1_point_ARIMA_JSdiv_geo), 4)  # 0.0073
round(mean(CoDa_MLFTS_1_point_ARIMA_JSdiv_geo), 4) # 0.0102
round(mean(RW_1_point_ARIMA_JSdiv_geo), 4)         # 0.0214
round(mean(RWD_1_point_ARIMA_JSdiv_geo), 4)        # 0.006

# 2nd population

round(mean(CoDa_FTS_2_point_ARIMA_err), 4)   # 58.3525
round(mean(CoDa_MFTS_2_point_ARIMA_err), 4)  # 48.9508
round(mean(CoDa_MLFTS_2_point_ARIMA_err), 4) # 46.1943
round(mean(RW_2_point_ARIMA_err), 4)         # 53.9087
round(mean(RWD_2_point_ARIMA_err), 4)        # 36.1416

round(mean(CoDa_FTS_2_point_ARIMA_KLdiv), 4)   # 0.0701
round(mean(CoDa_MFTS_2_point_ARIMA_KLdiv), 4)  # 0.0529
round(mean(CoDa_MLFTS_2_point_ARIMA_KLdiv), 4) # 0.0394
round(mean(RW_2_point_ARIMA_KLdiv), 4)         # 0.0853
round(mean(RWD_2_point_ARIMA_KLdiv), 4)        # 0.0183

round(mean(CoDa_FTS_2_point_ARIMA_JSdiv), 4)   # 0.0154
round(mean(CoDa_MFTS_2_point_ARIMA_JSdiv), 4)  #  0.013
round(mean(CoDa_MLFTS_2_point_ARIMA_JSdiv), 4) #  0.0096
round(mean(RW_2_point_ARIMA_JSdiv), 4)         # 0.0204
round(mean(RWD_2_point_ARIMA_JSdiv), 4)        # 0.0045

round(mean(CoDa_FTS_2_point_ARIMA_JSdiv_geo), 4)   # 0.0184
round(mean(CoDa_MFTS_2_point_ARIMA_JSdiv_geo), 4)  # 0.0133
round(mean(CoDa_MLFTS_2_point_ARIMA_JSdiv_geo), 4) # 0.0099
round(mean(RW_2_point_ARIMA_JSdiv_geo), 4)         # 0.0214
round(mean(RWD_2_point_ARIMA_JSdiv_geo), 4)        # 0.0046

#######
## ETS
#######

CoDa_FTS_1_point_ETS_err = CoDa_MFTS_1_point_ETS_err = CoDa_MLFTS_1_point_ETS_err =
    CoDa_FTS_1_point_ETS_KLdiv = CoDa_MFTS_1_point_ETS_KLdiv = CoDa_MLFTS_1_point_ETS_KLdiv =     
    CoDa_FTS_1_point_ETS_JSdiv = CoDa_MFTS_1_point_ETS_JSdiv = CoDa_MLFTS_1_point_ETS_JSdiv =     
    CoDa_FTS_1_point_ETS_JSdiv_geo = CoDa_MFTS_1_point_ETS_JSdiv_geo = CoDa_MLFTS_1_point_ETS_JSdiv_geo =     
    
    CoDa_FTS_2_point_ETS_err = CoDa_MFTS_2_point_ETS_err = CoDa_MLFTS_2_point_ETS_err = 
    CoDa_FTS_2_point_ETS_KLdiv = CoDa_MFTS_2_point_ETS_KLdiv = CoDa_MLFTS_2_point_ETS_KLdiv = 
    CoDa_FTS_2_point_ETS_JSdiv = CoDa_MFTS_2_point_ETS_JSdiv = CoDa_MLFTS_2_point_ETS_JSdiv = 
    CoDa_FTS_2_point_ETS_JSdiv_geo = CoDa_MFTS_2_point_ETS_JSdiv_geo = CoDa_MLFTS_2_point_ETS_JSdiv_geo = vector("numeric", 30)
for(iwk in 1:30)
{
    dum = CoDa_point(dat_1 = EW_female_pop, dat_2 = EW_male_pop, fh = iwk, fmethod = "ETS")
    
    ## 1st population
    
    # MAPE
    
    CoDa_FTS_1_point_ETS_err[iwk] = dum$CoDa_FTS_1_err
    CoDa_MFTS_1_point_ETS_err[iwk] = dum$CoDa_MFTS_1_err
    CoDa_MLFTS_1_point_ETS_err[iwk] = dum$CoDa_MLFTS_1_err
    
    
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

round(mean(CoDa_FTS_1_point_ETS_err), 4)   # 47.1566
round(mean(CoDa_MFTS_1_point_ETS_err), 4)  # 52.3170
round(mean(CoDa_MLFTS_1_point_ETS_err), 4) # 39.5333


round(mean(CoDa_FTS_1_point_ETS_KLdiv), 4)   # 0.0310
round(mean(CoDa_MFTS_1_point_ETS_KLdiv), 4)  # 0.0292
round(mean(CoDa_MLFTS_1_point_ETS_KLdiv), 4) # 0.0237

round(mean(CoDa_FTS_1_point_ETS_JSdiv), 4)   # 0.0075 
round(mean(CoDa_MFTS_1_point_ETS_JSdiv), 4)  # 0.0070
round(mean(CoDa_MLFTS_1_point_ETS_JSdiv), 4) # 0.0057

round(mean(CoDa_FTS_1_point_ETS_JSdiv_geo), 4)   # 0.0079
round(mean(CoDa_MFTS_1_point_ETS_JSdiv_geo), 4)  # 0.0074
round(mean(CoDa_MLFTS_1_point_ETS_JSdiv_geo), 4) # 0.0060

# 2nd population

round(mean(CoDa_FTS_2_point_ETS_err), 4)   # 55.3155
round(mean(CoDa_MFTS_2_point_ETS_err), 4)  # 48.7609
round(mean(CoDa_MLFTS_2_point_ETS_err), 4) # 44.0786


round(mean(CoDa_FTS_2_point_ETS_KLdiv), 4)   # 0.065
round(mean(CoDa_MFTS_2_point_ETS_KLdiv), 4)  # 0.0522
round(mean(CoDa_MLFTS_2_point_ETS_KLdiv), 4) # 0.0359

round(mean(CoDa_FTS_2_point_ETS_JSdiv), 4)   # 0.0143
round(mean(CoDa_MFTS_2_point_ETS_JSdiv), 4)  # 0.0128
round(mean(CoDa_MLFTS_2_point_ETS_JSdiv), 4) # 0.0087

round(mean(CoDa_FTS_2_point_ETS_JSdiv_geo), 4)   # 0.0017
round(mean(CoDa_MFTS_2_point_ETS_JSdiv_geo), 4)  # 0.0131
round(mean(CoDa_MLFTS_2_point_ETS_JSdiv_geo), 4) # 0.009

################
## RWF_no_drift
################

CoDa_FTS_1_point_RW_err = CoDa_MFTS_1_point_RW_err = CoDa_MLFTS_1_point_RW_err = 
    CoDa_FTS_1_point_RW_KLdiv = CoDa_MFTS_1_point_RW_KLdiv = CoDa_MLFTS_1_point_RW_KLdiv =     
    CoDa_FTS_1_point_RW_JSdiv = CoDa_MFTS_1_point_RW_JSdiv = CoDa_MLFTS_1_point_RW_JSdiv =     
    CoDa_FTS_1_point_RW_JSdiv_geo = CoDa_MFTS_1_point_RW_JSdiv_geo = CoDa_MLFTS_1_point_RW_JSdiv_geo =     
    
    CoDa_FTS_2_point_RW_err = CoDa_MFTS_2_point_RW_err = CoDa_MLFTS_2_point_RW_err = 
    CoDa_FTS_2_point_RW_KLdiv = CoDa_MFTS_2_point_RW_KLdiv = CoDa_MLFTS_2_point_RW_KLdiv = 
    CoDa_FTS_2_point_RW_JSdiv = CoDa_MFTS_2_point_RW_JSdiv = CoDa_MLFTS_2_point_RW_JSdiv = 
    CoDa_FTS_2_point_RW_JSdiv_geo = CoDa_MFTS_2_point_RW_JSdiv_geo = CoDa_MLFTS_2_point_RW_JSdiv_geo = vector("numeric", 30)
for(ik in 1:30)
{
    dum = CoDa_point(dat_1 = EW_female_pop, dat_2 = EW_male_pop, fh = ik, fmethod = "RWF_no_drift")
    
    ## 1st population
    
    # MAPE
    
    CoDa_FTS_1_point_RW_err[ik] = dum$CoDa_FTS_1_err
    CoDa_MFTS_1_point_RW_err[ik] = dum$CoDa_MFTS_1_err
    CoDa_MLFTS_1_point_RW_err[ik] = dum$CoDa_MLFTS_1_err
    
    
    
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

round(mean(CoDa_FTS_1_point_RW_err), 4)   # 77.9653
round(mean(CoDa_MFTS_1_point_RW_err), 4)  # 54.0564
round(mean(CoDa_MLFTS_1_point_RW_err), 4) # 64.9673


round(mean(CoDa_FTS_1_point_RW_KLdiv), 4)   # 0.1088
round(mean(CoDa_MFTS_1_point_RW_KLdiv), 4)  # 0.0548
round(mean(CoDa_MLFTS_1_point_RW_KLdiv), 4) # 0.0852

round(mean(CoDa_FTS_1_point_RW_JSdiv), 4)   # 0.0257
round(mean(CoDa_MFTS_1_point_RW_JSdiv), 4)  # 0.0131
round(mean(CoDa_MLFTS_1_point_RW_JSdiv), 4) # 0.0204

round(mean(CoDa_FTS_1_point_RW_JSdiv_geo), 4)   # 0.0275
round(mean(CoDa_MFTS_1_point_RW_JSdiv_geo), 4)  # 0.0139
round(mean(CoDa_MLFTS_1_point_RW_JSdiv_geo), 4) # 0.0214

# 2nd population

round(mean(CoDa_FTS_2_point_RW_err), 4)   # 72.6619
round(mean(CoDa_MFTS_2_point_RW_err), 4)  # 95.102
round(mean(CoDa_MLFTS_2_point_RW_err), 4) # 78.6318


round(mean(CoDa_FTS_2_point_RW_KLdiv), 4)   # 0.1251
round(mean(CoDa_MFTS_2_point_RW_KLdiv), 4)  # 0.1957
round(mean(CoDa_MLFTS_2_point_RW_KLdiv), 4) # 0.143

round(mean(CoDa_FTS_2_point_RW_JSdiv), 4)   # 0.0301
round(mean(CoDa_MFTS_2_point_RW_JSdiv), 4)  # 0.0464
round(mean(CoDa_MLFTS_2_point_RW_JSdiv), 4) # 0.0339

round(mean(CoDa_FTS_2_point_RW_JSdiv_geo), 4)   # 0.0312
round(mean(CoDa_MFTS_2_point_RW_JSdiv_geo), 4)  # 0.0487
round(mean(CoDa_MLFTS_2_point_RW_JSdiv_geo), 4) # 0.0359


#############
## RWF_drift
#############

CoDa_FTS_1_point_RWD_err = CoDa_MFTS_1_point_RWD_err = CoDa_MLFTS_1_point_RWD_err = 
    CoDa_FTS_1_point_RWD_KLdiv = CoDa_MFTS_1_point_RWD_KLdiv = CoDa_MLFTS_1_point_RWD_KLdiv =     
    CoDa_FTS_1_point_RWD_JSdiv = CoDa_MFTS_1_point_RWD_JSdiv = CoDa_MLFTS_1_point_RWD_JSdiv =     
    CoDa_FTS_1_point_RWD_JSdiv_geo = CoDa_MFTS_1_point_RWD_JSdiv_geo = CoDa_MLFTS_1_point_RWD_JSdiv_geo =     
    
    CoDa_FTS_2_point_RWD_err = CoDa_MFTS_2_point_RWD_err = CoDa_MLFTS_2_point_RWD_err = 
    CoDa_FTS_2_point_RWD_KLdiv = CoDa_MFTS_2_point_RWD_KLdiv = CoDa_MLFTS_2_point_RWD_KLdiv = 
    CoDa_FTS_2_point_RWD_JSdiv = CoDa_MFTS_2_point_RWD_JSdiv = CoDa_MLFTS_2_point_RWD_JSdiv = 
    CoDa_FTS_2_point_RWD_JSdiv_geo = CoDa_MFTS_2_point_RWD_JSdiv_geo = CoDa_MLFTS_2_point_RWD_JSdiv_geo = vector("numeric", 30)
for(iwk in 1:30)
{
    dum = CoDa_point(dat_1 = EW_female_pop, dat_2 = EW_male_pop, fh = iwk, fmethod = "RWF_drift")
    
    ## 1st population
    
    # MAPE
    
    CoDa_FTS_1_point_RWD_err[iwk] = dum$CoDa_FTS_1_err
    CoDa_MFTS_1_point_RWD_err[iwk] = dum$CoDa_MFTS_1_err
    CoDa_MLFTS_1_point_RWD_err[iwk] = dum$CoDa_MLFTS_1_err
    
   
    
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

round(mean(CoDa_FTS_1_point_RWD_err), 4)   # 47.6081
round(mean(CoDa_MFTS_1_point_RWD_err), 4)  # 51.9451
round(mean(CoDa_MLFTS_1_point_RWD_err), 4) # 38.4565


round(mean(CoDa_FTS_1_point_RWD_KLdiv), 4)   # 0.0323
round(mean(CoDa_MFTS_1_point_RWD_KLdiv), 4)  # 0.0288
round(mean(CoDa_MLFTS_1_point_RWD_KLdiv), 4) # 0.0224

round(mean(CoDa_FTS_1_point_RWD_JSdiv), 4)   # 0.0078
round(mean(CoDa_MFTS_1_point_RWD_JSdiv), 4)  # 0.0069
round(mean(CoDa_MLFTS_1_point_RWD_JSdiv), 4) # 0.0054

round(mean(CoDa_FTS_1_point_RWD_JSdiv_geo), 4)   # 0.0082
round(mean(CoDa_MFTS_1_point_RWD_JSdiv_geo), 4)  # 0.0073
round(mean(CoDa_MLFTS_1_point_RWD_JSdiv_geo), 4) # 0.0057

## 2nd population

round(mean(CoDa_FTS_2_point_RWD_err), 4)   # 43.131
round(mean(CoDa_MFTS_2_point_RWD_err), 4)  # 49.152
round(mean(CoDa_MLFTS_2_point_RWD_err), 4) # 44.353


round(mean(CoDa_FTS_2_point_RWD_KLdiv), 4)   # 0.032
round(mean(CoDa_MFTS_2_point_RWD_KLdiv), 4)  # 0.0535
round(mean(CoDa_MLFTS_2_point_RWD_KLdiv), 4) # 0.0352

round(mean(CoDa_FTS_2_point_RWD_JSdiv), 4)   # 0.0077
round(mean(CoDa_MFTS_2_point_RWD_JSdiv), 4)  # 0.0131
round(mean(CoDa_MLFTS_2_point_RWD_JSdiv), 4) # 0.0086

round(mean(CoDa_FTS_2_point_RWD_JSdiv_geo), 4)   # 0.0081
round(mean(CoDa_MFTS_2_point_RWD_JSdiv_geo), 4)  # 0.0134
round(mean(CoDa_MLFTS_2_point_RWD_JSdiv_geo), 4) # 0.0089


#######################
# product-ratio method
#######################
 
PR_1_mape = PR_1_KLdiv = PR_1_JSdiv = PR_1_JSdiv_geo =
    PR_2_err = PR_2_mape = PR_2_KLdiv = PR_2_JSdiv = PR_2_JSdiv_geo= vector("numeric", 30)

for(ik in 1:30)
{
    print(ik)
    dum = PR_eval(dat_1 = EW_female_pop, dat_2 = EW_male_pop, fh = ik,fmethod='ARIMA')
    
    ## 1st population
    
    # MAPE
    PR_1_mape[ik] = dum$PR_1_mape
    PR_1_KLdiv[ik] = dum$PR_1_KLdiv
    PR_1_JSdiv[ik] = dum$PR_1_JSdiv
    PR_1_JSdiv_geo[ik] =dum$PR_1_JSdiv_geo
    
    PR_2_mape[ik] = dum$PR_2_mape
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



##################
# Result summary
#
summary_result = data.frame(Modelling_method = c('FTS',rep('',3),'MFTS',rep('',3),'MLFTS',rep('',3),'RW','PR'),
                            Forecast_method = c(rep(c('ARIMA','ETS','RW','RWD'),length = 12),'',''),
                            MAPE1 = c(round(mean(CoDa_FTS_1_point_ARIMA_err), 4),round(mean(CoDa_FTS_1_point_ETS_err), 4),round(mean(CoDa_FTS_1_point_RW_err), 4),round(mean(CoDa_FTS_1_point_RWD_err), 4),
                                      round(mean(CoDa_MFTS_1_point_ARIMA_err), 4), round(mean(CoDa_MFTS_1_point_ETS_err), 4),round(mean(CoDa_MFTS_1_point_RW_err), 4),round(mean(CoDa_MFTS_1_point_RWD_err), 4),
                                      round(mean(CoDa_MLFTS_1_point_ARIMA_err), 4), round(mean(CoDa_MLFTS_1_point_ETS_err), 4),round(mean(CoDa_MLFTS_1_point_RW_err), 4),round(mean(CoDa_MLFTS_1_point_RWD_err), 4),
                                      round(mean(RW_1_point_ARIMA_err), 4),
                                      round(mean(PR_1_mape), 4) ),
                            KLD1 = c(round(mean(CoDa_FTS_1_point_ARIMA_KLdiv), 4),round(mean(CoDa_FTS_1_point_ETS_KLdiv), 4),round(mean(CoDa_FTS_1_point_RW_KLdiv), 4),round(mean(CoDa_FTS_1_point_RWD_KLdiv), 4),
                                     round(mean(CoDa_MFTS_1_point_ARIMA_KLdiv), 4), round(mean(CoDa_MFTS_1_point_ETS_KLdiv), 4),round(mean(CoDa_MFTS_1_point_RW_KLdiv), 4),round(mean(CoDa_MFTS_1_point_RWD_KLdiv), 4),
                                     round(mean(CoDa_MLFTS_1_point_ARIMA_KLdiv), 4), round(mean(CoDa_MLFTS_1_point_ETS_KLdiv), 4),round(mean(CoDa_MLFTS_1_point_RW_KLdiv), 4),round(mean(CoDa_MLFTS_1_point_RWD_KLdiv), 4),
                                     round(mean(RW_1_point_ARIMA_KLdiv), 4),
                                     round(mean(PR_1_KLdiv), 4) ),
                            JSD1 = c(round(mean(CoDa_FTS_1_point_ARIMA_JSdiv), 4),round(mean(CoDa_FTS_1_point_ETS_JSdiv), 4),round(mean(CoDa_FTS_1_point_RW_JSdiv), 4),round(mean(CoDa_FTS_1_point_RWD_JSdiv), 4),
                                     round(mean(CoDa_MFTS_1_point_ARIMA_JSdiv), 4), round(mean(CoDa_MFTS_1_point_ETS_JSdiv), 4),round(mean(CoDa_MFTS_1_point_RW_JSdiv), 4),round(mean(CoDa_MFTS_1_point_RWD_JSdiv), 4),
                                     round(mean(CoDa_MLFTS_1_point_ARIMA_JSdiv), 4), round(mean(CoDa_MLFTS_1_point_ETS_JSdiv), 4),round(mean(CoDa_MLFTS_1_point_RW_JSdiv), 4),round(mean(CoDa_MLFTS_1_point_RWD_JSdiv), 4),
                                     round(mean(RW_1_point_ARIMA_JSdiv), 4),
                                     round(mean(PR_2_JSdiv), 4) ),
                            JSDG1 = c(round(mean(CoDa_FTS_1_point_ARIMA_JSdiv_geo), 4),round(mean(CoDa_FTS_1_point_ETS_JSdiv_geo), 4),round(mean(CoDa_FTS_1_point_RW_JSdiv_geo), 4),round(mean(CoDa_FTS_1_point_RWD_JSdiv_geo), 4),
                                      round(mean(CoDa_MFTS_1_point_ARIMA_JSdiv_geo), 4), round(mean(CoDa_MFTS_1_point_ETS_JSdiv_geo), 4),round(mean(CoDa_MFTS_1_point_RW_JSdiv_geo), 4),round(mean(CoDa_MFTS_1_point_RWD_JSdiv_geo), 4),
                                      round(mean(CoDa_MLFTS_1_point_ARIMA_JSdiv_geo), 4), round(mean(CoDa_MLFTS_1_point_ETS_JSdiv_geo), 4),round(mean(CoDa_MLFTS_1_point_RW_JSdiv_geo), 4),round(mean(CoDa_MLFTS_1_point_RWD_JSdiv_geo), 4),
                                      round(mean(RW_1_point_ARIMA_JSdiv_geo), 4),
                                      round(mean(PR_2_JSdiv_geo), 4) ),
                            MAPE2 = c(round(mean(CoDa_FTS_2_point_ARIMA_err), 4),round(mean(CoDa_FTS_2_point_ETS_err), 4),round(mean(CoDa_FTS_2_point_RW_err), 4),round(mean(CoDa_FTS_2_point_RWD_err), 4),
                                      round(mean(CoDa_MFTS_2_point_ARIMA_err), 4), round(mean(CoDa_MFTS_2_point_ETS_err), 4),round(mean(CoDa_MFTS_2_point_RW_err), 4),round(mean(CoDa_MFTS_2_point_RWD_err), 4),
                                      round(mean(CoDa_MLFTS_2_point_ARIMA_err), 4), round(mean(CoDa_MLFTS_2_point_ETS_err), 4),round(mean(CoDa_MLFTS_2_point_RW_err), 4),round(mean(CoDa_MLFTS_2_point_RWD_err), 4),
                                      round(mean(RW_2_point_ARIMA_err), 4),
                                      round(mean(PR_2_mape), 4) ),
                            KLD2 = c(round(mean(CoDa_FTS_2_point_ARIMA_KLdiv), 4),round(mean(CoDa_FTS_2_point_ETS_KLdiv), 4),round(mean(CoDa_FTS_2_point_RW_KLdiv), 4),round(mean(CoDa_FTS_2_point_RWD_KLdiv), 4),
                                     round(mean(CoDa_MFTS_2_point_ARIMA_KLdiv), 4), round(mean(CoDa_MFTS_2_point_ETS_KLdiv), 4),round(mean(CoDa_MFTS_2_point_RW_KLdiv), 4),round(mean(CoDa_MFTS_2_point_RWD_KLdiv), 4),
                                     round(mean(CoDa_MLFTS_2_point_ARIMA_KLdiv), 4), round(mean(CoDa_MLFTS_2_point_ETS_KLdiv), 4),round(mean(CoDa_MLFTS_2_point_RW_KLdiv), 4),round(mean(CoDa_MLFTS_2_point_RWD_KLdiv), 4),
                                     round(mean(RW_2_point_ARIMA_KLdiv), 4),
                                     round(mean(PR_2_KLdiv), 4) ),
                            JSD2 = c(round(mean(CoDa_FTS_2_point_ARIMA_JSdiv), 4),round(mean(CoDa_FTS_2_point_ETS_JSdiv), 4),round(mean(CoDa_FTS_2_point_RW_JSdiv), 4),round(mean(CoDa_FTS_2_point_RWD_JSdiv), 4),
                                     round(mean(CoDa_MFTS_2_point_ARIMA_JSdiv), 4), round(mean(CoDa_MFTS_2_point_ETS_JSdiv), 4),round(mean(CoDa_MFTS_2_point_RW_JSdiv), 4),round(mean(CoDa_MFTS_2_point_RWD_JSdiv), 4),
                                     round(mean(CoDa_MLFTS_2_point_ARIMA_JSdiv), 4), round(mean(CoDa_MLFTS_2_point_ETS_JSdiv), 4),round(mean(CoDa_MLFTS_2_point_RW_JSdiv), 4),round(mean(CoDa_MLFTS_2_point_RWD_JSdiv), 4),
                                     round(mean(RW_2_point_ARIMA_JSdiv), 4),
                                     round(mean(PR_2_JSdiv), 4) ),
                            JSDG2 = c(round(mean(CoDa_FTS_2_point_ARIMA_JSdiv_geo), 4),round(mean(CoDa_FTS_2_point_ETS_JSdiv_geo), 4),round(mean(CoDa_FTS_2_point_RW_JSdiv_geo), 4),round(mean(CoDa_FTS_2_point_RWD_JSdiv_geo), 4),
                                      round(mean(CoDa_MFTS_2_point_ARIMA_JSdiv_geo), 4), round(mean(CoDa_MFTS_2_point_ETS_JSdiv_geo), 4),round(mean(CoDa_MFTS_2_point_RW_JSdiv_geo), 4),round(mean(CoDa_MFTS_2_point_RWD_JSdiv_geo), 4),
                                      round(mean(CoDa_MLFTS_2_point_ARIMA_JSdiv_geo), 4), round(mean(CoDa_MLFTS_2_point_ETS_JSdiv_geo), 4),round(mean(CoDa_MLFTS_2_point_RW_JSdiv_geo), 4),round(mean(CoDa_MLFTS_2_point_RWD_JSdiv_geo), 4),
                                      round(mean(RW_2_point_ARIMA_JSdiv_geo), 4),
                                      round(mean(PR_2_JSdiv_geo), 4) ))


##########
# Figures
##########

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



