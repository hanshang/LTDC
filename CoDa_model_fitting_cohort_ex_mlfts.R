#################
# load data
#################

female.int = as.matrix(read.table('uk_female_pop_complete.txt'))
male.int = as.matrix(read.table('uk_male_pop_complete.txt'))

#replace the 0 death count with 0.01
death0.f = which(is.na(female.int) | female.int==0)
death0.m = which(is.na(male.int) | male.int==0)

EW_female_pop = replace(female.int,death0.f,0.01) #total is 130 years
EW_male_pop = replace(male.int,death0.m,0.01)


#################
# load R package
#################
source("lib/CoDa_model_fitting_cohort.R")


R2_MLFTS_RWD = R_square_fit(dat_1 = t(EW_female_pop), dat_2 = t(EW_male_pop), fh = 30,
                            modeling_method = "MLFTS", forecasting_method = "RWF_drift")

#####################
# CoDa model fitting
#####################


plot(0:110, R2_MLFTS_RWD$alpha_x_1, xlab = "", ylab = expression(alpha[n](u)), type = "l", main = "Female")


plot(0:110, R2_MLFTS_RWD$alpha_x_2, xlab = "", ylab = "", type = "l", main = "Male")


plot(fts(0:110, t(R2_MLFTS_RWD$h_x_t_1)), xlab = "", ylab = expression(bold(beta)(u)))


plot(fts(0:110, t(R2_MLFTS_RWD$h_x_t_2)), xlab = "", ylab = "")


plot(fts(0:110, R2_MLFTS_RWD$recon_1), xlab = "", ylab = "Reconstructed life-table death count", ylim = c(0, 20000))


plot(fts(0:110, R2_MLFTS_RWD$recon_2), xlab = "", ylab = "", ylim = c(0, 20000))


plot(fts(0:110, R2_MLFTS_RWD$fore_count_1), xlab = "Age", ylab = "Forecast life-table death count")


plot(fts(0:110, R2_MLFTS_RWD$fore_count_2), xlab = "Age", ylab = "")



