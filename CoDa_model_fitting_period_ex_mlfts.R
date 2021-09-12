# An example of model fitting through a multilevel functional time series method with RWD forecasting method

source("CoDa_model_fitting.R")

R2_MLFTS_RWD = R_square_fit(dat_1 = t(EW_female_pop), dat_2 = t(EW_male_pop), fh = 30,
                            modeling_method = "MLFTS", forecasting_method = "RWF_drift")

#####################
# CoDa model fitting
#####################

# Fig_2a

plot(0:110, R2_MLFTS_RWD$alpha_x_1, xlab = "", ylab = expression(alpha[n](u)), type = "l", main = "Female")

# Fig_2b

plot(0:110, R2_MLFTS_RWD$alpha_x_2, xlab = "", ylab = "", type = "l", main = "Male")

# Fig_2c

plot(fts(0:110, t(R2_MLFTS_RWD$h_x_t_1)), xlab = "", ylab = expression(bold(beta)(u)))

# Fig_2d

plot(fts(0:110, t(R2_MLFTS_RWD$h_x_t_2)), xlab = "", ylab = "")

# Fig_2e

plot(fts(0:110, R2_MLFTS_RWD$recon_1), xlab = "", ylab = "Reconstructed life-table death count", ylim = c(0, 20000))

# Fig_2f

plot(fts(0:110, R2_MLFTS_RWD$recon_2), xlab = "", ylab = "", ylim = c(0, 20000))

# Fig_2g

plot(fts(0:110, R2_MLFTS_RWD$fore_count_1), xlab = "Age", ylab = "Forecast life-table death count")

# Fig_2h

plot(fts(0:110, R2_MLFTS_RWD$fore_count_2), xlab = "Age", ylab = "")

