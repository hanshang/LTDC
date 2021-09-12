####################################################################################################
# Raw data downloaded from HMD:
####################################################################################################
# EW_lt_female_death_period.txt: the raw female period life-table death counts
# EW_lt_male_death_period.txt: the raw male period life-table death counts
# Explanatory_analysis_qx.R: read the period life-table death counts
# uk_cohort_female.txt: the raw (incomplete) female cohort life-table death counts
# uk_cohort_male.txt: the raw (incomplete) male cohort life-table death counts

# the completed cohort data using PCLM:
# uk_female_pop_complete.txt: Completed female cohort life-table death counts for England and Wales
# uk_male_pop_complete.txt: Completed male cohort life-table death counts for England and Wales

########################################################################################################
# Functions:
########################################################################################################
# pclm_source.R: main code for completing the cohort using the penalized composite link model
# ER_GR.R: select the number of FPCA components
# product_ratio_func.R: main code for product-ratio method in CoDa space
# CoDa_naive_func.R: main code for random-walk with and without drift in CoDa space
# CoDa_model_fitting.R: main code for combining functional time series methods with CoDa
# CoDa_model_fitting_cohort.R: main code for combining functional time series methods with CoDa

################################################################################################################################
# Period results:
################################################################################################################################
# CoDa_model_fitting_period_ex_mlfts.R: an example of using the multilevel functional time series method
# CoDa_point_period.R: Point forecast evaluation and comparison
# CoDa_interval_fts_period.R: Interval forecast evaluation and comparison using the univariate functional time series method
# CoDa_interval_mfts_period.R: Interval forecast evaluation and comparison using the multivariate functional time series method
# CoDa_interval_mlfts_period.R: Interval forecast evaluation and comparison using the multilevel functional time series method

################################################################################################################################
# Cohort results:
################################################################################################################################
# CoDa_model_fitting_cohort_ex_mlfts.R: an example of using the multilevel functional time series method
# CoDa_point_cohort.R: Point forecast evaluation and comparison
# CoDa_interval_fts_cohort.R: Interval forecast evaluation and comparison using the univariate functional time series method
# CoDa_interval_mfts_cohort.R: Interval forecast evaluation and comparison using the multivariate functional time series method
# CoDa_interval_mlfts_cohort.R: Interval forecast evaluation and comparison using the multilevel functional time series method
