##################
# load R packages
##################

install.packages(c("ftsa", "psych", "meboot", "demography"))
require(ftsa)
require(psych)
require(meboot)
require(demography)
require(reldist)

##########################
# set a working directory
##########################

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

EW_female_pop = t(female_pop)
EW_male_pop = t(male_pop)
colnames(EW_female_pop) = colnames(EW_male_pop) = 1841:2018
rownames(EW_female_pop) = rownames(EW_male_pop) = 0:110 

################
# rainbow plots
################

# Fig_1a

plot(fts(0:110, EW_female_pop), xlab = "Age", ylab = "Life-table death counts",
     main = "England & Wales: female data (1841-2018)", ylim = c(0, 17500))
legend("topright", c("1841", "1876", "1911", "1946", "1981", "2016"), 
       col = c("#FF0000", "#FFF100", "#1CFF00", "#00FFD6", "#0037FF", "#BA00FF"),
       lty = rep(1, 6), cex = 0.8, ncol = 2)

# Fig_1b

plot(fts(0:110, EW_male_pop), xlab = "Age", ylab = "Life-table death counts",
     main = "England & Wales: male data (1841-2018)", ylim = c(0, 17500))

###########################
# compute Gini coefficient
###########################

# Fig_2a

plot(ts(apply(EW_female_pop, 2, gini), start = 1841, end = 2018),
     xlab = "Year", ylab = "Gini coefficient", main = "England & Wales: female data (1841-2018)")

# Fig_2b

plot(ts(apply(EW_male_pop, 2, gini), start = 1841, end = 2018),
     xlab = "Year", ylab = "Gini coefficient", main = "England & Wales: male data (1841-2018)")

