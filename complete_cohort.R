# Complete the cohort for England and Welash ('GERTENW') from 1841-1960 (the observed complete data is until 1925)

# R code article "Killing off cohorts: Forecasting mortality of non-extinct
# cohorts with the penalized composite link model" by Rizzi et al.
# Demo R code for Sweden. All other countries follow the same analysis. 

# load female and male death rates
load("Female_GBRTENW.RData" )
fdata = data
rm(data)

load("Male_GBRTENW.RData")
mdata = data
rm(data)

# load female and male life-tables
uk_female_cohort = read.table("uk_cohort_female.txt", header = TRUE, skip = 2)
uk_male_cohort = read.table("uk_cohort_male.txt", header = TRUE, skip = 2)


# load packages and functions
source("lib/pclm_source.r")


# Dowload packages
library(forecast)
library(vars)
library(MortalitySmooth)
library(ungroup)
library(plyr)
library(lmtest)
library(stsm)
library(KFKSDS)
library(KFAS)
library(ftsa)

options(scipen=999)

#################################
# CALCULATE LIFE TABLES #
#################################

# Female cohort life table
LT.fun_F<-function(mx, age, nx, ax, r){
  #Life table
  LT<-data.frame(Age=age,
                 mx=mx)
  #Adjust ax
  if(mx[1]< 0.01724){
    ax[1]<-0.1490 - 2.0867*LT$mx[LT$Age==0]
  }
  if(mx[1]>=0.0658){
    ax[1]<-0.31411
  }else{ax[1]<-0.0438 + 4.1075*LT$mx[LT$Age==0]}
  ax[length(ax)]<-1/LT$mx[length(LT$mx)]
  ax[ax==Inf]<-0
  LT$ax<-ax
  #qx : the n-year death prob for an individual aged x
  LT$qx<- LT$mx*nx/(1+(nx-ax)*LT$mx)
  LT$qx[LT$qx>1]<-1
  LT$qx[length(LT$qx)]<-1
  #px : the n-year survival prob for an individual aged x 
  LT$px<- 1-LT$qx
  #lx : the total number of survivals until age x
  LT$lx<-c(r, r*cumprod(LT$px[1:(length(LT$px)-1)]))
  #dx: the observed death counts at age x
  LT$dx<-LT$qx*LT$lx
  #LX: the total number of survivals until the mid-point (adjusted ax) of aged x
  LT$Lx<- LT$lx*nx - LT$dx*(nx-ax)
  #Tx
  LT$Tx<- rev(cumsum(LT$Lx[length(LT$Lx):1]))
  #ex
  LT$ex<- LT$Tx/LT$lx
  #outcomes: the table
  return(LT=LT)
}


years <- fdata$years[1]:1960
fmx <- fdata$f_mx
r<-100000 #radix of the life table
nx<-1
age <- 0:110
ax<-rep(nx/2, length(age))
# Create life tables for females for each calendar year
f.rate_list <- list()
for ( i in 1:dim(fmx)[2]){
    f.rate_list[[i]] <- fmx
}
f.LT <- apply(f.rate_list[[1]],2,LT.fun_F,
                age=age, nx=nx, ax=ax, r=r)
  

dx <- matrix(0, nrow=length(age), ncol=length(years))
for (i in 1:length(years)){
    dx[,i] <- f.LT[[i]]$dx
}
rownames(dx) <- age
colnames(dx) <- years
  
lx <- matrix(0, nrow=length(age), ncol=length(years))
for (i in 1:length(years)){
    lx[,i] <- f.LT[[i]]$lx
}
rownames(lx) <- age
colnames(lx) <- years
  
  
fdata$lf = list(f.LT = f.LT, dx = dx, lx = lx,
                 date = date())




# Male cohort life table
LT.fun_M<-function(mx, age, nx, ax, r){
  #Life table
  LT<-data.frame(Age=age,
                 mx=mx)
  #Adjust ax
  if(mx[1]< 0.0226){
    ax[1]<-0.14903 - 2.0367*LT$mx[LT$Age==0]
  }
  if(mx[1]>=0.0785){
    ax[1]<-0.2991
  }else{ax[1]<-0.0244 + 3.4994*LT$mx[LT$Age==0]}
  ax[length(ax)]<-1/LT$mx[length(LT$mx)]
  ax[ax==Inf]<-0
  LT$ax<-ax
  #qx
  LT$qx<- LT$mx*nx/(1+(nx-ax)*LT$mx)
  LT$qx[LT$qx>1]<-1
  LT$qx[length(LT$qx)]<-1
  #px
  LT$px<- 1-LT$qx
  #lx
  LT$lx<-c(r, r*cumprod(LT$px[1:(length(LT$px)-1)]))
  #dx
  LT$dx<-LT$qx*LT$lx
  #LX
  LT$Lx<- LT$lx*nx - LT$dx*(nx-ax)
  #Tx
  LT$Tx<- rev(cumsum(LT$Lx[length(LT$Lx):1]))
  #ex
  LT$ex<- LT$Tx/LT$lx
  #outcomes: the table
  return(LT=LT)
}



mmx <- mdata$m_mxc
r<-100000 #radix of the life table
nx<-1
age <- 0:110
ax<-rep(nx/2, length(age))
# Create life tables for females for each calendar year
m.rate_list <- list()
for ( i in 1:dim(mmx)[2]){
    m.rate_list[[i]] <- mmx
}
m.LT <- apply(m.rate_list[[1]],2,LT.fun_M,
                age=age, nx=nx, ax=ax, r=r)
  

dx <- matrix(0, nrow=length(age), ncol=length(years))
for (i in 1:length(years)){
    dx[,i] <- m.LT[[i]]$dx
}
rownames(dx) <- age
colnames(dx) <- years
  
lx <- matrix(0, nrow=length(age), ncol=length(years))
for (i in 1:length(years)){
    lx[,i] <- m.LT[[i]]$lx
}
rownames(lx) <- age
colnames(lx) <- years
  
  
mdata$lf = list(m.LT = m.LT, dx = dx, lx = lx,
                 date = date())


  
#######################
# STEP 3: FORECAST dx #
#######################

  


age <- 0:110
for (sex in c('M','F')){
  # Load files
  if(sex=="F") {data = fdata
  } else if(sex=="M") {data = mdata}
  dx <- data$lf$dx[,]
  
  # Define objects
  y <- NULL # sequence of age-specific death counts for each calendar year

  if(sex=="M"){
    Exp.all <- data$m_exc
    D.obs <- data$m_Dx
  } else if(sex=="F"){
    Exp.all <- data$f_exc
    D.obs <- data$f_Dx }
  
  # Define these values for each dataset
  c1 <- 1841:1906
  c2 <- 1907:1935
  c3 <- 1936:1960
  
  #check that the radix is 100,000 as required
  temp <- (dx[,colnames(dx) %in% 1900])
  temp[is.na(temp)] <- 0
  radix <- sum(temp)
  radix
  
  # Define time intervals
  dx.complete <- dx[,colnames(dx) %in% as.character(c1)] # define the observed completed cohorts
  dx.complete[is.na(dx.complete)] <- 0
  dx.na.0 <- dx[,colnames(dx) %in% as.character(c1)] # interval fitted with pclm only (named c1 ) (the observed death count)
  dx.na.1 <- dx[,colnames(dx) %in% as.character(c2)] # interval fitted with pclm only (named c2 ) 
  dx.na.2 <- dx[,colnames(dx) %in% as.character(c3)] # interval fitted by imposing mode (named c3) 
  
  dx.completed.0 <- matrix(NA,length(age),length(c1)) # define the complete cohort for c1 usinge the PCLM
  dx.completed.1 <- matrix(NA,length(age),length(c2))
  dx.completed.2 <- matrix(NA,length(age),length(c3))
  dx.completed.2.fine <- matrix(NA,length(age),length(c3))
  
  
  # ---------------------  Fit C1 ------------------------------- #
  for(t in 1:length((c1))){
    last.age.fi <- 110 # define last age in c1 interval
    dx.obs.temp.fi <- dx.na.0[1:length(0:(last.age.fi)),t]
    dx.obs.temp.fi[is.na(dx.obs.temp.fi)] <- 0
    y.fi <- c(dx.obs.temp.fi,(radix - sum(dx.obs.temp.fi)))
    y.fi <- ifelse(y.fi < 0,0,y.fi)
    names(y.fi) <- c(seq(0,length(0:(last.age.fi-1))), "open")
    # add extra bin with 0 counts
    y.fi <-c(y.fi,0)
    x=1:113
    mod.fi <- ungroup::pclm(x=0:112,y.fi,nlast=1,control = list(lambda=100))
    dx.completed.0[,t] <- c(mod.fi$fitted[1:110], sum(mod.fi$fitted[111:112]))
    print(t)
  }
  lamda.pr <- matrix(NA,length(c2),1)
  
  # -------------------------  Fit C2 ----------------------------- #
  for(t in 1:length((c2))){
    last.age <- 110-(c2[t]-max(c1)+1) # define last age in c2 interval
    dx.obs.temp <- dx.na.1[1:length(0:(last.age)),t]
    dx.obs.temp[is.na(dx.obs.temp)] <- 0
    y <- c(dx.obs.temp,(radix - sum(dx.obs.temp)))
    names(y) <- c(seq(0,length(0:(last.age-1))), "open")
    # add extra bin with 0 counts
    y <-c(y,0)
    # define the detailed grid where to estimate the underlying density
    m = 131
    x = 1:m # age in single year interval
    C0 <- make.c0(last.age = last.age)
    B = diag(m)
    # penalty parameters
    lambdas = 10^seq(2, 10, by = 0.5)
    pord = 2
    pord.2 = 3
    # containers for results
    R = matrix(0, m, 1)
    lamopt = rep(0, 1)
    s.e <- matrix(0, m, 1)
    # run the model
    aics = NULL
    for (lambda in lambdas) {
      mod  = pclm(y,  C0, B, lambda = lambda, deg = pord, show = F)
      aics = c(aics, mod$aic)
    }
    opt = which.min(aics)
    lamopt = lambdas[opt]
    mod  = pclm(y, C0[ , ], B, lambda = lamopt, deg = pord)
    mod.2  = pclm(y, C0[ , ], B, lambda = lamopt, deg = pord.2)
    R[, 1] = mod$gamma
    s.e[,1] = sqrt(diag(mod$H0))
    dx.completed.1[,t] <- c(R[-c(111:131)], sum(R[111:131]))
    lamda.pr[t,1] <- lamopt
    print(t)
  }
  
  
  # ----------- Calculate and forecast necessary time series  ------------ #
  # Define years
  years.all <- c(c1,c2,c3) # interval for observed and forecasts
  years.c1.c2 <- c(c1,c2) # c1 and c2 interval
  years.for <- c3
  
  M <- matrix(NA,length(years.c1.c2),2)
  dx.c1.c2 <- cbind(dx[,colnames(dx) %in% as.character(c1)],dx.completed.1)
  
  # Calculate mode and deaths at the mode
  for(u in 1:length(years.c1.c2)){
    M[u,2] <- max(dx.c1.c2[(5+40):111,u],na.rm=TRUE)
    M[u,1] <- which.max(dx.c1.c2[(5+40):111,u]) +(3+40)
  }
  
  # Calculate deaths after the mode
  part.up <- matrix(NA,length(years.c1.c2),1)
  for(u in 1:length(years.c1.c2)){
    temp <- ((M[u,1]) +1)
    part.up[u,1] <- radix - sum(dx.c1.c2[1:temp,u],na.rm=TRUE) # count number of deaths after the mode (mode not included)
  }
  part.re <- part.up/radix
  
  # deaths after the mode (with ARIMA model)
  
  M.1.for <- forecast::forecast(Arima(M[,1],order=c(0,1,1),include.drift = TRUE ),length(c3))
  
  M.2.for <- forecast::forecast(Arima(M[,2],order=c(1,1,1),include.drift = TRUE),length(c3))
  
  # Forecast number of deaths after the mode (with ARIMA model)
  part.re.for <- forecast::forecast(Arima(part.re,order=c(0,1,1),include.drift = T),length(c3))
  
  xreg = c(rep(0,times=19),1,rep(0,times=(nrow(M)-20)))
  
  
  auto.arima(M[,1],max.p=2,max.q=2,allowdrift=T,ic="aic")
  auto.arima(M[,2],max.p=2,max.q=2,allowdrift=T,ic="aic")
  
  
  # forecaste with LLT model
  # ----------------- Kalman Filter  ---------------- #
  # ------------ Mode ----------------- #
  mode.ts <- ts(M[,1],start=c1[1],end=max(c2))
  model_struc <- SSModel(mode.ts ~ + SSMtrend(degree = 2, Q = list(matrix(NA), matrix(0))), H = matrix(NA))
  fit_struc <- fitSSM(model_struc, inits = c(0, 0), method = "BFGS")$model
  out.struc <- KFS(fit_struc,filtering = "state", smoothing = "state")
  out.struc.res <- KFS(fit_struc,filtering = "state", smoothing = c("state", "mean", "disturbance"))
  pred.struc.1.res <- residuals(out.struc.res,type="pearson")
  fit.struc <- predict(out.struc$model,interval = "confidence",level = 0.95)
  pred.struc <- predict(fit_struc, n.ahead= length(c3))
  fit.struc.ts <- ts(fit.struc[,1],start=c1[1],end=max(c2))
  pred.struc.ts <- ts(pred.struc,start=min(c3),end=max(c3))
  
  # ------------ Deaths mode ----------------- #
  mode.ts.2 <- ts(M[,2],start=c1[1],end=max(c2))
  
  model_struc.2 <- SSModel(mode.ts.2 ~ + SSMtrend(degree = 2, Q = list(matrix(NA), matrix(0))), H = matrix(NA))
  fit_struc.2 <- fitSSM(model_struc.2, inits = c(0, 0), method = "BFGS")$model
  out.struc.2 <- KFS(fit_struc.2,filtering = "state", smoothing = "state")
  out.struc.2.res <- KFS(fit_struc.2,filtering = "state", smoothing = c("state", "mean", "disturbance"))
  pred.struc.2.res <- residuals(out.struc.2.res,type="pearson")
  fit.struc.2 <- predict(out.struc.2$model,interval = "confidence")
  pred.struc.2 <- predict(fit_struc.2, n.ahead= length(c3))
  fit.struc.ts.2 <- ts(fit.struc.2[,1],start=c1[1],end=max(c2))
  pred.struc.ts.2 <- ts(pred.struc.2,start=min(c3),end=max(c3))
  
  # ---------- Fraction -------------------- #
  mode.ts.3 <- ts(part.re,start=c1[1],end=max(c2))
  model_struc.3 <- SSModel(mode.ts.3 ~ + SSMtrend(degree = 2, Q = list(matrix(NA), matrix(0))), H = matrix(NA))
  fit_struc.3 <- fitSSM(model_struc.3, inits = c(0, 0), method = "BFGS")$model
  out.struc.3 <- KFS(fit_struc.3,filtering = "state")
  out.struc.3.res <- KFS(fit_struc.3,filtering = "state", smoothing = c("state", "mean", "disturbance"))
  pred.struc.3.res <- residuals(out.struc.3.res,type="pearson")
  fit.struc.3 <- predict(out.struc.3$model,interval = "confidence")
  pred.struc.3 <- predict(fit_struc.3, n.ahead= length(c3))
  fit.struc.ts.3 <- ts(fit.struc.3[,1],start=c1[1],end=max(c2))
  pred.struc.ts.3 <- ts(pred.struc.3,start=min(c3),end=max(c3))
  
  
  # Combine observed mode and forecast
  M.for <- rbind(M,cbind(pred.struc.ts,pred.struc.ts.2))
  share.for <- matrix(c(part.re,pred.struc.ts.3),length(years.all),1)
  rownames(M.for) <- years.all
  rownames(share.for) <- years.all
  colnames(M.for) <- c("M","M.dea")
  
  #------------ Fit C3 --------------------------#
  lamda.pr.T3 <- matrix(NA,length(c(c3)),1)
  # compute the group counts
  W=c(10^2,10^8,10^3,10^3)
  for(t in (length(c2)+1):(length(c(c2,c3))) ){
    year.temp<- min(c2)+t-1
    last.age.c3 <- max(age)-(year.temp-min(c2)+2)
    M.temp <- round(M.for[years.all %in% year.temp,1])
    M.temp.fine <- round_any((M.for[years.all %in% year.temp,1]), accuracy=0.5, f=round)
    M.dea.temp <- M.for[years.all %in% year.temp,2]
    share.temp <- share.for[years.all %in% year.temp,1]
    
    # define grid for half year age
    x.fine=c(0:length(0:(last.age.c3)),M.temp.fine,(M.temp.fine+1),120)
    y.c3 <- c(dx.na.2[1:length(0:(last.age.c3)),(t-length(c2))],((radix*(1-share.temp))-sum(dx.na.2[1:length(0:(last.age.c3)),(t-length(c2))]))-M.dea.temp,M.dea.temp,(radix*share.temp),0)
    names(y.c3) <- c(seq(0,length(0:(last.age.c3-1))), "open1","mode","open2","end")
    
    C3.fin.fine <- build_C_matrix(x=x.fine,y=y.c3,nlast=(131-max(x.fine)),out.step=0.5,type= "1D")$C
    C0.c3.fin <- make.c3(last.age.c3=last.age.c3,M.temp=M.temp)
    B = diag(m)
    basis = 1:(m*2)
    xl <- min(basis)
    xr <- max(basis)
    xmax <- xr + 0.01 * (xr - xl)
    xmin <- xl - 0.01 * (xr - xl)
    B.fine <-  MortSmooth_bbase(basis, xmin, xmax, ndx=floor(m/2), deg=3)
    
    # penalty parameters
    lambdas = 10^seq(1, 5, by = 0.25)
    pord = 2
    # containers for results
    R.c3 = matrix(0, m, 1)
    lamopt = rep(0, 1)
    s.e. <- matrix(0, m, 1)
    # run the model
    aics = NULL
    aics.fine = NULL
    for (lambda in lambdas) {
      mod.c3  = pclm.c3(y.c3,  C0.c3.fin, B, lambda = lambda, deg = pord, show = F,W=W,last.age.c3=last.age.c3) # grid is 1 age
      mod.c3.fine  = pclm.c3(y.c3 , C3.fin.fine, B.fine, lambda = lambda, deg = pord,W=W,last.age.c3=last.age.c3) #grid is 0.5 age
      aics = c(aics, mod.c3$aic)
      aics.fine = c(aics.fine, mod.c3.fine$aic)
    }
    
    opt = which.min(aics)
    opt.fine = which.min(aics.fine)
    lamopt = lambdas[opt]
    lamopt.fine = lambdas[opt.fine]
    mod.c3  = pclm.c3(y.c3 , C0.c3.fin[ , ], B, lambda = lamopt, deg = pord,W=W,last.age.c3=last.age.c3)
    mod.c3.fine  = pclm.c3(y.c3, C3.fin.fine[ , ], B.fine, lambda =  lamopt.fine, deg = pord,W=W,last.age.c3=last.age.c3)
    
    res.fine <- data.frame(gam=mod.c3.fine$gamma,bin=c(rep(age,each=2),rep(110,each=length(111:130)*2)))
    gam.sin <- aggregate(res.fine, by=list(res.fine$bin),FUN=sum, na.rm=TRUE)
    
    R.c3[, 1] = mod.c3$gamma
    s.e[,1] = sqrt(diag(mod.c3$H0))
    dx.completed.2[,(t-length(c2))] <- c(R.c3[-c(111:131)], sum(R.c3[111:131]))
    dx.completed.2.fine[,(t-length(c2))] <- gam.sin$gam
    range(0,dx[,years.all %in% year.temp],dx.completed.2[,1],na.rm = T)
    
    #compare the difference of using gird = 0.5 or 1 and use the grid = 0.5 in the end!
   
    lamda.pr.T3[t-(length(c2)+1),1] <- lamopt.fine
  }
  
  

  # completed death distributions
  fit.dx <- cbind(dx.complete,dx.completed.1,dx.completed.2.fine)
  fit.dx.fine <- cbind(dx.completed.0,dx.completed.1,dx.completed.2.fine)
  
  # Calculate life expectancy
  mx.fit.fine  <- ex.fit.fine <-  mx.fit <- ex.fit <- matrix(NA,length(0:110),length(years.all))
  dx.fit  <- fit.dx
  
  for(t in 1:length(years.all)){
    mx.fit[,t] <- mx_from_dx(dx.fit[,t])
    ex.fit[,t] <- lifetable.mx(x=(0:110), mx=mx.fit[,t], sex=sex)$ex
    mx.fit.fine[,t] <- mx_from_dx(fit.dx.fine[,t])
    ex.fit.fine[,t] <- lifetable.mx(x=(0:110), mx=mx.fit.fine[,t], sex=sex)$ex
  }
  
  colnames(fit.dx) <- c1[1]:1960
  colnames(fit.dx.fine) <- c1[1]:1960
  colnames(mx.fit) <- c1[1]:1960
  colnames(ex.fit) <- c1[1]:1960
  rownames(mx.fit) <- 0:110
  rownames(ex.fit) <- 0:110
  colnames(mx.fit.fine) <- c1[1]:1960
  colnames(ex.fit.fine) <- c1[1]:1960
  rownames(mx.fit.fine) <- 0:110
  rownames(ex.fit.fine) <- 0:110
  
  # Results storage
  data$fitted = list(fit.dx = fit.dx , fit.dx.fine = fit.dx.fine,
                     mx.fit = mx.fit, mx.fit.fine = mx.fit.fine,
                     ex.fit = ex.fit , ex.fit.fine = ex.fit.fine,
                     M.1.for = pred.struc.ts,
                     M.2.for = pred.struc.ts.2, part.re.for = pred.struc.ts.3,
                     M.for = M.for, share.for = share.for,
                     c1 = c1, c2 = c2, c3 = c3,
                     date = date())
  
  
  
  data$fitted$fit.dx
  
  if (sex == 'M') {data.M = data
  } else if (sex == 'F') {data.F = data}
  
  rm(data)
  
}


##########################################
# compelete the data and export the file #
# simple explanatory analysis            #
##########################################
year = 1841:1927
unobs_year = 1928:1960
age = 0:110
uk_female_dx = matrix(uk_female_cohort[,7], length(age), length(year))
uk_female_qx = matrix(as.numeric(as.matrix(uk_female_cohort[,4])), length(age), length(year))
uk_male_dx = matrix(uk_male_cohort[,7], length(age), length(year))
uk_male_qx = matrix(as.numeric(as.matrix(uk_male_cohort[,4])), length(age), length(year))
uk_male_ex = matrix(uk_male_cohort[,10], length(age), length(year))
uk_female_pop = uk_male_pop = matrix(NA, length(age), length(year))
for(ik in 1:length(year))
{
  start_pop_female = start_pop_male = 10^5
  for(ij in 1:length(age))
  {
    uk_female_pop[ij,ik] = uk_female_qx[ij,ik] * start_pop_female
    start_pop_female = start_pop_female - uk_female_pop[ij,ik]
      
    uk_male_pop[ij,ik] = uk_male_qx[ij,ik] * start_pop_male
    start_pop_male = start_pop_male - uk_male_pop[ij,ik]
  }
}
  
colnames(uk_female_pop) = colnames(uk_male_pop) = year
  
all(round(colSums(uk_female_pop, na.rm = TRUE), 8) == 10^5) # TRUE
all(round(colSums(uk_male_pop, na.rm = TRUE), 8) == 10^5)   # TRUE

uk_female_pop_complete = cbind(uk_female_pop, data.F$fitted$fit.dx.fine[,unobs_year-1840])
colSums(uk_female_pop_complete, na.rm = TRUE)
write.table(uk_female_pop_complete, file = "uk_female_pop_complete.txt", quote = FALSE, sep = " ")
  
uk_male_pop_complete = cbind(uk_male_pop, data.M$fitted$fit.dx.fine[,unobs_year-1840])
colSums(uk_male_pop_complete, na.rm = TRUE)
write.table(uk_male_pop_complete, file = "uk_male_pop_complete.txt", quote = FALSE, sep = " ")
  
uk_female_pop_complete = as.matrix(read.table("uk_female_pop_complete.txt", sep = " ", header = TRUE))
uk_male_pop_complete = as.matrix(read.table("uk_male_pop_complete.txt", sep = " ", header = TRUE))
colnames(uk_female_pop_complete) = colnames(uk_male_pop_complete) = 1841:1960
  



# visulization
plot(fts(0:110, uk_female_pop_complete), xlab = "Age", ylab = "Life-table death count",col='grey',
     main = "England & Wales: female data (1841-1960)", ylim = c(0, 18000))
lines(fts(0:110, (uk_female_pop_complete)[,(1:(1907-1841+1))]))
legend('topright',c('Observed cohort (1841 to 1907)','Completed cohort (1908 to 1960'),cex=0.9,text.col=c('red','dark grey'))


plot(fts(0:110, uk_male_pop_complete), xlab = "Age", ylab = "Life-table death count",col='grey',
     main = "England & Wales: male data (1841-1960)", ylim = c(0, 18000))
lines(fts(0:110, (uk_male_pop_complete)[,(1:(1907-1841+1))]))
# legend('topright',c('Observed cohort (1841 to 1907)','Completed cohort (1908 to 1960'),cex=0.9,text.col=c('red','dark grey'))

#the Gini coefficients
library(REAT)

plot(ts(apply(uk_female_pop_complete,2,gini),start=1841,end=1960),xlab='Year',ylab='Gini coefficient',
     main = "England & Wales: female data (1841-1960)")
lines(ts(apply(uk_female_pop_complete[,-(1:(1907-1841+1))],2,gini),start=1908,end=1960),col='blue')
legend('topleft',c('Observed cohort (1841 to 1907)','Completed cohort (1908 to 1960'),cex=0.9,col=c('black','blue'),lty=c(1,1))



plot(ts(apply(uk_male_pop_complete,2,gini),start=1841,end=1960),xlab='Year',ylab='Gini coefficient',
     main = "England & Wales: male data (1841-1960)")
lines(ts(apply(uk_male_pop_complete[,-(1:(1907-1841+1))],2,gini),start=1908,end=1960),col='blue')
legend('topleft',c('Observed cohort (1841 to 1907)','Completed cohort (1908 to 1960'),cex=0.9,col=c('black','blue'),lty=c(1,1))


