# R code article "Killing off cohorts: Forecasting mortality of non-extinct
# cohorts with the penalized composite link model" by Rizzi et al.
# SOURCE FUNCTIONS

# PCLM source function to forecast cohorts in T2
pclm <- function(y, C, X, lambda = 1, deg = 2, show = F){
  # Fit a PCLM (estimate b in ) E(y) = C %*% exp(X %*% b)
  # y = the vector of observed counts of length i
  # C = the composition matrix of dimension IxJ
  # X = the identity matrix of dimension JxJ; or B-spline basis
  # lambda = smoothing parameter
  # deg = order of differences of the components of b
  # show = indicates whether iteration details should be shown
  
  # Some preparations
  nx <- dim(X)[2]
  D <- diff(diag(nx), diff = deg)
  bstart <- log(sum(y) / nx);
  b <- rep(bstart, nx);
  
  # Perform the iterations
  for (it in 1:50) {
    b0 <- b
    eta <- X %*% b
    gam <- exp(eta)
    mu <- C %*% gam
    w <- c(1 / mu, rep(lambda, nx - deg)) 
    Gam <- gam %*% rep(1, nx)
    Q <- C %*% (Gam * X)
    z <- c(y - mu + Q %*% b, rep(0, nx - deg))
    Fit <- lsfit(rbind(Q, D), z, wt = w, intercept = F) 
    b <- Fit$coef
    
    db <- max(abs(b - b0))
    if (show)  cat(it, " ", db, "\n")
    if (db < 1e-6) break
    
  }
  
  # Regression diagnostics
  R <- t(Q) %*% diag(c(1 / mu)) %*% Q
  H <- solve(R + lambda * t(D) %*% D, R)
  H0 <- solve(R + lambda * t(D) %*% D) # variance-covariance matrix Bayesian approach
  H1 <- H0 %*% R %*% H0 # variance-covaraince matrix sandwitch estimator
  fit <- list()
  fit$trace <- sum(diag(H))
  ok <- y > 0
  fit$dev <- 2 * sum(y[ok] * log(y[ok] / mu[ok]))
  fit$gamma <- gam
  fit$aic <- fit$dev + 2 * fit$trace
  fit$bic <- fit$dev + log(length(y)) * fit$trace
  fit$mu <- mu
  fit$H0 <- H0
  fit$H1 <- H1
  fit$eta <- eta
  fit
}



# PCLM source function to forecast cohorts in T3
pclm.c3 <- function(y, C, X, lambda = 1, deg = 2, show = F,W=c(10,10^8,10,10^2),last.age.c3){
  # Fit a PCLM (estimate b in ) E(y) = C %*% exp(X %*% b)
  # y = the vector of observed counts of length i
  # C = the composition matrix of dimension IxJ
  # X = the identity matrix of dimension JxJ; or B-spline basis
  # lambda = smoothing parameter
  # deg = order of differences of the components of b
  # show = indicates whether iteration details should be shown
  
  # Some preparations
  nx <- dim(X)[2]
  D <- diff(diag(nx), diff = deg)
  D[1,] <- 0 # remove penalty for age 0  
  bstart <- log(sum(y) / nx);
  b <- rep(bstart, nx);
  W1 <- c(rep(1,length(0:last.age.c3)),W) # weights
  
  # Perform the iterations
  for (it in 1:60) {
    b0 <- b
    eta <- X %*% b
    gam <- exp(eta)
    mu <- C %*% gam
    w <- c((1/mu)*W1, rep(lambda, nx - deg)) 
    Gam <- gam %*% rep(1, nx)
    Q <- C %*% (Gam * X)
    z <- c(y - mu + Q %*% b, rep(0, nx - deg))
    Fit <- lsfit(rbind(Q, D), z, wt = w, intercept = F) 
    b <- Fit$coef
    db <- max(abs(b - b0))
    if (show)  cat(it, " ", db, "\n")
    if (db < 1e-8) break
    
  }
  
  # Regression diagnostics
  R <- t(Q) %*% diag(c(1 / mu)*W1) %*% Q 
  H <- solve(R + lambda * t(D) %*% D, R)
  H0 <- solve(R + lambda * t(D) %*% D)
  H1 <- H0 %*% R %*% H0 
  fit <- list()
  fit$trace <- sum(diag(H))
  ok <- y > 0
  fit$dev <- 2 * sum(y[ok] * log(y[ok] / mu[ok]))
  fit$gamma <- gam
  fit$aic <- fit$dev + 2 * fit$trace
  fit$bic <- fit$dev + log(length(y)) * fit$trace
  fit$mu <- mu
  fit$H0 <- H0
  fit$H1 <- H1
  fit$eta <- eta
  fit
}



# Compositional matrix for T2
make.c0 <- function(last.age){
  ilo_1 = seq(1, (length(0:(last.age))+1), by = 1)
  ilo_2 = 121
  ilo = c(ilo_1, ilo_2)
  ihi = c(seq(1,length(0:(last.age)),1),120,131)
  n = length(ihi)
  ihi[n] = 131
  # intervals lengths
  leng <- ihi-ilo+1
  # make C matrix 
  C0 = matrix(0, n, 131)
  for (i in 1:n) C0[i, ilo[i] : ihi[i]] =  1
  colnames(C0) <- 0:130
  return(C0)
} 

# Build C matrix
build_C_matrix <- function(x, y, nlast, offset, out.step, type) {
  # Build C matrix in the age direction
  nx <- length(x)
  gx <- seq(min(x), max(x) + nlast - out.step, by = out.step)
  gu <- c(diff(x), nlast)/out.step
  CA <- matrix(0, nrow = nx, ncol = sum(gu), dimnames = list(x, gx))
  xr <- c(x[-1], max(x) + nlast)
  for (j in 1:nx) CA[j, which(gx >= x[j] & gx < xr[j])] <- 1
  # Build C matrix in the year direction
  if (type == "1D") {
    ny <- length(y)
    CY <- NULL
    C  <- CA
  } else {
    ny <- ncol(y)
    CY <- diag(1, ncol = ny, nrow = ny) 
    C  <- CY %x% CA # Kronecker product
  }
  gy <- 1:ny
  # Output
  out <- as.list(environment())
  return(out)
}

# Compositional matrix for T3
make.c3 <- function(last.age.c3,M.temp){
  ilo_1.c3 = seq(1, (length(0:(last.age.c3))+1), by = 1)
  ilo_2.c3 = 121
  ilo.c3 = c(ilo_1.c3,(M.temp+2),ilo_2.c3)
  ihi.c3 = c(seq(1,length(0:(last.age.c3)),1),(M.temp),120,131)
  n.c3 = length(ihi.c3)
  ihi.c3[n.c3] = 131
  # intervals lengths
  leng.c3 <- ihi.c3-ilo.c3+1
  # make C matrix 
  C0.c3 = matrix(0, n.c3, 131)
  for (i in 1:n.c3) C0.c3[i, ilo.c3[i] : ihi.c3[i]] =  1
  C0.c3.fin <- rbind(C0.c3[1:(length(0:(last.age.c3))+1),],c(rep(0,M.temp),1,rep(0,(130-M.temp))),C0.c3[(length(0:(last.age.c3))+2):nrow(C0.c3),])
  colnames(C0.c3.fin) <- 0:130 
  return(C0.c3.fin)
}



# Function to compute mx from dx
mx_from_dx <- function(dx,ax=NULL){
  ## dimension of dx
  m <- length(dx)
  ## template vectors
  lx <- Lx <- mx <- rep(NA,m)
  ## set the radix of life table
  lx[1] <- sum(dx,na.rm = T)
  ## compute l(x+1)=l(x)-d(x) 
  for (i in 2:m){
    lx[i] <- lx[i-1] - dx[i-1]
  }
  ## set ax = 1/2
  if (is.null(ax)) ax <- rep(1/2,m)
  ## compute Lx = l(x+1) + ax*dx
  Lx[-m] <- lx[-1]+ax[-m]*dx[-m]
  ## compute mx
  mx <- dx/Lx
  ## return mx value
  return(mx)
}



# Function to calculate life tables
lifetable.mx <- function(x, mx, sex="F", ax=NULL){
  m <- length(x)
  n <- c(diff(x), NA)
  if(is.null(ax)){
    ax <- rep(0,m)
    if(x[1]!=0 | x[2]!=1){
      ax <- n/2
      ax[m] <- 1 / mx[m]
    }else{    
      if(sex=="F"){
        if(mx[1]>=0.107){
          ax[1] <- 0.350
        }else{
          ax[1] <- 0.053 + 2.800*mx[1]
        }
      }
      if(sex=="M"){
        if(mx[1]>=0.107){
          ax[1] <- 0.330
        }else{
          ax[1] <- 0.045 + 2.684*mx[1]
        }
      }
      ax[-1] <- n[-1]/2
      ax[m] <- 1 / mx[m]
    }
  }
  qx  <- n*mx / (1 + (n - ax) * mx)
  qx[m] <- 1
  px  <- 1-qx
  lx  <- cumprod(c(1,px))*100000
  dx  <- -diff(lx)
  Lx  <- n*lx[-1] + ax*dx
  lx <- lx[-(m+1)]
  Lx[m] <- lx[m]/mx[m]
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx/lx
  return.df <- data.frame(x, n, mx, ax, qx, px, lx, dx, Lx, Tx, ex)
  return(return.df)
}

