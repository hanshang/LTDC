###################################################
# Random-walk with and without drift in CoDa space
###################################################

naive_fun <- function(dat, fh, drift_term, level_ci)
{	
    n_year = nrow(dat)
    n_age = ncol(dat)
    
    # standardize life table death to sum to 1
    
    dat_center = sweep(dat, 1, apply(dat, 1, sum), "/")
    
    alpha_x = vector("numeric", n_age)
    for(ik in 1:n_age)
    {
      alpha_x[ik] = geometric.mean(dat_center[,ik])
    }
    
    f_x_t = matrix(NA, n_year, n_age)
    for(ik in 1:n_year)
    {
      f_x_t[ik,] = (dat[ik,]/alpha_x)/sum(dat[ik,]/alpha_x)
    }
    
    g_t = vector("numeric", n_year)
    h_x_t = matrix(NA, n_year, n_age)
    for(ik in 1:n_year)
    {
      g_t[ik] = geometric.mean(f_x_t[ik,])
      h_x_t[ik,] = log(f_x_t[ik,]/g_t[ik])
    }
    
    random_walk = matrix(NA, n_age, fh)
    random_walk_boot = array(NA, dim = c(n_age, 5000, fh))
    for(ik in 1:n_age)
    {
      dum = forecast.lagwalk(object = lagwalk(y = h_x_t[,ik], lag = 1, drift = drift_term), 
                             h = fh, level = level_ci, bootstrap = TRUE)
      random_walk[ik,] = dum$mean
      random_walk_boot[ik,,] = dum$sim
      rm(dum)
    }
    
    # back-transformation (point forecast)
    
    f_x_t_star_fore = d_x_t_star_fore = matrix(NA, n_age, fh)
    for(ik in 1:fh)
    {
      f_x_t_star_fore[,ik] = exp(random_walk[,ik])/sum(exp(random_walk[,ik]))
      d_x_t_star_fore[,ik] = (f_x_t_star_fore[,ik] * alpha_x)/sum((f_x_t_star_fore[,ik] * alpha_x))
    }
    
    # interval forecast
    
    f_x_t_star_fore_boot = d_x_t_star_fore_boot = array(NA, dim = c(n_age, 5000, fh))
    for(ik in 1:fh)
    {
      for(iw in 1:5000)
      {
        f_x_t_star_fore_boot[,iw,ik] = exp(random_walk_boot[,iw,ik])/sum(exp(random_walk_boot[,iw,ik]))
        d_x_t_star_fore_boot[,iw,ik] = (f_x_t_star_fore_boot[,iw,ik] * alpha_x)/sum(f_x_t_star_fore_boot[,iw,ik] * alpha_x)
      }
    }
    
    return(list(fore_count = d_x_t_star_fore * 10^5, 
                fore_lb = apply(d_x_t_star_fore_boot, c(1, 3), quantile, (100 - level_ci)/200) * 10^5,
                fore_ub = apply(d_x_t_star_fore_boot, c(1, 3), quantile, (100 + level_ci)/200) * 10^5))
}    

