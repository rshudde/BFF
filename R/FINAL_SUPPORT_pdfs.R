###################Prior and Posterior pdfs

########## One-sided t-test

######Prior

#density of the normal-moment random variable
dens_function_nm <- function(x, tau2, r) {
  dens = (abs(x))^(2*r) * exp(-x^2/(2*tau2))/((2*tau2)^(r+0.5)*gamma(r+0.5))
  return(dens)
}

#making sure it integrates to 1
dnm_1s <- function(x, tau2, r) {
  density = ifelse(x >= 0, dens_function_nm(x, tau2, r)*2, 0)
  return(density)
}

prior_plotting_func_1s <- function(x, tau2_w_r, r){
  prior_data <- data.frame(x = x, dens = dnm_1s(x = x,  tau2 = tau2_w_r, r = r))
  p <- ggplot(prior_data, aes(x = x, y = dens))+
    geom_line()+
    xlab(expression(lambda))+
    ylab("Prior Density")+
    theme_bw()
  print(p)
}

# t_stat = 1.5 ##Example
# n=100
# effect_size = seq(0, 1, by = 0.01)
# tau2 = get_one_sample_tau2(n = 100, w = effect_size, r =1 )
# df = n-1
# x = seq(0, 3, length.out = 1000)
#
# logbf <- t_val_r1(tau2 = tau2, t_stat = t_stat, df = df)
# omega <- effect_size[which.max(logbf)]
#
# tau2_w_r = n*omega^2/2
# mode_result_t1 <- optimize(function(x) dnm_1s(x = x, tau2 = n*omega^2/2 , r = 1), interval = c(0, 3), maximum = TRUE)

# prior_plotting_func_1s(x = x, tau2_w_r = tau2_w_r, r = 1)

######Posterior
#function for likelihood
dt_1s <- function(t, df, lambda){
  ifelse(t>0, 2*dt(x = t, df = df, ncp = lambda), 0)
}

#marginal under null
m0_1s_t <- function(x, df, lambda){
  m0_1s_t = dt(x, df = df, ncp = lambda)
  return(m0_1s_t)
}

#I1 in closed form expression of marginal under alternative (in supplemental material)
I1 <- function(tau2, r, df, t){
  c = 1/(2*((1+tau2)^(r+0.5)))
  gauss = Gauss2F1(a = (df+1)/2, b = (r + 0.5), c = 0.5, x = (tau2*(t^2))/((t^2 + df)*(tau2 + 1)))
  I1 = c*gauss
  return(I1)
}

#I2 in closed form expression of marginal under alternative (in supplemental material)
I2 <-  function(tau2, r, df, t){
  c = ((t*sqrt(tau2))*(gamma(df/2 + 1))*gamma(r+1))/((sqrt(t^2 + df))*((tau2+1)^(r+1))*gamma((df+1)/2)*gamma(r + 0.5))
  gauss = Gauss2F1(a = ((df/2)+1), b = (r+1), c = 1.5, x = (((t^2)*tau2)/((t^2 + df)*(1+tau2))))
  I2 = c*gauss
  return(I2)
}

#posterior pdf function
post_pdf_t_1s <- function(t, tau2_w_r, r, df, lambda){
  prior_1s_t = dnm_1s(x = t, tau2 = tau2_w_r, r = r)
  m0_1s_t_vals = m0_1s_t(x = t, df = df, lambda = lambda)
  I1_t1 = I1(tau2 = tau2_w_r, r = r, df = df, t = t)
  I2_t1 = I2(tau2 = tau2_w_r, r = r, df = df, t = t)
  m1_1s_t_vals = 2*m0_1s_t_vals*(I1_t1 + I2_t1)
  posterior_t_1s = (dt_1s(t = t, df = df, lambda = lambda) * prior_1s_t)/m1_1s_t_vals
  return(posterior_t_1s)
}


# t = seq(0, 3, length.out = 1000)
# t_stat = 1.5
# n=100
# effect_size = seq(0.01, 1, by = 0.01)
# tau2 = get_one_sample_tau2(n = 100, w = effect_size, r =1 )
# df = n-1
# r = 1
# lambda = 2
# logbf = t_val_r1(tau2 = tau2, t_stat = t_stat, df = df)
# omega = effect_size[which.max(logbf)]
# tau2_w_r = n*omega^2/2

# prior_1s_t = dnm_1s(x = t, tau2 = tau2_w_r, r = 1)
# posterior_t1_1s = post_pdf_t_1s(t = t, tau2_w_r = tau2_w_r, r = 1, df = df, lambda = 2)
#
# df = data.frame(x = t, prior = prior_1s_t, post = posterior_t1_1s)
# ggplot(df, aes(x = x))+
#   geom_line(aes(y = prior), color = "black")+
#   geom_line(aes(y = post), color = "blue")+
#   labs(x = expression(lambda))




########## Two-sided t-test

##Prior
dnm_2s <- function(x, tau2, r){
  density = dens_function_nm(x, tau2, r)
  return(density)
}

prior_plotting_func_2s <- function(x, tau2_w_r, r){
  prior_data <- data.frame(x = x, dens = dnm_2s(x = x,  tau2 = tau2_w_r, r = r))
  p <- ggplot(prior_data, aes(x = x, y = dens))+
    geom_line()+
    xlab(expression(lambda))+
    ylab("Prior Density")+
    theme_bw()
  print(p)
}


# n = 100
# t_stat = 1.5
# effect_size = seq(0.01, 1, by = 0.01)
# df = n - 2 ##from FINAL_FUNCTIONS_t_test.R
# n1 = n/2
# n2 = n/2
# x = seq(-3, 3, length.out = 1000)
# r=1
# tau2 = get_two_sample_tau2(n1 = n1,
#                            n2 = n2,
#                            w = effect_size,
#                            r = r)
#
# logbf <- t_val_r1(tau2 = tau2, t_stat = t_stat, df = df)
# omega <- effect_size[which.max(logbf)]
# tau2_w_r = n1*n2*(omega^2)/(2*(n1+n2))

##Posterior

dt_2s <- function(t, df, lambda){
  dt(x = t, df = df, ncp = lambda)
}

m0_2s_t <- function(x, df, lambda){
  m0_2s_t = dt(x, df = df, ncp = lambda)
  return(m0_2s_t)
}

#Requires I1 and I2 from above

#Function for pdf of 2 sided t test
post_pdf_t_2s <- function(t, tau2_w_r, r, df, lambda){
  prior_2s_t = dnm_2s(x = t, tau2 = tau2_w_r, r = r)
  m0_2s_t_vals = m0_2s_t(x = t, df = df, lambda = lambda)
  I1_t2 = I1(tau2 = tau2_w_r, r = r, df = df, t = t)
  I2_t2 = I2(tau2 = tau2_w_r, r = r, df = df, t = t)
  m1_2s_t_vals = m0_2s_t_vals*(I1_t2 + I2_t2)
  posterior_t2 = (dt_2s(t = t, df = df, lambda = lambda) * prior_2s_t)/m1_2s_t_vals
  return(posterior_t2)
}

#
# prior_2s_t = dnm_2s(x = t, tau2 = tau2_w_r, r = 1)
# posterior_t_2s = post_pdf_t_2s(t = t, tau2_w_r = tau2_w_r, r = 1, df = df, lambda = 2)
#
#
# df = data.frame(x = t, prior = prior_2s_t, post = posterior_t_2s)
# ggplot(df, aes(x = x))+
#   geom_line(aes(y = prior), color = "black")+
#   geom_line(aes(y = post), color = "blue")+
#   labs(x = expression(lambda))
