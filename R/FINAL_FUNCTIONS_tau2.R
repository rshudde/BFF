###########################################################################################################
###########functions to set tau2 - user does not interact with these ######################################
###########################################################################################################
get_one_sample_tau2 = function(n, w, r)
{
  to_return = n * w ^ 2 / (2 * r)
  return(to_return)
}

get_two_sample_tau2 = function(n1, n2, w, r)
{
  to_return = n1 * n2 * w ^ 2 / (2 * r *(n1 + n2))
  return(to_return)
}

get_count_tau2 = function(n, w, k, r)
{
  top = n * w ^ 2 * k
  bottom =2 * (k / 2 + r - 1)
  to_return = top / bottom
  return(to_return)
}

get_LRT_tau2 = function(n, w, k, r)
{
  top = n * k * w ^ 2
  bottom = 2 * (k / 2 + r - 1)
  to_return = top / bottom
  return(to_return)
}

get_linear_tau2 = function(n, w, k, r)
{
  top = n * k * w ^ 2
  bottom = 4 * (k / 2 + r - 1)
  to_return = top / bottom
  return(to_return)
}

# not currently using
get_corr_tau2 = function(n, w)
{
  top = (n-2) * w^2
  bottom = 2
  to_return = top / bottom
  return(to_return)
}

get_regression_tau2 = function(n, k, w, r)
{
  top = (n-k-1)*w^2
  bottom = (2* r)
  to_return = top / bottom
  return(to_return)
}

get_two_sample_invm_t_tau2 = function(n1, n2, w, nu)
{
  to_return = n1*n2*w^2*(nu+1)/(n1+n2) # M vector
  return(to_return)
}

get_one_sample_invm_t_tau2 = function(n, w, nu)
{
  to_return = n*w^2*(nu+1)/(2) # M vector
  return(to_return)
}

################## for fractional cases
# get_wbar = function(w, k)
# {
#   # a = 1/k
#   # b = sum(w^2)
#   # to_return = sqrt(a * b)
#
#   to_return = w^2
#
#   return(to_return)
# }
# get_tau_z_t_one_sample_frac = function(n, w, r)
# {
#   to_return = n * w^2 / (2*r)
#   return(to_return)
# }

# get_tau_z_t_two_sample_frac = function(n1, n2, w, r)
# {
#   N = n1 + n2
#   to_return = n1 * n2 * w^2 / (2*r*N)
#   return(to_return)
# }
#
# get_tau_poisson_frac = function(n, k, w, r)
# {
#
#   w_bar = get_wbar(w, k)
#   num  = n * k * w_bar^2
#   denom = 2 * (k/2 + r - 1)
#
#   to_return = num/denom
#   return(to_return)
# }
#
# get_tau_linear_frac = function(n, k, w, r)
# {
#
#   w_bar = get_wbar(w, k)
#   num  = n * k * w_bar^2
#   denom = 4 * (k/2 + r - 1)
#
#   to_return = num/denom
#   return(to_return)
# }
#
# get_tau_likelihood_frac = function(n, k, w, r)
# {
#   w_bar = get_wbar(w, k)
#   num  = n * k * w_bar^2
#   denom = 2 * (k/2 + r - 1)
#
#   to_return = num/denom
#   return(to_return)
# }
