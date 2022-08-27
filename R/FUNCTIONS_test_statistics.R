# functions to set tau2 - user does not interact with these
get_one_sample_tau2 = function(n, w)
{
  to_return = n * w^2 / 2
  return(to_return)
}

get_two_sample_tau2 = function(n1, n2, w)
{
  to_return = n1*n2*w^2 / (n1 + n2)
  return(to_return)
}

get_count_tau2 = function(n, w)
{
  to_return = n*w^2
  return(to_return)
}

get_LRT_tau2 = function(n, w)
{
  to_return = n*w^2
  return(to_return)
}

get_linear_tau2 = function(n, w)
{
  to_return = n*w^2 / 2
  return(to_return)
}

# user accessed functions
BFF_z_test = function(z_stat, effect_size, n, one_sample = TRUE, n1 = NULL, n2 = NULL, tau2 = NULL)
{
  # get correct tau^2 value
  if (is.null(tau2))
  {
    if (one_sample)
    {
      tau2 = get_one_sample_tau2(n=n, w=effect_size)
    } else {
      if (is.null(n)) n = n1
      tau2 = get_two_sample_tau2(n1=n, n2=n2, w=effect_size)
    }
  }

  term_one = (tau2 + 1)^(-3/2)
  term_two = 1 + tau2 * z_stat^2 / (tau2 + 1)
  term_three = exp(tau2 * z_stat^2/ (2 * (tau2 + 1)))

  to_return = term_one * term_two * term_three
  return(to_return)
}

BFF_t_test = function(t_stat, df, effect_size, n, one_sample = TRUE, n1 = NULL, n2 = NULL, tau2 = NULL)
{
  # get correct tau^2 value
  if (is.null(tau2))
  {
    if (one_sample)
    {
      tau2 = get_one_sample_tau2(n=n, w=effect_size)
    } else {
      if (is.null(n)) n = n1
      tau2 = get_two_sample_tau2(n1=n, n2=n2, w=effect_size)
    }
  }

  r = 1 + t_stat^2 / df
  s = 1 + t_stat^2 / (df*(1+ tau2))
  q = tau2 * (df + 1) / (df*(1+ tau2))

  BF = (tau2 + 1)^(-3/2) * (r/s)^((df+1)/2) * (1 + q*t_stat^2 / s)

  to_return = BF
  return(to_return)
}

BFF_chi2_test = function(chi_stat, df, effect_size, n, count = TRUE, tau2 = NULL)
{

  if (is.null(tau2))
  {
    if (count)
    {
      tau2 = get_count_tau2(n=n, w=effect_size)
    } else {
      tau2 = get_LRT_tau2(n=n, w=effect_size)
    }
  }

  BF = (tau2 + 1)^(-df/2 - 1) * (1 + tau2 * chi_stat / (df * (tau2 + 1))) * exp(tau2*chi_stat / (2*(tau2+1)))

  to_return = BF
  return(to_return)
}

BFF_F_test = function(f_stat, df1, df2, effect_size, n, tau2 = NULL)
{

  if (is.null(tau2)) tau2 = get_linear_tau2(n, effect_size)

  v = df2 * (tau2 + 1)
  term_one = (tau2 + 1)^(-df1/2 - 1)
  term_two = (1 + df1 * f_stat / df2) / (1 + df1 * f_stat / v)
  term_three = 1 + (df1 + df2) * tau2 * f_stat / (v * (1+df1 * f_stat/v))

  to_return = term_one * term_two * term_three
  return(to_return)
}
