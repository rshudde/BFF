## math support functions
double_factorial_even = function(n) {
  first = 2^(n/2)
  second = factorial(n/2)
  to_return = first * second
  return(to_return)
}

double_factorial_odd = function(n) {
  n1 = n+1
  numerator = factorial(n1)
  first = 2^(n1/2)
  second = factorial(n1/2)
  denomonator = first*second

  to_return = numerator / denomonator
  return(to_return)
}

# double factorial expression
double_factorial = function(n) {
  if (n %% 2 == 1) {
    to_return = double_factorial_odd(n)
  } else {
    to_return = double_factorial_even(n)
  }

  return(to_return)
}

# approximation of gamma for large n
gamma_approx = function(n)
{
  # if (n<7)
  # {
  #   const = gamma(n)
  # } else {
  #   const = n*log(n) - n
  # }

  return(exp(lgamma(n)))

  # return(const)
}


## functions used for the cases when r > 1
get_a = function(tau2, r) {
  to_return = (1 + tau2)^(-(r+1/2))
  return(to_return)
}

get_b = function(tau2, r, k) {
  to_return = (1 + tau2)^(-(r+k/2))
  return(to_return)
}

get_c = function(tau2, df, r) {
  numerator = gamma_approx(df/2+1) * gamma_approx(r + 1)
  denomonator = gamma_approx((df+1)/2) * gamma_approx(r + 1/2)
  to_return = numerator / denomonator
  return(to_return)
}

get_y_onesided_z_test = function(tau2, z) {
  to_return = sqrt(tau2) * z *(2 + 2*tau2)^(-1/2)
  return(to_return)
}

get_y_t_test = function(tau2, t, df) {
  inner = (df + t^2) * (1+tau2)
  to_return = sqrt(tau2) * t * inner^(-1/2)
  return(to_return)
}





