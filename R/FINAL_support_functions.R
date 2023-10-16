
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
sterling_gamma = function(n)
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
