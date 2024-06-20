
################# T functions if r is a fraction
log_T_frac = function(tau2, t, v, r)
{
  tp1 = 1 + tau2 # one plus tau^2
  tau = sqrt(tau2)
  c = 1 / (tp1 ^ (r + 1 / 2)) # c
  y = tau * t / sqrt((t ^ 2 + v) * tp1)

  a1 = (v + 1) / 2
  b1 = r + 1 / 2
  c1 = 1 / 2
  first_hypergeo = Gauss2F1(a1, b1, c1, y ^ 2)

  four = sterling_gamma(v / 2 + 1) * sterling_gamma(r + 1)
  five = sterling_gamma((v + 1) / 2) * sterling_gamma(r + 1 / 2)
  gamma_term = four / five

  aa = v / 2 + 1
  bb = r + 1
  cc = 3 / 2
  second_hypergeo = Gauss2F1(aa, bb, cc, y ^ 2)

  to_return = c * (first_hypergeo + 2*y * gamma_term * second_hypergeo)
  to_return = log(to_return)

  return(to_return)
}

log_T_frac_onesided = function(tau2, t, v, r)
{
  tp1 = 1 + tau2
  tau = sqrt(tau2)
  c = 1 / (tp1 ^ (r + 1 / 2)) # c
  y = tau * t / sqrt((t ^ 2 + v) * tp1)

  a1 = (v + 1) / 2
  b1 = r + 1 / 2
  c1 = 1 / 2
  first_hypergeo = Gauss2F1(a1, b1, c1, y ^ 2)

  four = sterling_gamma(v / 2 + 1) * sterling_gamma(r + 1)
  five = sterling_gamma((v + 1) / 2) * sterling_gamma(r + 1 / 2)
  gamma_term = four / five

  aa = v / 2 + 1
  bb = r + 1
  cc = 3 / 2
  second_hypergeo = Gauss2F1(aa, bb, cc, y ^ 2)

  to_return = c * (first_hypergeo + 2 * y * gamma_term * second_hypergeo)
  to_return = log(to_return)

  return(to_return)
}
