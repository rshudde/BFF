t_density = function(x, nu, mu = 0){
  first_term = sterling_approx((nu+1)/2)
  second_term = (pi * nu)^(1/2) * sterling_approx(nu/2)
  third_term = (1 + x^2 / nu)^(-(nu+1)/2)

  # non centrality terms
  non_central_terms = 1
  if (mu != 0) {
    exp_term = exp(-mu^2/2)
    A_term = hypergeom1F1((nu+1)/2, 1/2, mu^2*x^2/(2*(x^2+nu)))$f
    B_one = sqrt(2) * mu * x / sqrt(x^2 + nu)
    B_two= sterling_approx(nu/2 + 1) / sterling_approx((nu+1)/2)
    B_term = B_one * B_two * hypergeom1F1((nu+1)/2, 3/2, mu^2*x^2/(2*(x^2+nu)))$f
    non_central_terms = exp_term * (A_term + B_term)
  }

  print("here")
  print(first_term, second_term)
  density = first_term / second_term  # * third_term * non_central_terms
  return(density)
}

# Inverse Moment Prior
ginvmom = function(x, tau, nu, k = 1) {
  # x is lambda, k = 1 for z and t tests

  density = k * tau^(nu / 2) * abs(x)^(-nu - 1) * exp(-tau^k / x^(2 * k)) / gamma(nu / (2*k))
  # dens[!is.finite(dens)] = .Machine$double.eps
  return(density)
}
