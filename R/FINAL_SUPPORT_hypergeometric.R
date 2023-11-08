
# library(gsl)
Gauss2F1 <- function(a,b,c,x){} # at bottom of file


#hypergeom1F1: Confluent hypergeometric function
#In gajdosandrej/CharFunToolR: Numerical Computation Cumulative Distribution Function and Probability Density Function from Characteristic Function

hypergeom1F1 <- function(a, b, z, n) {
  ## CHECK THE INPUT PARAMETERS
  if(missing(n)) {
    n <- numeric()
  }

  if(length(n) == 0) {
    n <- 64
  }

  done <- FALSE
  transf <- FALSE
  szz <- dim(z)
  z <- c(z)
  sz <- length(z)
  f <- rep(NaN, sz)
  method <- -rep(1, sz)
  loops <- rep(0, sz)

  # 1F1(a,b,z) for special cases of the parameters a and b and any argument z
  if(b == 0 | b == -as.integer(abs(b))) {
    f <- Inf
    done <- TRUE
  } else if(a == b) {
    f <- exp(z)
    done <- TRUE
  } else if(a - b == 1) {
    f <- (1 + z / b) * exp(z)
    done <- TRUE
  } else if(a == 0) {
    f <- 1
    done <- TRUE
  } else if(a == -1) {
    f <- 1 - z / b
  } else if(a == 1 && b == 2) {
    f <- (exp(z) - 1) / z
  } else if(a == as.integer(a) && a < 0) {
    m <- -a
    cr <- 1
    f <- 1
    for(k in 1:m) {
      cr <- cr * (a + k - 1) / k / (b + k - 1) * z
      f <- f + cr
    }
    done <- TRUE
  }

  # 1F1(a,b,z) for other cases of the parameters a ,b and the argument z
  if(!done) {
    # If b < a  set 1F1(a,b,z) = exp(z)*1F1(b-a,b,-z)
    if(b < a){
      transf <- TRUE
      a <- b - a
      z <- -z
    }

    im_z <- Im(z)
    re_z <- Re(z)
    ind0 <- z == 0
    ind1 <- ((abs(z) < 10 & abs(im_z) < abs(re_z)) | abs(z) < 20 + abs(b) | a < 0) & !ind0
    ind2 <- (abs(z) >= 10) & (abs(im_z) >= abs(re_z)) & (a > 0) & (b > a)
    ind3 <- (!ind0 & !ind1 & !ind2)
    # z == 0 set 1F1(a,b,0) = 1
    if(any(ind0 != 0)) {
      f[ind0] <- 1
      method[ind0] <- 0
    }
    # 1F1(a,b,z) for small abs(z) or negative a:
    # abs(z) < 10 & abs(im_z) < abs(re_z) & ~ind0, OR
    # abs(z) < 20 + abs(b) & ~ind0, OR
    # a < 0
    if(any(ind1 != 0)) {
      chg <- 1
      crg <- 1
      chw <- 0
      zz <- z[ind1]
      for(j in 1:500) {
        crg <- crg * (a + j - 1) / (j * (b + j -1)) * zz
        chg <- chg + crg
        if(all(abs((chg - chw) / chg) < 1e-15)) {
          break
        }
        chw <- chg
      }
      method[ind1] <- 1
      loops[ind1] <- j
      f[ind1] <- chg
    }
    # 1F1(a,b,z) for large abs(z) such that abs(imag(z)) >= abs(re_z)
    # and a > b > 0 by using the steepest descent integration:
    # abs(z) >= 10 & abs(im_z) >= abs(re_z) & a > 0 & b > a
    if(any(ind2 != 0)) {
      x_w <- GaussLaguerre(n)
      x <- x_w$x
      w <- x_w$w
      gba <- log(gamma(b)) - (log(gamma(a)) + log(gamma(b - a)))
      rez <- re_z[ind2]
      imz <- im_z[ind2]
      ewa <- 1 / imz
      ewb <- exp(imz * 1i) / imz
      a1 <- a - 1
      ba1 <- b - a - 1
      r1 <- 0
      r2 <- 0
      for(j in 1:n) {
        x_i <- x[j] / imz * 1i
        aux1 <- rez * x_i + a1 * log(x_i) + ba1 * log(1 - x_i)
        aux2 <- rez * (1 + x_i) + a1 * log(1+x_i) + ba1 * log(-x_i)
        r1 <- r1 + w[j] * exp(gba + aux1)
        r2 <- r2 + w[j] * exp(gba + aux2)
      }
      method[ind2] <- 2
      loops[ind2] <- j
      f[ind2] <- (ewa * r1 - ewb * r2) * 1i
    }
    # 1F1(a,b,z) for (otherwise) large z by using asymptotic expansion
    if(any(ind3 != 0)) {
      zz <- z[ind3]
      g1 <- gamma(a)
      g2 <- gamma(b)
      ba <- b - a
      g3 <- gamma(ba)
      cs1 <- 1
      cs2 <- 1
      cr1 <- 1
      cr2 <- 1
      for(j in 1:500) {
        cr1 <- -cr1 * (a + j - 1) * (a - b + j) / (zz * j)
        cr2 <- cr2 * (b - a + j - 1) * (j - a) / (zz * j)
        cs1 <- cs1 + cr1
        cs2 <- cs2 + cr2
        if(all(abs(cr1+cr2) < 1e-15)) {
          break
        }
      }
      x <- Re(zz)
      y <- Im(zz)
      phi <- atan(y / x)
      phi[x == 0 & y > 0]  <- 0.5 * pi
      phi[x == 0 & y <= 0] <- -0.5 * pi
      ns <- rep(1,length(x))
      ns[phi > -1.5*pi & phi <= -0.5*pi] <- -1
      cfac <- exp(1i * pi * a * ns)
      cfac[y == 0] <- cos(pi * a)
      chg1 <- (g2 / g3) * zz^(-a) * cfac * cs1
      chg2 <- (g2 / g1) * exp(zz) * zz^(a-b) * cs2
      chg <- chg1 + chg2
      method[ind3] <- 3
      loops[ind3] <- j
      f[ind3] <- chg
    }
  }

  if(transf) {
    f <- exp(-z) * f
  }

  dim(f) <- szz
  dim(method) <- szz

  return(list("f" = f, "method" = method, "loops" = loops))

}


GaussLaguerre <- function(n, alpha) {
  # GaussLaguerre evaluates the Gauss-Laguerre Nodes and Weights on the interval (alpha,Inf).

  if(missing(alpha)) {
    alpha <- numeric()
  }

  if(length(alpha) == 0) {
    alpha <- 0
  }

  if(n == 64) {
    x <- c(2.241587414670593e-02, 1.181225120967662e-01, 2.903657440180303e-01,
           5.392862212279714e-01, 8.650370046481124e-01, 1.267814040775241e+00,
           1.747859626059435e+00, 2.305463739307505e+00, 2.940965156725248e+00,
           3.654752650207287e+00, 4.447266343313093e+00, 5.318999254496396e+00,
           6.270499046923656e+00, 7.302370002587399e+00, 8.415275239483027e+00,
           9.609939192796107e+00, 1.088715038388638e+01, 1.224776450424431e+01,
           1.369270784554751e+01, 1.522298111152473e+01, 1.683966365264873e+01,
           1.854391817085919e+01, 2.033699594873023e+01, 2.222024266595088e+01,
           2.419510487593325e+01, 2.626313722711848e+01, 2.842601052750102e+01,
           3.068552076752596e+01, 3.304359923643782e+01, 3.550232389114120e+01,
           3.806393216564646e+01, 4.073083544445863e+01, 4.350563546642153e+01,
           4.639114297861618e+01, 4.939039902562468e+01, 5.250669934134629e+01,
           5.574362241327837e+01, 5.910506191901708e+01, 6.259526440015138e+01,
           6.621887325124754e+01, 6.998098037714681e+01, 7.388718723248294e+01,
           7.794367743446311e+01, 8.215730377831930e+01, 8.653569334945649e+01,
           9.108737561313303e+01, 9.582194001552071e+01, 1.007502319695140e+02,
           1.058845994687999e+02, 1.112392075244396e+02, 1.168304450513065e+02,
           1.226774602685386e+02, 1.288028787692377e+02, 1.352337879495258e+02,
           1.420031214899315e+02, 1.491516659000494e+02, 1.567310751326712e+02,
           1.648086026551505e+02, 1.734749468364243e+02, 1.828582046914315e+02,
           1.931511360370729e+02, 2.046720284850595e+02, 2.180318519353285e+02,
           2.348095791713262e+02)

    w = c(5.625284233902887e-02, 1.190239873124205e-01, 1.574964038621475e-01,
          1.675470504157746e-01, 1.533528557792381e-01, 1.242210536093313e-01,
          9.034230098648389e-02, 5.947775576835545e-02, 3.562751890403607e-02,
          1.948041043116659e-02, 9.743594899382018e-03, 4.464310364166234e-03,
          1.875359581323119e-03, 7.226469815750032e-04, 2.554875328334960e-04,
          8.287143534397105e-05, 2.465686396788568e-05, 6.726713878829501e-06,
          1.681785369964073e-06, 3.850812981546759e-07, 8.068728040991898e-08,
          1.545723706757564e-08, 2.704480147613762e-09, 4.316775475431567e-10,
          6.277752541794292e-11, 8.306317376250609e-12, 9.984031787119531e-13,
          1.088353887008957e-13, 1.074017402970290e-14, 9.575737246084761e-16,
          7.697028063946171e-17, 5.564881054436309e-18, 3.609756216814263e-19,
          2.095095239662055e-20, 1.084792493734732e-21, 4.994712583627291e-23,
          2.037932077329677e-24, 7.340603648778086e-26, 2.324586950075985e-27,
          6.464528714804253e-29, 1.582906573670680e-30, 2.881154588676925e-32,
          5.412535994048359e-34, 2.241048506440640e-34, 5.283742844838896e-36,
          1.338202299148180e-34, 1.586340564468588e-35, 1.215192241351559e-34,
          1.217775335792122e-34, 1.673365556974291e-35, 2.735714461640009e-34,
          2.185380020634853e-34, 5.648495554594729e-35, 9.997398610925997e-36,
          1.500478177990158e-36, 1.416076744376295e-37, 5.444799396304293e-39,
          1.153008451969226e-40, 2.474260963687568e-42, 9.293338889710336e-45,
          2.974573897074668e-47, 1.941796748832940e-50, 5.776547415033449e-54,
          1.794991571658772e-58)

    return(list("x" = x, "w" = w))

  } else {
    idx <- 1:n
    a <- (2 * idx - 1) + alpha
    b <- sqrt(idx[1:(n-1)] * ((1:(n-1)) + alpha)) # bandSparse() Matrix package; sdiag() mgcv package
    CM <- diag(a) + as.matrix(matrix(0, n, n) + Matrix::bandSparse(n, n, c(1, -1), list(b, b)))

    eig <- eigen(CM)
    V <- eig$vectors
    L <- diag(eig$values)
    ind <- idxInOrigSeq(diag(L))
    x <- sort(diag(L))
    V <- t(Conj(V[,ind]))
    w <- gamma(alpha + 1) * V[,1] ^ 2

    return(list("x" = x, "w" = w))

  }
}

idxInOrigSeq <- function(orig_seq) {
  sort_seq <- sort(orig_seq)
  indices <- vector()
  used_indices <- vector()
  for(i in 1:length(orig_seq)) {
    idx <- 0
    for(j in 1:length(orig_seq)) {
      if(sort_seq[i] == orig_seq[j] && !j %in% used_indices) {
        idx <- j
        used_indices <- c(used_indices, j)
        break
      }
    }
    indices <- c(indices, idx)
  }
  return(indices)

}

Gauss2F1 <- function(a, b, c, x) {
  condition <- x >= 0 & x < 1
  to_return <- ifelse(condition,
                      gsl::hyperg_2F1(a, b, c, x),
                      gsl::hyperg_2F1(c - a, b, c, 1 - 1 / (1 - x)) / (1 - x)^b
  )
  return(to_return)
}
