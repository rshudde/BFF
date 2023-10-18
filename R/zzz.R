# # .onLoad <- function(libname, pkgname)
# # {
# #   library.dynam("BFF", pkgname, libname)
# # }
#
BFFStartupMessage <- function()
{
  # Startup message obtained as
  # > figlet -f slant BFF
  msg <- c(paste("version", utils::packageVersion("BFF")), "\nType 'citation(\"BFF\")' for citing this R package in publications")
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  # unlock .BFF variable allowing its modification
  # unlockBinding(".BFF", asNamespace("BFF"))
  # startup message
  msg <- BFFStartupMessage()
  if (!interactive())
    msg[1] <- paste("Package 'BFF' version", utils::packageVersion("BFF"), "for Bayesian hypothesis testing.")
  packageStartupMessage(msg)
  invisible()
}

# .onLoad = function(libname, pkgname)
# {
#   msg = BFFStartupMessage()
#   packageStartupMessage(msg)
#   invisible()
# }
