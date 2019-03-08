.onLoad <- function(libname, pkgname){
  library.dynam("barsN", package = pkgname, "aspline/R")
  source("barsN_Rwrapper")
}
# .onUnLoad <- function(libname, pkgname){
#   library.dynam.unload("barsN", package = c(pkgname), )
#   source("barsN_Rwrapper")
# }
