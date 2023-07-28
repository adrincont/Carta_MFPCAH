# =============================================================================#
# [Libraries                                                                 ] #
# [Comments] -----------------------------------------------------------------.#
# ============================================================================.#
# {Libraries} =============================================================.####
rm(list = ls())
install_load = function(package1, ...)  {   
  # convert arguments to vector
  packages = c(package1, ...)
  # start loop to determine if each package is installed
  for (package in packages) {
    # if package is installed locally, load
    #if (package %in% rownames(installed.packages()))
    if (TRUE)
      do.call('library', list(package))
    # if package is not installed locally, download, then load
    else {
      install.packages(package)
      do.call("library", list(package))
    }
  } 
}

rmvgamma = function(n,p,shape,rate=1/scale, rho=0, scale=1){
  res_i = matrix(data = rgamma(n*p,shape = shape,scale = scale),nrow = n,ncol = p)
  return(res_i)
  #stopifnot(length(rho) == 1,length(shape) == 1,length(rate) == 1,
  #          length(scale) == 1,length(n) == 1,length(p) == 1,
  #          rate > 0, shape > 0, rho >= 0, rho <= 1, scale > 0, n >= 0, p > 0)
  #n = round(n);p = ceiling(p)
  #theta = rate*rho/(1 - rho)
  #k = rnbinom(n,shape,rate/(rate + theta))
  #return(matrix(rgamma(p*n,shape + k,rate + theta),n))
}

list_lib = c(
  # Operaciones con datos funcionales
  "fda.usc",
  "fda",
  "MFPCA",
  "FactoMineR",
  # Graficas
  "ggplot2",
  "dplyr",
  "GGally",
  "patchwork",
  # Tratamiento de datos
  "R.matlab",
  "gridExtra",
  "readxl",
  "tidyr",
  "sp",
  "ggdensity",
  "foreach",
  "mvtnorm"
  )
# Load libraries and install if they do not exist
install_load(list_lib)
#theme_set(theme_bw())
#tab_lib = installed.packages()
#tab_lib = tab_lib[rownames(tab_lib) %in% c("comoOdeCpp",list_lib),]
#write.csv(tab_lib,file = "info_library.csv")
# {End of code} ===========================================================.####