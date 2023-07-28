# [Funcion para generar direcciones hibridas Xi = (Psi,Theta)] ===============.#
Xi_gen = function(n_fun,n_vec,n_grid,conf) {
  # n_fun: Tamaño de la parte funcional
  # n_vec: Tamaño de la parte vectorial
  # n_grid: numero de puntos en la grilla
  # return(): multiFunData class base de tamaño n_vec
  # Generacion direccion vectorial Theta
  Mat_aux = matrix(0.3,nrow = n_vec,ncol = n_vec)
  diag(Mat_aux) = rep(1,n_vec)
  Theta = t(eigen(Mat_aux)$vectors)
  # Generacion direccion funcional Psi
  arg_grill = (1:n_grid - 1)/(n_grid - 1)
  base_multi = simMultiFunData(type = "weighted",
                                argvals = rep(list(list(arg_grill)),n_fun),
                                M = rep(n_vec,n_fun), eFunType = conf, 
                                eValType = "linear", N = 1,ignoreDeg = rep(1,n_fun))$trueFuns
  Xi = mfh_data(fun_data = base_multi,vec_data = Theta)
  Xi = (1/sqrt(2))*Xi  
  return(Xi)
}
# [Fin del codigo] ===========================================================.#