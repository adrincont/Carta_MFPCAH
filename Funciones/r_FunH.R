# [Simulacion de observaciones funcionales hibridas] =========================.#
# Ruido hibrido
epsilon_H = function(Z_j,sd){
  mean_i = rep(0,length(Z_j@funData)) 
  ei = mvrnorm(n = length(Z_j@funData[[1]]@argvals[[1]]),mu = mean_i,Sigma = sd)
  vec = Z_j@vecData*0
  fun_list = list()
  for (i in 1:length(Z_j@funData)) {
    Ei = matrix(ei[,i],nrow = dim(Z_j@funData[[i]]@X)[1],ncol = dim(Z_j@funData[[i]]@X)[2],byrow = TRUE)
    fun_list[[i]] = funData(Z_j@funData[[i]]@argvals[[1]],Ei)
  }
  fun = multiFunData(fun_list)
  vec = Z_j@vecData*0
  ei = mfh_data(fun_data = fun,vec_data = vec)
  return(ei)
}
# distribución t estandarizada
rt_modified = function(N, nu, mu = 0, standard_dev){
  x1 = rt(N, nu) # 1
  x2 = x1/sqrt(nu/(nu - 2)) # 2
  x3 = x2 * standard_dev # 3
  x4 = x3 + mu # 4
  return(x4)
}
r_scores = function(N,Sigma_k,nu,dist) {
  var_i = diag(Sigma_k)
  matriz_t = matrix(nrow = N, ncol = length(var_i))
  if (dist == "norm") {
    matriz_t = MASS::mvrnorm(n = N,mu = rep(0,length(var_i)),Sigma = Sigma_k)
    return(matriz_t)
  }else {
    for (i in 1:length(var_i)) {
      if (dist == "t") {
        matriz_t[,i] = rt_modified(N = N,nu = nu,standard_dev = sqrt(var_i[i]))
      }
      if (dist == "gamma") {
        matriz_t[,i] = rgamma(n = N,shape = nu,scale = sqrt(var_i[i]/nu))
        matriz_t[,i] = matriz_t[,i] - mean(matriz_t[,i])
      }
    }
    return(matriz_t)
  }
}
# Generación de las observaciones
r_FunH = function(n_fun,n_vec,n_grid,N,conf, eValType = "exponential", 
                  dist_i = 'norm'){
  N = N + 10
  # Generacion de direccion principal Xi = (Psi,Theta)
  n_vec_aux = max(c(n_fun*2,10,n_vec))
  Xi = Xi_gen(n_fun = n_fun,n_vec = n_vec_aux,n_grid = n_grid,conf = conf)
  # Generacion de las observaciones
  ## Matriz de covarianza de los scores
  Sigma_k = matrix(0,nrow = n_vec_aux,ncol = n_vec_aux)
  diag(Sigma_k) = eVal(n_vec_aux, eValType)
  if (dist_i == 't') { nu_i = 5 }
  if (dist_i != 't') { nu_i = 3 }
  scor_i = r_scores(N = N,Sigma_k = Sigma_k,nu = nu_i,dist = dist_i)
  dir_i = Xi
  ## Generación de observaciones
  Z = Inv_PCA_H(dir_j = dir_i,scor_j = scor_i)
  if (n_vec == 1) {
    Z@vecData = matrix(c(Z@vecData[,1:n_vec]),ncol = 1)
  }else{
    Z@vecData = Z@vecData[,1:n_vec]
  }
  return(select_Mfun(Z,11:N))
}
# Generación de medias OC para datos híbridos
Z_oc_funtion = function(Z_k,delta_k,fun_list_k,type_OC = "H") {
  # OC parte vectorial
  set.seed(1997) # Fijamos semilla para localizacion de 1 en vector media OC
  if (type_OC == "H" | type_OC == "V") {
    vec = matrix(Z_k@vecData[1,]*0 + round(runif(ncol(Z_k@vecData)))*delta_k,nrow = 1)
  }else{
    vec = matrix(Z_k@vecData[1,]*0,nrow = 1)
  }
  # OC parte funcional
  if (type_OC == "H" | type_OC == "F") {
    fun_list = list()
    for (i in 1:length(Z_k@funData)) {
      Ei = matrix(fun_list_k[[i]](Z_k@funData[[i]]@argvals[[1]])*delta_k,nrow = 1)
      fun_list[[i]] = funData(Z_k@funData[[i]]@argvals[[1]],Ei)
    }
    fun = multiFunData(fun_list)
  }else{
    fun_list = list()
    for (i in 1:length(Z_k@funData)) {
      Ei = matrix(fun_list_k[[i]](Z_k@funData[[i]]@argvals[[1]])*0,nrow = 1)
      fun_list[[i]] = funData(Z_k@funData[[i]]@argvals[[1]],Ei)
    }
    fun = multiFunData(fun_list)
  }
  # OC híbrido
  Z_shiff = mfh_data(fun_data = fun,vec_data = vec)
  return(Z_shiff)
}
Z_tau_oc = function(Z_k,tau,mu_i){
  n_k = nrow(Z_k@vecData)
  Z_k1 = select_Mfun(Z_k,1:tau)
  Z_k2 = select_Mfun(Z_k,(tau + 1):n_k) + mu_i
  Z_k = append_Mfun(Z_k1,Z_k2)
}
# [Fin del código] ===========================================================.#