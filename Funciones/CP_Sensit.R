# [Calculo de componentes sensibles (Codigo 2)] ==============================.#
# {Calculo de percentiles usando kernel} ----------------------------------.####
quant_dens = function(x, prob){
  res_hdr = get_hdr_1d(x, method = method_kde_1d(kernel = "gaussian",bw = "nrd0"), 
                       probs = prob, n = 2000)
  #res_hdr = get_hdr_1d(x, method = 'kde',probs = prob,n = 2000)
  return(min(res_hdr$df_est$x[cumsum(res_hdr$df_est$fhat_discretized) >= prob]))
}
# {Calculos auxiliares de la estadística T2} ------------------------------.####
ID_sens = function(i,score_dat_j,CPs_Fj,type){
  score_i = matrix(score_dat_j[i,],nrow = 1)
  if (type == 'H') {
    lambda_j = CPs_Fj$PCA$MFPCAH$mfdaH$values
  }
  if (type == 'F') {
    lambda_j = CPs_Fj$PCA$values
  }
  if (type == 'V') {
    M = dim(CPs_Fj$Tm)[2]
    lambda_j = c(CPs_Fj$PCA$eig[,1])[1:M]
  }
  Tm_aux = T2_m(score_i,lambda_j, vect = "cp")
  Tm_new = Tm_aux(i = 1,cp = 1:length(lambda_j))
  Rm_new = Tm_new/apply(CPs_Fj$Tm,2,mean)
  Rm_max = mean(sort(Rm_new/CPs_Fj$limits,decreasing = TRUE)[1:2])
  #Rm_max = max(Rm_new/CPs_Fj$limits)
  if (Rm_max > CPs_Fj$Cl) {
    S = c(1:length(lambda_j))[Rm_new > CPs_Fj$limits]
  }else{
    S = c(1:length(lambda_j))[Rm_new > max(Rm_new)*10]
  }
  result = list(R_max = Rm_max, S = S)
  return(result)
}
S_sens = function(j,list_j,CPs_Fj, score_dat_j, scor_train_j,alpha_j ,type) {
  if (type == 'H') {
    lambda_j = CPs_Fj$PCA$MFPCAH$mfdaH$values
  }
  if (type == 'F') {
    lambda_j = CPs_Fj$PCA$values
  }
  if (type == 'V') {
    J = ncol(CPs_Fj$PCA$ind$coord)
    lambda_j = CPs_Fj$PCA$eig[1:J,1]
  }
  score_j = matrix(score_dat_j[j,],nrow = 1)
  Tm_aux = T2_m(score_j,lambda_j, vect = "i")
  S = list_j$S
  if (length(S) != 0) {
    Tm_new = Tm_aux(i = 1,cp = S)
    delta_i = T_delta(scor_train_j,S,lambda_j,delta = alpha_j)
  }else{
    S = 1:length(lambda_j)
    Tm_new = Tm_aux(i = 1,cp = S)
    delta_i = T_delta(scor_train_j,S,lambda_j,delta = alpha_j)
  }
  return(c(T2 = Tm_new, LS = delta_i))
}
T2_m = function(data_j,lambda_j,vect = NULL) {
  T2_mi = function(data_w,lambda_w,cp_w) {
    scor_i = data_w[cp_w]
    if (length(cp_w) == 1) {
      Lmbda_i = 1/lambda_w[cp_w]
      return(scor_i^2*Lmbda_i)
    }else{
      Lmbda_i = diag(1/lambda_w[cp_w])
      return(t(scor_i) %*% Lmbda_i %*% scor_i) 
    }
  }
  fun_aux = function(i,cp){
    return(T2_mi(data_w = data_j[i,],lambda_w = lambda_j,cp_w = cp))
  }
  if (is.null(vect)) {
    return(T2_mi)
  } else {
    return(Vectorize(fun_aux,vect))
  }
}
# {Calculo de limites kernel para T2} -------------------------------------.####
T_delta = function(scor_data_j,S_j,lambda_j,delta) {
  Tm_aux = T2_m(scor_data_j,lambda_j,vect = "i")
  T2_deff = Tm_aux(1:nrow(scor_data_j),cp = S_j)
  delta_j = quant_dens(T2_deff,prob = delta) # Limite de control estandar
  return(delta_j)
}
# {Calculo de componentes sensibles} --------------------------------------.####
## Para clase hibridos .....................................................####
CP_Sensit_H = function(data_i,M = NULL, delta, omega, wind, alpha_rm){
  # {Antes del monitoreo}
  n = dim(data_i@vecData)[1]
  n_a = round(n*(1 - wind))
  id_a = sample(1:n,n_a)
  id_b = 1:n;id_b = id_b[!(id_b %in% id_a)]
  data_a = select_Mfun(data_i,id_a)
  data_b = select_Mfun(data_i,id_b)
  # 1) Estandarizamos los datos
  scalar_a = data_a
  standar_obj = NULL
  scalar_b = data_b
  # 2) Hacemos PCA con datos de A
  PCA_a = MFHPCA(scalar_a,delta = delta,omega = omega,scale = FALSE,M = M)
  # 3) Calculamos los puntajes de proyecciones sobre B
  scores_data_b = pryc_pca_H(data_j = scalar_b,PCA_j = PCA_a)
  # 4) Calculamos las estadísticas T2 para cada componente principal
  lambda_i = PCA_a$MFPCAH$mfdaH$values
  Tm_aux = T2_m(scores_data_b,lambda_i,vect = "i")
  Tm = list();Rm = list()
  for (i in 1:ncol(scores_data_b)) {
    Tm[[paste0("Tm_",i)]] = Tm_aux(1:nrow(scores_data_b),cp = i)
    Rm[[paste0("Rm_",i)]] = Tm[[paste0("Tm_",i)]]/mean(Tm[[paste0("Tm_",i)]])
  }
  Tm = do.call("cbind", Tm)
  Rm = do.call("cbind", Rm)
  # 5) Calculamos limites de los componentes sensibles quantil
  Clm = apply(Rm,2,function(x){quant_dens(x,prob = alpha_rm)})
  CL_inv = matrix(1/Clm,ncol = dim(Rm)[2],nrow = dim(Rm)[1],byrow = TRUE)
  #Cl = quant_dens(apply(Rm*CL_inv,1,function(x){max(x)}),prob = alpha_rm)
  Cl = quant_dens(apply(Rm*CL_inv,1,function(x){mean(sort(x,decreasing = TRUE)[1:2])}),prob = alpha_rm)
  result_i = list(PCA = PCA_a,standar_obj = standar_obj,Tm = Tm,Rm = Rm,limits = Clm,Cl = Cl)
  return(result_i)
}
## Para clase vectorial ....................................................####
CP_Sensit_V = function(data_i, delta, wind, alpha_rm){
  data_i = as.data.frame(data_i)
  # {Antes del monitoreo}
  n = dim(data_i)[1]
  n_a = round(n*(1 - wind))
  id_a = sample(1:n,n_a)
  id_b = 1:n;id_b = id_b[!(id_b %in% id_a)]
  data_a = data_i[id_a,]
  data_b = data_i[id_b,]
  # 1) Estandarizamos los datos
  scalar_a = data_a
  standar_obj = NULL
  scalar_b = data_b
  # 2) Hacemos PCA con datos de A
  cat('Calculating PCA of the vectorial part...','\n')
  p_vec = ncol(data_a)
  pca_r = PCA(data_a, scale.unit = FALSE, ncp = p_vec, graph = FALSE)
  J = c(1:p_vec)[cumsum(pca_r$eig[,1]/sum(pca_r$eig[,1])) > delta][1]
  pca_r = PCA(data_a, scale.unit = FALSE, ncp = J, graph = FALSE)
  # 3) Calculamos los puntajes de proyecciones sobre B
  scores_data_b = predict.PCA(pca_r,scalar_b)$coord
  # 4) Calculamos las estadísticas T2 para cada componente principal
  lambda_i = pca_r$eig[1:J,1]
  Tm_aux = T2_m(scores_data_b,lambda_i,vect = "i")
  Tm = list();Rm = list()
  for (i in 1:ncol(scores_data_b)) {
    Tm[[paste0("Tm_",i)]] = Tm_aux(1:nrow(scores_data_b),cp = i)
    Rm[[paste0("Rm_",i)]] = Tm[[paste0("Tm_",i)]]/mean(Tm[[paste0("Tm_",i)]])
  }
  Tm = do.call("cbind", Tm)
  Rm = do.call("cbind", Rm)
  # 5) Calculamos limites de los componentes sensibles quantil 0.99
  Clm = apply(Rm,2,function(x){quant_dens(x,prob = alpha_rm)})
  CL_inv = matrix(1/Clm,ncol = dim(Rm)[2],nrow = dim(Rm)[1],byrow = TRUE)
  #Cl = quant_dens(apply(Rm*CL_inv,1,function(x){max(x)}),prob = alpha_rm)
  Cl = quant_dens(apply(Rm*CL_inv,1,function(x){mean(sort(x,decreasing = TRUE)[1:2])}),prob = alpha_rm)
  # 6) Resultados
  result_i = list(PCA = pca_r,standar_obj = standar_obj,Tm = Tm,Rm = Rm,limits = Clm,Cl = Cl)
  return(result_i)
}
## Para clase funcional ....................................................####
CP_Sensit_F = function(data_i, delta, wind, alpha_rm){
  # {Antes del monitoreo}
  n = dim(data_i[[1]]@X)[1]
  n_a = round(n*(1 - wind))
  id_a = sample(1:n,n_a)
  id_b = 1:n;id_b = id_b[!(id_b %in% id_a)]
  data_a =  data_i[id_a]
  data_b = data_i[id_b]
  p_fun = length(data_i)
  # 1) Estandarizamos los datos
  scalar_a = data_a
  standar_obj = NULL
  scalar_b = data_b
  # 2) Hacemos PCA con datos de A
  cat('Calculating MFPCA of the multivariate functional part...','\n')
  uniExpansions_r = list()
  for (i in 1:p_fun) {
    uniExpansions_r[[i]] = list(type = 'uFPCA',pve = delta) # PCA univariados
  }
  mfda_r  = MFPCA(scalar_a, M = 10,uniExpansions = uniExpansions_r)
  aport_mfda_r = summary(mfda_r)
  L = c(1:10)[cumsum(mfda_r$values/sum(mfda_r$values)) > delta][1]
  if (is.na(L)) {
    mfda_r = MFPCA(scalar_a, M = 20,uniExpansions = uniExpansions_r)
    L = c(1:20)[cumsum(mfda_r$values/sum(mfda_r$values)) > delta][1]
  }
  mfda_r = MFPCA(scalar_a ,M = L, uniExpansions = uniExpansions_r)
  PCA_a = mfda_r
  # 3) Calculamos los puntajes de proyecciones sobre B
  scores_data_b = pryc_pca_F(data_j = scalar_b,PCA_j = PCA_a)
  # 4) Calculamos las estadísticas T2 para cada componente principal
  lambda_i = PCA_a$values[1:L]
  Tm_aux = T2_m(scores_data_b,lambda_i,vect = "i")
  Tm = list();Rm = list()
  for (i in 1:ncol(scores_data_b)) {
    Tm[[paste0("Tm_",i)]] = Tm_aux(1:nrow(scores_data_b),cp = i)
    Rm[[paste0("Rm_",i)]] = Tm[[paste0("Tm_",i)]]/mean(Tm[[paste0("Tm_",i)]])
  }
  Tm = do.call("cbind", Tm)
  Rm = do.call("cbind", Rm)
  # 5) Calculamos limites de los componentes sensibles quantil 0.99
  Clm = apply(Rm,2,function(x){quant_dens(x,prob = alpha_rm)})
  CL_inv = matrix(1/Clm,ncol = dim(Rm)[2],nrow = dim(Rm)[1],byrow = TRUE)
  #Cl = quant_dens(apply(Rm*CL_inv,1,function(x){max(x)}),prob = alpha_rm)
  Cl = quant_dens(apply(Rm*CL_inv,1,function(x){mean(sort(x,decreasing = TRUE)[1:2])}),prob = alpha_rm)
  result_i = list(PCA = PCA_a,standar_obj = standar_obj,Tm = Tm,Rm = Rm,limits = Clm,Cl = Cl)
  return(result_i)
}
# [Fin del codigo] ========================================================.####