# [Carta de control con residual] ============================================.#
# Definir clase Cntrl_QT2 -------------------------------------------------.####
check_Cntrl_QT2 = function(object){
  errors = character()
  if (length(errors) == 0) TRUE else errors
}
setClass("Cntrl_QT2", representation(density = "list", levels = "list", 
                                     standar_obj = "list", PCA = "list",
                                     M = "numeric"),
         validity = check_Cntrl_QT2)
# Funcion para el calculo de parametros en fase 1 -------------------------.####
Cntrl_QT2 = function(Z_i, delta, omega, ARL_0){
  alpha_i = 1 - 1/(ARL_0 + 1)
  ### 1) Calculo de componentes principales
  N = dim(Z_i@vecData)[1]
  PCA_a = MFHPCA(Z_i,delta = delta,omega = omega,scale = FALSE)
  M = length(PCA_a$MFPCAH$mfdaH$values)
  ### 3) Calculo de las estadisticas T2
  lambda_i = PCA_a$MFPCAH$mfdaH$values[1:M]
  scor_data = PCA_a$MFPCAH$mfdaH$scores
  Tm_aux_1 = T2_m(scor_data,lambda_i,vect = "i")
  T2_deff = Tm_aux_1(1:nrow(scor_data),cp = 1:M)
  ### 4) Calculo de la estadistica Q
  M_dir = 1:M
  dir_i = select_Mfun(PCA_a$MFPCAH$dir,M_dir)
  if (length(M_dir) == 1) {
    scor_i = matrix(PCA_a$MFPCAH$mfdaH$scores[,M_dir],ncol = 1)
  } else{
    scor_i = PCA_a$MFPCAH$mfdaH$scores[,M_dir] 
  }
  recost_i = Inv_PCA_H(dir_j = dir_i,scor_j = scor_i)
  residuals_i = Z_i + (-1*recost_i)
  Q2_deff = c(sapply(list(1:N),function(i,res_i){
    return(scalarProduct_mfh(select_Mfun(res_i,i),select_Mfun(res_i,i)))
    },res_i = residuals_i))
  # Densidad para T2_deff
  dens_T2_deff = density(T2_deff)
  lev_1 = quant_dens(T2_deff,prob = 0.99)
  lev_2 = quant_dens(T2_deff,prob = alpha_i)
  levels_T2 = c("0.95" = lev_1,"lev" = lev_2)
  # Densidad para Q2_deff
  dens_Q2_deff = density(Q2_deff)
  lev_1 = quant_dens(Q2_deff,prob = 0.99)
  lev_2 = quant_dens(Q2_deff,prob = alpha_i)
  levels_Q2 = c("0.95" = lev_1,"lev" = lev_2)
  # Densidad conjunta para T2_deff,Q2_deff
  dens2d = MASS::kde2d(T2_deff,Q2_deff,n = 5000,lims = c(-5,max(T2_deff)*(1 + 0.5),-5,max(Q2_deff)*(1 + 0.5)))
  res_hdr = get_hdr(data.frame(x = T2_deff,y = Q2_deff), method = "kde",probs = c(0.99,alpha_i),hdr_membership  = FALSE)
  levels_K2 = c("0.95" = res_hdr$breaks[2],"lev" = res_hdr$breaks[1])
  # Configurar resultados
  density_rs = list(T2 = dens_T2_deff,Q2 = dens_Q2_deff,Q2_T2 = dens2d)
  levels_rs = list(T2 = levels_T2,Q2 = levels_Q2,Q2_T2 = levels_K2)
  PCA_rs = PCA_a
  class(PCA_rs) = "list"
  standar_obj = list(NULL)
  result = new("Cntrl_QT2", density = density_rs, levels = levels_rs, standar_obj = standar_obj,PCA = PCA_rs,M = M) 
  return(result)
}
# Metodos para la clase ---------------------------------------------------.####
print_Cntrl_k2 = function(object){
  cat('_________________________________________________________________\n')
  cat('Limits for T2 ___________________________________________________\n')
  print(object@levels$T2)
  cat('Limits for Q2 ___________________________________________________\n')
  print(object@levels$Q2)
  cat('Limits for K2 ___________________________________________________\n')
  print(object@levels$Q2_T2)
  cat('_________________________________________________________________\n')
}
setMethod("show",signature = "Cntrl_QT2",function(object){print_Cntrl_k2(object)})
# Seguimiento en linea para nuevos datos ----------------------------------.####
CntrK2_point = function(Z_i,object){
  PCA_a = object@PCA
  i_lev = 2
  standar_obj = NULL
  w = NULL
  lambda_i = PCA_a$MFPCAH$mfdaH$values
  M = object@M
  N = nrow(Z_i@vecData)
  ### 1) estandarizamos y proyectamos con la información de fase 1
  standar_dat = Z_i
  score_dat = pryc_pca_H(data_j = standar_dat,PCA_j = PCA_a)
  ### 2) Calculamos T2
  Tm_aux = T2_m(score_dat,lambda_i,vect = "i")
  T2_deff = Tm_aux(1:nrow(score_dat),cp = 1:M)
  ### 2) Calculamos Q2
  n_dir = length(lambda_i)
  M_dir = 1:M
  dir_i = select_Mfun(PCA_a$MFPCAH$dir,M_dir)
  if (length(M_dir) == 1) {
    scor_i = matrix(score_dat[,M_dir],ncol = 1)
  } else{
    scor_i = score_dat[,M_dir] 
  }
  recost_i = Inv_PCA_H(dir_j = dir_i,scor_j = scor_i)
  residuals_i = Z_i + (-1*recost_i)
  Q2_deff = c(sapply(list(1:N),function(i,res_i){
    return(scalarProduct_mfh(select_Mfun(res_i,i),select_Mfun(res_i,i)))
  },res_i = residuals_i))
  df_rs = data.frame(T2 = T2_deff,Q2 = Q2_deff,Cnt_T2 = 0,Cnt_Q2 = 0,Cnt_K2 = 0)
  # Para la carta T2
  df_rs$Cnt_T2[df_rs$T2 > object@levels$T2[i_lev]] = 1
  # Para la carta Q2
  df_rs$Cnt_Q2[df_rs$Q2 > object@levels$Q2[i_lev]] = 1
  # Carta para el conjunto
  puntos = data.frame(x = df_rs$T2,y = df_rs$Q2)
  ls = contourLines(object@density$Q2_T2, level = object@levels$Q2_T2[i_lev])
  aux = rep(0,length(df_rs$Cnt_K2))
  for (cn_i in 1:length(ls)) {
    aux = aux + point.in.polygon(puntos$x, puntos$y, ls[[cn_i]]$x, ls[[cn_i]]$y)
  }
  df_rs$Cnt_K2 = abs(aux - 1)
  return(df_rs)
}
# [Fin del codigo] ===========================================================.#