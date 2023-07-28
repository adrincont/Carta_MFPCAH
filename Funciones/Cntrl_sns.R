# [Calculo de componentes sensibles (Codigo 1)] ==============================.#
# Definir clase Cntrl_senns -----------------------------------------------.####
check_Cntrl_senns = function(object){
  errors = character()
  if (length(errors) == 0) TRUE else errors
}
setClass("Cntrl_senns", representation(CPs = "list", delta_st = "numeric", 
                                       scor_data = "matrix",alpha = "numeric"),
         validity = check_Cntrl_senns)
# Funcion para el calculo de parámetros en Fase I -------------------------.####
## Para clase híbridos .....................................................####
Cntrl_sns_H = function(Z_i, M = NULL, delta, omega, wind, ARL_0){
  alpha_i = 1 - 1/(ARL_0 + 1)
  ### 1) Calculo de componentes principales sensibles 
  CPs_F1_j = CP_Sensit_H(Z_i,M = M,delta = delta,omega = omega,wind = wind,alpha_rm = alpha_i)
  ### 2) Rutina para el calculo de los limites de acuerdo a las componentes
  lambda_i = CPs_F1_j$PCA$MFPCAH$mfdaH$values
  scor_data = pryc_pca_H(data_j = Z_i,PCA_j = CPs_F1_j$PCA)
  delta_deff = T_delta(scor_data,1:length(lambda_i),lambda_i,delta = alpha_i)
  result = new("Cntrl_senns", CPs = CPs_F1_j,delta_st = delta_deff,scor_data = scor_data,
               alpha = alpha_i)
  return(result)
}
## Para clase vectorial ....................................................####
Cntrl_sns_V = function(Z_i, delta, wind, ARL_0){
  alpha_i = 1 - 1/(ARL_0 + 1)
  Z_i = as.data.frame(Z_i)
  ### 1) Calculo de componentes principales sensibles 
  CPs_F1_j = CP_Sensit_V(Z_i,delta = delta,wind = wind,alpha_rm = alpha_i)
  J = ncol(CPs_F1_j$PCA$ind$coord)
  ### 2) Rutina para el calculo de los limites de acuerdo a las componentes
  lambda_i = CPs_F1_j$PCA$eig[1:J,1]
  scor_data = predict.PCA(CPs_F1_j$PCA,Z_i)$coord
  delta_deff = T_delta(scor_data,1:length(lambda_i),lambda_i,delta = alpha_i)
  result = new("Cntrl_senns", CPs = CPs_F1_j,delta_st = delta_deff,scor_data = scor_data,
               alpha = alpha_i)
  return(result)
}
## Para clase funcional ....................................................####
Cntrl_sns_F = function(Z_i, delta, wind, ARL_0){
  alpha_i = 1 - 1/(ARL_0 + 1)
  ### 1) Calculo de componentes principales sensibles 
  CPs_F1_j = CP_Sensit_F(Z_i,delta = delta,wind = wind,alpha_rm = alpha_i)
  L = dim(CPs_F1_j$PCA$scores)[2]
  ### 2) Rutina para el calculo de los limites de acuerdo a las componentes
  lambda_i = CPs_F1_j$PCA$values[1:L]
  scor_data = pryc_pca_F(data_j = Z_i,PCA_j = CPs_F1_j$PCA)
  delta_deff = T_delta(scor_data,1:length(lambda_i),lambda_i,delta = alpha_i)
  result = new("Cntrl_senns", CPs = CPs_F1_j,delta_st = delta_deff,scor_data = scor_data,
               alpha = alpha_i)
  return(result)
}
# Resumen de l salida -----------------------------------------------------.####
print_Cntrl = function(CPs){
  cat('Summarize for PCA _______________________________________________\n')
  PCA_i = summary(CPs@CPs$PCA)
  cat('_________________________________________________________________\n')
  cat('Limits for main components ______________________________________\n')
  print(CPs@CPs$limits)
  cat('_________________________________________________________________\n')
}
setMethod("show",signature = "Cntrl_senns",
          function(object){print_Cntrl(object)})
# Grafico de limites sensibles --------------------------------------------.####
autoplot.Cntrl_senns = function(object){
  Rm = object@CPs$Rm
  Tab_lim = as_tibble(t(object@CPs$limits)) %>% gather()
  plot_i = as_tibble(Rm) %>% gather() %>% ggplot(aes(value, fill = key, colour = key)) +
    geom_density(alpha = 0.1) + 
    geom_vline(data = Tab_lim,aes(xintercept = value, colour = key),
               linetype = "dashed")
  return(plot_i)
}
# Seguimiento en linea para nuevos datos ----------------------------------.####
## Para clase híbrida ......................................................####
Cntr_point = function(Z_i,Cntr_j, type){
  scor_train = Cntr_j@scor_data
  if (type == 'H') {
    alpha_i = Cntr_j@alpha
    CPs_F1_j = Cntr_j@CPs
    lambda_i = Cntr_j@CPs$PCA$MFPCAH$mfdaH$values
    standar_dat = Z_i
    score_dat = pryc_pca_H(data_j = standar_dat,PCA_j = CPs_F1_j$PCA)
  }else if (type == 'F') {
    alpha_i = Cntr_j@alpha
    CPs_F1_j = Cntr_j@CPs
    L = dim(CPs_F1_j$PCA$scores)[2]
    lambda_i = CPs_F1_j$PCA$values[1:L]
    standar_dat = Z_i
    score_dat = pryc_pca_F(data_j = standar_dat,PCA_j = CPs_F1_j$PCA)
  }else if (type == 'V') {
    alpha_i = Cntr_j@alpha
    CPs_F1_j = Cntr_j@CPs
    J = ncol(CPs_F1_j$PCA$ind$coord)
    lambda_i = CPs_F1_j$PCA$eig[1:J,1]
    standar_dat = Z_i
    score_dat = predict.PCA(CPs_F1_j$PCA,standar_dat)$coord
  } else {stop('LETRA NO VALIDA')}
  N = nrow(score_dat)
  ### 3) Calculamos Tm,Rm y determinamos los componentes sensibles
  Cal_m = sapply(c(1:N),ID_sens,score_dat_j = score_dat, 
                 CPs_Fj = CPs_F1_j,type = type,simplify = FALSE)
  R_max = sapply(Cal_m, function(x){x$R_max})
  Args_i = list(CPs_Fj = CPs_F1_j,score_dat_j = score_dat, 
                scor_train_j = scor_train, type = type, alpha_j = alpha_i)
  T_calc_i = mapply(S_sens, as.list(1:N), Cal_m, MoreArgs = Args_i,SIMPLIFY = TRUE)
  df_control = as.data.frame(cbind(Rm2 = c(R_max), t(T_calc_i)))
  df_control$Cnt = as.numeric(unlist(df_control$T2) > unlist(df_control$LS))
  df_control$i = 1:nrow(df_control)
  return(as_tibble(df_control))
}
# Carta de control híbridos por separados ---------------------------------.####
Cntrl_sns_HS = function(Z_i, delta, wind, ARL_0){
  ARL_0 = ARL_0*2
  Fun_i = Z_i@funData
  Vect_i = as.data.frame(Z_i@vecData)
  Cntr_F1_F = Cntrl_sns_F(Z_i = Fun_i,delta = delta,wind = wind,ARL_0 = ARL_0)
  Cntr_F1_V = Cntrl_sns_V(Z_i = Vect_i,delta = delta,wind = wind,ARL_0 = ARL_0)
  res_i = list('Fun' = Cntr_F1_F, 'Vect' = Cntr_F1_V)
  return(res_i)
}
Cntr_point_HS = function(Z_i,Cntr_j){
  Fun_i_New = Z_i@funData
  Vect_New = as.data.frame(Z_i@vecData)
  Cntr_F1_F = Cntr_j$Fun
  Cntr_F1_V = Cntr_j$Vect
  Tab_C1_F = Cntr_point(Fun_i_New,Cntr_F1_F,type = 'F')
  Tab_C1_V = Cntr_point(Vect_New,Cntr_F1_V,type = 'V')
  colnames(Tab_C1_F) = paste0(colnames(Tab_C1_F),'_F')
  colnames(Tab_C1_V) = paste0(colnames(Tab_C1_V),'_V')
  Tab_C1 = cbind(Tab_C1_F,Tab_C1_V)
  Tab_C1$Cnt = Tab_C1$Cnt_F + Tab_C1$Cnt_V
  return(Tab_C1)
}
# [Fin del codigo] ===========================================================.#