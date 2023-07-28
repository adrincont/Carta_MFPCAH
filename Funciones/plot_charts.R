# [Grafica para las cartas de control] =======================================.#
# [Grafica para carta de componentes sensibles]----------------------------.####
plot_chart_sns = function(Cntr,Tab_i = NULL, type) {
  alpha_i = Cntr@alpha
  CPs_Fj = Cntr@CPs
  if (type == 'H') {
    lambda_j = CPs_Fj$PCA$MFPCAH$mfdaH$values
  }
  if (type == 'F') {
    lambda_j = CPs_Fj$PCA$values
  }
  if (type == 'V') {
    lambda_j = c(CPs_Fj$PCA$eig[,1])
  }
  scor_train = Cntr@scor_data
  score_dat = scor_train
  N = nrow(score_dat)
  ### 3) Calculamos Tm,Rm y determinamos los componentes sensibles
  Cal_m = sapply(c(1:N),ID_sens,score_dat_j = score_dat, 
                 CPs_Fj = CPs_Fj,type = type,simplify = FALSE)
  R_max = sapply(Cal_m, function(x){x$R_max})
  Args_i = list(CPs_Fj = CPs_Fj,score_dat_j = score_dat, 
                scor_train_j = scor_train, type = type, alpha_j = alpha_i)
  T_calc_i = mapply(S_sens, as.list(1:N), Cal_m, MoreArgs = Args_i,SIMPLIFY = TRUE)
  df_control = as.data.frame(cbind(Rm2 = c(R_max), t(T_calc_i)))
  df_control$Cnt = as.numeric(unlist(df_control$T2) > unlist(df_control$LS))
  df_control$i = 1:nrow(df_control)
  df_control$i = df_control$i - max(df_control$i)
  rownames(df_control) = df_control$i
  plot_i1 = df_control %>% ggplot(aes(i,T2)) + geom_point(color = "gray")  + 
    geom_line(color = "gray") + theme_bw()
  plot_i2 = df_control %>% ggplot(aes(i,Rm2)) + geom_point(color = "gray") +
    geom_line(color = "gray")
  if (!is.null(Tab_i)) {
    plot_i1 = plot_i1 + geom_point(data = Tab_i, aes(i,T2,color = as.factor(Cnt)))  + 
      geom_line(data = Tab_i, aes(i,T2)) + geom_line(data = Tab_i,aes(y = LS),color = "green") +
      geom_vline(xintercept = max(df_control$i),color = 'black') +
      scale_color_manual(values = c("blue","red"), name = "",labels = c("In","Out")) + theme_bw() +
      ylab(expression(paste(T[s]^2)))
    Tab_i$CL = CPs_Fj$Cl
    plot_i2 = plot_i2 + geom_point(data = Tab_i, aes(i,Rm2,color = as.factor(Cnt))) +
      geom_line(data = Tab_i, aes(i,Rm2)) + geom_line(data = Tab_i,aes(y = CL),color = "green") + 
      geom_vline(xintercept = max(df_control$i),color = 'black') +
      scale_color_manual(values = c("blue","red"), name = "", labels = c("In","Out")) + theme_bw() +
      ylab(expression(paste(R[m]^2)))
  }
  return(list(plot_i1,plot_i2))
}
# [Grafica para carta de residuales]---------------------------------------.####
plot_chart_QT = function(Z_i,Cntr,Tab_i = NULL){
  i_lev = 2
  score_dat = Cntr@PCA$MFPCAH$mfdaH$scores
  lambda_i = Cntr@PCA$MFPCAH$mfdaH$values
  M = Cntr@M
  N = dim(Z_i@vecData)[1]
  PCA_a = Cntr@PCA
  Tm_aux = T2_m(score_dat,lambda_i,vect = "i")
  T2_deff = Tm_aux(1:nrow(score_dat),cp = 1:M)
  Tm_aux = T2_m(score_dat,rep(1,length(lambda_i)),vect = "i")
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

  df_rs = data.frame(T2 = T2_deff,Q2 = Q2_deff,Cnt_T2 = 0,Cnt_Q2 = 0,Cnt_K2 = 0)
  # Para la carta T2
  df_rs$Cnt_T2[df_rs$T2 > Cntr@levels$T2[i_lev]] = 1
  # Para la carta Q2
  df_rs$Cnt_Q2[df_rs$Q2 > Cntr@levels$Q2[i_lev]] = 1
  # Carta para el conjunto
  puntos = data.frame(x = df_rs$T2,y = df_rs$Q2)
  ls = contourLines(Cntr@density$Q2_T2, level = Cntr@levels$Q2_T2[i_lev])
  aux = rep(0,length(df_rs$Cnt_K2))
  for (cn_i in 1:length(ls)) {
    aux = aux + point.in.polygon(puntos$x, puntos$y, ls[[cn_i]]$x, ls[[cn_i]]$y)
  }
  df_rs$Cnt_K2 = abs(aux - 1)
  
  ls1 = contourLines(Cntr@density$Q2_T2, level = Cntr@levels$Q2_T2[1])
  ls2 = contourLines(Cntr@density$Q2_T2, level = Cntr@levels$Q2_T2[2])
  T2_df = data.frame(nombres = names(Cntr@levels$T2), valores = Cntr@levels$T2)
  Q2_df = data.frame(nombres = names(Cntr@levels$Q2), valores = Cntr@levels$Q2)
  p = df_rs %>% ggplot(aes(T2,Q2)) + geom_point(aes(shape = as.factor(Cnt_K2)))
  for (i in 1:length(ls1)) {
    p = p + geom_path(data = data.frame(x = ls1[[i]]$x,y = ls1[[i]]$y),aes(x,y,color = "0.95"))
  }
  for (j in 1:length(ls2)) {
    p = p + geom_path(data = data.frame(x = ls2[[j]]$x,y = ls2[[j]]$y),aes(x,y,color = "lev"))
  }
  p = p + 
      geom_vline(data = T2_df,aes(xintercept = valores,color = nombres)) +
      geom_hline(data = Q2_df,aes(yintercept = valores,color = nombres))
  if (!is.null(Tab_i)) {
    p = p + geom_point(data = Tab_i, aes(x = T2,y = Q2,color = as.factor(Cnt_K2)),shape = 2) 
  }
  p = p + scale_color_manual(values = c("blue","black","red","green"), name = "") + theme_bw()
  return(p)
}
# [Fin del código] ===========================================================.#