#.============================================================================.#
#                            SIMULACIONES ARTICULO                             #
#.============================================================================.#
# Autor: Andrey Duvan Rincon Torres    Fecha: 27/07/2023
# [Librerias] ================================================================.#
rm(list =  ls())
source("library.R")
# [Funciones] ================================================================.#
source("Funciones/multiFunDataH.R") # Clase de objetos multiFunDataH
source("Funciones/multiFunDataH_plot.R") # Gráficos para la clase multiFunDataH
source("Funciones/Xi_gen.R") # Generar direcciones híbridas
source("Funciones/r_FunH.R") # Generar observaciones híbridas
source("Funciones/scale_FunDataH.R") # Estandarizar híbridos
source("Funciones/MFHPCA.R") # PCA híbrido
source("Funciones/CP_Sensit.R") # Componentes sensibles
source("Funciones/Cntrl_sns.R") # Carta de control con componentes sensibles
source("Funciones/Cntrl_QT2.R") # Carta de control residual
source("Funciones/plot_charts.R") # Grafico para las cartas de control
RL_fun = function(RL_j,x){length(RL_j[RL_j <= x])/length(RL_j)}
RL_fun = Vectorize(RL_fun,'x')
# Configuración proceso híbrido OC
fun_list_delta_1 = list(f1 = function(x){3*x + x^2},
                        f2 = function(x){3*x^2 + x})
fun_list_delta_2 = list(f1 = Vectorize(function(x) {
  if (x >= 1/4 & x <= 3/4) {return(sin(4*pi*x))} else {return(0)}},'x'),
  f2 = Vectorize(function(x) {
    if (x >= 1/4 & x <= 3/4) {return(cos(4*pi*x))} else {return(0)}},'x')
)
fun_list_delta_3 = list(f1 = function(x){3*exp(-x)},
                        f2 = function(x){sin(4*pi*x)})
# Niveles de cambio en media a considerar
list_delta = list('I' = 0.1, 'II' = 0.3, 'III' = 0.8, 'IV' = 2)
# [ARL0 diferentes distribuciones y tamaños muestra (Tabla 1)] ============.####
Celd = 11;Tab_i = 1;P_i = 2;N_iter = 1050
# Hiperparametros fijos ---------------------------------------------------.####
ARL_0 = 200;delta = 0.99;omega = 0.95;wind = 0.25;n_grill = 100
# Hiperparametros variables -----------------------------------------------.####
n_fun = 5;n_vec = 10;n_data = 1500;dist_i = 'norm'
conf_1 = c("PolyHigh","FourierLin","Wiener")
conf_2 = c("PolyHigh","FourierLin","Wiener","PolyHigh","FourierLin")
conf_3 = c("PolyHigh","FourierLin","Wiener","PolyHigh","FourierLin",
           "PolyHigh","FourierLin","Wiener","PolyHigh","FourierLin")
conf = conf_2 # Se selecciona la configuracion de la parte funcional
# Simulación --------------------------------------------------------------.####
ic = 1 + N_iter*(P_i - 1)
results = list()
results$Data = list();results$Data$F1 = list();results$Data$F2 = list()
results$Char = list();results$ARL_tab = list();ARL_list = list()
for (i in 1:N_iter) {
  # Calculo de la carta
  Z = r_FunH(n_fun = n_fun,n_vec = n_vec,n_grid = n_grill,N = n_data,conf = conf,dist_i = dist_i)
  results$Data$F1[[i]] = Z
  Cntr_F1 = Cntrl_sns_H(Z_i = Z,delta = delta,wind = wind, ARL_0 = ARL_0,omega = omega)
  results$Char[[i]] = Cntr_F1
  Z_new_i = r_FunH(n_fun = n_fun,n_vec = n_vec,n_grid = n_grill,N = 1000,conf = conf,dist_i = dist_i)
  results$Data$F2[[i]] = Z_new_i
  results$ARL[[i]] = Cntr_point(Z_new_i,Cntr_F1,type = 'H')
  ARL_list = append(ARL_list,list(which(results$ARL[[i]]$Cnt == 1)))
  # Objetos de las listas
  if (ic %% 50 == 0) {
    name_i = paste0('Resultados/C',Celd,'_T',Tab_i,'_S_',ic/50,".rds")
    name_i2 = paste0('Resultados/C',Celd,'_T',Tab_i,'_RL_P',P_i,".rds")
    cat("Guardando ",name_i,"....\n")
    saveRDS(results, file = name_i);saveRDS(ARL_list, file = name_i2)
    results$Data = list();results$Data$F1 = list();results$Data$F2 = list()
    results$Char = list();results$ARL_tab = list();results$ARL = list()
  }
  # Resultados por iteracion
  if (ic %% 50 == 0) {
    RL_i = sapply(ARL_list,function(x){x[1]})
    RL_i[is.na(RL_i)] = 1050
    ARL = mean(RL_i); SDRL = sd(RL_i)
    lab_i = c(10,30,50,70,100,150,200,250,300,500,1000)
    FAR = RL_fun(RL_i,lab_i)
    names(FAR) = lab_i
    cat('|Resultados (',i,')------------------------------------------------|\n')
    cat('ARL0 = ',round(ARL,3),'| ','SDRL = ',round(SDRL,3),'|\n')
    cat('FAR \n');print(round(FAR,3));cat('\n')
    cat('|------------------------------------------------------------------|\n') 
  }
  ic = ic + 1
  gc()
}
# [ARL1 diferentes distribuciones y tamaños muestra (Tabla 2)] ============.####
Celd = 11;Tab_i = 1;P_i = 2;N_iter = 1050
# Hiperparametros fijos ---------------------------------------------------.####
ARL_0 = 200;delta = 0.99;omega = 0.95;wind = 0.25
n_grill = 100; conf = c("FourierLin","Wiener");n_fun = 2;n_vec = 5
fun_list_delta = fun_list_delta_3
# Hiperparametros variables -----------------------------------------------.####
n_data = 1500;dist_i = 'norm'
# Simulación --------------------------------------------------------------.####
ic = 1 + N_iter*(P_i - 1)
results = list()
results$Data = list();results$Data$F1 = list();results$Data$F2 = list()
results$Char = list();results$ARL_tab = list();ARL_list = list()
for (name_i in names(list_delta)) {
  results$Data$F2[[name_i]] = list()
  results$ARL[[name_i]] = list()
  ARL_list[[name_i]] = list()
}
for (i in 1:N_iter) {
  # Calculo de la carta
  Z = r_FunH(n_fun = n_fun,n_vec = n_vec,n_grid = n_grill,N = n_data,conf = conf,dist_i = dist_i)
  results$Data$F1[[i]] = Z
  Cntr_F1 = Cntrl_sns_H(Z_i = Z,delta = delta,wind = wind, ARL_0 = ARL_0,omega = omega)
  results$Char[[i]] = Cntr_F1
  Z_new_i = r_FunH(n_fun = n_fun,n_vec = n_vec,n_grid = n_grill,N = 400,conf = conf,dist_i = dist_i)
  # fuera de control
  for (name_i in names(list_delta)) {
    Z_oc_i = Z_oc_funtion(Z_new_i,delta_k = list_delta[[name_i]],fun_list_delta)
    Z_oc = Z_tau_oc(Z_new_i,tau = 25,mu_i = Z_oc_i)
    results$Data$F2[[name_i]][[i]] = Z_oc
    results$ARL[[name_i]][[i]] = Cntr_point(Z_oc,Cntr_F1,type = 'H')
    ARL_list[[name_i]] = append(ARL_list[[name_i]],list(which(results$ARL[[name_i]][[i]]$Cnt == 1)))
    if (ic %% 50 == 0) {
      name_j = paste0('Resultados/C',name_i,'_',Celd,'_T',Tab_i,'_S_',ic/50,".rds")
      name_j2 = paste0('Resultados/C',name_i,'_',Celd,'_T',Tab_i,'_RL_P',P_i,".rds")
      cat("Guardando ",name_j,"....\n")
      saveRDS(results, file = name_j);saveRDS(ARL_list[[name_i]], file = name_j2)
      results$Data = list();results$Data$F1 = list();results$Data$F2 = list()
      results$Char = list();results$ARL_tab = list();results$ARL = list()
      for (w_i in names(list_delta)) {
        results$Data$F2[[w_i]] = list()
        results$ARL[[w_i]] = list()
      }
    }
    if (ic %% 50 == 0) {
      RL_i = sapply(ARL_list[[name_i]],function(x){x[1]})
      RL_i[is.na(RL_i)] = 450
      ARL = mean(RL_i); SDRL = sd(RL_i)
      lab_i = c(10,30,50,70,100,150,200,250,300,500,1000)
      FAR = RL_fun(RL_i,lab_i)
      names(FAR) = lab_i
      cat('|Resultados (',i,'_',name_i,')-------------------------------------|\n')
      cat('ARL0 = ',round(ARL,3),'| ','SDRL = ',round(SDRL,3),'|\n')
      cat('FAR \n');print(round(FAR,3));cat('\n')
      cat('|------------------------------------------------------------------|\n') 
    }
  }
  ic = ic + 1
  gc()
}
# [Gráficos de los procesos fuera de control] ..............................####
Z_new_i = r_FunH(n_fun = 2,n_vec = 5,n_grid = 100,
                 N = 500,conf =  c("FourierLin","FourierLin"),dist_i = 'norm')
# Cambio 1
Z_oc_i = Z_oc_funtion(Z_new_i,delta_k = 0.8,fun_list_delta_1)
Z_oc = Z_tau_oc(Z_new_i,tau = 1,mu_i = Z_oc_i)
par(mfrow = c(3,3), mar = c(4,5,2,2))
fun_dat_h = Z_new_i@funData[[1]]
roahd::fbplot(Data = roahd::fData(grid = fun_dat_h@argvals[[1]],
                                  values = fun_dat_h@X),
              Fvalue = 1.5,ylab = 'IC (I)') 
fun_dat_h = Z_oc@funData[[1]]
roahd::fbplot(Data = roahd::fData(grid = fun_dat_h@argvals[[1]],
                                  values = fun_dat_h@X),Fvalue = 1.5,
              ylab = expression(paste("OC I (",3*x + x^2,")"))) 
fun_dat_h = Z_oc@funData[[2]]
roahd::fbplot(Data = roahd::fData(grid = fun_dat_h@argvals[[1]],
                                  values = fun_dat_h@X),Fvalue = 1.5,
              ylab = expression(paste("OC (",3*x^2 + x,")")))
# Cambio 2
Z_oc_i = Z_oc_funtion(Z_new_i,delta_k = 0.8,fun_list_delta_2)
Z_oc = Z_tau_oc(Z_new_i,tau = 1,mu_i = Z_oc_i)
fun_dat_h = Z_new_i@funData[[1]]
roahd::fbplot(Data = roahd::fData(grid = fun_dat_h@argvals[[1]],
                                  values = fun_dat_h@X),
              Fvalue = 1.5,ylab = 'IC (II)') 
fun_dat_h = Z_oc@funData[[1]]
roahd::fbplot(Data = roahd::fData(grid = fun_dat_h@argvals[[1]],
                                  values = fun_dat_h@X),Fvalue = 1.5,
              ylab = expression(paste("OC (",sin(4*pi*x)," [1/4,3/4])"))) 
fun_dat_h = Z_oc@funData[[2]]
roahd::fbplot(Data = roahd::fData(grid = fun_dat_h@argvals[[1]],
                                  values = fun_dat_h@X),Fvalue = 1.5,
              ylab = expression(paste("OC (",cos(4*pi*x)," [1/4,3/4])"))) 
# Cambio 3
Z_oc_i = Z_oc_funtion(Z_new_i,delta_k = 0.8,fun_list_delta_3)
Z_oc = Z_tau_oc(Z_new_i,tau = 1,mu_i = Z_oc_i)
fun_dat_h = Z_new_i@funData[[1]]
roahd::fbplot(Data = roahd::fData(grid = fun_dat_h@argvals[[1]],
                                  values = fun_dat_h@X),
              Fvalue = 1.5,ylab = 'IC (III)') 
fun_dat_h = Z_oc@funData[[1]]
roahd::fbplot(Data = roahd::fData(grid = fun_dat_h@argvals[[1]],
                                  values = fun_dat_h@X),Fvalue = 1.5,
              ylab = expression(paste("OC I (",exp(-x),")"))) 
fun_dat_h = Z_oc@funData[[2]]
roahd::fbplot(Data = roahd::fData(grid = fun_dat_h@argvals[[1]],
                                  values = fun_dat_h@X),Fvalue = 1.5,
              ylab = expression(paste("OC II (",sin(4*pi*x),")")))
# [Gráficos de cartas fuera de control] ....................................####
Z = r_FunH(n_fun = 2,n_vec = 5,n_grid = 100,
           N = 1000,conf =  c("FourierLin","FourierLin"),dist_i = 'norm')
Z_new_i = r_FunH(n_fun = 2,n_vec = 5,n_grid = 100,
                 N = 200,conf =  c("FourierLin","FourierLin"),dist_i = 'norm')
Cntr_F1 = Cntrl_sns_H(Z_i = Z,delta = 0.99,wind = 0.3, ARL_0 = 200,omega = 0.90)
Z_oc_i = Z_oc_funtion(Z_new_i,delta_k = 1.5,fun_list_delta_1)
Z_oc = Z_tau_oc(Z_new_i,tau = 50,mu_i = Z_oc_i)
Tab_i = Cntr_point(Z_oc,Cntr_F1,type = 'H')
plot_chart = plot_chart_sns(Cntr_F1,Tab_i,type = 'H')
p1 = plot_chart[[1]] + xlim(-100,100)
p2 = plot_chart[[2]] + xlim(-100,100) 
p1/p2
# [Simulación para ARL1 en propuesta P1 (propuesta híbrida)] ==============.####
Celd = 1;Tab_i = 3;P_i = 1;N_iter = 550
# Hiperparametros fijos ---------------------------------------------------.####
ARL_0 = 200;delta = 0.99;dist_i = 'norm'
n_grill = 100; conf = c("FourierLin","Wiener");n_fun = 2;n_vec = 5
n_data = 1500;omega = 0.95;wind = 0.25
# Hiperparametros variables -----------------------------------------------.####
fun_list_delta = fun_list_delta_3;type_OC = "H"
# Simulación --------------------------------------------------------------.####
ic = 1 + N_iter*(P_i - 1)
results = list()
results$Data = list();results$Data$F1 = list();results$Data$F2 = list()
results$Char = list();results$ARL_tab = list();ARL_list = list()
for (name_i in names(list_delta)) {
  results$Data$F2[[name_i]] = list()
  results$ARL[[name_i]] = list()
  ARL_list[[name_i]] = list()
}
for (i in 1:N_iter) {
  # Calculo de la carta
  Z = r_FunH(n_fun = n_fun,n_vec = n_vec,n_grid = n_grill,N = n_data,conf = conf,dist_i = dist_i)
  results$Data$F1[[i]] = Z
  Cntr_F1 = Cntrl_sns_H(Z_i = Z,delta = delta,wind = wind, ARL_0 = ARL_0,omega = omega)
  results$Char[[i]] = Cntr_F1
  Z_new_i = r_FunH(n_fun = n_fun,n_vec = n_vec,n_grid = n_grill,N = 300,conf = conf,dist_i = dist_i)
  # fuera de control
  for (name_i in names(list_delta)) {
    Z_oc_i = Z_oc_funtion(Z_new_i,delta_k = list_delta[[name_i]],fun_list_delta,type_OC = type_OC)
    Z_oc = Z_tau_oc(Z_new_i,tau = 25,mu_i = Z_oc_i)
    results$Data$F2[[name_i]][[i]] = Z_oc
    results$ARL[[name_i]][[i]] = Cntr_point(Z_oc,Cntr_F1,type = 'H')
    ARL_list[[name_i]] = append(ARL_list[[name_i]],list(which(results$ARL[[name_i]][[i]]$Cnt == 1)))
    if (ic %% 50 == 0) {
      name_j = paste0('Resultados/C',name_i,'_',Celd,'_T',Tab_i,'_S_',ic/50,".rds")
      name_j2 = paste0('Resultados/C',name_i,'_',Celd,'_T',Tab_i,'_RL_P',P_i,".rds")
      cat("Guardando ",name_j,"....\n")
      saveRDS(results, file = name_j);saveRDS(ARL_list[[name_i]], file = name_j2)
      results$Data = list();results$Data$F1 = list();results$Data$F2 = list()
      results$Char = list();results$ARL_tab = list();results$ARL = list()
      for (w_i in names(list_delta)) {
        results$Data$F2[[w_i]] = list()
        results$ARL[[w_i]] = list()
      }
    }
    if (ic %% 50 == 0) {
      RL_i = sapply(ARL_list[[name_i]],function(x){x[1]})
      RL_i[is.na(RL_i)] = 320
      ARL = mean(RL_i); SDRL = sd(RL_i)
      lab_i = c(10,30,50,70,100,150,200,250,300,500,1000)
      FAR = RL_fun(RL_i,lab_i)
      names(FAR) = lab_i
      cat('|Resultados (',i,'_',name_i,')-------------------------------------|\n')
      cat('ARL0 = ',round(ARL,3),'| ','SDRL = ',round(SDRL,3),'|\n')
      cat('FAR \n');print(round(FAR,3));cat('\n')
      cat('|------------------------------------------------------------------|\n') 
    }
  }
  ic = ic + 1
  gc()
}
# [Simulación para ARL1 en propuesta P1 (propuesta individual)] ===========.####
Celd = 1;Tab_i = 3;P_i = 1;N_iter = 550
# Hiperparametros fijos ---------------------------------------------------.####
ARL_0 = 200;delta = 0.99;dist_i = 'norm'
n_grill = 100; conf = c("FourierLin","Wiener");n_fun = 2;n_vec = 5
n_data = 1500;omega = 0.95;wind = 0.25
# Hiperparametros variables -----------------------------------------------.####
fun_list_delta = fun_list_delta_3;type_OC = "H"
# Simulación --------------------------------------------------------------.####
ic = 1 + N_iter*(P_i - 1)
results = list()
results$Data = list();results$Data$F1 = list();results$Data$F2 = list()
results$Char = list();results$ARL_tab = list();ARL_list = list()
for (name_i in names(list_delta)) {
  results$Data$F2[[name_i]] = list()
  results$ARL[[name_i]] = list()
  ARL_list[[name_i]] = list()
}
for (i in 1:N_iter) {
  # Calculo de la carta
  Z = r_FunH(n_fun = n_fun,n_vec = n_vec,n_grid = n_grill,N = n_data,conf = conf,dist_i = dist_i)
  results$Data$F1[[i]] = Z
  Cntr_F1 = Cntrl_sns_HS(Z_i = Z,delta = delta,wind = wind,ARL_0 = ARL_0)
  results$Char[[i]] = Cntr_F1
  Z_new_i = r_FunH(n_fun = n_fun,n_vec = n_vec,n_grid = n_grill,N = 300,conf = conf,dist_i = dist_i)
  # fuera de control
  for (name_i in names(list_delta)) {
    Z_oc_i = Z_oc_funtion(Z_new_i,delta_k = list_delta[[name_i]],fun_list_delta,type_OC = type_OC)
    Z_oc = Z_tau_oc(Z_new_i,tau = 25,mu_i = Z_oc_i)
    results$Data$F2[[name_i]][[i]] = Z_oc
    results$ARL[[name_i]][[i]] = Cntr_point_HS(Z_oc,Cntr_F1)
    ARL_list[[name_i]] = append(ARL_list[[name_i]],list(which(results$ARL[[name_i]][[i]]$Cnt != 0)))
    if (ic %% 50 == 0) {
      name_j = paste0('Resultados/C',name_i,'_',Celd,'_T',Tab_i,'_S_',ic/50,".rds")
      name_j2 = paste0('Resultados/C',name_i,'_',Celd,'_T',Tab_i,'_RL_P',P_i,".rds")
      cat("Guardando ",name_j,"....\n")
      saveRDS(results, file = name_j);saveRDS(ARL_list[[name_i]], file = name_j2)
      results$Data = list();results$Data$F1 = list();results$Data$F2 = list()
      results$Char = list();results$ARL_tab = list();results$ARL = list()
      for (w_i in names(list_delta)) {
        results$Data$F2[[w_i]] = list()
        results$ARL[[w_i]] = list()
      }
    }
    if (ic %% 50 == 0) {
      RL_i = sapply(ARL_list[[name_i]],function(x){x[1]})
      RL_i[is.na(RL_i)] = 320
      ARL = mean(RL_i); SDRL = sd(RL_i)
      lab_i = c(10,30,50,70,100,150,200,250,300,500,1000)
      FAR = RL_fun(RL_i,lab_i)
      names(FAR) = lab_i
      cat('|Resultados (',i,'_',name_i,')-------------------------------------|\n')
      cat('ARL0 = ',round(ARL,3),'| ','SDRL = ',round(SDRL,3),'|\n')
      cat('FAR \n');print(round(FAR,3));cat('\n')
      cat('|------------------------------------------------------------------|\n') 
    }
  }
  ic = ic + 1
  gc()
}
# [Tabla 4] ===============================================================.####
Celd = 1;Tab_i = 4;P_i = 1;N_iter = 550
# Hiperparametros fijos ---------------------------------------------------.####
ARL_0 = 200;delta = 0.99;dist_i = 'norm'
n_grill = 100; conf = c("FourierLin","Wiener");n_fun = 2;n_vec = 5
fun_list_delta = fun_list_delta_3;n_data = 1500
# Hiperparametros variables -----------------------------------------------.####
omega = 0.95;wind = 0.25
# Simulación --------------------------------------------------------------.####
ic = 1 + N_iter*(P_i - 1)
results = list()
results$Data = list();results$Data$F1 = list();results$Data$F2 = list()
results$Char = list();results$ARL_tab = list();ARL_list = list()
for (name_i in names(list_delta)) {
  results$Data$F2[[name_i]] = list()
  results$ARL[[name_i]] = list()
  ARL_list[[name_i]] = list()
}
for (i in 1:N_iter) {
  # Calculo de la carta
  Z = r_FunH(n_fun = n_fun,n_vec = n_vec,n_grid = n_grill,N = n_data,conf = conf,dist_i = dist_i)
  results$Data$F1[[i]] = Z
  Cntr_F1 = Cntrl_sns_H(Z_i = Z,delta = delta,wind = wind, ARL_0 = ARL_0,omega = omega)
  results$Char[[i]] = Cntr_F1
  Z_new_i = r_FunH(n_fun = n_fun,n_vec = n_vec,n_grid = n_grill,N = 300,conf = conf,dist_i = dist_i)
  # fuera de control
  for (name_i in names(list_delta)) {
    Z_oc_i = Z_oc_funtion(Z_new_i,delta_k = list_delta[[name_i]],fun_list_delta)
    Z_oc = Z_tau_oc(Z_new_i,tau = 25,mu_i = Z_oc_i)
    results$Data$F2[[name_i]][[i]] = Z_oc
    results$ARL[[name_i]][[i]] = Cntr_point(Z_oc,Cntr_F1,type = 'H')
    ARL_list[[name_i]] = append(ARL_list[[name_i]],list(which(results$ARL[[name_i]][[i]]$Cnt == 1)))
    if (ic %% 50 == 0) {
      name_j = paste0('Resultados/C',name_i,'_',Celd,'_T',Tab_i,'_S_',ic/50,".rds")
      name_j2 = paste0('Resultados/C',name_i,'_',Celd,'_T',Tab_i,'_RL_P',P_i,".rds")
      cat("Guardando ",name_j,"....\n")
      saveRDS(results, file = name_j);saveRDS(ARL_list[[name_i]], file = name_j2)
      results$Data = list();results$Data$F1 = list();results$Data$F2 = list()
      results$Char = list();results$ARL_tab = list();results$ARL = list()
      for (w_i in names(list_delta)) {
        results$Data$F2[[w_i]] = list()
        results$ARL[[w_i]] = list()
      }
    }
    if (ic %% 50 == 0) {
      RL_i = sapply(ARL_list[[name_i]],function(x){x[1]})
      RL_i[is.na(RL_i)] = 320
      ARL = mean(RL_i); SDRL = sd(RL_i)
      lab_i = c(10,30,50,70,100,150,200,250,300,500,1000)
      FAR = RL_fun(RL_i,lab_i)
      names(FAR) = lab_i
      cat('|Resultados (',i,'_',name_i,')-------------------------------------|\n')
      cat('ARL0 = ',round(ARL,3),'| ','SDRL = ',round(SDRL,3),'|\n')
      cat('FAR \n');print(round(FAR,3));cat('\n')
      cat('|------------------------------------------------------------------|\n') 
    }
  }
  ic = ic + 1
  gc()
}
# [Fin del código]=========================================================.####