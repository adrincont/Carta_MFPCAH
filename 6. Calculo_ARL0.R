#.============================================================================.#
#                      CARTA DE CONTROL 1 (CPS hibridos)                       #
#.============================================================================.#
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
# [Simulación de proceso] =================================================.####
N_iter = 1000;ARL_0 = 200;delta = 0.99;omega = 0.90;wind = 0.25
n_fun = 3;n_vec = 5;n_grill = 100;n_data = 800
conf = sample(c("PolyHigh","FourierLin","Wiener"),size = 10,replace = TRUE)
# [Calculo de los ARL C1] ===================================================.####
results = list()
results$Data = list();results$Data$F1 = list();results$Data$F2 = list()
results$Char = list();results$ARL_tab = list();ARL_list = list()
ic = 0
for (i in 1:N_iter) {
  # Calculo de la carta
  Z = r_FunH(n_fun = n_fun,n_vec = n_vec,n_grid = n_grill,N = n_data,conf = conf)
  results$Data$F1[[i]] = Z
  Cntr_F1 = Cntrl_sns_H(Z_i = Z,delta = delta,wind = wind, ARL_0 = ARL_0,omega = omega)
  results$Char[[i]] = Cntr_F1
  Z_new_i = r_FunH(n_fun = n_fun,n_vec = n_vec,n_grid = n_grill,N = 1000,conf = conf)
  results$Data$F2[[i]] = Z_new_i
  results$ARL[[i]] = Cntr_point(Z_new_i,Cntr_F1,type = 'H')
  ARL_list = append(ARL_list,list(which(results$ARL[[i]]$Cnt == 1)))
  # Objetos de las listas
  if (ic %% 100 == 0) {
    name_i = paste0('Resultados_C1/C1_T1_S_',ic/50,".rds")
    name_i2 = paste0('Resultados_C1/C1_T1_RL',".rds")
    cat("Guardando ",name_i,"....\n")
    saveRDS(results, file = name_i);saveRDS(ARL_list, file = name_i2)
    results$Data = list();results$Data$F1 = list();results$Data$F2 = list()
    results$Char = list();results$ARL_tab = list()
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
# [Calculo de los ARL C2] ===================================================.####
conf = c("PolyHigh","FourierLin")
# [Proceso fuera de control ] =============================================.####
fun_list_delta = list(f1 = function(x){3*x + x^2},
                      f2 = function(x){x + 3*x^2},
                      f3 = function(x){rep(0,length(x))})
Z_oc_i = Z_oc_funtion(Z_new,delta_k = 1,fun_list_delta)
Z_oc = Z_tau_oc(Z_new,tau = 30,mu_i = Z_oc_i)
autoplot(Z_oc)
# [Fin del código]===========================================================.#