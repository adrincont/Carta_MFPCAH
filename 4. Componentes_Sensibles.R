#.============================================================================.#
#                    Componentes principales sensibles                         #
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
# [Simulación de proceso] =================================================.####
conf = c("PolyHigh","FourierLin","Wiener")
Z_a = r_FunH(n_fun = 3,n_vec = 5,n_grid = 100,N = 800,conf = conf)
Z_b = r_FunH(n_fun = 3,n_vec = 5,n_grid = 100,N = 200,conf = conf)
# [Calculo de competentes sensibles] ======================================.####
## {Paso 1 limites de componentes sensibles} ------------------------------.#### 
Comp_sens_i = CP_Sensit_H(Z_a,delta = 0.99,omega = 0.90,wind = 0.25,alpha_rm = 0.99)
## {Paso 2 proyectar nuevas observaciones} --------------------------------.####
score_dat_i = pryc_pca_H(Z_b,Comp_sens_i$PCA)
scor_train = pryc_pca_H(Z_a,Comp_sens_i$PCA)
N = nrow(score_dat_i)
## {Paso 3 calculamos los Rm, Tm, CL} -------------------------------------.####
Cal_m = sapply(c(1:N),ID_sens,score_dat_j = score_dat_i, 
               CPs_Fj = Comp_sens_i,simplify = FALSE, type = 'H')
R_max = sapply(Cal_m, function(x){x$R_max})
## {Paso 4 calculamos los T_delta con los sensibles} ----------------------.####
Args_i = list(CPs_Fj = Comp_sens_i,score_dat_j = score_dat_i, scor_train_j = scor_train,alpha_j = 1/200)
T_calc_i = mapply(S_sens, as.list(1:N), Cal_m, MoreArgs = Args_i,type = 'H')
## {Paso 5 resultados fonales} --------------------------------------------.####
resultados_i = as.data.frame(cbind(R_max = c(R_max), t(T_calc_i)))
resultados_i$Cnt = as.numeric(unlist(resultados_i$T2) > unlist(resultados_i$LS))
# [Fin del codigo] ========================================================.####