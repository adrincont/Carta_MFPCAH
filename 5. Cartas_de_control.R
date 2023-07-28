#.============================================================================.#
#         Carta para datos vectoriales y funcionales individuales              #
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
source("Funciones/scale_FunDataH.R") # Estandarizar hibridos
source("Funciones/MFHPCA.R") # PCA híbrido
source("Funciones/CP_Sensit.R") # Componentes sensibles
source("Funciones/Cntrl_sns.R") # Carta de control con componentes sensibles
source("Funciones/plot_charts.R") # Gráficos carta de control
source("Funciones/Cntrl_QT2.R") # Gráficos carta de control
# [Simulacion de proceso] =================================================.####
conf = c("PolyHigh","FourierLin","Wiener","PolyHigh","FourierLin")
Z = r_FunH(n_fun = 3,n_vec = 4,n_grid = 100,N = 650,conf = conf)
# [Carta de control de componentes sensibles (híbridos)]===================.####
## {Fase 1 limites de componentes sensibles} ------------------------------.####
Cntr_F1 = Cntrl_sns_H(Z,delta = 0.99,omega = 0.90,wind = 0.25, ARL_0 = 200)
autoplot(Cntr_F1)
## {Fase 2 de seguimiento} ------------------------------------------------.####
Z_new = r_FunH(n_fun = 3,n_vec = 4,n_grid = 100,N = 1000,conf = conf)
Tab_C1 = Cntr_point(Z_new,Cntr_F1,type = 'H')
which(Tab_C1$Cnt == 1)
plot_chart_sns(Cntr_F1,Tab_C1,type = 'H')
# [Carta de control de componentes sensibles (funcional)]==================.####
## {Fase 1 limites de componentes sensibles} ------------------------------.#### 
Fun_i = Z@funData
Cntr_F1_F = Cntrl_sns_F(Fun_i,delta = 0.90,wind = 0.25,ARL_0 = 200)
autoplot(Cntr_F1_F)
## {Fase 2 de seguimiento} ------------------------------------------------.####
Fun_i_New = Z_new@funData
Tab_C1_F = Cntr_point(Fun_i_New,Cntr_F1_F,type = 'F')
plot_chart_sns(Cntr_F1_F,Tab_C1_F,type = 'F')
# [Carta de control de componentes sensibles (vectorial)]==================.####
## {Fase 1 limites de componentes sensibles} ------------------------------.#### 
Vect = as.data.frame(Z@vecData)
Cntr_F1_V = Cntrl_sns_V(Vect,delta = 0.90,wind = 0.25,ARL_0 = 200)
autoplot(Cntr_F1_V)
## {Fase 2 de seguimiento} ------------------------------------------------.####
Vect_New = as.data.frame(Z_new@vecData)
Tab_C1_V = Cntr_point(Vect_New,Cntr_F1_V,type = 'V')
plot_chart_sns(Cntr_F1_V,Tab_C1_V,type = 'V')
# [Carta de control de componentes sensibles (hibridos separados)]=========.####
## {Fase 1 limites de componentes sensibles} ------------------------------.####
Cntr_F1_HS = Cntrl_sns_HS(Z,delta = 0.90,wind = 0.25,ARL_0 = 200)
## {Fase 2 de seguimiento} ------------------------------------------------.####
Tab_C1_HS = Cntr_point_HS(Z_new,Cntr_F1_HS)
which(Tab_C1_HS$Cnt != 0)[1]
# [Fin del código] ===========================================================.#