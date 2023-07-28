#.============================================================================.#
#                 Simulación de observaciones híbridas                         #
#.============================================================================.#
# Autor: Andrey Duvan Rincon Torres    Fecha: 27/07/2023
# [Librerias] ================================================================.#
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
# [Codigo] ===================================================================.#
n_fun = 3;n_vec = 4;n_grid = 200;N = 100
conf = c("Wiener","PolyHigh","FourierLin")
Z = r_FunH(n_fun,n_vec,n_grid,N,conf = conf,dist_i = 'norm')
autoplot(Z)
# Correlacion entre diferentes puntos de tiempo de un mismo perfil
plot(Z@funData[[1]]@X[,100],Z@funData[[1]]@X[,80])
# Correlación entre mismo tiempo de diferentes perfiles
plot(Z@funData[[1]]@X[,161],Z@funData[[2]]@X[,161])
# Correlaciopn entre vectores
plot(Z@vecData[,1],Z@funData[[1]]@X[,1])
plot(Z@vecData[,1],Z@funData[[1]]@X[,65])
plot(Z@vecData[,2],Z@funData[[1]]@X[,20])
plot(Z@vecData[,3],Z@funData[[1]]@X[,10])
plot(Z@vecData[,4],Z@funData[[1]]@X[,90])
# [Fin del codigo] ========================================================.####