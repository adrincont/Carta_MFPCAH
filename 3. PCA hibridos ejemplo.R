# ============================================================================ #
#                           MFDA a datos híbridos                             .#
#.============================================================================.#
# Autor: Andrey Duvan Rincon Torres    Fecha: 27/07/2023
# [Librerias] ================================================================.#
source("library.R")
# [Funciones] ================================================================.#
source("Funciones/multiFunDataH.R") # Clase de objetos multiFunDataH
source("Funciones/multiFunDataH_plot.R") # Gráficos para la clase multiFunDataH
source("Funciones/Xi_gen.R") # Generar direcciones híbridas
source("Funciones/r_FunH.R") # Generar observaciones híbridas
source("Funciones/scale_FunDataH.R") # Estandarizar hibridos
source("Funciones/MFHPCA.R") # PCA híbrido
# {Simulacion de proceso} =================================================.####
conf = c("PolyHigh","FourierLin","Wiener")
Z = r_FunH(n_fun = 3,n_vec = 5,n_grid = 100,N = 800,conf = conf,dist_i = 'gamma')
autoplot(Z)
# {Calculo del PCA híbrido} ===============================================.####
PCAH = MFHPCA(Z,delta = 0.99,omega = 0.99,scale = FALSE)
# Resultados
summary(PCAH)
autoplot(PCAH,'var')
autoplot(PCAH,'eigen')
autoplot(PCAH,'ind')
# {Reconstruccion de los datos} ===========================================.####
Z_real = Z
scors_dat = pryc_pca_H(data_j = Z,PCA_j = PCAH)
Z_hat = Inv_PCA_H(dir_j = PCAH$MFPCAH$mfdaH$dir,scor_j = scors_dat)
autoplot(Z_real,1)$fun/autoplot(Z_hat,1)$fun
autoplot(Z_real,2)$vec;autoplot(Z_hat,2)$vec
Z_diff = Z_real + (-1)*Z_hat
autoplot(Z_diff)
# {Reconstruccion nuevos datos} ===========================================.####
Z_new = r_FunH(n_fun = 3,n_vec = 5,n_grid = 100,N = 500,conf = conf,dist_i = 'gamma')
scor_j = pryc_pca_H(data_j = Z_new,PCA_j = PCAH)
Z_hat = Inv_PCA_H(dir_j = PCAH$MFPCAH$mfdaH$dir,scor_j = scor_j)
autoplot(Z_new,1)$fun/autoplot(Z_hat,1)$fun
Z_diff = Z_new + (-1)*Z_hat
autoplot(Z_diff)
# [Fin del codigo] ========================================================.####