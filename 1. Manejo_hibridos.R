#.============================================================================.#
#                           Manejo de objetos hibridos                         #
#.============================================================================.#
# Autor: Andrey Duvan Rincon Torres    Fecha: 27/07/2023
# [Librerias] ================================================================.#
source("library.R")
# [Funciones] ================================================================.#
source("Funciones/multiFunDataH.R") # Clase de objetos multiFunDataH
source("Funciones/multiFunDataH_plot.R") # Graficos para la clase multiFunDataH
# [Codigo] ===================================================================.#
## {datos maiz} -----------------------------------------------------------.####
data_2 = readMat('Datos/Datos1.mat')
dat_2 = list()
## Parte funcional
X_1 = data_2$m5spec$data
X_2 = data_2$mp5spec$data
X_3 = data_2$mp6spec$data
dimnames(X_1)[[2]] = seq(1100,2498,2)
dimnames(X_2)[[2]] = seq(1100,2498,2)
dimnames(X_3)[[2]] = seq(1100,2498,2)
fun_1 = funData(seq(1100,2498,2),X_1)
fun_2 = funData(seq(1100,2498,2),X_2)
fun_3 = funData(seq(1100,2498,2),X_3)
dat_2$fmult = multiFunData(list(fun_1,fun_2,fun_3))
## Parte vectorial
dat_2$vect = data_2$propvals$data
colnames(dat_2$vect) = c("humedad","aceite","proteína","almidon")
## Crear objeto hibrido
data_h = mfh_data(fun_data = dat_2$fmult,vec_data = dat_2$vect)
## Grafico de un objeto híbrido
autoplot(data_h)
## {Aceites de oliva} -----------------------------------------------------.####
data_3 = read_xlsx("Datos/Datos2.xlsx")
dat_3 = list()
## Parte funcional
X1 = data_3[,-c(1:5)]
args_i = as.numeric(colnames(X1))
X1 = as.matrix(X1)
fun_0 = fdata(X1,argvals = args_i)
fun_1 = fdata.deriv(fun_0)
fun_0 = funData(fun_0$argvals,fun_0$data)
fun_1 = funData(fun_1$argvals,fun_1$data)
dat_3$fmult = multiFunData(list(fun_0,fun_1))
## Parte vectorial
dat_3$vect = data_3[,c(2:5)]
## Crear objeto hibrido
data_h2 = mfh_data(fun_data = dat_3$fmult,vec_data = as.matrix(dat_3$vect))
## Grafico de un objeto híbrido
autoplot(data_h2)
# [Fin del codigo] ========================================================.####