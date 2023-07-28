#.============================================================================.#
#                         Aplicaciones carta hibrida                           #
#.============================================================================.#
# Autor: Andrey Duvan Rincon Torres    Fecha: 27/07/2023
# [Librerias] ================================================================.#
source("library.R")
library('roahd')
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
## {datos maiz} -----------------------------------------------------------.####
### {datos} ................................................................####
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
## Análisis exploratorio datos
par(mfrow = c(3,1), mar = c(4,4,2,2))
fun_dat_h = data_h@funData[[1]]
roahd::fbplot(Data = roahd::fData(grid = fun_dat_h@argvals[[1]],
                                  values = fun_dat_h@X),
              Fvalue = 1.5,xlab = 'i (nm)',ylab = 'f1 (NIR M5)') 
fun_dat_h = data_h@funData[[2]]
roahd::fbplot(Data = roahd::fData(grid = fun_dat_h@argvals[[1]],
                                  values = fun_dat_h@X),
              Fvalue = 1.5,xlab = 'i (nm)',ylab = 'f2 (NIR MP5)') 
fun_dat_h = data_h@funData[[2]]
roahd::fbplot(Data = roahd::fData(grid = fun_dat_h@argvals[[1]],
                                  values = fun_dat_h@X),
              Fvalue = 1.5,xlab = 'i (nm)',ylab = 'f3 (NIR MP6)')
autoplot(data_h)$vec + theme_bw()
### {Fase I y Fase II} .....................................................####
set.seed(123)
index_I = sample(1:80,70)
index_II = c(1:80)[!(c(1:80) %in% index_I)]
data_I = select_Mfun(data_h,index_I)
data_II = select_Mfun(data_h,index_II)
### {Fase I } ..............................................................####
scale_data_h = scale_FunDataH_Train(data_I)
data_h_sc = scale_FunDataH(data_I,scale_obj = scale_data_h,center = TRUE,w = 1)
delta = 0.99;omega = 0.95;wind = 0.25;ARL_0 = 250
Cntr_F1 = Cntrl_sns_H(Z_i = data_h_sc,delta = delta,wind = wind, ARL_0 = ARL_0,
                      omega = omega)
# PCA hibrido
summary(Cntr_F1@CPs$PCA)
autoplot(Cntr_F1@CPs$PCA,type = 'eigen')$fun + plot_layout(ncol = 1)
autoplot(Cntr_F1@CPs$PCA) + theme_bw()
### {Fase II } .............................................................####
data_h_II_sc = scale_FunDataH(data_II,scale_obj = scale_data_h,center = TRUE,w = 1)
Cntr_h = Cntr_point(data_h_II_sc,Cntr_F1,type = 'H')
tab_2 = plot_chart_sns(Cntr_F1,Cntr_h,type = 'H')
p2 = autoplot(Cntr_F1) + theme_bw() + ggtitle('Aportes a componentes')
p2_1 = tab_2[[1]] + geom_segment(data = NULL,aes(x = -70,y = Cntr_F1@delta_st,
                                                 xend = 1, 
                                                 yend = Cntr_F1@delta_st),
                                 color = "green")
p2_2 = tab_2[[2]] + geom_segment(data = NULL,aes(x = -70,y = Cntr_F1@CPs$Cl,
                                                 xend = 1, 
                                                 yend = Cntr_F1@CPs$Cl),
                                 color = "green")
p2 + (p2_1/p2_2) 
## {datos aceite de oliva} ------------------------------------------------.####
### {datos} ................................................................####
data_3 = read_xlsx("Datos/Datos_3.xlsx")
dat_3 = list()
## Parte funcional
X1 = data_3[,-c(1:5)]
args_i = as.numeric(colnames(X1))
X1 = as.matrix(X1)
tempbasis = create.bspline.basis(c(min(args_i),max(args_i)),60)
fun_smoot = smooth.basis(argvals = args_i,y = t(X1),
                         fdParobj = tempbasis,dfscale = 1)
args_smoot = seq(min(args_i),max(args_i),len = 200)
X_smoot = eval.fd(evalarg = args_smoot,fdobj = fun_smoot$fd)
fun_0 = fdata(t(X_smoot),argvals = args_smoot)
fun_1 = fdata.deriv(fun_0)
fun_0 = funData(fun_0$argvals,fun_0$data)
fun_1 = funData(fun_1$argvals,fun_1$data)
dat_3$fmult = multiFunData(list(fun_0,fun_1))
## Parte vectorial
dat_3$vect = data_3[,c(2:5)]
## Crear objeto hibrido
data_h2 = mfh_data(fun_data = dat_3$fmult,vec_data = as.matrix(dat_3$vect))
## Análisis exploratorio datos
par(mfrow = c(2,1), mar = c(4,4,2,2))
fun_dat_h = data_h2@funData[[1]]
roahd::fbplot(Data = roahd::fData(grid = fun_dat_h@argvals[[1]],
                                  values = fun_dat_h@X),
              Fvalue = 1.5,xlab = 'i (nm)',ylab = 'f1') 
fun_dat_h = data_h2@funData[[2]]
roahd::fbplot(Data = roahd::fData(grid = fun_dat_h@argvals[[1]],
                                  values = fun_dat_h@X),
              Fvalue = 1.5,xlab = 'i (nm)',ylab = 'f2') 
autoplot(data_h2)$vec + theme_bw()
### {Fase I y Fase II} .....................................................####
set.seed(123)
index2_I = sample(1:123,100)
index2_II = c(1:123)[!(c(1:123) %in% index2_I)]
data2_I = select_Mfun(data_h2,index2_I)
data2_II = select_Mfun(data_h2,index2_II)
### {Fase I } ..............................................................####
scale_data_h2 = scale_FunDataH_Train(data2_I)
data_h2_sc = scale_FunDataH(data2_I,scale_obj = scale_data_h2,center = TRUE,w = 1)
delta = 0.99;omega = 0.95;wind = 0.25;ARL_0 = 200
Cntr_F2 = Cntrl_sns_H(Z_i = data_h2_sc,delta = delta,wind = wind, ARL_0 = ARL_0,
                      omega = omega)
# PCA hibrido
summary(Cntr_F2@CPs$PCA)
autoplot(Cntr_F2@CPs$PCA,type = 'eigen')$fun + plot_layout(ncol = 1)
autoplot(Cntr_F2@CPs$PCA) + theme_bw()
### {Fase II } .............................................................####
data_h2_II_sc = scale_FunDataH(data2_II,scale_obj = scale_data_h2,
                               center = TRUE,w = 1)
Cntr_h2 = Cntr_point(data_h2_II_sc,Cntr_F2,type = 'H')

tab_2 = plot_chart_sns(Cntr_F2,Cntr_h2,type = 'H')
p2 = autoplot(Cntr_F2) + theme_bw() + ggtitle('Aportes a componentes')
p2_1 = tab_2[[1]] + geom_segment(data = NULL,aes(x = -100,y = Cntr_F2@delta_st,
                                                 xend = 1, 
                                                 yend = Cntr_F2@delta_st),
                                 color = "green")
p2_2 = tab_2[[2]] + geom_segment(data = NULL,aes(x = -100,y = Cntr_F2@CPs$Cl,
                                                 xend = 1, 
                                                 yend = Cntr_F2@CPs$Cl),
                                 color = "green")
p2 + (p2_1/p2_2) 
# [Fin del codigo] ===========================================================.#