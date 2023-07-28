#.============================================================================.#
#       Construcción de tablas con salidas archivo simulaciones                #
#.============================================================================.#
# [Librerías] ================================================================.#
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
# [Lectura partes ARL0] ===================================================.####
Celd = 36;Tab_i = 1;ARL_list = list()
for (j in 1:4) {
  name_i = paste0("~/GitHub/MFPCAH/Codigo/Resultados/C",Celd,"/C",
                  Celd,"_T",Tab_i,'_RL_P',j,".rds")
  print(name_i)
  res_i = readRDS(name_i)
  ARL_list = append(ARL_list,res_i)
}
RL_i = sapply(ARL_list,function(x){x[1]})
RL_i[is.na(RL_i)] = 1050
ARL = mean(RL_i); SDRL = sd(RL_i)
lab_i = c(10,30,50,70,100,150,200,250,300,500,1000)
FAR = RL_fun(RL_i,lab_i)
names(FAR) = lab_i
cat('|Resultados -------------------------------------------------------|\n')
cat('ARL0 = ',round(ARL,3),'| ','SDRL = ',round(SDRL,3),'|\n')
cat('FAR \n');print(round(FAR,3));cat('\n')
cat('|------------------------------------------------------------------|\n')
name_i2 = paste0('Resultados/C',Celd,"/C",Celd,'_T',Tab_i,'_RL',".rds")
saveRDS(ARL_list, file = name_i2)
# [Lectura partes ARL1] ===================================================.####
Celd = 3;Tab_i = 3;ARL_list = list()
list_ch = c('I','II','III','IV')
for (name_i in list_ch) {
  ARL_list[[name_i]] = list()
  for (j in 1:4) {
    name_j2 = paste0('Resultados/C',Celd,'/C',name_i,'_',Celd,'_T',Tab_i,'_RL_P',j,".rds")
    print(name_j2)
    res_i = readRDS(name_j2)
    ARL_list[[name_i]] = append(ARL_list[[name_i]],res_i)
  }
}
for (name_i in list_ch) {
  ARL_list_i = ARL_list[[name_i]]
  RL_i = sapply(ARL_list_i,function(x){x[1]})
  RL_i[is.na(RL_i)] = 1050
  ARL = mean(RL_i); SDRL = sd(RL_i)
  lab_i = c(10,30,50,70,100,150,200,250,300,500,1000)
  FAR = RL_fun(RL_i,lab_i)
  names(FAR) = lab_i
  cat('|Resultados (',name_i,')-------------------------------------|\n')
  cat('ARL0 = ',round(ARL,3),'| ','SDRL = ',round(SDRL,3),'|\n')
  cat('FAR \n');print(round(FAR,3));cat('\n')
  cat('|------------------------------------------------------------------|\n')
}
for (name_i in list_ch) {
  ARL_list_i = ARL_list[[name_i]]
  name_j2 = paste0('Resultados/C',Celd,'/C',name_i,'_',Celd,'_T',Tab_i,'_RL.rds')
  print(name_j2)
  saveRDS(ARL_list_i, file = name_j2)
}
# [Fin del código]=========================================================.####